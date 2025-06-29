/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <CCA/Components/MPM/ConstitutiveModel/J2PlasticModels/UCNH.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/MPMEquationOfStateFactory.h>
#include <CCA/Components/MPM/GradientComputer/DisplacementGradientComputer.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Exceptions/ConvergenceFailure.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include<CCA/Components/MPM/Core/MPMLabel.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/Gaussian.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/Math/Weibull.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <fstream>
#include <iostream>

using namespace Uintah;
using namespace Vaango;

//#define Comer
#undef Comer

// Constructors //
//////////////////
UCNH::UCNH(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
  , ImplicitCM()
{
  d_useModifiedEOS = false;
  ps->require("bulk_modulus", d_initialData.Bulk);
  ps->require("shear_modulus", d_initialData.tauDev);
  ps->get("useModifiedEOS", d_useModifiedEOS);
  d_8or27 = Mflag->d_8or27;

  // Plasticity
  ps->getWithDefault("usePlasticity", d_usePlasticity, false);
  if (d_usePlasticity) {
    ps->getWithDefault("alpha", d_initialData.Alpha, 0.0);
    ps->require("yield_stress", d_initialData.FlowStress);
    ps->require("hardening_modulus", d_initialData.K);

    getYieldStressDistribution(ps);

    pPlasticStrain_label = VarLabel::create(
      "p.plasticStrain", ParticleVariable<double>::getTypeDescription());
    pPlasticStrain_label_preReloc = VarLabel::create(
      "p.plasticStrain+", ParticleVariable<double>::getTypeDescription());
    pYieldStress_label = VarLabel::create(
      "p.yieldStress", ParticleVariable<double>::getTypeDescription());
    pYieldStress_label_preReloc = VarLabel::create(
      "p.yieldStress+", ParticleVariable<double>::getTypeDescription());
  } // End Plasticity

  // Initial stress
  // Fix: Need to make it more general.  Add gravity turn-on option and
  //      read from file option etc.
  ps->getWithDefault("useInitialStress", d_useInitialStress, false);
  d_init_pressure = 0.0;
  if (d_useInitialStress) {
    ps->getWithDefault("initial_pressure", d_init_pressure, 0.0);
  }

  // Equation of state factory for pressure (default is DefaultHyperEOS)
  d_eos = MPMEquationOfStateFactory::create(ps);
  d_eos->setBulkModulus(d_initialData.Bulk);
  if (!d_eos) {
     std::ostringstream desc;
    desc << "An error occured in the MPM EquationOfStateFactory that has \n"
         << " slipped through the existing bullet proofing. Please check and "
            "correct."
         << "\n";
    throw ParameterNotFound(desc.str(), __FILE__, __LINE__);
  }

  // Universal Labels
  bElBarLabel = VarLabel::create(
    "p.bElBar", ParticleVariable<Matrix3>::getTypeDescription());
  bElBarLabel_preReloc = VarLabel::create(
    "p.bElBar+", ParticleVariable<Matrix3>::getTypeDescription());
  pDeformRateLabel = VarLabel::create(
    "p.deformRate", ParticleVariable<Matrix3>::getTypeDescription());
  pDeformRateLabel_preReloc = VarLabel::create(
    "p.deformRate+", ParticleVariable<Matrix3>::getTypeDescription());
}

UCNH::UCNH(ProblemSpecP& ps, MPMFlags* Mflag, bool plas, [[maybe_unused]] bool dam)
  : ConstitutiveModel(Mflag)
  , ImplicitCM()
{
  d_useModifiedEOS = false;
  ps->require("bulk_modulus", d_initialData.Bulk);
  ps->require("shear_modulus", d_initialData.tauDev);
  ps->get("useModifiedEOS", d_useModifiedEOS);
  d_8or27 = Mflag->d_8or27;

  // Plasticity
  ps->getWithDefault("usePlasticity", d_usePlasticity, plas);
  if (d_usePlasticity) {
    ps->getWithDefault("alpha", d_initialData.Alpha, 0.0);
    ps->require("yield_stress", d_initialData.FlowStress);
    ps->require("hardening_modulus", d_initialData.K);

    getYieldStressDistribution(ps);

    pPlasticStrain_label = VarLabel::create(
      "p.plasticStrain", ParticleVariable<double>::getTypeDescription());
    pPlasticStrain_label_preReloc = VarLabel::create(
      "p.plasticStrain+", ParticleVariable<double>::getTypeDescription());
    pYieldStress_label = VarLabel::create(
      "p.yieldStress", ParticleVariable<double>::getTypeDescription());
    pYieldStress_label_preReloc = VarLabel::create(
      "p.yieldStress+", ParticleVariable<double>::getTypeDescription());
  } // End Plasticity

  // Initial stress
  // Fix: Need to make it more general.  Add gravity turn-on option and
  //      read from file option etc.
  ps->getWithDefault("useInitialStress", d_useInitialStress, false);
  d_init_pressure = 0.0;
  if (d_useInitialStress) {
    ps->getWithDefault("initial_pressure", d_init_pressure, 0.0);
  }

  // Equation of state factory for pressure
  d_eos = MPMEquationOfStateFactory::create(ps);
  d_eos->setBulkModulus(d_initialData.Bulk);
  if (!d_eos) {
     std::ostringstream desc;
    desc << "An error occured in the MPM EquationOfStateFactory that has \n"
         << " slipped through the existing bullet proofing. Please check and "
            "correct."
         << "\n";
    throw ParameterNotFound(desc.str(), __FILE__, __LINE__);
  }

  // Universal Labels
  bElBarLabel = VarLabel::create(
    "p.bElBar", ParticleVariable<Matrix3>::getTypeDescription());
  bElBarLabel_preReloc = VarLabel::create(
    "p.bElBar+", ParticleVariable<Matrix3>::getTypeDescription());
  pDeformRateLabel = VarLabel::create(
    "p.deformRate", ParticleVariable<Matrix3>::getTypeDescription());
  pDeformRateLabel_preReloc = VarLabel::create(
    "p.deformRate+", ParticleVariable<Matrix3>::getTypeDescription());
}

UCNH::UCNH(const UCNH* cm)
  : ConstitutiveModel(cm)
  , ImplicitCM(cm)
{
  d_useModifiedEOS     = cm->d_useModifiedEOS;
  d_initialData.Bulk   = cm->d_initialData.Bulk;
  d_initialData.tauDev = cm->d_initialData.tauDev;

  // Plasticity Setup
  d_usePlasticity = cm->d_usePlasticity;
  if (d_usePlasticity) {
    d_initialData.FlowStress = cm->d_initialData.FlowStress;
    d_initialData.K          = cm->d_initialData.K;
    d_initialData.Alpha      = cm->d_initialData.Alpha;

    setYieldStressDistribution(cm);

    pPlasticStrain_label = VarLabel::create(
      "p.plasticStrain", ParticleVariable<double>::getTypeDescription());
    pPlasticStrain_label_preReloc = VarLabel::create(
      "p.plasticStrain+", ParticleVariable<double>::getTypeDescription());
    pYieldStress_label = VarLabel::create(
      "p.yieldStress", ParticleVariable<double>::getTypeDescription());
    pYieldStress_label_preReloc = VarLabel::create(
      "p.yieldStress+", ParticleVariable<double>::getTypeDescription());
  } // End Plasticity Setup

  // Initial stress
  d_useInitialStress = cm->d_useInitialStress;
  d_init_pressure    = cm->d_init_pressure;

  // EOS from factory
  d_eos = MPMEquationOfStateFactory::createCopy(cm->d_eos.get());
  d_eos->setBulkModulus(d_initialData.Bulk);

  // Universal Labels
  bElBarLabel = VarLabel::create(
    "p.bElBar", ParticleVariable<Matrix3>::getTypeDescription());
  bElBarLabel_preReloc = VarLabel::create(
    "p.bElBar+", ParticleVariable<Matrix3>::getTypeDescription());
  pDeformRateLabel = VarLabel::create(
    "p.deformRate", ParticleVariable<Matrix3>::getTypeDescription());
  pDeformRateLabel_preReloc = VarLabel::create(
    "p.deformRate+", ParticleVariable<Matrix3>::getTypeDescription());
}

void
UCNH::getYieldStressDistribution(ProblemSpecP& ps)
{
  d_yield.dist  = "constant";
  d_yield.range = 0.0; // yield stress = FlowStress +- range
  d_yield.seed  = 0;   // seed for distribution generator
  ps->getWithDefault("yield_distrib", d_yield.dist, "constant");
  //"constant", "uniform", "weibull" or "gauss" not implemented
  if (d_yield.dist == "uniform") {
    ps->require("yield_range", d_yield.range);
    ps->getWithDefault("yield_seed", d_yield.seed, 0);
  }
}

void
UCNH::setYieldStressDistribution(const UCNH* cm)
{
  d_yield.dist  = cm->d_yield.dist;
  d_yield.range = cm->d_yield.range;
  d_yield.seed  = cm->d_yield.seed;
}

void
UCNH::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "UCNH");
  }

  cm_ps->appendElement("bulk_modulus", d_initialData.Bulk);
  cm_ps->appendElement("shear_modulus", d_initialData.tauDev);
  cm_ps->appendElement("useModifiedEOS", d_useModifiedEOS);
  cm_ps->appendElement("usePlasticity", d_usePlasticity);

  // Plasticity
  if (d_usePlasticity) {
    cm_ps->appendElement("yield_stress", d_initialData.FlowStress);
    cm_ps->appendElement("hardening_modulus", d_initialData.K);
    cm_ps->appendElement("alpha", d_initialData.Alpha);
    cm_ps->appendElement("yield_distrib", d_yield.dist);
    if (d_yield.dist == "uniform") {
      cm_ps->appendElement("yield_range", d_yield.range);
      cm_ps->appendElement("yield_seed", d_yield.seed);
    }
  }

  cm_ps->appendElement("useInitialStress", d_useInitialStress);
  if (d_useInitialStress) {
    cm_ps->appendElement("initial_pressure", d_init_pressure);
  }

  // EOS from factory
  d_eos->outputProblemSpec(cm_ps);
}

std::unique_ptr<ConstitutiveModel>
UCNH::clone()
{
  return std::make_unique<UCNH>(this);
}

UCNH::~UCNH()
{
  // Plasticity Deletes
  if (d_usePlasticity) {
    VarLabel::destroy(pPlasticStrain_label);
    VarLabel::destroy(pPlasticStrain_label_preReloc);
    VarLabel::destroy(pYieldStress_label);
    VarLabel::destroy(pYieldStress_label_preReloc);
  }

  // Universal Deletes
  VarLabel::destroy(bElBarLabel);
  VarLabel::destroy(bElBarLabel_preReloc);
  VarLabel::destroy(pDeformRateLabel);
  VarLabel::destroy(pDeformRateLabel_preReloc);
}

void
UCNH::addParticleState(std::vector<const VarLabel*>& from,
                       std::vector<const VarLabel*>& to)
{
  // Add the local particle state data for this constitutive model.
  // Plasticity
  if (d_usePlasticity) {
    from.push_back(pPlasticStrain_label);
    from.push_back(pYieldStress_label);
    to.push_back(pPlasticStrain_label_preReloc);
    to.push_back(pYieldStress_label_preReloc);
  }

  // Universal
  from.push_back(bElBarLabel);
  to.push_back(bElBarLabel_preReloc);
  if (flag->d_integrator != MPMFlags::Implicit) {
    from.push_back(pDeformRateLabel);
    to.push_back(pDeformRateLabel_preReloc);
  }
}

/************************************
 * Initialize model
 */
void
UCNH::addInitialComputesAndRequires(Task* task,
                                    const MPMMaterial* matl,
                                    [[maybe_unused]] const PatchSet* patches) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  // Plasticity
  if (d_usePlasticity) {
    task->computes(pPlasticStrain_label, matlset);
    task->computes(pYieldStress_label, matlset);
  }

  // Universal
  task->computes(bElBarLabel, matlset);
  task->computes(pDeformRateLabel, matlset);
}

void
UCNH::initializeCMData(const Patch* patch,
                       const MPMMaterial* matl,
                       DataWarehouse* new_dw)
{
  // Put stuff in here to initialize each particle's
  // constitutive model parameters and deformationMeasure
  Matrix3 Identity;
  Identity.Identity();
  Matrix3 zero(0.0);

  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  if (flag->d_integrator == MPMFlags::Implicit)
    initSharedDataForImplicit(patch, matl, new_dw);
  else {
    // initSharedDataForExplicit(patch, matl, new_dw);
    ParticleVariable<double> pdTdt;
    ParticleVariable<Matrix3> pStress;
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel, pset);
    new_dw->allocateAndPut(pStress, lb->pStressLabel, pset);

    ParticleSubset::iterator iter = pset->begin();
    // Initial stress option
    if (!d_useInitialStress) {
      for (; iter != pset->end(); iter++) {
        particleIndex idx = *iter;
        pdTdt[idx]        = 0.0;
        pStress[idx]      = zero;
      }
    } else {
      ParticleVariable<double> pVolume;
      ParticleVariable<Matrix3> pDefGrad;
      new_dw->getModifiable(pVolume, lb->pVolumeLabel, pset);
      new_dw->getModifiable(pDefGrad, lb->pDefGradLabel, pset);

      double p = d_init_pressure;
      Matrix3 sigInit(p, 0.0, 0.0, 0.0, p, 0.0, 0.0, 0.0, p);
      double rho_orig = matl->getInitialDensity();
      double rho_cur  = d_eos->computeDensity(rho_orig, p);
      double diag     = cbrt(rho_cur / rho_orig);
      for (; iter != pset->end(); iter++) {
        particleIndex idx = *iter;
        pdTdt[idx]        = 0.0;
        pDefGrad[idx]     = Matrix3(diag, 0., 0., 0., diag, 0., 0., 0., diag);
        pStress[idx]      = sigInit;
        pVolume[idx] *= pDefGrad[idx].Determinant();
      }
    }
  }

  ParticleSubset::iterator iterUniv = pset->begin();
  ParticleSubset::iterator iterPlas = pset->begin();

  // Plasticity
  if (d_usePlasticity) {
    ParticleVariable<double> pPlasticStrain, pYieldStress;

    new_dw->allocateAndPut(pPlasticStrain, pPlasticStrain_label, pset);
    new_dw->allocateAndPut(pYieldStress, pYieldStress_label, pset);

    // std::cerr << "d_usePlasticity = " << d_usePlasticity << " dist = " <<
    // d_yield.dist
    //     << " range = " << d_yield.range << "\n";
    if (d_yield.dist == "uniform") {
      // Initialize a random number generator
      // Make the seed differ for each patch, otherwise each patch gets the
      // same set of random #s.
      int patchID              = patch->getID();
      int patch_div_32         = patchID / 32;
      patchID                  = patchID % 32;
      unsigned int unique_seed = ((d_yield.seed + patch_div_32 + 1) << patchID);
      auto randGen             = scinew MusilRNG(unique_seed);
      // std::cout << "   seed = " << unique_seed << " first rand = " << (*randGen)()
      // << "\n";
      for (; iterPlas != pset->end(); iterPlas++) {
        pPlasticStrain[*iterPlas] = d_initialData.Alpha;
        double rand               = (*randGen)();
        pYieldStress[*iterPlas] =
          d_initialData.FlowStress + (2 * rand - 1) * d_yield.range;
      }
      delete randGen;
    } else {
      for (; iterPlas != pset->end(); iterPlas++) {
        pPlasticStrain[*iterPlas] = d_initialData.Alpha;
        pYieldStress[*iterPlas]   = d_initialData.FlowStress;
      }
    }
  }

  // Universal
  ParticleVariable<Matrix3> pstress, bElBar, pDeformRate;

  new_dw->allocateAndPut(pDeformRate, pDeformRateLabel, pset);
  new_dw->allocateAndPut(bElBar, bElBarLabel, pset);

  for (; iterUniv != pset->end(); iterUniv++) {
    bElBar[*iterUniv]      = Identity;
    pDeformRate[*iterUniv] = zero;
  }

  // If not implicit, compute timestep
  if (!(flag->d_integrator == MPMFlags::Implicit)) {
    // End by computing the stable timestep
    computeStableTimestep(patch, matl, new_dw);
  }
}

void
UCNH::computeStableTimestep(const Patch* patch,
                            const MPMMaterial* matl,
                            DataWarehouse* new_dw)
{
  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of computeStressTensor
  Vector dx = patch->dCell();
  int dwi   = matl->getDWIndex();
  // Retrieve the array of constitutive parameters
  ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
  constParticleVariable<double> pMass, pVolume;
  constParticleVariable<Vector> pVelocity;

  new_dw->get(pMass, lb->pMassLabel, pset);
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pVelocity, lb->pVelocityLabel, pset);

  double c_dil = 0.0;
  Vector WaveSpeed(1.e-12, 1.e-12, 1.e-12);

  double mu   = d_initialData.tauDev;
  double bulk = d_initialData.Bulk;

  for (int idx : *pset) {
    // Compute wave speed at each particle, store the maximum
    Vector pVelocity_idx = pVelocity[idx];
    if (pMass[idx] > 0) {
      c_dil = sqrt((bulk + 4. * mu / 3.) * pVolume[idx] / pMass[idx]);
    } else {
      c_dil         = 0.0;
      pVelocity_idx = Vector(0.0, 0.0, 0.0);
    }
    WaveSpeed = Vector(Max(c_dil + std::abs(pVelocity_idx.x()), WaveSpeed.x()),
                       Max(c_dil + std::abs(pVelocity_idx.y()), WaveSpeed.y()),
                       Max(c_dil + std::abs(pVelocity_idx.z()), WaveSpeed.z()));
  }

  WaveSpeed       = dx / WaveSpeed;
  double delT_new = WaveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

/*****************************************
 * Explicit stress update
 */
void
UCNH::addComputesAndRequires(Task* task,
                             const MPMMaterial* matl,
                             const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  if (flag->d_integrator == MPMFlags::Implicit) {
    bool reset = flag->d_doGridReset;
    addSharedCRForImplicit(task, matlset, reset);
  } else {
    addSharedCRForExplicit(task, matlset, patches);
  }

  // Other constitutive model and input dependent computes and requires
  Ghost::GhostType gnone = Ghost::None;

  // Plasticity
  if (d_usePlasticity) {
    task->needs(Task::OldDW, pPlasticStrain_label, matlset, gnone);
    task->needs(Task::OldDW, pYieldStress_label, matlset, gnone);
    task->computes(pPlasticStrain_label_preReloc, matlset);
    task->computes(pYieldStress_label_preReloc, matlset);
  }

  // for pParticleID
  task->needs(Task::OldDW, lb->pParticleIDLabel, matlset, gnone);

  if (flag->d_withColor) {
    task->needs(Task::OldDW, lb->pColorLabel, Ghost::None);
  }

  // Universal
  task->needs(Task::OldDW, bElBarLabel, matlset, gnone);
  task->computes(bElBarLabel_preReloc, matlset);
  task->computes(pDeformRateLabel_preReloc, matlset);
}

void
UCNH::computeStressTensor(const PatchSubset* patches,
                          const MPMMaterial* matl,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw)
{
  // Constants
  double onethird = (1.0 / 3.0), sqtwthds = sqrt(2.0 / 3.0);
  Matrix3 Identity;
  Identity.Identity();

  // Grab initial data
  double shear    = d_initialData.tauDev;
  double bulk     = d_initialData.Bulk;
  double rho_orig = matl->getInitialDensity();
  double flow     = 0.0;
  double K        = 0.0;

  // Get delT
  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));
  // Ghost::GhostType  gac   = Ghost::AroundCells;

  // Normal patch loop
  for (int pp = 0; pp < patches->size(); pp++) {
    const Patch* patch = patches->get(pp);

    // Temporary and "get" variables
    double delgamma = 0.0, fTrial = 0.0, IEl = 0.0, J = 0.0;
    double muBar = 0.0, p = 0.0, sTnorm = 0.0, U = 0.0, W = 0.0;
    double se    = 0.0; // Strain energy placeholder
    double c_dil = 0.0; // Speed of sound
    Matrix3 pBBar_new(0.0), bEB_new(0.0), bElBarTrial(0.0), pDefGradInc(0.0);
    Matrix3 displacementGradient(0.0), fBar(0.0), defGrad(0.0), normal(0.0);
    Matrix3 tauDev(0.0), tauDevTrial(0.0);
    Vector WaveSpeed(1.e-12, 1.e-12, 1.e-12);

    // Get particle info and patch info
    int dwi              = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);
    Vector dx            = patch->dCell();
    // double time = d_mat_manager->getElapsedTime();

    // Get Interpolator
    auto interpolator = flag->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());
    std::vector<double> S(interpolator->size());

    // Particle and grid data universal to model type
    // Old data containers
    constParticleVariable<Short27> pgCode;
    constParticleVariable<double> pMass;
    constParticleVariable<double> pPlasticStrain_old, pYieldStress_old, pcolor;
    constParticleVariable<long64> pParticleID;
    constParticleVariable<Point> px;
    constParticleVariable<Matrix3> pDefGrad, bElBar, pVelGrad;
    constParticleVariable<Matrix3> pSize;
    constParticleVariable<Vector> pVelocity;

    constParticleVariable<double> pVolume_new;
    constParticleVariable<Matrix3> pDefGrad_new;
    constParticleVariable<Matrix3> pVelGrad_new;

    // New data containers
    ParticleVariable<double> pPlasticStrain, pYieldStress;
    ParticleVariable<double> pdTdt, p_q;
    ParticleVariable<Matrix3> pDeformRate;
    ParticleVariable<Matrix3> pStress, bElBar_new;

    // Particle and grid data
    // double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};

    // Plasticity gets
    if (d_usePlasticity) {
      old_dw->get(pPlasticStrain_old, pPlasticStrain_label, pset);
      old_dw->get(pYieldStress_old, pYieldStress_label, pset);
      new_dw->allocateAndPut(
        pPlasticStrain, pPlasticStrain_label_preReloc, pset);
      new_dw->allocateAndPut(pYieldStress, pYieldStress_label_preReloc, pset);

      pPlasticStrain.copyData(pPlasticStrain_old);
      pYieldStress.copyData(pYieldStress_old);

      // Copy initial data
      flow = d_initialData.FlowStress;
      K    = d_initialData.K;
    }

    // Universal Gets
    old_dw->get(px, lb->pXLabel, pset);
    old_dw->get(pMass, lb->pMassLabel, pset);
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);
    old_dw->get(pSize, lb->pSizeLabel, pset);
    old_dw->get(bElBar, bElBarLabel, pset);
    old_dw->get(pVelGrad, lb->pVelGradLabel, pset);
    old_dw->get(pDefGrad, lb->pDefGradLabel, pset);

    new_dw->get(pVolume_new, lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pVelGrad_new, lb->pVelGradLabel_preReloc, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    // Universal Allocations
    new_dw->allocateAndPut(bElBar_new, bElBarLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress, lb->pStressLabel_preReloc, pset);
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);
    new_dw->allocateAndPut(pDeformRate, pDeformRateLabel_preReloc, pset);

    if (flag->d_withColor) {
      old_dw->get(pcolor, lb->pColorLabel, pset);
    }

    ParticleSubset::iterator iter = pset->begin();
    for (; iter != pset->end(); iter++) {
      particleIndex idx = *iter;

      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;

#ifdef Comer
      // gcd change to set shear = pcolor for each particle
      if (flag->d_with_color) {
        shear = pcolor[idx];
      }
#endif

      pDeformRate[idx] =
        (pVelGrad_new[idx] + pVelGrad_new[idx].Transpose()) * 0.5;
      defGrad     = pDefGrad_new[idx];
      pDefGradInc = pDefGrad_new[idx] * pDefGrad[idx].Inverse();

      // 1) Get the volumetric part of the deformation
      // 2) Compute the deformed volume and new density
      J              = defGrad.Determinant();
      double rho_cur = rho_orig / J;

      // Get the volume preserving part of the deformation gradient increment
      //      fBar = pDefGradInc*pow(Jinc, -onethird);
      double Jinc = pDefGradInc.Determinant();
      fBar        = pDefGradInc / cbrt(Jinc);

      // Compute the trial elastic part of the volume preserving
      // part of the left Cauchy-Green deformation tensor
      bElBarTrial = fBar * bElBar[idx] * fBar.Transpose();
      if (!d_usePlasticity) {
        double cubeRootJ       = cbrt(J);
        double Jtothetwothirds = cubeRootJ * cubeRootJ;
        bElBarTrial =
          pDefGrad_new[idx] * pDefGrad_new[idx].Transpose() / Jtothetwothirds;
      }
      IEl   = onethird * bElBarTrial.Trace();
      muBar = IEl * shear;

      // tauDevTrial is equal to the shear modulus times dev(bElBar)
      // Compute ||tauDevTrial||
      tauDevTrial = (bElBarTrial - Identity * IEl) * shear;
      sTnorm      = tauDevTrial.Norm();

      // Check for plastic loading
      double alpha;
      if (d_usePlasticity) {
        flow   = pYieldStress[idx];
        alpha  = pPlasticStrain[idx];
        fTrial = sTnorm - sqtwthds * (K * alpha + flow);
      }
      if (d_usePlasticity && (fTrial > 0.0)) {
        // plastic
        // Compute increment of slip in the direction of flow
        delgamma = (fTrial / (2.0 * muBar)) / (1.0 + (K / (3.0 * muBar)));
        normal   = tauDevTrial / sTnorm;

        // The actual shear stress
        tauDev = tauDevTrial - normal * 2.0 * muBar * delgamma;

        // Deal with history variables
        pPlasticStrain[idx] = alpha + sqtwthds * delgamma;
        bElBar_new[idx]     = tauDev / shear + Identity * IEl;
      } else {
        // The actual shear stress
        tauDev          = tauDevTrial;
        bElBar_new[idx] = bElBarTrial;
      }

      // get the hydrostatic part of the stress
      p = d_eos->computePressure(rho_orig, rho_cur);
      // p = 0.5*bulk*(J - 1.0/J);

      // compute the total stress (volumetric + deviatoric)
      pStress[idx] = Identity * p + tauDev / J;
      if (std::isnan(pStress[idx].Norm())) {
        std::cerr << "particle = " << idx << " velGrad = " << pVelGrad[idx] << "\n";
        std::cerr << " stress = " << pStress[idx] << "\n";
        std::cerr << " pMass = " << pMass[idx] << " pvol = " << pVolume_new[idx]
             << "\n";
        std::cerr << " delgamma = " << delgamma << " ftrial = " << fTrial
             << " mubar = " << muBar << " K = " << K << "\n";
        std::cerr << " fBar = " << fBar << " Jinc = " << Jinc
             << " pDefGradInc = " << pDefGradInc << "\n";
        std::cerr << " IEl = " << IEl << " bElBarTrial = " << bElBarTrial
             << " shear = " << shear << "\n";
        std::cerr << " normal = " << normal << " tauDevTrial = " << tauDevTrial
             << " sTnorm = " << sTnorm << "\n";
        std::cerr << " p = " << p << " rho_orig = " << rho_orig
             << " rho_cur = " << rho_cur << "\n";
        throw InvalidValue(
          "**UCNH ERROR**: Nan in stress value", __FILE__, __LINE__);
      }

      // Compute the strain energy for non-localized particles
      U    = d_eos->computeStrainEnergy(rho_orig, rho_cur);
      bulk = d_eos->computeBulkModulus(rho_orig, rho_cur);
      //  U = .5*bulk*(.5*(J*J - 1.0) - log(J));

      W        = .5 * shear * (bElBar_new[idx].Trace() - 3.0);
      double e = (U + W) * pVolume_new[idx] / J;
      se += e;

      // Compute the local sound speed (uniaxial strain, p-wave modulus)
      c_dil = sqrt((bulk + 4. * shear / 3.) / rho_cur);

      // Compute wave speed at each particle, store the maximum
      Vector pvel = pVelocity[idx];
      WaveSpeed   = Vector(Max(c_dil + std::abs(pvel.x()), WaveSpeed.x()),
                         Max(c_dil + std::abs(pvel.y()), WaveSpeed.y()),
                         Max(c_dil + std::abs(pvel.z()), WaveSpeed.z()));

      // Compute artificial viscosity term
      if (flag->d_artificialViscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;
        double c_bulk = sqrt(bulk / rho_cur);
        p_q[idx]      = artificialBulkViscosity(
          pDeformRate[idx].Trace(), c_bulk, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.;
      }
    } // end loop over particles

    WaveSpeed       = dx / WaveSpeed;
    double delT_new = WaveSpeed.minComponent();

    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(se), lb->StrainEnergyLabel);
    }

    // delete interpolator;
    // std::cout << "End compute stress." << "\n";
  }
}

/**********************************************
 * Implicit stres update
 */
void
UCNH::addComputesAndRequires(Task* task,
                             const MPMMaterial* matl,
                             [[maybe_unused]] const PatchSet* patches,
                             [[maybe_unused]] const bool recurse,
                             const bool SchedParent) const
{
  const MaterialSubset* matlset = matl->thisMaterial();

  if (flag->d_integrator == MPMFlags::Implicit) {
    bool reset = flag->d_doGridReset;
    addSharedCRForImplicit(task, matlset, reset, true, SchedParent);
  }

  Ghost::GhostType gnone = Ghost::None;
  if (d_usePlasticity) {
    if (SchedParent) {
      task->needs(Task::ParentOldDW, pPlasticStrain_label, matlset, gnone);
    } else {
      task->needs(Task::OldDW, pPlasticStrain_label, matlset, gnone);
    }
  }

  if (SchedParent) {
    task->needs(Task::ParentOldDW, bElBarLabel, matlset, gnone);
  } else {
    task->needs(Task::OldDW, bElBarLabel, matlset, gnone);
  }
}
void
UCNH::computeStressTensorImplicit(const PatchSubset* patches,
                                  const MPMMaterial* matl,
                                  [[maybe_unused]] DataWarehouse* old_dw,
                                  DataWarehouse* new_dw,
                                  Solver* solver,
                                  const bool)

{
  // Constants
  int dwi         = matl->getDWIndex();
  double onethird = (1.0 / 3.0);
  double shear    = d_initialData.tauDev;
  double bulk     = d_initialData.Bulk;
  double rho_orig = matl->getInitialDensity();

  // Ghost::GhostType gac = Ghost::AroundCells;
  Matrix3 Identity;
  Identity.Identity();
  DataWarehouse* parent_old_dw =
    new_dw->getOtherDataWarehouse(Task::ParentOldDW);

  // Particle and grid variables
  constParticleVariable<double> pVol, pMass, pVolumeold;
  constParticleVariable<Point> px;
  constParticleVariable<Matrix3> pSize;
  constParticleVariable<Matrix3> pDefGrad, pBeBar;
  constParticleVariable<double> pVolume_new;
  constParticleVariable<Matrix3> pDefGrad_new;
  ParticleVariable<Matrix3> pBeBar_new, pStress;

  // Local variables
  Matrix3 tauDev(0.0), pDefGradInc(0.0), pDispGrad(0.0), pRelDefGradBar(0.0);
  double D[6][6];
  double B[6][24];
  double Bnl[3][24];
  // Unused because not using computeStiffnessMatrix() as in CNHPDamage
  //     double Kmatrix[24][24];
  int dof[24];
  // Unused because each 8 and 27 option have their owndouble v[576];

  IntVector lowIndex = IntVector(0, 0, 0), highIndex = IntVector(0, 0, 0);
  for (int pp = 0; pp < patches->size(); pp++) {
    const Patch* patch = patches->get(pp);

    if (d_8or27 == 8) {
      lowIndex  = patch->getNodeLowIndex();
      highIndex = patch->getNodeHighIndex() + IntVector(1, 1, 1);
    } else if (d_8or27 == 27) {
      lowIndex  = patch->getExtraNodeLowIndex();
      highIndex = patch->getExtraNodeHighIndex() + IntVector(1, 1, 1);
    }

    Array3<int> l2g(lowIndex, highIndex);
    solver->copyL2G(l2g, patch);

    Vector dx      = patch->dCell();
    double oodx[3] = { 1. / dx.x(), 1. / dx.y(), 1. / dx.z() };

    ParticleSubset* pset = parent_old_dw->getParticleSubset(dwi, patch);
    parent_old_dw->get(px, lb->pXLabel, pset);
    parent_old_dw->get(pSize, lb->pSizeLabel, pset);
    parent_old_dw->get(pMass, lb->pMassLabel, pset);
    parent_old_dw->get(pVolumeold, lb->pVolumeLabel, pset);
    parent_old_dw->get(pDefGrad, lb->pDefGradLabel, pset);
    parent_old_dw->get(pBeBar, bElBarLabel, pset);

    new_dw->get(pVolume_new, lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    new_dw->allocateAndPut(pStress, lb->pStressLabel_preReloc, pset);
    new_dw->allocateTemporary(pBeBar_new, pset);

    double volold = 0.0, volnew = 0.0;
    if (matl->getIsRigid()) { // Rigid test
      for (auto idx : *pset) {
        pStress[idx] = Matrix3(0.0);
      }
    } else {
      auto interpolator = flag->d_interpolator->clone(patch);
      std::vector<IntVector> ni(interpolator->size());
      std::vector<Vector> d_S(interpolator->size());

      for (auto idx : *pset) {

        // Compute the displacement gradient and B matrices
        if (d_usePlasticity) {

          interpolator->findCellAndShapeDerivatives(
            px[idx], ni, d_S, pSize[idx], pDefGrad[idx]);

          DisplacementGradientComputer dg(flag);
          dg.computeBmats(ni, d_S, oodx, l2g, B, Bnl, dof);
        }

        // Compute the deformation gradient increment using the pDispGrad
        // Update the deformation gradient tensor to its time n+1 value.
        double J    = pDefGrad_new[idx].Determinant();
        pDefGradInc = pDefGrad_new[idx] * pDefGrad[idx].Inverse();
        if (d_usePlasticity) {
          // Compute BeBar
          pRelDefGradBar = pDefGradInc / cbrt(pDefGradInc.Determinant());

          pBeBar_new[idx] =
            pRelDefGradBar * pBeBar[idx] * pRelDefGradBar.Transpose();
        } else {
          Matrix3 bElBar_new = pDefGrad_new[idx] *
                               pDefGrad_new[idx].Transpose() *
                               pow(J, -(2. / 3.));
          pBeBar_new[idx] = bElBar_new;
        }

        // Update the particle volume
        volold = (pMass[idx] / rho_orig);

        // tauDev is equal to the shear modulus times dev(bElBar)
        double mubar   = onethird * pBeBar_new[idx].Trace() * shear;
        Matrix3 shrTrl = (pBeBar_new[idx] * shear - Identity * mubar);

        // get the hydrostatic part of the stress
        double p = bulk * log(J) / J;

        // compute the total stress (volumetric + deviatoric)
        pStress[idx] = Identity * p + shrTrl / J;

        // Compute the tangent stiffness matrix
        computeTangentStiffnessMatrix(shrTrl, mubar, J, bulk, D);

        double sig[3][3];
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            sig[i][j] = pStress[idx](i, j);
          }
        }

        int nDOF = 3 * d_8or27;

        if (d_8or27 == 8) {

          double B[6][24];
          double Bnl[3][24];
          int dof[24];
          double v[24 * 24];
          double kmat[24][24];
          double kgeo[24][24];

          // Fill in the B and Bnl matrices and the dof vector
          interpolator->findCellAndShapeDerivatives(
            px[idx], ni, d_S, pSize[idx], pDefGrad[idx]);
          loadBMats(l2g, dof, B, Bnl, d_S, ni, oodx);
          // kmat = B.transpose()*D*B*volold
          BtDB(B, D, kmat);
          // kgeo = Bnl.transpose*sig*Bnl*volnew;
          BnltDBnl(Bnl, sig, kgeo);

          for (int I = 0; I < nDOF; I++) {
            for (int J = 0; J < nDOF; J++) {
              v[nDOF * I + J] = kmat[I][J] * volold + kgeo[I][J] * volnew;
            }
          }
          solver->fillMatrix(nDOF, dof, nDOF, dof, v);
        } else {
          double B[6][81];
          double Bnl[3][81];
          int dof[81];
          double v[81 * 81];
          double kmat[81][81];
          double kgeo[81][81];

          // the code that computes kmat doesn't yet know that D is symmetric
          D[1][0] = D[0][1];
          D[2][0] = D[0][2];
          D[3][0] = D[0][3];
          D[4][0] = D[0][4];
          D[5][0] = D[0][5];
          D[1][1] = D[1][1];
          D[2][1] = D[1][2];
          D[3][1] = D[1][3];
          D[4][1] = D[1][4];
          D[1][2] = D[2][1];
          D[2][2] = D[2][2];
          D[3][2] = D[2][3];
          D[1][3] = D[3][1];
          D[2][3] = D[3][2];
          D[4][3] = D[3][4];

          // Fill in the B and Bnl matrices and the dof vector
          interpolator->findCellAndShapeDerivatives(
            px[idx], ni, d_S, pSize[idx], pDefGrad[idx]);
          loadBMatsGIMP(l2g, dof, B, Bnl, d_S, ni, oodx);
          // kmat = B.transpose()*D*B*volold
          BtDBGIMP(B, D, kmat);
          // kgeo = Bnl.transpose*sig*Bnl*volnew;
          BnltDBnlGIMP(Bnl, sig, kgeo);

          for (int I = 0; I < nDOF; I++) {
            for (int J = 0; J < nDOF; J++) {
              v[nDOF * I + J] = kmat[I][J] * volold + kgeo[I][J] * volnew;
            }
          }
          solver->fillMatrix(nDOF, dof, nDOF, dof, v);
        }
      }
      // delete interpolator;
    } // end rigid
  }   // end of loop over particles

  solver->flushMatrix();
}

void
UCNH::allocateCMDataAddRequires(Task* task,
                                const MPMMaterial* matl,
                                const PatchSet* patches,
                                MPMLabel* lb) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  Ghost::GhostType gnone        = Ghost::None;

  // Allocate the variables shared by all constitutive models
  // for the particle convert operation
  // This method is defined in the ConstitutiveModel base class.

  // Add requires local to this model
  // Plasticity
  if (d_usePlasticity) {
    task->needs(Task::NewDW, pPlasticStrain_label_preReloc, matlset, gnone);
    task->needs(Task::NewDW, pYieldStress_label_preReloc, matlset, gnone);
  }

  // Universal
  task->needs(Task::NewDW, bElBarLabel_preReloc, matlset, gnone);
  if (flag->d_integrator != MPMFlags::Implicit) { // non implicit
    addSharedRForConvertExplicit(task, matlset, patches);
    task->needs(Task::NewDW, pDeformRateLabel_preReloc, matlset, gnone);
  } else { // Implicit only stuff
    task->needs(Task::NewDW, lb->pStressLabel_preReloc, matlset, gnone);
  }
}

void
UCNH::allocateCMDataAdd(DataWarehouse* new_dw,
                        ParticleSubset* addset,
                        ParticleLabelVariableMap* newState,
                        ParticleSubset* delset,
                        [[maybe_unused]] DataWarehouse* old_dw)
{
  // Copy the data common to all constitutive models from the particle to be
  // deleted to the particle to be added.
  // This method is defined in the ConstitutiveModel base class.

  if (flag->d_integrator != MPMFlags::Implicit) {
    copyDelToAddSetForConvertExplicit(new_dw, delset, addset, newState);
  } else { // Implicit
    ParticleVariable<Matrix3> pstress;
    new_dw->allocateTemporary(pstress, addset);

    constParticleVariable<Matrix3> o_stress;
    new_dw->get(o_stress, lb->pStressLabel_preReloc, delset);

    ParticleSubset::iterator o, n = addset->begin();
    for (o = delset->begin(); o != delset->end(); o++, n++) {
      pstress[*n] = o_stress[*o];
    }

    (*newState)[lb->pStressLabel] = pstress.clone();
  }

  // Copy the data local to this constitutive model from the particles to
  // be deleted to the particles to be added
  ParticleSubset::iterator nPlas = addset->begin();
  ParticleSubset::iterator nUniv = addset->begin();

  // Plasticity
  if (d_usePlasticity) {
    ParticleVariable<double> pPlasticStrain;
    ParticleVariable<double> pYieldStress;
    new_dw->allocateTemporary(pPlasticStrain, addset);
    new_dw->allocateTemporary(pYieldStress, addset);

    constParticleVariable<double> o_pPlasticStrain;
    constParticleVariable<double> o_pYieldStress;
    new_dw->get(o_pPlasticStrain, pPlasticStrain_label_preReloc, delset);
    new_dw->get(o_pYieldStress, pYieldStress_label_preReloc, delset);

    ParticleSubset::iterator o;
    for (o = delset->begin(); o != delset->end(); o++, nPlas++) {
      pPlasticStrain[*nPlas] = o_pPlasticStrain[*o];
      pYieldStress[*nPlas]   = o_pYieldStress[*o];
    }

    (*newState)[pPlasticStrain_label] = pPlasticStrain.clone();
    (*newState)[pYieldStress_label]   = pYieldStress.clone();
  } // End Plasticity

  // Universal
  ParticleVariable<Matrix3> bElBar;
  ParticleVariable<Matrix3> pDeformRate;
  new_dw->allocateTemporary(pDeformRate, addset);
  new_dw->allocateTemporary(bElBar, addset);

  constParticleVariable<Matrix3> o_bElBar;
  constParticleVariable<Matrix3> o_pDeformRate;
  constParticleVariable<Matrix3> o_pVelGrad;
  new_dw->get(o_bElBar, bElBarLabel_preReloc, delset);
  new_dw->get(o_pDeformRate, pDeformRateLabel_preReloc, delset);

  ParticleSubset::iterator o;
  for (o = delset->begin(); o != delset->end(); o++, nUniv++) {
    bElBar[*nUniv]      = o_bElBar[*o];
    pDeformRate[*nUniv] = o_pDeformRate[*o];
  }

  (*newState)[bElBarLabel]      = bElBar.clone();
  (*newState)[pDeformRateLabel] = pDeformRate.clone();
}

void
UCNH::carryForward(const PatchSubset* patches,
                   const MPMMaterial* matl,
                   DataWarehouse* old_dw,
                   DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch   = patches->get(p);
    int dwi              = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

    // Carry forward the data common to all constitutive models
    // when using RigidMPM.
    // This method is defined in the ConstitutiveModel base class.
    carryForwardSharedData(pset, old_dw, new_dw, matl);

    // Carry forward the data local to this constitutive model
    // Universal
    ParticleVariable<Matrix3> bElBar_new;
    constParticleVariable<Matrix3> bElBar;
    old_dw->get(bElBar, bElBarLabel, pset);
    new_dw->allocateAndPut(bElBar_new, bElBarLabel_preReloc, pset);
    for (int idx : *pset) {
      bElBar_new[idx] = bElBar[idx];
    }
    if (flag->d_integrator != MPMFlags::Implicit) {

      ParticleVariable<Matrix3> pDeformRate;
      new_dw->allocateAndPut(pDeformRate, pDeformRateLabel_preReloc, pset);
      pDeformRate.copyData(bElBar);
    }

    // Plasticity
    if (d_usePlasticity) {
      ParticleVariable<double> pPlasticStrain, pYieldStress;
      constParticleVariable<double> pPlasticStrain_old, pYieldStress_old;
      old_dw->get(pPlasticStrain_old, pPlasticStrain_label, pset);
      old_dw->get(pYieldStress_old, pYieldStress_label, pset);
      new_dw->allocateAndPut(
        pPlasticStrain, pPlasticStrain_label_preReloc, pset);
      new_dw->allocateAndPut(pYieldStress, pYieldStress_label_preReloc, pset);
      pPlasticStrain.copyData(pPlasticStrain_old);
      pYieldStress.copyData(pYieldStress_old);
    } // End Plasticity

    new_dw->put(delt_vartype(1.e10), lb->delTLabel, patch->getLevel());
    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(0.), lb->StrainEnergyLabel);
    }
  } // End Particle Loop
}

void
UCNH::computeStressTensorImplicit(const PatchSubset* patches,
                                  const MPMMaterial* matl,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw)
{
  // Constants
  double onethird = (1.0 / 3.0);
  double sqtwthds = sqrt(2.0 / 3.0);
  Matrix3 Identity;
  Identity.Identity();
  // Ghost::GhostType gac = Ghost::AroundCells;

  // double rho_orig    = matl->getInitialDensity();
  double shear       = d_initialData.tauDev;
  double bulk        = d_initialData.Bulk;
  double flowStress  = d_initialData.FlowStress;
  double hardModulus = d_initialData.K;
  double se          = 0.0;

  int dwi = matl->getDWIndex();

  // Particle and grid data
  constParticleVariable<double> pMass, pPlasticStrain;
  constParticleVariable<long64> pParticleID;
  constParticleVariable<Point> pX;
  constParticleVariable<Matrix3> pSize;
  constParticleVariable<Matrix3> pDefGrad, pBeBar;
  constParticleVariable<double> pVolume_new;
  constParticleVariable<Matrix3> pDefGrad_new;
  ParticleVariable<double> pdTdt, pPlasticStrain_new;
  ParticleVariable<Matrix3> pBeBar_new, pStress_new;

  // Local variables
  Matrix3 dispGrad(0.0), tauDev(0.0), defGradInc(0.0);
  Matrix3 beBarTrial(0.0), tauDevTrial(0.0), normal(0.0), relDefGradBar(0.0);
  Matrix3 defGrad(0.0);

  // Loop thru patches
  for (int pp = 0; pp < patches->size(); pp++) {
    const Patch* patch = patches->get(pp);

    // Get particle info
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);
    // Loop thru particles
    ParticleSubset::iterator iter = pset->begin();

    // Initialize patch variables
    se = 0.0;

    // Get patch info
    // Vector dx = patch->dCell();
    // Unused    double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};

    // Plastic gets and allocates
    if (d_usePlasticity) {
      old_dw->get(pPlasticStrain, pPlasticStrain_label, pset);
      new_dw->allocateAndPut(
        pPlasticStrain_new, pPlasticStrain_label_preReloc, pset);
    }

    // Universal gets and allocates
    old_dw->get(pMass, lb->pMassLabel, pset);
    old_dw->get(pX, lb->pXLabel, pset);
    old_dw->get(pSize, lb->pSizeLabel, pset);
    old_dw->get(pDefGrad, lb->pDefGradLabel, pset);
    old_dw->get(pBeBar, bElBarLabel, pset);

    new_dw->get(pVolume_new, lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    // Allocate space for updated particle variables
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(pBeBar_new, bElBarLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);

    if (matl->getIsRigid()) {
      for (iter = pset->begin(); iter != pset->end(); iter++) {
        particleIndex idx = *iter;
        // Assign zero internal heating by default - modify if necessary.
        pdTdt[idx]       = 0.0;
        pStress_new[idx] = Matrix3(0.0);
      }
    } else { /*if(!matl->getIsRigid()) */
      // Compute the displacement gradient and the deformation gradient
      auto interpolator = flag->d_interpolator->clone(patch);
      std::vector<IntVector> ni(interpolator->size());
      std::vector<Vector> d_S(interpolator->size());

      // Unused because no "active stress carried over from CNHImplicit
      // double time = d_mat_manager->getElapsedTime();

      for (iter = pset->begin(); iter != pset->end(); iter++) {
        particleIndex idx = *iter;

        // Assign zero internal heating by default - modify if necessary.
        pdTdt[idx] = 0.0;

        defGradInc  = pDefGrad_new[idx] * pDefGrad[idx].Inverse();
        double Jinc = defGradInc.Determinant();

        // Update the deformation gradient tensor to its time n+1 value.
        double J = pDefGrad_new[idx].Determinant();

        if (d_usePlasticity) {

          // Compute trial BeBar
          relDefGradBar = defGradInc / cbrt(Jinc);

          // Compute the trial elastic part of the volume preserving
          // part of the left Cauchy-Green deformation tensor
          beBarTrial = relDefGradBar * pBeBar[idx] * relDefGradBar.Transpose();
        } else {
          beBarTrial = pDefGrad_new[idx] * pDefGrad_new[idx].Transpose() *
                       pow(J, -(2. / 3.));
        }

        // Compute the deformed volume
        // double rho_cur   = rho_orig/J;

        double IEl   = onethird * beBarTrial.Trace();
        double muBar = IEl * shear;

        // tauDevTrial is equal to the shear modulus times dev(bElBar)
        // Compute ||tauDevTrial||
        tauDevTrial   = (beBarTrial - Identity * IEl) * shear;
        double sTnorm = tauDevTrial.Norm();

        // get the hydrostatic part of the stress
        double p = bulk * log(J) / J;

        // Check for plastic loading
        double alpha = 0.0;
        if (d_usePlasticity) {
          alpha = pPlasticStrain[idx];
          p     = 0.5 * bulk * (J - 1.0 / J);
        }
        double fTrial = sTnorm - sqtwthds * (hardModulus * alpha + flowStress);

        if (d_usePlasticity && (fTrial > 0.0)) {
          // plastic
          // Compute increment of slip in the direction of flow
          double delgamma =
            (fTrial / (2.0 * muBar)) / (1.0 + (hardModulus / (3.0 * muBar)));
          normal = tauDevTrial / sTnorm;

          // The actual shear stress
          tauDev = tauDevTrial - normal * 2.0 * muBar * delgamma;

          // Deal with history variables
          pPlasticStrain_new[idx] = alpha + sqtwthds * delgamma;
          pBeBar_new[idx]         = tauDev / shear + Identity * IEl;
        } else {

          // The actual shear stress
          tauDev          = tauDevTrial;
          pBeBar_new[idx] = beBarTrial;

          // carry forward in implicit
          if (d_usePlasticity) {
            pPlasticStrain_new[idx] = alpha;
          }
        }

        // compute the total stress (volumetric + deviatoric)
        pStress_new[idx] = Identity * p + tauDev / J;

        // Compute the strain energy for non-localized particles
        double U = .5 * bulk * (.5 * (J * J - 1.0) - log(J));
        double W = .5 * shear * (pBeBar_new[idx].Trace() - 3.0);
        double e = (U + W) * pVolume_new[idx] / J;
        se += e;

      } // end loop over particles
      if (flag->d_reductionVars->accStrainEnergy ||
          flag->d_reductionVars->strainEnergy) {
        new_dw->put(sum_vartype(se), lb->StrainEnergyLabel);
      }
      // delete interpolator;
    } // End rigid else
  }   // End Patch For Loop
}

/*! Compute tangent stiffness matrix */
void
UCNH::computeTangentStiffnessMatrix(const Matrix3& sigdev,
                                    const double& mubar,
                                    const double& J,
                                    const double& bulk,
                                    double D[6][6])
{
  double twth  = 2.0 / 3.0;
  double frth  = 2.0 * twth;
  double coef1 = bulk;
  double coef2 = 2. * bulk * log(J);

  for (int ii = 0; ii < 6; ++ii) {
    for (int jj = 0; jj < 6; ++jj) {
      D[ii][jj] = 0.0;
    }
  }
  D[0][0] = coef1 - coef2 + mubar * frth - frth * sigdev(0, 0);
  D[0][1] = coef1 - mubar * twth - twth * (sigdev(0, 0) + sigdev(1, 1));
  D[0][2] = coef1 - mubar * twth - twth * (sigdev(0, 0) + sigdev(2, 2));
  D[0][3] = -twth * (sigdev(0, 1));
  D[0][4] = -twth * (sigdev(0, 2));
  D[0][5] = -twth * (sigdev(1, 2));
  D[1][1] = coef1 - coef2 + mubar * frth - frth * sigdev(1, 1);
  D[1][2] = coef1 - mubar * twth - twth * (sigdev(1, 1) + sigdev(2, 2));
  D[1][3] = D[0][3];
  D[1][4] = D[0][4];
  D[1][5] = D[0][5];
  D[2][2] = coef1 - coef2 + mubar * frth - frth * sigdev(2, 2);
  D[2][3] = D[0][3];
  D[2][4] = D[0][4];
  D[2][5] = D[0][5];
  D[3][3] = -.5 * coef2 + mubar;
  D[4][4] = D[3][3];
  D[5][5] = D[3][3];
}

/*! Compute K matrix */
void
UCNH::computeStiffnessMatrix(const double B[6][24],
                             const double Bnl[3][24],
                             const double D[6][6],
                             const Matrix3& sig,
                             const double& vol_old,
                             const double& vol_new,
                             double Kmatrix[24][24])
{

  // Kmat = B.transpose()*D*B*volold
  double Kmat[24][24];
  BtDB(B, D, Kmat);

  // Kgeo = Bnl.transpose*sig*Bnl*volnew;
  double Kgeo[24][24];
  BnlTSigBnl(sig, Bnl, Kgeo);

  /*
   std::cout.setf(ios::scientific,ios::floatfield);
   std::cout.precision(10);
   std::cout << "Kmat = " << "\n";
   for(int kk = 0; kk < 24; kk++) {
   for (int ll = 0; ll < 24; ++ll) {
   std::cout << Kmat[ll][kk] << " " ;
   }
   std::cout << "\n";
   }
   std::cout << "Kgeo = " << "\n";
   for(int kk = 0; kk < 24; kk++) {
   for (int ll = 0; ll < 24; ++ll) {
   std::cout << Kgeo[ll][kk] << " " ;
   }
   std::cout << "\n";
   }
   */

  for (int ii = 0; ii < 24; ii++) {
    for (int jj = 0; jj < 24; jj++) {
      Kmatrix[ii][jj] = Kmat[ii][jj] * vol_old + Kgeo[ii][jj] * vol_new;
    }
  }
}

void
UCNH::BnlTSigBnl(const Matrix3& sig,
                 const double Bnl[3][24],
                 double Kgeo[24][24]) const
{
  double t1, t10, t11, t12, t13, t14, t15, t16, t17;
  double t18, t19, t2, t20, t21, t22, t23, t24, t25;
  double t26, t27, t28, t29, t3, t30, t31, t32, t33;
  double t34, t35, t36, t37, t38, t39, t4, t40, t41;
  double t42, t43, t44, t45, t46, t47, t48, t49, t5;
  double t50, t51, t52, t53, t54, t55, t56, t57, t58;
  double t59, t6, t60, t61, t62, t63, t64, t65, t66;
  double t67, t68, t69, t7, t70, t71, t72, t73, t74;
  double t75, t77, t78, t8, t81, t85, t88, t9, t90;
  double t79, t82, t83, t86, t87, t89;

  t1  = Bnl[0][0] * sig(0, 0);
  t4  = Bnl[0][0] * sig(0, 0);
  t2  = Bnl[0][0] * sig(0, 1);
  t3  = Bnl[0][0] * sig(0, 2);
  t5  = Bnl[1][1] * sig(1, 1);
  t8  = Bnl[1][1] * sig(1, 1);
  t6  = Bnl[1][1] * sig(1, 2);
  t7  = Bnl[1][1] * sig(0, 1);
  t9  = Bnl[2][2] * sig(2, 2);
  t12 = Bnl[2][2] * sig(2, 2);
  t10 = Bnl[2][2] * sig(0, 2);
  t11 = Bnl[2][2] * sig(1, 2);
  t13 = Bnl[0][3] * sig(0, 0);
  t16 = Bnl[0][3] * sig(0, 0);
  t14 = Bnl[0][3] * sig(0, 1);
  t15 = Bnl[0][3] * sig(0, 2);
  t17 = Bnl[1][4] * sig(1, 1);
  t20 = Bnl[1][4] * sig(1, 1);
  t18 = Bnl[1][4] * sig(1, 2);
  t19 = Bnl[1][4] * sig(0, 1);
  t21 = Bnl[2][5] * sig(2, 2);
  t22 = Bnl[2][5] * sig(0, 2);
  t23 = Bnl[2][5] * sig(1, 2);
  t24 = Bnl[2][5] * sig(2, 2);
  t25 = Bnl[0][6] * sig(0, 0);
  t26 = Bnl[0][6] * sig(0, 1);
  t27 = Bnl[0][6] * sig(0, 2);
  t28 = Bnl[0][6] * sig(0, 0);
  t29 = Bnl[1][7] * sig(1, 1);
  t30 = Bnl[1][7] * sig(1, 2);
  t31 = Bnl[1][7] * sig(0, 1);
  t32 = Bnl[1][7] * sig(1, 1);
  t33 = Bnl[2][8] * sig(2, 2);
  t34 = Bnl[2][8] * sig(0, 2);
  t35 = Bnl[2][8] * sig(1, 2);
  t36 = Bnl[2][8] * sig(2, 2);
  t37 = Bnl[0][9] * sig(0, 0);
  t38 = Bnl[0][9] * sig(0, 1);
  t39 = Bnl[0][9] * sig(0, 2);
  t40 = Bnl[0][9] * sig(0, 0);
  t41 = Bnl[1][10] * sig(1, 1);
  t42 = Bnl[1][10] * sig(1, 2);
  t43 = Bnl[1][10] * sig(0, 1);
  t44 = Bnl[1][10] * sig(1, 1);
  t45 = Bnl[2][11] * sig(2, 2);
  t46 = Bnl[2][11] * sig(0, 2);
  t47 = Bnl[2][11] * sig(1, 2);
  t48 = Bnl[2][11] * sig(2, 2);
  t49 = Bnl[0][12] * sig(0, 0);
  t50 = Bnl[0][12] * sig(0, 1);
  t51 = Bnl[0][12] * sig(0, 2);
  t52 = Bnl[0][12] * sig(0, 0);
  t53 = Bnl[1][13] * sig(1, 1);
  t54 = Bnl[1][13] * sig(1, 2);
  t55 = Bnl[1][13] * sig(0, 1);
  t56 = Bnl[1][13] * sig(1, 1);
  t57 = Bnl[2][14] * sig(2, 2);
  t58 = Bnl[2][14] * sig(0, 2);
  t59 = Bnl[2][14] * sig(1, 2);
  t60 = Bnl[2][14] * sig(2, 2);
  t61 = Bnl[0][15] * sig(0, 0);
  t62 = Bnl[0][15] * sig(0, 1);
  t63 = Bnl[0][15] * sig(0, 2);
  t64 = Bnl[0][15] * sig(0, 0);
  t65 = Bnl[1][16] * sig(1, 1);
  t66 = Bnl[1][16] * sig(1, 2);
  t67 = Bnl[1][16] * sig(0, 1);
  t68 = Bnl[1][16] * sig(1, 1);
  t69 = Bnl[2][17] * sig(2, 2);
  t70 = Bnl[2][17] * sig(0, 2);
  t71 = Bnl[2][17] * sig(1, 2);
  t72 = Bnl[2][17] * sig(2, 2);
  t73 = Bnl[0][18] * sig(0, 0);
  t74 = Bnl[0][18] * sig(0, 1);
  t75 = Bnl[0][18] * sig(0, 2);
  t77 = Bnl[1][19] * sig(1, 1);
  t78 = Bnl[1][19] * sig(1, 2);
  t79 = Bnl[1][19] * sig(0, 1);
  t81 = Bnl[2][20] * sig(2, 2);
  t82 = Bnl[2][20] * sig(0, 2);
  t83 = Bnl[2][20] * sig(1, 2);
  t85 = Bnl[0][21] * sig(0, 0);
  t86 = Bnl[0][21] * sig(0, 1);
  t87 = Bnl[0][21] * sig(0, 2);
  t88 = Bnl[1][22] * sig(1, 1);
  t89 = Bnl[1][22] * sig(1, 2);
  t90 = Bnl[2][23] * sig(2, 2);

  Kgeo[0][0]   = t1 * Bnl[0][0];
  Kgeo[0][1]   = t2 * Bnl[1][1];
  Kgeo[0][2]   = t3 * Bnl[2][2];
  Kgeo[0][3]   = t4 * Bnl[0][3];
  Kgeo[0][4]   = t2 * Bnl[1][4];
  Kgeo[0][5]   = t3 * Bnl[2][5];
  Kgeo[0][6]   = t4 * Bnl[0][6];
  Kgeo[0][7]   = t2 * Bnl[1][7];
  Kgeo[0][8]   = t3 * Bnl[2][8];
  Kgeo[0][9]   = t4 * Bnl[0][9];
  Kgeo[0][10]  = t2 * Bnl[1][10];
  Kgeo[0][11]  = t3 * Bnl[2][11];
  Kgeo[0][12]  = t4 * Bnl[0][12];
  Kgeo[0][13]  = t2 * Bnl[1][13];
  Kgeo[0][14]  = t3 * Bnl[2][14];
  Kgeo[0][15]  = t4 * Bnl[0][15];
  Kgeo[0][16]  = t2 * Bnl[1][16];
  Kgeo[0][17]  = t3 * Bnl[2][17];
  Kgeo[0][18]  = t4 * Bnl[0][18];
  Kgeo[0][19]  = t2 * Bnl[1][19];
  Kgeo[0][20]  = t3 * Bnl[2][20];
  Kgeo[0][21]  = t4 * Bnl[0][21];
  Kgeo[0][22]  = t2 * Bnl[1][22];
  Kgeo[0][23]  = t3 * Bnl[2][23];
  Kgeo[1][0]   = Kgeo[0][1];
  Kgeo[1][1]   = t5 * Bnl[1][1];
  Kgeo[1][2]   = t6 * Bnl[2][2];
  Kgeo[1][3]   = t7 * Bnl[0][3];
  Kgeo[1][4]   = Bnl[1][4] * t8;
  Kgeo[1][5]   = t6 * Bnl[2][5];
  Kgeo[1][6]   = t7 * Bnl[0][6];
  Kgeo[1][7]   = Bnl[1][7] * t8;
  Kgeo[1][8]   = t6 * Bnl[2][8];
  Kgeo[1][9]   = t7 * Bnl[0][9];
  Kgeo[1][10]  = Bnl[1][10] * t8;
  Kgeo[1][11]  = t6 * Bnl[2][11];
  Kgeo[1][12]  = t7 * Bnl[0][12];
  Kgeo[1][13]  = Bnl[1][13] * t8;
  Kgeo[1][14]  = t6 * Bnl[2][14];
  Kgeo[1][15]  = t7 * Bnl[0][15];
  Kgeo[1][16]  = Bnl[1][16] * t8;
  Kgeo[1][17]  = t6 * Bnl[2][17];
  Kgeo[1][18]  = t7 * Bnl[0][18];
  Kgeo[1][19]  = Bnl[1][19] * t8;
  Kgeo[1][20]  = t6 * Bnl[2][20];
  Kgeo[1][21]  = t7 * Bnl[0][21];
  Kgeo[1][22]  = Bnl[1][22] * t8;
  Kgeo[1][23]  = t6 * Bnl[2][23];
  Kgeo[2][0]   = Kgeo[0][2];
  Kgeo[2][1]   = Kgeo[1][2];
  Kgeo[2][2]   = t9 * Bnl[2][2];
  Kgeo[2][3]   = t10 * Bnl[0][3];
  Kgeo[2][4]   = Bnl[1][4] * t11;
  Kgeo[2][5]   = t12 * Bnl[2][5];
  Kgeo[2][6]   = t10 * Bnl[0][6];
  Kgeo[2][7]   = Bnl[1][7] * t11;
  Kgeo[2][8]   = t12 * Bnl[2][8];
  Kgeo[2][9]   = t10 * Bnl[0][9];
  Kgeo[2][10]  = Bnl[1][10] * t11;
  Kgeo[2][11]  = t12 * Bnl[2][11];
  Kgeo[2][12]  = t10 * Bnl[0][12];
  Kgeo[2][13]  = Bnl[1][13] * t11;
  Kgeo[2][14]  = t12 * Bnl[2][14];
  Kgeo[2][15]  = t10 * Bnl[0][15];
  Kgeo[2][16]  = Bnl[1][16] * t11;
  Kgeo[2][17]  = t12 * Bnl[2][17];
  Kgeo[2][18]  = t10 * Bnl[0][18];
  Kgeo[2][19]  = t11 * Bnl[1][19];
  Kgeo[2][20]  = t12 * Bnl[2][20];
  Kgeo[2][21]  = t10 * Bnl[0][21];
  Kgeo[2][22]  = t11 * Bnl[1][22];
  Kgeo[2][23]  = t12 * Bnl[2][23];
  Kgeo[3][0]   = Kgeo[0][3];
  Kgeo[3][1]   = Kgeo[1][3];
  Kgeo[3][2]   = Kgeo[2][3];
  Kgeo[3][3]   = t13 * Bnl[0][3];
  Kgeo[3][4]   = t14 * Bnl[1][4];
  Kgeo[3][5]   = Bnl[2][5] * t15;
  Kgeo[3][6]   = t16 * Bnl[0][6];
  Kgeo[3][7]   = t14 * Bnl[1][7];
  Kgeo[3][8]   = Bnl[2][8] * t15;
  Kgeo[3][9]   = t16 * Bnl[0][9];
  Kgeo[3][10]  = t14 * Bnl[1][10];
  Kgeo[3][11]  = Bnl[2][11] * t15;
  Kgeo[3][12]  = t16 * Bnl[0][12];
  Kgeo[3][13]  = t14 * Bnl[1][13];
  Kgeo[3][14]  = Bnl[2][14] * t15;
  Kgeo[3][15]  = t16 * Bnl[0][15];
  Kgeo[3][16]  = t14 * Bnl[1][16];
  Kgeo[3][17]  = Bnl[2][17] * t15;
  Kgeo[3][18]  = t16 * Bnl[0][18];
  Kgeo[3][19]  = t14 * Bnl[1][19];
  Kgeo[3][20]  = Bnl[2][20] * t15;
  Kgeo[3][21]  = t16 * Bnl[0][21];
  Kgeo[3][22]  = t14 * Bnl[1][22];
  Kgeo[3][23]  = Bnl[2][23] * t15;
  Kgeo[4][0]   = Kgeo[0][4];
  Kgeo[4][1]   = Kgeo[1][4];
  Kgeo[4][2]   = Kgeo[2][4];
  Kgeo[4][3]   = Kgeo[3][4];
  Kgeo[4][4]   = t17 * Bnl[1][4];
  Kgeo[4][5]   = t18 * Bnl[2][5];
  Kgeo[4][6]   = t19 * Bnl[0][6];
  Kgeo[4][7]   = t20 * Bnl[1][7];
  Kgeo[4][8]   = t18 * Bnl[2][8];
  Kgeo[4][9]   = t19 * Bnl[0][9];
  Kgeo[4][10]  = t20 * Bnl[1][10];
  Kgeo[4][11]  = t18 * Bnl[2][11];
  Kgeo[4][12]  = t19 * Bnl[0][12];
  Kgeo[4][13]  = t20 * Bnl[1][13];
  Kgeo[4][14]  = t18 * Bnl[2][14];
  Kgeo[4][15]  = t19 * Bnl[0][15];
  Kgeo[4][16]  = t20 * Bnl[1][16];
  Kgeo[4][17]  = t18 * Bnl[2][17];
  Kgeo[4][18]  = t19 * Bnl[0][18];
  Kgeo[4][19]  = t20 * Bnl[1][19];
  Kgeo[4][20]  = t18 * Bnl[2][20];
  Kgeo[4][21]  = t19 * Bnl[0][21];
  Kgeo[4][22]  = t20 * Bnl[1][22];
  Kgeo[4][23]  = t18 * Bnl[2][23];
  Kgeo[5][0]   = Kgeo[0][5];
  Kgeo[5][1]   = Kgeo[1][5];
  Kgeo[5][2]   = Kgeo[2][5];
  Kgeo[5][3]   = Kgeo[3][5];
  Kgeo[5][4]   = Kgeo[4][5];
  Kgeo[5][5]   = t21 * Bnl[2][5];
  Kgeo[5][6]   = t22 * Bnl[0][6];
  Kgeo[5][7]   = t23 * Bnl[1][7];
  Kgeo[5][8]   = t24 * Bnl[2][8];
  Kgeo[5][9]   = t22 * Bnl[0][9];
  Kgeo[5][10]  = t23 * Bnl[1][10];
  Kgeo[5][11]  = t24 * Bnl[2][11];
  Kgeo[5][12]  = t22 * Bnl[0][12];
  Kgeo[5][13]  = t23 * Bnl[1][13];
  Kgeo[5][14]  = t24 * Bnl[2][14];
  Kgeo[5][15]  = t22 * Bnl[0][15];
  Kgeo[5][16]  = t23 * Bnl[1][16];
  Kgeo[5][17]  = t24 * Bnl[2][17];
  Kgeo[5][18]  = t22 * Bnl[0][18];
  Kgeo[5][19]  = t23 * Bnl[1][19];
  Kgeo[5][20]  = t24 * Bnl[2][20];
  Kgeo[5][21]  = t22 * Bnl[0][21];
  Kgeo[5][22]  = t23 * Bnl[1][22];
  Kgeo[5][23]  = t24 * Bnl[2][23];
  Kgeo[6][0]   = Kgeo[0][6];
  Kgeo[6][1]   = Kgeo[1][6];
  Kgeo[6][2]   = Kgeo[2][6];
  Kgeo[6][3]   = Kgeo[3][6];
  Kgeo[6][4]   = Kgeo[4][6];
  Kgeo[6][5]   = Kgeo[5][6];
  Kgeo[6][6]   = t25 * Bnl[0][6];
  Kgeo[6][7]   = t26 * Bnl[1][7];
  Kgeo[6][8]   = t27 * Bnl[2][8];
  Kgeo[6][9]   = t28 * Bnl[0][9];
  Kgeo[6][10]  = t26 * Bnl[1][10];
  Kgeo[6][11]  = t27 * Bnl[2][11];
  Kgeo[6][12]  = t28 * Bnl[0][12];
  Kgeo[6][13]  = t26 * Bnl[1][13];
  Kgeo[6][14]  = t27 * Bnl[2][14];
  Kgeo[6][15]  = t28 * Bnl[0][15];
  Kgeo[6][16]  = t26 * Bnl[1][16];
  Kgeo[6][17]  = t27 * Bnl[2][17];
  Kgeo[6][18]  = t28 * Bnl[0][18];
  Kgeo[6][19]  = t26 * Bnl[1][19];
  Kgeo[6][20]  = t27 * Bnl[2][20];
  Kgeo[6][21]  = t28 * Bnl[0][21];
  Kgeo[6][22]  = t26 * Bnl[1][22];
  Kgeo[6][23]  = t27 * Bnl[2][23];
  Kgeo[7][0]   = Kgeo[0][7];
  Kgeo[7][1]   = Kgeo[1][7];
  Kgeo[7][2]   = Kgeo[2][7];
  Kgeo[7][3]   = Kgeo[3][7];
  Kgeo[7][4]   = Kgeo[4][7];
  Kgeo[7][5]   = Kgeo[5][7];
  Kgeo[7][6]   = Kgeo[6][7];
  Kgeo[7][7]   = t29 * Bnl[1][7];
  Kgeo[7][8]   = t30 * Bnl[2][8];
  Kgeo[7][9]   = t31 * Bnl[0][9];
  Kgeo[7][10]  = t32 * Bnl[1][10];
  Kgeo[7][11]  = t30 * Bnl[2][11];
  Kgeo[7][12]  = t31 * Bnl[0][12];
  Kgeo[7][13]  = t32 * Bnl[1][13];
  Kgeo[7][14]  = t30 * Bnl[2][14];
  Kgeo[7][15]  = t31 * Bnl[0][15];
  Kgeo[7][16]  = t32 * Bnl[1][16];
  Kgeo[7][17]  = t30 * Bnl[2][17];
  Kgeo[7][18]  = t31 * Bnl[0][18];
  Kgeo[7][19]  = t32 * Bnl[1][19];
  Kgeo[7][20]  = t30 * Bnl[2][20];
  Kgeo[7][21]  = t31 * Bnl[0][21];
  Kgeo[7][22]  = t32 * Bnl[1][22];
  Kgeo[7][23]  = t30 * Bnl[2][23];
  Kgeo[8][0]   = Kgeo[0][8];
  Kgeo[8][1]   = Kgeo[1][8];
  Kgeo[8][2]   = Kgeo[2][8];
  Kgeo[8][3]   = Kgeo[3][8];
  Kgeo[8][4]   = Kgeo[4][8];
  Kgeo[8][5]   = Kgeo[5][8];
  Kgeo[8][6]   = Kgeo[6][8];
  Kgeo[8][7]   = Kgeo[7][8];
  Kgeo[8][8]   = t33 * Bnl[2][8];
  Kgeo[8][9]   = t34 * Bnl[0][9];
  Kgeo[8][10]  = t35 * Bnl[1][10];
  Kgeo[8][11]  = t36 * Bnl[2][11];
  Kgeo[8][12]  = t34 * Bnl[0][12];
  Kgeo[8][13]  = t35 * Bnl[1][13];
  Kgeo[8][14]  = t36 * Bnl[2][14];
  Kgeo[8][15]  = t34 * Bnl[0][15];
  Kgeo[8][16]  = t35 * Bnl[1][16];
  Kgeo[8][17]  = t36 * Bnl[2][17];
  Kgeo[8][18]  = t34 * Bnl[0][18];
  Kgeo[8][19]  = t35 * Bnl[1][19];
  Kgeo[8][20]  = t36 * Bnl[2][20];
  Kgeo[8][21]  = t34 * Bnl[0][21];
  Kgeo[8][22]  = t35 * Bnl[1][22];
  Kgeo[8][23]  = t36 * Bnl[2][23];
  Kgeo[9][0]   = Kgeo[0][9];
  Kgeo[9][1]   = Kgeo[1][9];
  Kgeo[9][2]   = Kgeo[2][9];
  Kgeo[9][3]   = Kgeo[3][9];
  Kgeo[9][4]   = Kgeo[4][9];
  Kgeo[9][5]   = Kgeo[5][9];
  Kgeo[9][6]   = Kgeo[6][9];
  Kgeo[9][7]   = Kgeo[7][9];
  Kgeo[9][8]   = Kgeo[8][9];
  Kgeo[9][9]   = t37 * Bnl[0][9];
  Kgeo[9][10]  = t38 * Bnl[1][10];
  Kgeo[9][11]  = t39 * Bnl[2][11];
  Kgeo[9][12]  = t40 * Bnl[0][12];
  Kgeo[9][13]  = t38 * Bnl[1][13];
  Kgeo[9][14]  = t39 * Bnl[2][14];
  Kgeo[9][15]  = t40 * Bnl[0][15];
  Kgeo[9][16]  = t38 * Bnl[1][16];
  Kgeo[9][17]  = t39 * Bnl[2][17];
  Kgeo[9][18]  = t40 * Bnl[0][18];
  Kgeo[9][19]  = t38 * Bnl[1][19];
  Kgeo[9][20]  = t39 * Bnl[2][20];
  Kgeo[9][21]  = t40 * Bnl[0][21];
  Kgeo[9][22]  = t38 * Bnl[1][22];
  Kgeo[9][23]  = t39 * Bnl[2][23];
  Kgeo[10][0]  = Kgeo[0][10];
  Kgeo[10][1]  = Kgeo[1][10];
  Kgeo[10][2]  = Kgeo[2][10];
  Kgeo[10][3]  = Kgeo[3][10];
  Kgeo[10][4]  = Kgeo[4][10];
  Kgeo[10][5]  = Kgeo[5][10];
  Kgeo[10][6]  = Kgeo[6][10];
  Kgeo[10][7]  = Kgeo[7][10];
  Kgeo[10][8]  = Kgeo[8][10];
  Kgeo[10][9]  = Kgeo[9][10];
  Kgeo[10][10] = t41 * Bnl[1][10];
  Kgeo[10][11] = t42 * Bnl[2][11];
  Kgeo[10][12] = t43 * Bnl[0][12];
  Kgeo[10][13] = t44 * Bnl[1][13];
  Kgeo[10][14] = t42 * Bnl[2][14];
  Kgeo[10][15] = t43 * Bnl[0][15];
  Kgeo[10][16] = t44 * Bnl[1][16];
  Kgeo[10][17] = t42 * Bnl[2][17];
  Kgeo[10][18] = t43 * Bnl[0][18];
  Kgeo[10][19] = t44 * Bnl[1][19];
  Kgeo[10][20] = t42 * Bnl[2][20];
  Kgeo[10][21] = t43 * Bnl[0][21];
  Kgeo[10][22] = t44 * Bnl[1][22];
  Kgeo[10][23] = t42 * Bnl[2][23];
  Kgeo[11][0]  = Kgeo[0][11];
  Kgeo[11][1]  = Kgeo[1][11];
  Kgeo[11][2]  = Kgeo[2][11];
  Kgeo[11][3]  = Kgeo[3][11];
  Kgeo[11][4]  = Kgeo[4][11];
  Kgeo[11][5]  = Kgeo[5][11];
  Kgeo[11][6]  = Kgeo[6][11];
  Kgeo[11][7]  = Kgeo[7][11];
  Kgeo[11][8]  = Kgeo[8][11];
  Kgeo[11][9]  = Kgeo[9][11];
  Kgeo[11][10] = Kgeo[10][11];
  Kgeo[11][11] = t45 * Bnl[2][11];
  Kgeo[11][12] = t46 * Bnl[0][12];
  Kgeo[11][13] = t47 * Bnl[1][13];
  Kgeo[11][14] = t48 * Bnl[2][14];
  Kgeo[11][15] = t46 * Bnl[0][15];
  Kgeo[11][16] = t47 * Bnl[1][16];
  Kgeo[11][17] = t48 * Bnl[2][17];
  Kgeo[11][18] = t46 * Bnl[0][18];
  Kgeo[11][19] = t47 * Bnl[1][19];
  Kgeo[11][20] = t48 * Bnl[2][20];
  Kgeo[11][21] = t46 * Bnl[0][21];
  Kgeo[11][22] = t47 * Bnl[1][22];
  Kgeo[11][23] = t48 * Bnl[2][23];
  Kgeo[12][0]  = Kgeo[0][12];
  Kgeo[12][1]  = Kgeo[1][12];
  Kgeo[12][2]  = Kgeo[2][12];
  Kgeo[12][3]  = Kgeo[3][12];
  Kgeo[12][4]  = Kgeo[4][12];
  Kgeo[12][5]  = Kgeo[5][12];
  Kgeo[12][6]  = Kgeo[6][12];
  Kgeo[12][7]  = Kgeo[7][12];
  Kgeo[12][8]  = Kgeo[8][12];
  Kgeo[12][9]  = Kgeo[9][12];
  Kgeo[12][10] = Kgeo[10][12];
  Kgeo[12][11] = Kgeo[11][12];
  Kgeo[12][12] = t49 * Bnl[0][12];
  Kgeo[12][13] = t50 * Bnl[1][13];
  Kgeo[12][14] = t51 * Bnl[2][14];
  Kgeo[12][15] = t52 * Bnl[0][15];
  Kgeo[12][16] = t50 * Bnl[1][16];
  Kgeo[12][17] = t51 * Bnl[2][17];
  Kgeo[12][18] = t52 * Bnl[0][18];
  Kgeo[12][19] = t50 * Bnl[1][19];
  Kgeo[12][20] = t51 * Bnl[2][20];
  Kgeo[12][21] = t52 * Bnl[0][21];
  Kgeo[12][22] = t50 * Bnl[1][22];
  Kgeo[12][23] = t51 * Bnl[2][23];
  Kgeo[13][0]  = Kgeo[0][13];
  Kgeo[13][1]  = Kgeo[1][13];
  Kgeo[13][2]  = Kgeo[2][13];
  Kgeo[13][3]  = Kgeo[3][13];
  Kgeo[13][4]  = Kgeo[4][13];
  Kgeo[13][5]  = Kgeo[5][13];
  Kgeo[13][6]  = Kgeo[6][13];
  Kgeo[13][7]  = Kgeo[7][13];
  Kgeo[13][8]  = Kgeo[8][13];
  Kgeo[13][9]  = Kgeo[9][13];
  Kgeo[13][10] = Kgeo[10][13];
  Kgeo[13][11] = Kgeo[11][13];
  Kgeo[13][12] = Kgeo[12][13];
  Kgeo[13][13] = t53 * Bnl[1][13];
  Kgeo[13][14] = t54 * Bnl[2][14];
  Kgeo[13][15] = t55 * Bnl[0][15];
  Kgeo[13][16] = t56 * Bnl[1][16];
  Kgeo[13][17] = t54 * Bnl[2][17];
  Kgeo[13][18] = t55 * Bnl[0][18];
  Kgeo[13][19] = t56 * Bnl[1][19];
  Kgeo[13][20] = t54 * Bnl[2][20];
  Kgeo[13][21] = t55 * Bnl[0][21];
  Kgeo[13][22] = t56 * Bnl[1][22];
  Kgeo[13][23] = t54 * Bnl[2][23];
  Kgeo[14][0]  = Kgeo[0][14];
  Kgeo[14][1]  = Kgeo[1][14];
  Kgeo[14][2]  = Kgeo[2][14];
  Kgeo[14][3]  = Kgeo[3][14];
  Kgeo[14][4]  = Kgeo[4][14];
  Kgeo[14][5]  = Kgeo[5][14];
  Kgeo[14][6]  = Kgeo[6][14];
  Kgeo[14][7]  = Kgeo[7][14];
  Kgeo[14][8]  = Kgeo[8][14];
  Kgeo[14][9]  = Kgeo[9][14];
  Kgeo[14][10] = Kgeo[10][14];
  Kgeo[14][11] = Kgeo[11][14];
  Kgeo[14][12] = Kgeo[12][14];
  Kgeo[14][13] = Kgeo[13][14];
  Kgeo[14][14] = t57 * Bnl[2][14];
  Kgeo[14][15] = t58 * Bnl[0][15];
  Kgeo[14][16] = t59 * Bnl[1][16];
  Kgeo[14][17] = t60 * Bnl[2][17];
  Kgeo[14][18] = t58 * Bnl[0][18];
  Kgeo[14][19] = t59 * Bnl[1][19];
  Kgeo[14][20] = t60 * Bnl[2][20];
  Kgeo[14][21] = t58 * Bnl[0][21];
  Kgeo[14][22] = t59 * Bnl[1][22];
  Kgeo[14][23] = t60 * Bnl[2][23];
  Kgeo[15][0]  = Kgeo[0][15];
  Kgeo[15][1]  = Kgeo[1][15];
  Kgeo[15][2]  = Kgeo[2][15];
  Kgeo[15][3]  = Kgeo[3][15];
  Kgeo[15][4]  = Kgeo[4][15];
  Kgeo[15][5]  = Kgeo[5][15];
  Kgeo[15][6]  = Kgeo[6][15];
  Kgeo[15][7]  = Kgeo[7][15];
  Kgeo[15][8]  = Kgeo[8][15];
  Kgeo[15][9]  = Kgeo[9][15];
  Kgeo[15][10] = Kgeo[10][15];
  Kgeo[15][11] = Kgeo[11][15];
  Kgeo[15][12] = Kgeo[12][15];
  Kgeo[15][13] = Kgeo[13][15];
  Kgeo[15][14] = Kgeo[14][15];
  Kgeo[15][15] = t61 * Bnl[0][15];
  Kgeo[15][16] = t62 * Bnl[1][16];
  Kgeo[15][17] = t63 * Bnl[2][17];
  Kgeo[15][18] = t64 * Bnl[0][18];
  Kgeo[15][19] = t62 * Bnl[1][19];
  Kgeo[15][20] = t63 * Bnl[2][20];
  Kgeo[15][21] = t64 * Bnl[0][21];
  Kgeo[15][22] = t62 * Bnl[1][22];
  Kgeo[15][23] = t63 * Bnl[2][23];
  Kgeo[16][0]  = Kgeo[0][16];
  Kgeo[16][1]  = Kgeo[1][16];
  Kgeo[16][2]  = Kgeo[2][16];
  Kgeo[16][3]  = Kgeo[3][16];
  Kgeo[16][4]  = Kgeo[4][16];
  Kgeo[16][5]  = Kgeo[5][16];
  Kgeo[16][6]  = Kgeo[6][16];
  Kgeo[16][7]  = Kgeo[7][16];
  Kgeo[16][8]  = Kgeo[8][16];
  Kgeo[16][9]  = Kgeo[9][16];
  Kgeo[16][10] = Kgeo[10][16];
  Kgeo[16][11] = Kgeo[11][16];
  Kgeo[16][12] = Kgeo[12][16];
  Kgeo[16][13] = Kgeo[13][16];
  Kgeo[16][14] = Kgeo[14][16];
  Kgeo[16][15] = Kgeo[15][16];
  Kgeo[16][16] = t65 * Bnl[1][16];
  Kgeo[16][17] = t66 * Bnl[2][17];
  Kgeo[16][18] = t67 * Bnl[0][18];
  Kgeo[16][19] = t68 * Bnl[1][19];
  Kgeo[16][20] = t66 * Bnl[2][20];
  Kgeo[16][21] = t67 * Bnl[0][21];
  Kgeo[16][22] = t68 * Bnl[1][22];
  Kgeo[16][23] = t66 * Bnl[2][23];
  Kgeo[17][0]  = Kgeo[0][17];
  Kgeo[17][1]  = Kgeo[1][17];
  Kgeo[17][2]  = Kgeo[2][17];
  Kgeo[17][3]  = Kgeo[3][17];
  Kgeo[17][4]  = Kgeo[4][17];
  Kgeo[17][5]  = Kgeo[5][17];
  Kgeo[17][6]  = Kgeo[6][17];
  Kgeo[17][7]  = Kgeo[7][17];
  Kgeo[17][8]  = Kgeo[8][17];
  Kgeo[17][9]  = Kgeo[9][17];
  Kgeo[17][10] = Kgeo[10][17];
  Kgeo[17][11] = Kgeo[11][17];
  Kgeo[17][12] = Kgeo[12][17];
  Kgeo[17][13] = Kgeo[13][17];
  Kgeo[17][14] = Kgeo[14][17];
  Kgeo[17][15] = Kgeo[15][17];
  Kgeo[17][16] = Kgeo[16][17];
  Kgeo[17][17] = t69 * Bnl[2][17];
  Kgeo[17][18] = t70 * Bnl[0][18];
  Kgeo[17][19] = t71 * Bnl[1][19];
  Kgeo[17][20] = t72 * Bnl[2][20];
  Kgeo[17][21] = t70 * Bnl[0][21];
  Kgeo[17][22] = t71 * Bnl[1][22];
  Kgeo[17][23] = t72 * Bnl[2][23];
  Kgeo[18][0]  = Kgeo[0][18];
  Kgeo[18][1]  = Kgeo[1][18];
  Kgeo[18][2]  = Kgeo[2][18];
  Kgeo[18][3]  = Kgeo[3][18];
  Kgeo[18][4]  = Kgeo[4][18];
  Kgeo[18][5]  = Kgeo[5][18];
  Kgeo[18][6]  = Kgeo[6][18];
  Kgeo[18][7]  = Kgeo[7][18];
  Kgeo[18][8]  = Kgeo[8][18];
  Kgeo[18][9]  = Kgeo[9][18];
  Kgeo[18][10] = Kgeo[10][18];
  Kgeo[18][11] = Kgeo[11][18];
  Kgeo[18][12] = Kgeo[12][18];
  Kgeo[18][13] = Kgeo[13][18];
  Kgeo[18][14] = Kgeo[14][18];
  Kgeo[18][15] = Kgeo[15][18];
  Kgeo[18][16] = Kgeo[16][18];
  Kgeo[18][17] = Kgeo[17][18];
  Kgeo[18][18] = t73 * Bnl[0][18];
  Kgeo[18][19] = t74 * Bnl[1][19];
  Kgeo[18][20] = t75 * Bnl[2][20];
  Kgeo[18][21] = t73 * Bnl[0][21];
  Kgeo[18][22] = t74 * Bnl[1][22];
  Kgeo[18][23] = t75 * Bnl[2][23];
  Kgeo[19][0]  = Kgeo[0][19];
  Kgeo[19][1]  = Kgeo[1][19];
  Kgeo[19][2]  = Kgeo[2][19];
  Kgeo[19][3]  = Kgeo[3][19];
  Kgeo[19][4]  = Kgeo[4][19];
  Kgeo[19][5]  = Kgeo[5][19];
  Kgeo[19][6]  = Kgeo[6][19];
  Kgeo[19][7]  = Kgeo[7][19];
  Kgeo[19][8]  = Kgeo[8][19];
  Kgeo[19][9]  = Kgeo[9][19];
  Kgeo[19][10] = Kgeo[10][19];
  Kgeo[19][11] = Kgeo[11][19];
  Kgeo[19][12] = Kgeo[12][19];
  Kgeo[19][13] = Kgeo[13][19];
  Kgeo[19][14] = Kgeo[14][19];
  Kgeo[19][15] = Kgeo[15][19];
  Kgeo[19][16] = Kgeo[16][19];
  Kgeo[19][17] = Kgeo[17][19];
  Kgeo[19][18] = Kgeo[18][19];
  Kgeo[19][19] = t77 * Bnl[1][19];
  Kgeo[19][20] = t78 * Bnl[2][20];
  Kgeo[19][21] = t79 * Bnl[0][21];
  Kgeo[19][22] = t77 * Bnl[1][22];
  Kgeo[19][23] = t78 * Bnl[2][23];
  Kgeo[20][0]  = Kgeo[0][20];
  Kgeo[20][1]  = Kgeo[1][20];
  Kgeo[20][2]  = Kgeo[2][20];
  Kgeo[20][3]  = Kgeo[3][20];
  Kgeo[20][4]  = Kgeo[4][20];
  Kgeo[20][5]  = Kgeo[5][20];
  Kgeo[20][6]  = Kgeo[6][20];
  Kgeo[20][7]  = Kgeo[7][20];
  Kgeo[20][8]  = Kgeo[8][20];
  Kgeo[20][9]  = Kgeo[9][20];
  Kgeo[20][10] = Kgeo[10][20];
  Kgeo[20][11] = Kgeo[11][20];
  Kgeo[20][12] = Kgeo[12][20];
  Kgeo[20][13] = Kgeo[13][20];
  Kgeo[20][14] = Kgeo[14][20];
  Kgeo[20][15] = Kgeo[15][20];
  Kgeo[20][16] = Kgeo[16][20];
  Kgeo[20][17] = Kgeo[17][20];
  Kgeo[20][18] = Kgeo[18][20];
  Kgeo[20][19] = Kgeo[19][20];
  Kgeo[20][20] = t81 * Bnl[2][20];
  Kgeo[20][21] = t82 * Bnl[0][21];
  Kgeo[20][22] = t83 * Bnl[1][22];
  Kgeo[20][23] = t81 * Bnl[2][23];
  Kgeo[21][0]  = Kgeo[0][21];
  Kgeo[21][1]  = Kgeo[1][21];
  Kgeo[21][2]  = Kgeo[2][21];
  Kgeo[21][3]  = Kgeo[3][21];
  Kgeo[21][4]  = Kgeo[4][21];
  Kgeo[21][5]  = Kgeo[5][21];
  Kgeo[21][6]  = Kgeo[6][21];
  Kgeo[21][7]  = Kgeo[7][21];
  Kgeo[21][8]  = Kgeo[8][21];
  Kgeo[21][9]  = Kgeo[9][21];
  Kgeo[21][10] = Kgeo[10][21];
  Kgeo[21][11] = Kgeo[11][21];
  Kgeo[21][12] = Kgeo[12][21];
  Kgeo[21][13] = Kgeo[13][21];
  Kgeo[21][14] = Kgeo[14][21];
  Kgeo[21][15] = Kgeo[15][21];
  Kgeo[21][16] = Kgeo[16][21];
  Kgeo[21][17] = Kgeo[17][21];
  Kgeo[21][18] = Kgeo[18][21];
  Kgeo[21][19] = Kgeo[19][21];
  Kgeo[21][20] = Kgeo[20][21];
  Kgeo[21][21] = t85 * Bnl[0][21];
  Kgeo[21][22] = t86 * Bnl[1][22];
  Kgeo[21][23] = t87 * Bnl[2][23];
  Kgeo[22][0]  = Kgeo[0][22];
  Kgeo[22][1]  = Kgeo[1][22];
  Kgeo[22][2]  = Kgeo[2][22];
  Kgeo[22][3]  = Kgeo[3][22];
  Kgeo[22][4]  = Kgeo[4][22];
  Kgeo[22][5]  = Kgeo[5][22];
  Kgeo[22][6]  = Kgeo[6][22];
  Kgeo[22][7]  = Kgeo[7][22];
  Kgeo[22][8]  = Kgeo[8][22];
  Kgeo[22][9]  = Kgeo[9][22];
  Kgeo[22][10] = Kgeo[10][22];
  Kgeo[22][11] = Kgeo[11][22];
  Kgeo[22][12] = Kgeo[12][22];
  Kgeo[22][13] = Kgeo[13][22];
  Kgeo[22][14] = Kgeo[14][22];
  Kgeo[22][15] = Kgeo[15][22];
  Kgeo[22][16] = Kgeo[16][22];
  Kgeo[22][17] = Kgeo[17][22];
  Kgeo[22][18] = Kgeo[18][22];
  Kgeo[22][19] = Kgeo[19][22];
  Kgeo[22][20] = Kgeo[20][22];
  Kgeo[22][21] = Kgeo[21][22];
  Kgeo[22][22] = t88 * Bnl[1][22];
  Kgeo[22][23] = t89 * Bnl[2][23];
  Kgeo[23][0]  = Kgeo[0][23];
  Kgeo[23][1]  = Kgeo[1][23];
  Kgeo[23][2]  = Kgeo[2][23];
  Kgeo[23][3]  = Kgeo[3][23];
  Kgeo[23][4]  = Kgeo[4][23];
  Kgeo[23][5]  = Kgeo[5][23];
  Kgeo[23][6]  = Kgeo[6][23];
  Kgeo[23][7]  = Kgeo[7][23];
  Kgeo[23][8]  = Kgeo[8][23];
  Kgeo[23][9]  = Kgeo[9][23];
  Kgeo[23][10] = Kgeo[10][23];
  Kgeo[23][11] = Kgeo[11][23];
  Kgeo[23][12] = Kgeo[12][23];
  Kgeo[23][13] = Kgeo[13][23];
  Kgeo[23][14] = Kgeo[14][23];
  Kgeo[23][15] = Kgeo[15][23];
  Kgeo[23][16] = Kgeo[16][23];
  Kgeo[23][17] = Kgeo[17][23];
  Kgeo[23][18] = Kgeo[18][23];
  Kgeo[23][19] = Kgeo[19][23];
  Kgeo[23][20] = Kgeo[20][23];
  Kgeo[23][21] = Kgeo[21][23];
  Kgeo[23][22] = Kgeo[22][23];
  Kgeo[23][23] = t90 * Bnl[2][23];
}

/*********************
 * For copying damage parameter to be used by SerialMPM or ImpMPM
 */
void
UCNH::addRequiresDamageParameter([[maybe_unused]] Task* task,
                                 [[maybe_unused]] const MPMMaterial* matl,
                                 [[maybe_unused]] const PatchSet* patches) const
{
}

void
UCNH::getDamageParameter([[maybe_unused]] const Patch* patch,
                         [[maybe_unused]] ParticleVariable<int>& damage,
                         [[maybe_unused]] int dwi,
                         [[maybe_unused]] DataWarehouse* old_dw,
                         [[maybe_unused]] DataWarehouse* new_dw)
{
}

/***********************************
 * For MPMICE only
 */
void
UCNH::computePressEOSCM(const double rho_cur,
                        double& pressure,
                        const double p_ref,
                        double& dp_drho,
                        double& cSquared,
                        const MPMMaterial* matl,
                        [[maybe_unused]] double temperature)
{
  double bulk     = d_initialData.Bulk;
  double rho_orig = matl->getInitialDensity();

  if (d_useModifiedEOS && rho_cur < rho_orig) {

    double A                = p_ref; // MODIFIED EOS
    double n                = bulk / p_ref;
    double rho_rat_to_the_n = pow(rho_cur / rho_orig, n);
    pressure                = A * rho_rat_to_the_n;
    dp_drho                 = (bulk / rho_cur) * rho_rat_to_the_n;
    cSquared                = dp_drho; // speed of sound squared

  } else { // STANDARD EOS

    double p = 0.0;
    d_eos->computePressure(rho_orig, rho_cur, p, dp_drho, cSquared);
    pressure = -p + p_ref;
    dp_drho  = -dp_drho;

    // double p_g = .5*bulk*(rho_cur/rho_orig - rho_orig/rho_cur);
    // pressure   = p_ref + p_g;
    // dp_drho    = .5*bulk*(rho_orig/(rho_cur*rho_cur) + 1./rho_orig);
    // cSquared   = bulk/rho_cur;  // speed of sound squared
  }
}

/***********************************
 * For MPMICE only
 */
// The "CM" versions use the pressure-volume relationship of the CNH model
double
UCNH::computeRhoMicroCM(double pressure,
                        const double p_ref,
                        const MPMMaterial* matl,
                        [[maybe_unused]] double temperature,
                        double rho_guess)
{
  /*
  if (std::isnan(pressure)) {
    std::cout << " Pressure = " << pressure << std::"\n";
  }
  */
  double rho_orig = matl->getInitialDensity();
  double bulk     = d_initialData.Bulk;

  double p_gauge = pressure - p_ref;
  double rho_cur = -1.0;
  bool error     = false;

  if (d_useModifiedEOS && p_gauge < 0.0) {

    double A = p_ref; // MODIFIED EOS
    double n = p_ref / bulk;
    rho_cur  = rho_orig * pow(pressure / A, n);

  } else { // STANDARD EOS

    try {
      rho_cur = d_eos->computeDensity(rho_orig, -p_gauge);
    } catch (ConvergenceFailure& e) {
      std::cout << e.message() << "\n";
      error = true;
    }
    if (error || rho_cur < 0.0 || std::isnan(rho_cur)) {
       std::ostringstream desc;
      desc << "rho_cur = " << rho_cur << " pressure = " << -p_gauge
           << " p_ref = " << p_ref << " 1/sp_vol_CC = " << rho_guess << "\n";
      throw InvalidValue(desc.str(), __FILE__, __LINE__);
    }

    // double p_g_over_bulk = p_gauge/bulk;
    // rho_cur=rho_orig*(p_g_over_bulk + sqrt(p_g_over_bulk*p_g_over_bulk +1.));
  }
  return rho_cur;
}

/***********************************
 * For MPMICE only
 */
double
UCNH::getCompressibility()
{
  return 1.0 / d_initialData.Bulk;
}

