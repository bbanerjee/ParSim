/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#include <CCA/Components/MPM/ConstitutiveModel/RockSoilModels/CamClay.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>

#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/MPMEquationOfStateFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/InternalVariableModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_CamClay.h>
#include <CCA/Components/MPM/ConstitutiveModel/ShearModulusModels/ShearModulusModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldConditionFactory.h>

#include <CCA/Ports/DataWarehouse.h>

#include <Core/Exceptions/ConvergenceFailure.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MPMInterpolators/LinearInterpolator.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/ParticleSubset.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/FastMatrix.h>
#include <Core/Math/Gaussian.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/Math/SymmMatrix3.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Util/DOUT.hpp>
#include <Core/Util/DebugStream.h>

#include <cmath>
#include <iostream>

namespace {

Uintah::Dout g_cc_dbg("CamClay", "CamClay", "general info on CamClay", false);
Uintah::Dout g_cc_conv("CamClayConv",
                       "CamClay",
                       "convergence of CamClay",
                       false);
} // namespace

using namespace Uintah;
using namespace Vaango;

//__________________________________
//  To turn on debug flags use, e.g.,
//  csh/tcsh : setenv SCI_DEBUG
//  "CamClay:+,CamClayDefGrad:+,CamClayStrain:+".....
//  bash     : export SCI_DEBUG="CamClay:+,CamClayDefGrad:+,CamClayConv:+" )
//  default is OFF

static DebugStream cout_CC_F("CamClayDefGrad", false);
static DebugStream cout_CC_Eps("CamClayStrain", false);

const double CamClay::sqrtThreeTwo = std::sqrt(1.5);
const double CamClay::sqrtTwoThird = 1.0 / CamClay::sqrtThreeTwo;

CamClay::CamClay(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  d_eos = MPMEquationOfStateFactory::create(ps);
  if (!d_eos) {
    std::ostringstream desc;
    desc << "**ERROR** Internal error while creating "
            "CamClay->MPMEquationOfStatFactory."
         << std::endl;
    throw InternalError(desc.str(), __FILE__, __LINE__);
  }

  d_shear = Vaango::ShearModulusModelFactory::create(ps, d_eos);
  if (!d_shear) {
    std::ostringstream desc;
    desc << "**ERROR** Internal error while creating "
            "CamClay->ShearModulusModelFactory."
         << std::endl;
    throw InternalError(desc.str(), __FILE__, __LINE__);
  }

  ProblemSpecP intvar_ps = ps->findBlock("internal_variable_model");
  if (!intvar_ps) {
    std::ostringstream err;
    err << "**ERROR** Please add an 'internal_variable_model' tag to the\n"
        << " 'camclay' block in the input .ups file.  The\n"
        << " default type is 'borja_consolidation_pressure'.\n";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }
  d_intvar = std::make_shared<Vaango::IntVar_BorjaPressure>(intvar_ps, d_shear);
  if (!d_intvar) {
    std::ostringstream err;
    err << "**ERROR** An error occured while creating the internal variable \n"
        << " model for CamClay. Please file a bug report.\n";
    throw InternalError(err.str(), __FILE__, __LINE__);
  }

  d_yield = Vaango::YieldConditionFactory::create(ps, d_intvar.get());
  if (!d_yield) {
    std::ostringstream desc;
    desc << "**ERROR** Internal error while creating "
            "CamClay->YieldConditionFactory."
         << std::endl;
    throw InternalError(desc.str(), __FILE__, __LINE__);
  }

  initializeLocalMPMLabels();
}

CamClay::CamClay(const CamClay* cm)
  : ConstitutiveModel(cm)
{
  d_eos    = MPMEquationOfStateFactory::createCopy(cm->d_eos);
  d_shear  = Vaango::ShearModulusModelFactory::createCopy(cm->d_shear);
  d_yield  = Vaango::YieldConditionFactory::createCopy(cm->d_yield);
  d_intvar = std::make_shared<Vaango::IntVar_BorjaPressure>(cm->d_intvar.get());

  initializeLocalMPMLabels();
}

CamClay::~CamClay()
{
  // Destructor
  VarLabel::destroy(pStrainLabel);
  VarLabel::destroy(pElasticStrainLabel);
  VarLabel::destroy(pDeltaGammaLabel);

  VarLabel::destroy(pStrainLabel_preReloc);
  VarLabel::destroy(pElasticStrainLabel_preReloc);
  VarLabel::destroy(pDeltaGammaLabel_preReloc);

  delete d_eos;
  delete d_shear;
  delete d_yield;
}

void
CamClay::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "camclay");
  }

  d_eos->outputProblemSpec(cm_ps);
  d_shear->outputProblemSpec(cm_ps);
  d_yield->outputProblemSpec(cm_ps);
  d_intvar->outputProblemSpec(cm_ps);
}

std::unique_ptr<ConstitutiveModel>
CamClay::clone()
{
  return std::make_unique<CamClay>(*this);
}

void
CamClay::initializeLocalMPMLabels()
{
  pStrainLabel = VarLabel::create(
    "p.strain", ParticleVariable<Matrix3>::getTypeDescription());
  pElasticStrainLabel = VarLabel::create(
    "p.elasticStrain", ParticleVariable<Matrix3>::getTypeDescription());
  pDeltaGammaLabel = VarLabel::create(
    "p.deltaGamma", ParticleVariable<double>::getTypeDescription());

  pStrainLabel_preReloc = VarLabel::create(
    "p.strain+", ParticleVariable<Matrix3>::getTypeDescription());
  pElasticStrainLabel_preReloc = VarLabel::create(
    "p.elasticStrain+", ParticleVariable<Matrix3>::getTypeDescription());
  pDeltaGammaLabel_preReloc = VarLabel::create(
    "p.deltaGamma+", ParticleVariable<double>::getTypeDescription());
}

void
CamClay::addParticleState(std::vector<const VarLabel*>& from,
                          std::vector<const VarLabel*>& to)
{
  // Add the local particle state data for this constitutive model.
  from.push_back(pStrainLabel);
  from.push_back(pElasticStrainLabel);
  from.push_back(pDeltaGammaLabel);

  to.push_back(pStrainLabel_preReloc);
  to.push_back(pElasticStrainLabel_preReloc);
  to.push_back(pDeltaGammaLabel_preReloc);

  // Add the particle state for the internal variable models
  d_intvar->addParticleState(from, to);
}

void
CamClay::addInitialComputesAndRequires(Task* task,
                                       const MPMMaterial* matl,
                                       const PatchSet* patch) const
{
  const MaterialSubset* matlset = matl->thisMaterial();

  task->computes(pStrainLabel, matlset);
  task->computes(pElasticStrainLabel, matlset);
  task->computes(pDeltaGammaLabel, matlset);

  // Add internal evolution variables computed by internal variable model
  d_intvar->addInitialComputesAndRequires(task, matl, patch);
}

void
CamClay::initializeCMData(const Patch* patch,
                          const MPMMaterial* matl,
                          DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);
  computeStableTimestep(patch, matl, new_dw);

  // Put stuff in here to initialize each particle's
  // constitutive model parameters
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  ParticleVariable<Matrix3> pStrain, pElasticStrain;
  ParticleVariable<double> pDeltaGamma;

  new_dw->allocateAndPut(pStrain, pStrainLabel, pset);
  new_dw->allocateAndPut(pElasticStrain, pElasticStrainLabel, pset);
  new_dw->allocateAndPut(pDeltaGamma, pDeltaGammaLabel, pset);

  for (int& iter : *pset) {

    pStrain[iter]        = CamClay::zero;
    pElasticStrain[iter] = CamClay::zero;
    pDeltaGamma[iter]    = 0.0;
  }

  // Initialize the data for the internal variable model
  d_intvar->initializeInternalVariable(pset, new_dw);
}

void
CamClay::computeStableTimestep(const Patch* patch,
                               const MPMMaterial* matl,
                               DataWarehouse* new_dw)
{
  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of computeStressTensor
  Vector dx     = patch->dCell();
  int matlindex = matl->getDWIndex();

  // Retrieve the array of constitutive parameters
  ParticleSubset* pset = new_dw->getParticleSubset(matlindex, patch);

  constParticleVariable<double> pMass, pVol_new;
  constParticleVariable<Vector> pVelocity;

  new_dw->get(pMass, lb->pMassLabel, pset);
  new_dw->get(pVol_new, lb->pVolumeLabel, pset);
  new_dw->get(pVelocity, lb->pVelocityLabel, pset);

  double c_dil = 0.0;
  Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);

  double shear = d_shear->computeInitialShearModulus();
  double bulk  = d_eos->computeInitialBulkModulus();

  ParticleSubset::iterator iter = pset->begin();
  for (; iter != pset->end(); iter++) {
    particleIndex idx = *iter;

    // Compute wave speed at each particle, store the maximum
    Vector pVelocity_idx = pVelocity[idx];
    if (pMass[idx] > 0) {
      // ** WARNING ** assuming incrementally linear elastic
      //               this is the volumetric wave speed
      c_dil = sqrt((bulk + 4.0 * shear / 3.0) * pVol_new[idx] / pMass[idx]);
    } else {
      c_dil         = 0.0;
      pVelocity_idx = Vector(0.0, 0.0, 0.0);
    }
    waveSpeed = Vector(Max(c_dil + std::abs(pVelocity_idx.x()), waveSpeed.x()),
                       Max(c_dil + std::abs(pVelocity_idx.y()), waveSpeed.y()),
                       Max(c_dil + std::abs(pVelocity_idx.z()), waveSpeed.z()));
  }

  waveSpeed       = dx / waveSpeed;
  double delT_new = waveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

void
CamClay::addComputesAndRequires(Task* task,
                                const MPMMaterial* matl,
                                const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  Ghost::GhostType gnone        = Ghost::None;
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForHypoExplicit(task, matlset, patches);

  // Other constitutive model and input dependent computes and requires
  task->requires(Task::OldDW, pStrainLabel, matlset, gnone);
  task->requires(Task::OldDW, pElasticStrainLabel, matlset, gnone);
  task->requires(Task::OldDW, pDeltaGammaLabel, matlset, gnone);

  task->computes(pStrainLabel_preReloc, matlset);
  task->computes(pElasticStrainLabel_preReloc, matlset);
  task->computes(pDeltaGammaLabel_preReloc, matlset);

  // Add internal evolution variables computed by internal variable model
  d_intvar->addComputesAndRequires(task, matl, patches);
}

void
CamClay::computeStressTensor(const PatchSubset* patches,
                             const MPMMaterial* matl,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw)
{

  DOUTALL(g_cc_dbg, "Compute stress tensor");

  // Constants
  // double sqrtThreeTwo  = sqrt(1.5);
  // double sqrtTwoThird  = 1.0 / sqrtThreeTwo;
  Ghost::GhostType gac = Ghost::AroundCells;
  Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);
  double rho_0             = matl->getInitialDensity();
  double totalStrainEnergy = 0.0;

  // Strain variables  (declared later)
  // Matrix3 strain(0.0);                  // Total strain
  // double strain_v = 0.0;                // Volumeric strain (eps_v)
  // Matrix3 strain_dev(0.0);              // Deviatoric strain (e)
  // double strain_dev_norm = 0.0;         // ||e||
  // double strain_s = 0.0;                // eps_s = sqrt(2/3) ||e||

  // Matrix3 strain_elast_tr(0.0);         // Trial elastic strain
  // double strain_elast_v_tr(0.0);        // Trial volumetric elastic strain
  // Matrix3 strain_elast_devtr(0.0);      // Trial deviatoric elastic strain
  // double strain_elast_devtr_norm = 0.0; // ||ee||
  // double strain_elast_s_tr = 0.0;       // epse_s = sqrt(2/3) ||ee||

  // double strain_elast_v_n = 0.0;        // last volumetric elastic strain
  // Matrix3 strain_elast_dev_n(0.0);      // last devaitoric elastic strain
  // double strain_elast_dev_n_norm = 0.0;
  // double strain_elast_s_n = 0.0;

  // Plasticity related variables
  // Matrix3 nn(0.0);                    // Plastic flow direction n = ee/||ee||

  // Loop thru patches
  for (int patchIndex = 0; patchIndex < patches->size(); patchIndex++) {
    const Patch* patch = patches->get(patchIndex);

    auto interpolator = flag->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());
    std::vector<double> S(interpolator->size());

    // Get grid size
    Vector dx = patch->dCell();
    // double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};

    // Get the set of particles
    int dwi              = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

    // GET GLOBAL DATA

    // Get the deformation gradient (F) and velocity gradient (L)
    constParticleVariable<Matrix3> pDefGrad_old, pDefGrad_new;
    constParticleVariable<Matrix3> pVelGrad_old, pVelGrad_new;
    old_dw->get(pDefGrad_old, lb->pDefGradLabel, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);
    old_dw->get(pVelGrad_old, lb->pVelGradLabel, pset);
    new_dw->get(pVelGrad_new, lb->pVelGradLabel_preReloc, pset);

    // Get the particle location, particle size, particle mass, particle volume
    constParticleVariable<Point> px;
    constParticleVariable<Matrix3> pSize;
    constParticleVariable<double> pMass, pVol_new;
    old_dw->get(px, lb->pXLabel, pset);
    old_dw->get(pSize, lb->pSizeLabel, pset);
    old_dw->get(pMass, lb->pMassLabel, pset);
    new_dw->get(pVol_new, lb->pVolumeLabel_preReloc, pset);

    // Get the velocity from the grid and particle velocity
    constParticleVariable<Vector> pVelocity;
    constNCVariable<Vector> gVelocity;
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);
    new_dw->get(gVelocity, lb->gVelocityStarLabel, dwi, patch, gac, NGN);

    // Get the particle stress
    constParticleVariable<Matrix3> pStress_old;
    old_dw->get(pStress_old, lb->pStressLabel, pset);

    // Get the time increment (delT)
    delt_vartype delT_orig;
    old_dw->get(delT_orig, lb->delTLabel, getLevel(patches));

    // GET LOCAL DATA
    constParticleVariable<Matrix3> pStrain_old, pElasticStrain_old;
    constParticleVariable<double> pDeltaGamma_old;
    old_dw->get(pStrain_old, pStrainLabel, pset);
    old_dw->get(pElasticStrain_old, pElasticStrainLabel, pset);
    old_dw->get(pDeltaGamma_old, pDeltaGammaLabel, pset);

    // Create and allocate arrays for storing the updated information
    // GLOBAL
    ParticleVariable<Matrix3> pStress_new;
    ParticleVariable<double> pdTdt, p_q;

    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);

    // LOCAL
    ParticleVariable<Matrix3> pStrain_new, pElasticStrain_new;
    ParticleVariable<double> pDeltaGamma_new;
    new_dw->allocateAndPut(pStrain_new, pStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pElasticStrain_new, pElasticStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pDeltaGamma_new, pDeltaGammaLabel_preReloc, pset);

    // Get the internal variable and allocate space for the updated internal
    // variables
    constParticleVariable<double> pPc;
    ParticleVariable<double> pPc_new;
    d_intvar->getInternalVariable(pset, old_dw, pPc);
    d_intvar->allocateAndPutInternalVariable(pset, new_dw, pPc_new);

    // Loop thru particles
    for (auto idx : *pset) {

      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;

      //-----------------------------------------------------------------------
      // Stage 1: Compute rate of deformation
      //-----------------------------------------------------------------------
      Matrix3 defGrad_new   = pDefGrad_new[idx];
      double J_new          = defGrad_new.Determinant();
      Matrix3 velGrad_new   = pVelGrad_new[idx];
      Matrix3 rateOfDef_new = (velGrad_new + velGrad_new.Transpose()) * 0.5;

      // Calculate the current mass density and deformed volume
      double rho_cur = rho_0 / J_new;

      // Compute polar decompositions of F_old and F_new (F = RU)
      Matrix3 rightStretch_old;
      rightStretch_old.Identity();
      Matrix3 rotation_old;
      rotation_old.Identity();
      pDefGrad_old[idx].polarDecompositionRMB(rightStretch_old, rotation_old);
      Matrix3 rightStretch_new;
      rightStretch_new.Identity();
      Matrix3 rotation_new;
      rotation_new.Identity();
      defGrad_new.polarDecompositionRMB(rightStretch_new, rotation_new);

      // Unrotate the spatial rate of deformation tensor and elastic strain
      rateOfDef_new =
        (rotation_old.Transpose()) * (rateOfDef_new * rotation_old);
      Matrix3 elasticStrain_old =
        (rotation_old.Transpose()) * (pElasticStrain_old[idx] * rotation_old);

      DOUT(g_cc_dbg, "F = " << defGrad_new)
      DOUT(g_cc_dbg, "Eps_e = " << elasticStrain_old)

      // Set up substepping (default is 1 step)
      double c_dil{ 0.0 };
      double bulk{ 0.0 };
      double delT      = delT_orig;
      int num_substeps = 1;
      for (int step = 0; step < num_substeps; step++) {

        //-----------------------------------------------------------------------
        // Stage 2: Compute trial state
        //-----------------------------------------------------------------------
        auto [strain_elast_tr, strain_elast_devtr, state] =
          computeTrialState(idx,
                            matl,
                            delT,
                            rateOfDef_new,
                            elasticStrain_old,
                            pVol_new[idx],
                            pMass[idx],
                            pPc[idx],
                            rho_0,
                            rho_cur);

        //-----------------------------------------------------------------------
        // Stage 3: Elastic-plastic stress update
        //-----------------------------------------------------------------------

        // Compute plastic flow direction (n = ee/||ee||)
        // Magic to deal with small strains
        double small             = 1.0e-12;
        double strain_elast_v_tr = state.epse_v_tr;
        double strain_elast_s_tr = state.epse_s_tr;
        double oo_strain_elast_s_tr =
          (strain_elast_s_tr > small) ? 1.0 / strain_elast_s_tr : 1.0;
        Matrix3 nn = strain_elast_devtr * (sqrtTwoThird * oo_strain_elast_s_tr);

        // Calculate yield function
        auto [ftrial, status_trial] = d_yield->evalYieldCondition(&state);

        DOUT(g_cc_dbg, "ftrial = " << ftrial)

        small = 1.0e-8; // **WARNING** Should not be hard coded (use d_tol)
        if (ftrial > small) { // Plastic loading

          auto [stress_new,
                eps_new,
                epse_new,
                delgamma_new,
                pc_new,
                status,
                msg] = plasticUpdate(idx,
                                     matl,
                                     ftrial,
                                     elasticStrain_old,
                                     strain_elast_tr,
                                     strain_elast_v_tr,
                                     strain_elast_s_tr,
                                     nn,
                                     pDeltaGamma_old[idx],
                                     pPc[idx],
                                     state);

          if (status == SolveStatus::CONVERGENCE_FAILURE) {
            std::cout << msg << std::endl;
            std::cout << "Increasing number of substeps\n";
            step = 0;
            num_substeps *= 2;
            delT = delT_orig / num_substeps;
            continue;
            if (num_substeps > 64) {
              throw ConvergenceFailure(
                msg, num_substeps, delT, delT_orig, __FILE__, __LINE__);
            }
          } else {
            pStress_new[idx]        = stress_new;
            pStrain_new[idx]        = eps_new;
            pElasticStrain_new[idx] = epse_new;
            pDeltaGamma_new[idx]    = delgamma_new;
            pPc_new[idx]            = pc_new;
          }
        } else { // Elastic range
          auto [stress_new, eps_new, epse_new, delgamma_new, pc_new] =
            elasticUpdate(
              strain_elast_tr, nn, pDeltaGamma_old[idx], pPc[idx], state);

          pStress_new[idx]        = stress_new;
          pStrain_new[idx]        = eps_new;
          pElasticStrain_new[idx] = epse_new;
          pDeltaGamma_new[idx]    = delgamma_new;
          pPc_new[idx]            = pc_new;
        }

        // Compute the strain energy
        double W_vol      = d_eos->computeStrainEnergy(&state);
        double W_dev      = d_shear->computeStrainEnergy(&state);
        totalStrainEnergy = (W_vol + W_dev) * pVol_new[idx];

        // compute the local sound wave speed
        bulk  = d_eos->computeBulkModulus(rho_0, rho_cur);
        c_dil = std::sqrt((bulk + 4.0 * state.shearModulus / 3.0) / rho_cur);
      }

      //-----------------------------------------------------------------------
      // Stage 4:
      //-----------------------------------------------------------------------
      // Rotate back to spatial configuration
      pStress_new[idx] =
        rotation_new * (pStress_new[idx] * rotation_new.Transpose());
      pStrain_new[idx] =
        rotation_new * (pStrain_new[idx] * rotation_new.Transpose());
      pElasticStrain_new[idx] =
        rotation_new * (pElasticStrain_new[idx] * rotation_new.Transpose());

      // Compute wave speed at each particle, store the maximum
      Vector pVel = pVelocity[idx];
      waveSpeed   = Vector(Max(c_dil + std::abs(pVel.x()), waveSpeed.x()),
                         Max(c_dil + std::abs(pVel.y()), waveSpeed.y()),
                         Max(c_dil + std::abs(pVel.z()), waveSpeed.z()));

      // Compute artificial viscosity term
      if (flag->d_artificialViscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;
        double c_bulk = sqrt(bulk / rho_cur);
        double Dkk    = rateOfDef_new.Trace();
        p_q[idx]      = artificialBulkViscosity(Dkk, c_bulk, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.;
      }

      // DOUTALL(g_cc_dbg,
      //         "Compute stress tensor: idx: " << idx << " stress: "
      //                                        << pStress_new[idx]);
    } // end loop over particles

    waveSpeed       = dx / waveSpeed;
    double delT_new = waveSpeed.minComponent();

    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(totalStrainEnergy), lb->StrainEnergyLabel);
    }
  }

  DOUTALL(g_cc_dbg, "Compute stress tensor done: " << getpid());
}

std::tuple<Matrix3, Matrix3, Matrix3, double, double>
CamClay::elasticUpdate(const Matrix3& strain_elast_tr,
                       const Matrix3& nn,
                       double pDeltaGamma_old,
                       double pPc_old,
                       const ModelState_CamClay& state)
{
  // update stress from trial elastic strain
  auto stress = CamClay::one * state.p + nn * (CamClay::sqrtTwoThird * state.q);

  // update elastic strain from trial value
  auto total_strain   = strain_elast_tr;
  auto elastic_strain = strain_elast_tr;

  // update delta gamma (plastic strain increament)
  auto delta_gamma = pDeltaGamma_old;

  // Update consolidation pressure
  auto pc = pPc_old;

  return std::make_tuple(stress, total_strain, elastic_strain, delta_gamma, pc);
}

std::tuple<Matrix3,
           Matrix3,
           Matrix3,
           double,
           double,
           CamClay::SolveStatus,
           std::string>
CamClay::plasticUpdate(particleIndex idx,
                       const MPMMaterial* matl,
                       double ftrial,
                       const Matrix3& elasticStrain_old,
                       const Matrix3& strain_elast_tr,
                       double strain_elast_v_tr,
                       double strain_elast_s_tr,
                       const Matrix3& nn,
                       double pDeltaGamma_old,
                       double pPc_old,
                       ModelState_CamClay& state)
{
  // Calc volumetric and deviatoric elastic strains at beginninging of
  // timestep (t_n)
  double strain_elast_v_n = elasticStrain_old.Trace();
  Matrix3 strain_elast_dev_n =
    elasticStrain_old - CamClay::one * (strain_elast_v_n / 3.0);
  double strain_elast_dev_n_norm = strain_elast_dev_n.Norm();
  double strain_elast_s_n = CamClay::sqrtTwoThird * strain_elast_dev_n_norm;

  if (cout_CC_Eps.active()) {
    cout_CC_Eps << "CamClay: idx = " << idx
                << " t_n: eps_v_e = " << strain_elast_v_n
                << " eps_s_e = " << strain_elast_s_n << std::endl;
  }

  double strain_elast_v = strain_elast_v_n;
  double strain_elast_s = strain_elast_s_n;
  double delgamma       = pDeltaGamma_old;

  auto [status, iters, tol, val, msg] = doNewtonSolve(idx,
                                                      matl,
                                                      ftrial,
                                                      strain_elast_tr,
                                                      strain_elast_v_tr,
                                                      strain_elast_s_tr,
                                                      elasticStrain_old,
                                                      delgamma,
                                                      strain_elast_v,
                                                      strain_elast_s,
                                                      state);
  if (status == SolveStatus::INVALID_VALUE) {
    throw InvalidValue(msg, __FILE__, __LINE__);
  } else if (status == SolveStatus::CONVERGENCE_FAILURE) {
    return std::make_tuple(
      CamClay::zero, CamClay::zero, CamClay::zero, 0.0, 0.0, status, msg);
    // throw ConvergenceFailure(msg, iters, tol, val, __FILE__, __LINE__);
  }

  // update stress
  auto stress = CamClay::one * state.p + nn * (CamClay::sqrtTwoThird * state.q);

  // update elastic strain
  auto total_strain   = strain_elast_tr;
  auto elastic_strain = nn * (CamClay::sqrtThreeTwo * strain_elast_s) +
                        CamClay::one * (strain_elast_v / 3.0);

  // update delta gamma (plastic strain)
  auto delta_gamma = pDeltaGamma_old + delgamma;

  // Update consolidation pressure
  auto pc = state.p_c;

  return std::make_tuple(
    stress, total_strain, elastic_strain, delta_gamma, pc, status, msg);
}

std::tuple<Matrix3, Matrix3, ModelState_CamClay>
CamClay::computeTrialState(particleIndex idx,
                           const MPMMaterial* matl,
                           double delT,
                           const Matrix3& rateOfDef_new,
                           const Matrix3& elasticStrain_old,
                           double pVol_new,
                           double pMass,
                           double pPc,
                           double rho_0,
                           double rho_cur)
{
  // Compute strain increment from rotationally corrected rate of
  // deformation
  // (Forward Euler)
  Matrix3 strainInc = rateOfDef_new * delT;

  // Calculate the total strain
  //   Volumetric strain &  Deviatoric strain
  // Matrix3 strain = pStrain_old[idx] + strainInc;
  // double strain_v = strain.Trace();
  // Matrix3 strain_dev = strain - CamClay::one*(strain_v/3.0);

  // Trial elastic strain
  //   Volumetric elastic strain &  Deviatoric elastic strain
  Matrix3 strain_elast_tr  = elasticStrain_old + strainInc;
  double strain_elast_v_tr = strain_elast_tr.Trace();
  Matrix3 strain_elast_devtr =
    strain_elast_tr - CamClay::one * (strain_elast_v_tr / 3.0);
  double strain_elast_devtr_norm = strain_elast_devtr.Norm();
  double strain_elast_s_tr = CamClay::sqrtTwoThird * strain_elast_devtr_norm;

  // Set up the ModelState (for t_n)
  Vaango::ModelState_CamClay state;
  state.density                  = rho_cur;
  state.initialDensity           = rho_0;
  state.volume                   = pVol_new;
  state.initialVolume            = pMass / rho_0;
  state.elasticStrainTensor      = strain_elast_tr;
  state.epse_v                   = strain_elast_v_tr;
  state.epse_s                   = strain_elast_s_tr;
  state.elasticStrainTensorTrial = strain_elast_tr;
  state.epse_v_tr                = strain_elast_v_tr;
  state.epse_s_tr                = strain_elast_s_tr;
  state.p_c0                     = pPc;
  state.p_c                      = state.p_c0;

  // Compute mu and q
  double mu          = d_shear->computeShearModulus(&state);
  state.shearModulus = mu;
  double q           = d_shear->computeQ(&state);
  state.q            = q;

  if (std::isnan(q)) {
    std::ostringstream desc;
    desc << "idx = " << idx << " epse_v = " << state.epse_v
         << " epse_s = " << state.epse_s << " q = " << q << std::endl;
    throw InvalidValue(desc.str(), __FILE__, __LINE__);
  }

  // Compute p and bulk modulus
  double p =
    d_eos->computePressure(matl, &state, CamClay::zero, CamClay::zero, 0.0);
  state.p = p;

  if (std::isnan(p)) {
    std::ostringstream desc;
    desc << "idx = " << idx << " epse_v = " << state.epse_v
         << " epse_s = " << state.epse_s << " p = " << p << std::endl;
    throw InvalidValue(desc.str(), __FILE__, __LINE__);
  }

  // Get internal state variable (p_c)
  double pc_n = d_intvar->computeInternalVariable("dummy", nullptr, &state);
  state.p_c   = pc_n;

  DOUT(g_cc_dbg,
       "p = " << state.p << " q = " << state.q << " p_c = " << state.p_c);

  return std::make_tuple(strain_elast_tr, strain_elast_devtr, state);
}

std::tuple<CamClay::SolveStatus, int, double, double, std::string>
CamClay::doNewtonSolve(particleIndex idx,
                       const MPMMaterial* matl,
                       double ftrial,
                       const Matrix3& strain_elast_tr,
                       double strain_elast_v_tr,
                       double strain_elast_s_tr,
                       const Matrix3& elasticStrain_old,
                       double& delgamma,
                       double& strain_elast_v,
                       double& strain_elast_s,
                       Vaango::ModelState_CamClay& state)
{
  double p_old  = state.p;
  double q_old  = state.q;
  double pc_old = state.p_c;
  double f_old  = ftrial;

  // Newton iteration constants
  double tolr    = 1.0e-4; // 1e-4
  double tolf    = 1.0e-8; // 1.0e-8;
  int iter_break = 20;

  double fyield           = ftrial;
  double strain_elast_v_k = strain_elast_v;
  double strain_elast_s_k = strain_elast_s;
  double delgamma_k       = delgamma;

  // Derivatives
  double dfdp = d_yield->df_dp(&state);
  double dfdq = d_yield->df_dq(&state);

  // Residual
  double rv = strain_elast_v_k - strain_elast_v_tr + delgamma_k * dfdp;
  double rs = strain_elast_s_k - strain_elast_s_tr + delgamma_k * dfdq;
  double rf = fyield;

  // Set up Newton iterations
  double rtolv = 1.0;
  double rtols = 1.0;
  double rtolf = 1.0;
  int klocal   = 0;

  double mu{ 0.0 };
  double p{ 0.0 };
  double q{ 0.0 };
  double pc{ 0.0 };

  // Do Newton iterations
  state.elasticStrainTensor      = elasticStrain_old;
  state.elasticStrainTensorTrial = strain_elast_tr;
  double fmax                    = 1.0;
  while ((rtolv > tolr) || (rtols > tolr) || (rtolf > tolf)) {
    klocal++;

    auto [deldelgamma, delvoldev] =
      computeDeltaGammaIncrement(state, delgamma_k, dfdp, dfdq, rv, rs, rf);

    if (std::abs(delvoldev[0]) > 10.0 * std::abs(strain_elast_v_k) &&
        std::abs(delvoldev[0]) > 1.0e-4) {
      std::ostringstream desc;
      desc << "**WARNING** Increment of elastic volumetric strain too large.\n";
      desc << "delvoldev[0] = " << delvoldev[0]
           << " strain_elast_v_k = " << strain_elast_v_k << "\n";
      return std::make_tuple(SolveStatus::CONVERGENCE_FAILURE,
                             klocal,
                             delvoldev[0],
                             strain_elast_v_k,
                             desc.str());
    }

    // Allow for line search step (not fully implemented)
    double delgamma_local{ 0.0 };
    bool f_residuals_equal{ false };
    bool do_line_search = false;
    do {
      // update
      strain_elast_v = strain_elast_v_k + delvoldev[0];
      strain_elast_s = strain_elast_s_k + delvoldev[1];
      if (strain_elast_s < 0.0) {
        strain_elast_s = strain_elast_s_k;
      }
      delgamma_local = delgamma_k + deldelgamma;
      // if (delgamma_local < 0.0) delgamma_local = delgamma_k;

      state.epse_v = strain_elast_v;
      state.epse_s = strain_elast_s;

      mu = d_shear->computeShearModulus(&state);
      q  = d_shear->computeQ(&state);
      p =
        d_eos->computePressure(matl, &state, CamClay::zero, CamClay::zero, 0.0);
      pc = d_intvar->computeInternalVariable("dummy", nullptr, &state);

      if (std::isnan(p) || std::isnan(q) || std::isnan(pc)) {
        std::ostringstream desc;
        desc << "idx = " << idx << " k = " << klocal
             << " epse_v = " << state.epse_v << " epse_s = " << state.epse_s
             << " p = " << p << " q = " << q << " pc = " << pc
             << " f = " << fyield << std::endl;
        desc << " rf = " << rf << " rv = " << rv << " rs = " << rs << std::endl;
        desc << " deldelgamma = " << deldelgamma << " mu = " << mu << std::endl;
        desc << __FILE__ << " : " << __LINE__ << std::endl;
        return std::make_tuple(
          SolveStatus::INVALID_VALUE, iter_break, rtolf, tolf, desc.str());
      }

      // Update actual state
      state.shearModulus = mu;
      state.q            = q;
      state.p            = p;
      state.p_c          = pc;

      // compute updated derivatives
      dfdp = d_yield->df_dp(&state);
      dfdq = d_yield->df_dq(&state);

      // compute updated yield condition
      auto [fyield_inner, status_inner] = d_yield->evalYieldCondition(&state);
      fyield                            = fyield_inner;

#ifdef CATCH_NOT_FINITE
      if (!std::isfinite(fyield)) {
        std::ostringstream desc;
        desc << "mu = " << state.shearModulus << " p = " << state.p
             << " q = " << state.q << " p_c = " << state.p_c
             << " dfdp = " << dfdp << " dfdq = " << dfdq
             << " fyield = " << fyield << "\n";
        throw InvalidValue(desc.str(), __FILE__, __LINE__);
      }
#endif

      // Calculate max value of f
      auto fmax = d_yield->evalYieldConditionMax(&state);

      // save old residuals
      double rf_old = rf;
      // double rv_old = rv;
      // double rs_old = rs;

      // update residual
      rv = strain_elast_v - strain_elast_v_tr + delgamma_local * dfdp;
      rs = strain_elast_s - strain_elast_s_tr + delgamma_local * dfdq;
      rf = fyield;

      if (idx == 1) {
        DOUTALL(g_cc_conv,
                "idx = " << idx << " k = " << klocal << " rv = " << rv
                         << " rs = " << rs << " rf = " << rf
                         << " rf_old = " << rf_old << " fmax = " << fmax);
        DOUTALL(g_cc_conv,
                " rtolv = " << rtolv << " rtols = " << rtols
                            << " rtolf = " << rtolf << " tolr = " << tolr
                            << " tolf = " << tolf);
        DOUTALL(g_cc_conv,
                " pqpc = [" << p << " " << q << " " << pc << "]"
                            << " pqpc_old = [" << p_old << " " << q_old << " "
                            << pc_old << "]"
                            << " fold = " << f_old);
        DOUTALL(g_cc_conv,
                " epsv = " << strain_elast_v << " epss = " << strain_elast_s
                           << " f = " << fyield);
      }

      // std::cout << "rf = " << rf << " rf_old = " << rf_old << "\n";
      if (std::abs(rf / rf_old) < 1.0e-2) {
        do_line_search = false;
      } else {
        if ((std::abs(rf) > std::abs(rf_old)) ||
            (rf < 0.0 && std::abs(rf - rf_old) > fmax)) {
          // std::cout << "idx = " << idx << " rf = " << rf
          //           << " rf_old = " << rf_old << std::endl;
          do_line_search = true;
          delvoldev[0] *= 0.5;
          delvoldev[1] *= 0.5;
          deldelgamma *= 0.5;
        } else {
          do_line_search = false;
        }
      }

      if (std::abs(rf - rf_old) < tolr || std::abs(rf + rf_old) < tolr) {
        f_residuals_equal = true;
        do_line_search    = false;
      }

    } while (do_line_search);

    // Get ready for next iteration
    strain_elast_v_k = strain_elast_v;
    strain_elast_s_k = strain_elast_s;
    delgamma_k       = delgamma_local;

    // calculate tolerances
    rtolv = std::abs(delvoldev[0]);
    rtols = std::abs(delvoldev[1]);
    rtolf = std::abs(fyield / fmax);

    // Check max iters
    if (klocal == iter_break) {
      std::ostringstream desc;
      desc << "**ERROR** Newton iterations did not converge" << std::endl
           << " idx = " << idx << " rtolv = " << rtolv << "(tolr = " << tolr
           << ")"
           << " rtols = " << rtols << "(tolr = " << tolr << ")"
           << " rtolf = " << rtolf << "(tolf = " << tolf << ")"
           << " klocal = " << klocal << std::endl;
      // desc << " p_old = " << p << " q_old = " << q << " pc_old = " <<
      // pc_n << " f_old = " << ftrial << std::endl;
      desc << " p = " << p << " q = " << q << " pc = " << pc
           << " f = " << fyield << "(fmax = " << fmax << ")" << std::endl
           << " eps_v_e = " << strain_elast_v << " eps_s_e = " << strain_elast_s
           << std::endl;
#if 0
      desc << "L = " << velGrad_new << std::endl;
      desc << "F_old = " << pDefGrad_old[idx] << std::endl;
      desc << "F_new = " << defGrad_new << std::endl;
      desc << "J_old = " << pDefGrad_old[idx].Determinant() << std::endl;
      desc << "J_new = " << J_new << std::endl;
#endif
      desc << __FILE__ << " : " << __LINE__ << std::endl;
      return std::make_tuple(
        SolveStatus::CONVERGENCE_FAILURE, iter_break, rtolf, tolf, desc.str());
    }

#if 0
    if (cout_CC_Eps.active()) {
      if (idx == 14811) {
        Matrix3 eps_e = nn * (sqrtThreeTwo * strain_elast_s) +
                        CamClay::one * (strain_elast_v / 3.0);
        Matrix3 eps_p  = strain_elast_tr - eps_e;
        double eps_v_p = eps_p.Trace();
        cout_CC_Eps << "idx = " << idx << " k = " << klocal << std::endl
                    << " eps_v_e = " << strain_elast_v
                    << " eps_v_p = " << eps_v_p
                    << " eps_s_e = " << strain_elast_s << std::endl
                    << " f_n+1 = " << fyield << " pqpc = [" << p << " " << q
                    << " " << pc << "]"
                    << " rtolf = " << rtolf << " tolf = " << tolf << std::endl;
      }
    }
#endif

    delgamma = delgamma_local;
    if (f_residuals_equal) {
      break;
    }

  } // End of Newton-Raphson while

  if ((delgamma < 0.0) && (std::abs(delgamma) > 1.0e-10)) {
    std::ostringstream desc;
    desc << "**ERROR** delgamma less than 0.0 in local converged solution."
         << std::endl;
    desc << __FILE__ << " : " << __LINE__ << std::endl;
    return std::make_tuple(
      SolveStatus::CONVERGENCE_FAILURE, klocal, rtolf, delgamma, desc.str());
  }

  return std::make_tuple(SolveStatus::OK, klocal, rtolf, delgamma, "Converged");
}

std::tuple<double, std::vector<double>>
CamClay::computeDeltaGammaIncrement(Vaango::ModelState_CamClay& state,
                                    double delgamma_k,
                                    double dfdp,
                                    double dfdq,
                                    double rv,
                                    double rs,
                                    double rf)
{
  // Compute needed derivatives
  double dpdepsev = d_eos->computeDpDepse_v(&state);
  double dpdepses = d_eos->computeDpDepse_s(&state);
  double dqdepsev = dpdepses;
  double dqdepses = d_shear->computeDqDepse_s(&state);
  double dpcdepsev =
    d_intvar->computeVolStrainDerivOfInternalVariable("dummy", &state);

  // Compute derivatives of residuals
  double dr1_dx1 = 1.0 + delgamma_k * (2.0 * dpdepsev - dpcdepsev);
  double dr1_dx2 = 2.0 * delgamma_k * dpdepses;
  double dr1_dx3 = dfdp;

  double d2fdqdepsv = d_yield->d2f_dq_depsVol(&state, d_eos, d_shear);
  double d2fdqdepss = d_yield->d2f_dq_depsDev(&state, d_eos, d_shear);
  double dr2_dx1    = delgamma_k * d2fdqdepsv;
  double dr2_dx2    = 1.0 + delgamma_k * d2fdqdepss;
  double dr2_dx3    = dfdq;

  double dr3_dx1 = dfdq * dqdepsev + dfdp * dpdepsev - state.p * dpcdepsev;
  double dr3_dx2 = dfdq * dqdepses + dfdp * dpdepses;

#ifdef COMPARE_BORJA_TAMAGNINI
  double M            = 1.05;
  double alpha        = 60.0;
  double eps_v0_e     = 0.0;
  double kappa_tilde  = 0.018;
  double lambda_tilde = 0.13;
  double p0           = -9.0;
  Matrix3 Amat_BT     = computeBorjaTamagniniAmatrix(lambda_tilde,
                                                 kappa_tilde,
                                                 strain_elast_v_k,
                                                 eps_v0_e,
                                                 p,
                                                 pc,
                                                 q,
                                                 M,
                                                 mu,
                                                 p0,
                                                 alpha,
                                                 strain_elast_s_k,
                                                 delgamma_k);
  std::cout << "Amat_BT = " << Amat_BT << "\n";
#endif

  FastMatrix A_MAT(2, 2), inv_A_MAT(2, 2);
  A_MAT(0, 0) = dr1_dx1;
  A_MAT(0, 1) = dr1_dx2;
  A_MAT(1, 0) = dr2_dx1;
  A_MAT(1, 1) = dr2_dx2;

  inv_A_MAT.destructiveInvert(A_MAT);

  std::vector<double> B_MAT(2), C_MAT(2), AinvB(2), rvs_vec(2), Ainvrvs(2);
  B_MAT[0] = dr1_dx3;
  B_MAT[1] = dr2_dx3;

  C_MAT[0] = dr3_dx1;
  C_MAT[1] = dr3_dx2;

  rvs_vec[0] = rv;
  rvs_vec[1] = rs;

  inv_A_MAT.multiply(B_MAT, AinvB);

  inv_A_MAT.multiply(rvs_vec, Ainvrvs);

  double denom       = C_MAT[0] * AinvB[0] + C_MAT[1] * AinvB[1];
  double deldelgamma = 0.0;
  if (std::abs(denom) > 1e-20) {
    deldelgamma = (-C_MAT[0] * Ainvrvs[0] - C_MAT[1] * Ainvrvs[1] + rf) / denom;
  }

  std::vector<double> delvoldev(2);
  delvoldev[0] = -Ainvrvs[0] - AinvB[0] * deldelgamma;
  delvoldev[1] = -Ainvrvs[1] - AinvB[1] * deldelgamma;

  if (std::isnan(deldelgamma)) {
    std::ostringstream desc;
    desc << "A_MAT = ";
    A_MAT.print(desc);
    desc << "A_MAT_inv = ";
    inv_A_MAT.print(desc);
    desc << "B_MAT = " << B_MAT[0] << " " << B_MAT[1] << "\n";
    desc << "C_MAT = " << C_MAT[0] << " " << C_MAT[1] << "\n";
    desc << "rv = " << rv << " rs = " << rs << " rf = " << rf << "\n";
    desc << "deldelgamma = " << deldelgamma << " delvoldev = " << delvoldev[0]
         << " , " << delvoldev[1] << std::endl;
    desc << "delvoldev[0] = " << -Ainvrvs[0] << " - " << AinvB[0] * deldelgamma
         << "\n";
    desc << "state.q = " << state.q << " state.p " << state.p
         << " state.pc = " << state.p_c << "\n";
    throw InvalidValue(desc.str(), __FILE__, __LINE__);
  }

  return std::make_tuple(deldelgamma, delvoldev);
}

Matrix3
CamClay::computeBorjaTamagniniAmatrix(double lambda_tilde,
                                      double kappa_tilde,
                                      double eps_v_e,
                                      double eps_v0_e,
                                      double p,
                                      double p_c,
                                      double q,
                                      double M,
                                      double mu_e,
                                      double p0,
                                      double alpha,
                                      double eps_s_e,
                                      double deltaPhi)
{
  // Already converted internally
  double lambda_hat = lambda_tilde / (1.0 - lambda_tilde);
  double kappa_hat  = kappa_tilde / (1.0 - kappa_tilde);
  double theta      = 1.0 / (lambda_hat - kappa_hat);

  double Omega = -(eps_v_e - eps_v0_e) / kappa_hat;

  double del_p_f     = 2.0 * p - p_c;
  double del_q_f     = 2.0 * q / (M * M);
  double del_pc_f    = -p;
  double del2_p_pc_f = -1.0;
  double del2_q_pc_f = 0.0;
  double D11_e       = -p / kappa_hat;
  double D22_e       = 3.0 * mu_e;
  double mu_vol      = p0 * alpha * std::exp(Omega);
  double D12_e       = 3.0 * mu_vol * eps_s_e / kappa_hat;
  double D21_e       = D12_e;
  // double H11         = 2.0;
  // double H22         = 2.0 / (M * M);
  // double H12         = 0.0;
  // double H21         = H12;
  double G11 = 2.0 * D11_e;
  double G12 = 2.0 * D12_e;
  double G21 = 2.0 * D21_e / (M * M);
  double G22 = 2.0 * D22_e / (M * M);
  double K_p = theta * p_c;

  Matrix3 A_mat;
  A_mat(0, 0) = 1.0 + deltaPhi * (G11 + K_p * del2_p_pc_f);
  A_mat(0, 1) = deltaPhi * G12;
  A_mat(0, 2) = del_p_f;

  A_mat(1, 0) = deltaPhi * (G21 + K_p * del2_q_pc_f);
  A_mat(1, 1) = 1 + deltaPhi * G22;
  A_mat(1, 2) = del_q_f;

  A_mat(2, 0) = D11_e * del_p_f + D21_e * del_q_f + K_p * del_pc_f;
  A_mat(2, 1) = D12_e * del_p_f + D22_e * del_q_f;
  A_mat(2, 2) = 0.0;

  return A_mat;
}

void
CamClay::carryForward(const PatchSubset* patches,
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
    constParticleVariable<Matrix3> pStrain, pElasticStrain;
    constParticleVariable<double> pDeltaGamma;

    old_dw->get(pStrain, pStrainLabel, pset);
    old_dw->get(pElasticStrain, pElasticStrainLabel, pset);
    old_dw->get(pDeltaGamma, pDeltaGammaLabel, pset);

    ParticleVariable<Matrix3> pStrain_new, pElasticStrain_new;
    ParticleVariable<double> pDeltaGamma_new;

    new_dw->allocateAndPut(pStrain_new, pStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pElasticStrain_new, pElasticStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pDeltaGamma_new, pDeltaGammaLabel_preReloc, pset);

    // Get and copy the internal variables
    constParticleVariable<double> pPc;
    d_intvar->getInternalVariable(pset, old_dw, pPc);
    d_intvar->allocateAndPutRigid(pset, new_dw, pPc);

    for (int idx : *pset) {
      pStrain_new[idx]        = pStrain[idx];
      pElasticStrain_new[idx] = pElasticStrain[idx];
      pDeltaGamma_new[idx]    = pDeltaGamma[idx];
    }

    new_dw->put(delt_vartype(1.e10), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(0.), lb->StrainEnergyLabel);
    }
  }
}

void
CamClay::allocateCMDataAddRequires(Task* task,
                                   const MPMMaterial* matl,
                                   const PatchSet* patch,
                                   MPMLabel* lb) const
{
  const MaterialSubset* matlset = matl->thisMaterial();

  // Allocate the variables shared by all constitutive models
  // for the particle convert operation
  // This method is defined in the ConstitutiveModel base class.
  addSharedRForConvertExplicit(task, matlset, patch);

  // Add requires local to this model
  Ghost::GhostType gnone = Ghost::None;
  task->requires(Task::NewDW, pStrainLabel_preReloc, matlset, gnone);
  task->requires(Task::NewDW, pElasticStrainLabel_preReloc, matlset, gnone);
  task->requires(Task::NewDW, pDeltaGammaLabel_preReloc, matlset, gnone);
  d_intvar->allocateCMDataAddRequires(task, matl, patch, lb);
}

void
CamClay::allocateCMDataAdd(DataWarehouse* new_dw,
                           ParticleSubset* addset,
                           ParticleLabelVariableMap* newState,
                           ParticleSubset* delset,
                           DataWarehouse* old_dw)
{
  // Copy the data common to all constitutive models from the particle to be
  // deleted to the particle to be added.
  // This method is defined in the ConstitutiveModel base class.
  copyDelToAddSetForConvertExplicit(new_dw, delset, addset, newState);

  // Copy the data local to this constitutive model from the particles to
  // be deleted to the particles to be added
  ParticleSubset::iterator n, o;

  ParticleVariable<Matrix3> pStrain, pElasticStrain;
  ParticleVariable<double> pDeltaGamma;

  constParticleVariable<Matrix3> o_Strain, o_ElasticStrain;
  constParticleVariable<double> o_DeltaGamma;

  new_dw->allocateTemporary(pStrain, addset);
  new_dw->allocateTemporary(pElasticStrain, addset);
  new_dw->allocateTemporary(pDeltaGamma, addset);

  new_dw->get(o_Strain, pStrainLabel_preReloc, delset);
  new_dw->get(o_ElasticStrain, pElasticStrainLabel_preReloc, delset);
  new_dw->get(o_DeltaGamma, pDeltaGammaLabel_preReloc, delset);

  n = addset->begin();
  for (o = delset->begin(); o != delset->end(); o++, n++) {
    pStrain[*n]        = o_Strain[*o];
    pElasticStrain[*n] = o_ElasticStrain[*o];
    pDeltaGamma[*n]    = o_DeltaGamma[*o];
  }

  (*newState)[pStrainLabel]        = pStrain.clone();
  (*newState)[pElasticStrainLabel] = pElasticStrain.clone();
  (*newState)[pDeltaGammaLabel]    = pDeltaGamma.clone();

  // Initialize the data for the internal variable model
  d_intvar->allocateCMDataAdd(new_dw, addset, newState, delset, old_dw);
}

double
CamClay::computeRhoMicroCM(double pressure,
                           const double p_ref,
                           const MPMMaterial* matl,
                           double temperature,
                           double rho_guess)
{
  double rho_orig = matl->getInitialDensity();
  pressure -= p_ref;
  double rho_cur = d_eos->computeDensity(rho_orig, pressure);

  if (std::isnan(rho_cur)) {
    std::ostringstream desc;
    desc << "rho_cur = " << rho_cur << " pressure = " << pressure
         << " p_ref = " << p_ref << " rho_orig = " << rho_orig << std::endl;
    throw InvalidValue(desc.str(), __FILE__, __LINE__);
  }

  return rho_cur;
}

void
CamClay::computePressEOSCM(double rho_cur,
                           double& pressure,
                           double p_ref,
                           double& dp_drho,
                           double& csquared,
                           const MPMMaterial* matl,
                           double temperature)
{
  double rho_orig = matl->getInitialDensity();
  d_eos->computePressure(rho_orig, rho_cur, pressure, dp_drho, csquared);
  pressure += p_ref;

  if (std::isnan(pressure)) {
    std::ostringstream desc;
    desc << "rho_cur = " << rho_cur << " pressure = " << pressure
         << " p_ref = " << p_ref << " dp_drho = " << dp_drho << std::endl;
    throw InvalidValue(desc.str(), __FILE__, __LINE__);
  }
}

double
CamClay::getCompressibility()
{
  return 1.0 / d_eos->initialBulkModulus();
}

void
CamClay::scheduleCheckNeedAddMPMMaterial(Task* task,
                                         const MPMMaterial*,
                                         const PatchSet*) const
{
  task->computes(lb->NeedAddMPMMaterialLabel);
}

void
CamClay::checkNeedAddMPMMaterial(const PatchSubset* patches,
                                 const MPMMaterial* matl,
                                 DataWarehouse*,
                                 DataWarehouse* new_dw)
{
  DOUTALL(g_cc_dbg,
          getpid() << "checkNeedAddMPMMaterial: In : Matl = " << matl
                   << " id = " << matl->getDWIndex()
                   << " patch = " << (patches->get(0))->getID());

  double need_add = 0.;
  new_dw->put(sum_vartype(need_add), lb->NeedAddMPMMaterialLabel);
}
