/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2022 Parresia Research Limited, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/HypoElastic.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Labels/MPMLabel.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/Math/Short27.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <fstream>
#include <iostream>

using namespace Uintah;

HypoElastic::HypoElastic(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  ps->require("G", d_modelParam.G);
  ps->require("K", d_modelParam.K);

  // Thermal expansion coefficient
  d_modelParam.alpha = 0.0;
  ps->get("alpha", d_modelParam.alpha); // for thermal stress
}

HypoElastic::HypoElastic(const HypoElastic* cm)
  : ConstitutiveModel(cm)
{
  d_modelParam.G = cm->d_modelParam.G;
  d_modelParam.K = cm->d_modelParam.K;
  d_modelParam.alpha = cm->d_modelParam.alpha; // for thermal stress
}

void
HypoElastic::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "hypo_elastic");
  }

  cm_ps->appendElement("G", d_modelParam.G);
  cm_ps->appendElement("K", d_modelParam.K);
  cm_ps->appendElement("alpha", d_modelParam.alpha);
}

HypoElastic*
HypoElastic::clone()
{
  return scinew HypoElastic(*this);
}

void
HypoElastic::addParticleState(std::vector<const VarLabel*>& from,
                              std::vector<const VarLabel*>& to)
{
  // This is an INCREMENTAL model. Needs polar decomp R to be saved.
  from.push_back(lb->pPolarDecompRLabel);
  to.push_back(lb->pPolarDecompRLabel_preReloc);
}

void
HypoElastic::initializeCMData(const Patch* patch, const MPMMaterial* matl,
                              DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);
  computeStableTimestep(patch, matl, new_dw);
}

void
HypoElastic::computeStableTimestep(const Patch* patch, const MPMMaterial* matl,
                                   DataWarehouse* new_dw)
{
  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of computeStressTensor
  Vector dx = patch->dCell();
  int matID = matl->getDWIndex();
  ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);
  constParticleVariable<double> pMass, pVolume;
  constParticleVariable<Vector> pVelocity;

  new_dw->get(pMass, lb->pMassLabel, pset);
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pVelocity, lb->pVelocityLabel, pset);

  double G = d_modelParam.G;
  double bulk = d_modelParam.K;

  Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);
  for (int idx : *pset) {
    // Compute wave speed at each particle, store the maximum
    double c_dil = sqrt((bulk + 4. * G / 3.) * pVolume[idx] / pMass[idx]);
    Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
    waveSpeed = Max(velMax, waveSpeed);
  }
  waveSpeed = dx / waveSpeed;
  double delT_new = waveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());

  //std::cout << "Hypoelastic init: delT = " << delT_new << "\n";
}

void
HypoElastic::addComputesAndRequires(Task* task, const MPMMaterial* matl,
                                    const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addComputesAndRequiresForRotatedExplicit(task, matlset, patches);

  Ghost::GhostType gnone = Ghost::None;
  // for thermal stress
  task->requires(Task::OldDW, lb->pTempPreviousLabel, matlset, gnone);
}

void
HypoElastic::computeStressTensor(const PatchSubset* patches,
                                 const MPMMaterial* matl, DataWarehouse* old_dw,
                                 DataWarehouse* new_dw)
{
  double onethird = (1.0 / 3.0);
  Matrix3 Identity, zero(0.), One(1.);
  Identity.Identity();

  double rho_orig = matl->getInitialDensity();
  double mu = d_modelParam.G;
  double lambda = d_modelParam.K - 2.0 / 3.0 * d_modelParam.G;
  double M_dil = lambda + 2.0 * mu;
  double alpha = d_modelParam.alpha; // for thermal stress

  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    Vector dx = patch->dCell();
    int matID = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    constParticleVariable<Point> pX;
    constParticleVariable<double> pMass, pVolume, pTemperature, pTempPrevious;
    constParticleVariable<Vector> pVelocity;
    constParticleVariable<Matrix3> pSize, pDefGrad;
    old_dw->get(pX, lb->pXLabel, pset);
    old_dw->get(pMass, lb->pMassLabel, pset);
    old_dw->get(pVolume, lb->pVolumeLabel, pset);
    old_dw->get(pTempPrevious, lb->pTempPreviousLabel, pset);
    old_dw->get(pTemperature, lb->pTemperatureLabel, pset);
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);
    old_dw->get(pSize, lb->pSizeLabel, pset);
    old_dw->get(pDefGrad, lb->pDefGradLabel, pset);

    constParticleVariable<double> pVolume_new;
    constParticleVariable<Matrix3> pDefRate_mid, pDefGrad_new, pStress_old;
    new_dw->get(pVolume_new, lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pDefRate_mid, lb->pDeformRateMidLabel, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);
    new_dw->get(pStress_old, lb->pStressUnrotatedLabel, pset);

    ParticleVariable<double> pdTdt, p_q;
    ParticleVariable<Matrix3> pStress_new;
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);

    CCVariable<double> vol_0_CC, dvol_CC;
    CCVariable<int> ppc_CC;
    new_dw->allocateTemporary(vol_0_CC, patch);
    new_dw->allocateTemporary(dvol_CC, patch);
    new_dw->allocateTemporary(ppc_CC, patch);

    vol_0_CC.initialize(0.);
    dvol_CC.initialize(0.);
    ppc_CC.initialize(0);

    double strainEnergy = 0.0;
    Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);

    for (int idx : *pset) {
      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;

      // Rate of particle temperature change for thermal stress
      double pTempRate = (pTemperature[idx] - pTempPrevious[idx]) / delT;

      // Calculate rate of deformation D, and deviatoric rate DPrime,
      // including effect of thermal strain
      Matrix3 D = pDefRate_mid[idx] - Identity * (alpha * pTempRate);

      IntVector cell_index;
      patch->findCell(pX[idx], cell_index);

      vol_0_CC[cell_index] += pVolume[idx];
      dvol_CC[cell_index] += D.Trace() * pVolume[idx];
      ppc_CC[cell_index]++;
    }

    double press_stab = 0.;
    if (flag->d_doPressureStabilization) {
      press_stab = 1.;
      for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        dvol_CC[c] /= vol_0_CC[c];
      }
    }

    for (int idx : *pset) {

      double pTempRate = (pTemperature[idx] - pTempPrevious[idx]) / delT;

      // Calculate rate of deformation D, and deviatoric rate DPrime,
      // including effect of thermal strain
      IntVector cell_index;
      patch->findCell(pX[idx], cell_index);

      Matrix3 D = pDefRate_mid[idx] - Identity * (alpha * pTempRate);
      double trD = D.Trace();

      // Alter D to stabilize the pressure in each cell
      D += Identity * (onethird * press_stab * (dvol_CC[cell_index] - trD));
      trD = D.Trace();

      // This is the (updated) Cauchy stress
      pStress_new[idx] = pStress_old[idx] + 
        (Identity * (lambda * trD) + D * (2.0 * mu)) * delT;

      // Compute the local sound speed
      double J = pDefGrad_new[idx].Determinant();
      double rho_cur = rho_orig / J;
      double c_dil = std::sqrt(M_dil / rho_cur);

      // Compute the strain energy for all the particles
      Matrix3 avgStress = (pStress_new[idx] + pStress_old[idx]) * .5;
      double rateOfWork = computeRateOfWork(avgStress, pDefRate_mid[idx]);
      strainEnergy += (rateOfWork * pVolume_new[idx] * delT);

      // Compute wave speed at each particle, store the maximum
      Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
      waveSpeed     = Max(velMax, waveSpeed);

      // Compute artificial viscosity term
      if (flag->d_artificialViscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;
        double c_bulk = sqrt(d_modelParam.K / rho_cur);
        p_q[idx] = artificialBulkViscosity(trD, c_bulk, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.;
      }
    } // end loop over particles

    waveSpeed = dx / waveSpeed;
    double delT_new = waveSpeed.minComponent();
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());

    // std::cout << "Hypoelastic: delT = " << delT << "\n";

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(strainEnergy), lb->StrainEnergyLabel);
    }
  }
}

void
HypoElastic::carryForward(const PatchSubset* patches, const MPMMaterial* matl,
                          DataWarehouse* old_dw, DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    int matID = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    // Carry forward the data common to all constitutive models
    // when using RigidMPM.
    // This method is defined in the ConstitutiveModel base class.
    carryForwardSharedData(pset, old_dw, new_dw, matl);

    // Carry forward the data local to this constitutive model
    new_dw->put(delt_vartype(1.e10), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(0.), lb->StrainEnergyLabel);
    }
  }
}

void
HypoElastic::addComputesAndRequires(Task*, const MPMMaterial*, const PatchSet*,
                                    const bool, const bool) const
{
}

void
HypoElastic::allocateCMDataAddRequires(Task* task, const MPMMaterial* matl,
                                       const PatchSet* patches,
                                       MPMLabel* lb) const
{
  const MaterialSubset* matlset = matl->thisMaterial();

  // Allocate the variables shared by all constitutive models
  // for the particle convert operation
  // This method is defined in the ConstitutiveModel base class.
  addSharedRForConvertExplicit(task, matlset, patches);
}

void
HypoElastic::allocateCMDataAdd(DataWarehouse* new_dw, ParticleSubset* addset,
                               ParticleLabelVariableMap* newState,
                               ParticleSubset* delset, DataWarehouse*)
{
  // Copy the data common to all constitutive models from the particle to be
  // deleted to the particle to be added.
  // This method is defined in the ConstitutiveModel base class.
  copyDelToAddSetForConvertExplicit(new_dw, delset, addset, newState);
}

double
HypoElastic::computeRhoMicroCM(double pressure, const double p_ref,
                               const MPMMaterial* matl, double temperature,
                               double rho_guess)
{
  double rho_orig = matl->getInitialDensity();
  // double p_ref=101325.0;
  double p_gauge = pressure - p_ref;
  double rho_cur;
  // double G = d_modelParam.G;
  double bulk = d_modelParam.K;

  rho_cur = rho_orig / std::max(1.0e-12, (1 - p_gauge / bulk));

  return rho_cur;
}

void
HypoElastic::computePressEOSCM(double rho_cur, double& pressure, double p_ref,
                               double& dp_drho, double& csquared,
                               const MPMMaterial* matl, double temperature)
{

  // double G = d_modelParam.G;
  double bulk = d_modelParam.K;
  double rho_orig = matl->getInitialDensity();

  double p_g = bulk * (1.0 - rho_orig / rho_cur);
  pressure = p_ref + p_g;
  dp_drho = bulk * rho_orig / (rho_cur * rho_cur);
  csquared = bulk / rho_cur; // speed of sound squared
}

double
HypoElastic::getCompressibility()
{
  return 1.0 / d_modelParam.K;
}

