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

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/HypoElasticFortran.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/Constants.h>
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
#include<CCA/Components/MPM/Core/MPMLabel.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <Core/Malloc/Allocator.h>
#include <Core/Math/MinMax.h>

#include <sci_defs/uintah_defs.h>

#include <fstream>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////
// The following functions are found in fortran/*.F

extern "C" {

#if defined(FORTRAN_UNDERSCORE_END)
#define HOOKECHK hookechk_
#define HOOKE_INCREMENTAL hooke_incremental_
#elif defined(FORTRAN_UNDERSCORE_LINUX)
#define HOOKECHK hookechk_
#define HOOKE_INCREMENTAL hooke_incremental__
#else // NONE
#define HOOKECHK hookechk
#define HOOKE_INCREMENTAL hooke_incremental
#endif

void HOOKECHK(double UI[], double UJ[], double UK[]);
void HOOKE_INCREMENTAL(int& nblk, int& ninsv, double& dt, double UI[],
                       double stress[], double D[], double svarg[],
                       double& USM);
}

// End fortran functions.
////////////////////////////////////////////////////////////////////////////////

using std::cerr;
using namespace Uintah;

HypoElasticFortran::HypoElasticFortran(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  ps->require("G", d_modelParam.G);
  ps->require("K", d_modelParam.K);

  double UI[2];
  UI[0] = d_modelParam.K;
  UI[1] = d_modelParam.G;
  HOOKECHK(UI, UI, UI);
}

HypoElasticFortran::HypoElasticFortran(const HypoElasticFortran* cm)
  : ConstitutiveModel(cm)
{
  d_modelParam.G = cm->d_modelParam.G;
  d_modelParam.K = cm->d_modelParam.K;
}

std::unique_ptr<ConstitutiveModel>
HypoElasticFortran::clone()
{
  return std::make_unique<HypoElasticFortran>(*this);
}

void
HypoElasticFortran::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "hypo_elastic_fortran");
  }

  cm_ps->appendElement("G", d_modelParam.G);
  cm_ps->appendElement("K", d_modelParam.K);
}

void
HypoElasticFortran::addParticleState(std::vector<const VarLabel*>& from,
                                     std::vector<const VarLabel*>& to)
{
  // This is an INCREMENTAL model. Needs polar decomp R to be saved.
  from.push_back(lb->pPolarDecompRLabel);
  to.push_back(lb->pPolarDecompRLabel_preReloc);
}

void
HypoElasticFortran::initializeCMData(const Patch* patch,
                                     const MPMMaterial* matl,
                                     DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);
  computeStableTimestep(patch, matl, new_dw);
}

void
HypoElasticFortran::computeStableTimestep(const Patch* patch,
                                          const MPMMaterial* matl,
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
    double c_dil = std::sqrt((bulk + 4. * G / 3.) * pVolume[idx] / pMass[idx]);
    Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
    waveSpeed = Max(velMax, waveSpeed);
  }
  waveSpeed = dx / waveSpeed;
  double delT_new = waveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

void
HypoElasticFortran::addComputesAndRequires(Task* task, const MPMMaterial* matl,
                                           const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addComputesAndRequiresForRotatedExplicit(task, matlset, patches);
}

void
HypoElasticFortran::computeStressTensor(const PatchSubset* patches,
                                        const MPMMaterial* matl,
                                        DataWarehouse* old_dw,
                                        DataWarehouse* new_dw)
{
  double rho_orig = matl->getInitialDensity();
  int matID = matl->getDWIndex();

  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    Vector dx = patch->dCell();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    constParticleVariable<double> pVolume_new;
    constParticleVariable<Vector> pVelocity;
    constParticleVariable<Matrix3> pStress_old, pDefRate_mid, pDefGrad_new;
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);
    new_dw->get(pVolume_new, lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pStress_old, lb->pStressUnrotatedLabel, pset);
    new_dw->get(pDefRate_mid, lb->pDeformRateMidLabel, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    ParticleVariable<double> pdTdt, p_q;
    ParticleVariable<Matrix3> pStress_new;
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);

    double UI[2];
    UI[0] = d_modelParam.K;
    UI[1] = d_modelParam.G;

    double strainEnergy = 0.0;
    Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);

    for (int idx : *pset) {
      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;

      double sigarg[6];
      sigarg[0] = pStress_old[idx](0, 0);
      sigarg[1] = pStress_old[idx](1, 1);
      sigarg[2] = pStress_old[idx](2, 2);
      sigarg[3] = pStress_old[idx](0, 1);
      sigarg[4] = pStress_old[idx](1, 2);
      sigarg[5] = pStress_old[idx](2, 0);
      double Darray[6];
      Darray[0] = pDefRate_mid[idx](0, 0);
      Darray[1] = pDefRate_mid[idx](1, 1);
      Darray[2] = pDefRate_mid[idx](2, 2);
      Darray[3] = pDefRate_mid[idx](0, 1);
      Darray[4] = pDefRate_mid[idx](1, 2);
      Darray[5] = pDefRate_mid[idx](2, 0);

      double svarg[1];
      double USM = 9e99;
      double dt = delT;
      int nblk = 1;
      int ninsv = 1;
      HOOKE_INCREMENTAL(nblk, ninsv, dt, UI, sigarg, Darray, svarg, USM);

      pStress_new[idx](0, 0) = sigarg[0];
      pStress_new[idx](1, 1) = sigarg[1];
      pStress_new[idx](2, 2) = sigarg[2];
      pStress_new[idx](0, 1) = sigarg[3];
      pStress_new[idx](1, 0) = sigarg[3];
      pStress_new[idx](2, 1) = sigarg[4];
      pStress_new[idx](1, 2) = sigarg[4];
      pStress_new[idx](2, 0) = sigarg[5];
      pStress_new[idx](0, 2) = sigarg[5];

      // Compute the strain energy for all the particles
      Matrix3 avgStress = (pStress_new[idx] + pStress_old[idx]) * .5;
      double rateOfWork = computeRateOfWork(avgStress, pDefRate_mid[idx]);
      strainEnergy += (rateOfWork * pVolume_new[idx] * delT);

      // Compute wave speed at each particle, store the maximum
      double J = pDefGrad_new[idx].Determinant();
      double rho_cur = rho_orig / J;
      double c_dil = std::sqrt(USM / rho_cur);
      Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
      waveSpeed = Max(velMax, waveSpeed);

      // Compute artificial viscosity term
      if (flag->d_artificialViscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;
        double c_bulk = sqrt(UI[0] / rho_cur);
        p_q[idx] = artificialBulkViscosity(pDefRate_mid[idx].Trace(), 
                                           c_bulk, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.;
      }
    } // end loop over particles

    waveSpeed = dx / waveSpeed;
    double delT_new = waveSpeed.minComponent();
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(strainEnergy), lb->StrainEnergyLabel);
    }
  }
}

void
HypoElasticFortran::addComputesAndRequires(Task*, const MPMMaterial*,
                                           const PatchSet*, const bool,
                                           const bool) const
{
}

void
HypoElasticFortran::allocateCMDataAddRequires(Task* task,
                                              const MPMMaterial* matl,
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
HypoElasticFortran::allocateCMDataAdd(DataWarehouse* new_dw,
                                      ParticleSubset* addset,
                                      ParticleLabelVariableMap* newState,
                                      ParticleSubset* delset, DataWarehouse*)
{
  // Copy the data common to all constitutive models from the particle to be
  // deleted to the particle to be added.
  // This method is defined in the ConstitutiveModel base class.
  copyDelToAddSetForConvertExplicit(new_dw, delset, addset, newState);
}

void
HypoElasticFortran::carryForward(const PatchSubset* patches,
                                 const MPMMaterial* matl, DataWarehouse* old_dw,
                                 DataWarehouse* new_dw)
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

double
HypoElasticFortran::computeRhoMicroCM(double pressure, const double p_ref,
                                      const MPMMaterial* matl,
                                      double temperature, double rho_guess)
{
  double rho_orig = matl->getInitialDensity();
  // double p_ref=101325.0;
  double p_gauge = pressure - p_ref;
  double rho_cur;
  // double G = d_modelParam.G;
  double bulk = d_modelParam.K;

  rho_cur = rho_orig / (1 - p_gauge / bulk);

  return rho_cur;

#if 0
  std::cout << "NO VERSION OF computeRhoMicroCM EXISTS YET FOR HypoElasticFortran"
       << std::endl;
#endif
}

void
HypoElasticFortran::computePressEOSCM(double rho_cur, double& pressure,
                                      double p_ref, double& dp_drho,
                                      double& tmp, const MPMMaterial* matl,
                                      double temperature)
{
  // double G = d_modelParam.G;
  double bulk = d_modelParam.K;
  double rho_orig = matl->getInitialDensity();

  double p_g = bulk * (1.0 - rho_orig / rho_cur);
  pressure = p_ref + p_g;
  dp_drho = bulk * rho_orig / (rho_cur * rho_cur);
  tmp = bulk / rho_cur; // speed of sound squared

#if 0
  std::cout << "NO VERSION OF computePressEOSCM EXISTS YET FOR HypoElasticFortran"
       << std::endl;
#endif
}

double
HypoElasticFortran::getCompressibility()
{
  return 1.0 / d_modelParam.K;
}
