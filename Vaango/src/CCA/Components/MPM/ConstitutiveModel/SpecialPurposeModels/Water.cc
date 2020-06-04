/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/SpecialPurposeModels/Water.h>
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
#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <iostream>

using namespace Uintah;

Water::Water(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  d_useModifiedEOS = false;
  ps->require("bulk_modulus", d_initialData.d_Bulk);
  ps->require("viscosity", d_initialData.d_Viscosity);
  ps->require("gamma", d_initialData.d_Gamma);
}

Water::Water(const Water* cm)
  : ConstitutiveModel(cm)
{
  d_useModifiedEOS = cm->d_useModifiedEOS;
  d_initialData.d_Bulk = cm->d_initialData.d_Bulk;
  d_initialData.d_Viscosity = cm->d_initialData.d_Viscosity;
  d_initialData.d_Gamma = cm->d_initialData.d_Gamma;
}

void
Water::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "water");
  }

  cm_ps->appendElement("bulk_modulus", d_initialData.d_Bulk);
  cm_ps->appendElement("viscosity", d_initialData.d_Viscosity);
  cm_ps->appendElement("gamma", d_initialData.d_Gamma);
}

Water*
Water::clone()
{
  return scinew Water(*this);
}

void
Water::initializeCMData(const Patch* patch, const MPMMaterial* matl,
                        DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);

  computeStableTimestep(patch, matl, new_dw);
}

void
Water::computeStableTimestep(const Patch* patch, const MPMMaterial* matl,
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

  double c_dil = 0.0;
  Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);

  double bulk = d_initialData.d_Bulk;
  for (int idx : *pset) {
    // Compute wave speed at each particle, store the maximum
    c_dil = sqrt((bulk)*pVolume[idx] / pMass[idx]);
    Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
    waveSpeed = Max(velMax, waveSpeed);
  }
  waveSpeed = dx / waveSpeed;
  double delT_new = waveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

void
Water::addParticleState(std::vector<const VarLabel*>&,
                        std::vector<const VarLabel*>&)
{
  // Add the local particle state data for this constitutive model.
}

void
Water::addComputesAndRequires(Task* task, const MPMMaterial* matl,
                              const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForExplicit(task, matlset, patches);
}

void
Water::computeStressTensor(const PatchSubset* patches, const MPMMaterial* matl,
                           DataWarehouse* old_dw, DataWarehouse* new_dw)
{
  double oneThird = (1.0 / 3.0);
  Matrix3 Identity;
  Identity.Identity();

  double viscosity = d_initialData.d_Viscosity;
  double bulk = d_initialData.d_Bulk;
  double gamma = d_initialData.d_Gamma;
  double rho_orig = matl->getInitialDensity();
  int matID = matl->getDWIndex();

  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  if (!flag->d_doGridReset) {
    std::cerr << "The water model doesn't work without resetting the grid" << endl;
  }

  for (int pp = 0; pp < patches->size(); pp++) {
    const Patch* patch = patches->get(pp);
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);
    Vector dx = patch->dCell();

    constParticleVariable<Vector> pVelocity;
    constParticleVariable<Matrix3> pDefGrad_new, pVelGrad;
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);
    new_dw->get(pVelGrad, lb->pVelGradLabel_preReloc, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    ParticleVariable<double> pdTdt, p_q;
    ParticleVariable<Matrix3> pStress;
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress, lb->pStressLabel_preReloc, pset);

    double strainEnergy = 0.;
    Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);

    for (int idx : *pset) {

      pdTdt[idx] = 0.0;

      // Viscous part of the stress
      Matrix3 D = (pVelGrad[idx] + pVelGrad[idx].Transpose()) * 0.5;
      Matrix3 DPrime = D - Identity * oneThird * D.Trace();
      Matrix3 shearStress = DPrime * (2. * viscosity);

      // Hydrostatic part of the stress
      double J = pDefGrad_new[idx].Determinant();
      double jtotheminusgamma = std::pow(J, -gamma);
      double p = bulk * (jtotheminusgamma - 1.0);

      // Compute the total stress (volumetric + deviatoric)
      pStress[idx] = Identity * (-p) + shearStress;

      // Update wave speed
      double rho_cur = rho_orig / J;
      double c_dil = std::sqrt((gamma * jtotheminusgamma * bulk) / rho_cur);
      Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
      waveSpeed = Max(velMax, waveSpeed);

      // Compute artificial viscosity term
      if (flag->d_artificialViscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;
        double c_bulk = sqrt(bulk / rho_cur);
        Matrix3 D = (pVelGrad[idx] + pVelGrad[idx].Transpose()) * 0.5;
        p_q[idx] = artificialBulkViscosity(D.Trace(), c_bulk, rho_cur, dx_ave);
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
  } // end loop over patches
}

void
Water::addComputesAndRequires(Task*, const MPMMaterial*, const PatchSet*,
                              const bool, const bool) const
{
}

void
Water::allocateCMDataAddRequires(Task* task, const MPMMaterial* matl,
                                 const PatchSet* patches, MPMLabel*) const
{
  const MaterialSubset* matlset = matl->thisMaterial();

  // Allocate the variables shared by all constitutive models
  // for the particle convert operation
  // This method is defined in the ConstitutiveModel base class.
  addSharedRForConvertExplicit(task, matlset, patches);
}

void
Water::allocateCMDataAdd(
  DataWarehouse* new_dw, ParticleSubset* addset,
  ParticleLabelVariableMap* newState,
  ParticleSubset* delset, DataWarehouse*)
{
  // Copy the data common to all constitutive models from the particle to be
  // deleted to the particle to be added.
  // This method is defined in the ConstitutiveModel base class.
  copyDelToAddSetForConvertExplicit(new_dw, delset, addset, newState);
}

void
Water::carryForward(const PatchSubset* patches, const MPMMaterial* matl,
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

// The "CM" versions use the pressure-volume relationship of the CNH model
double
Water::computeRhoMicroCM(double pressure, const double p_ref,
                         const MPMMaterial* matl, double temperature,
                         double rho_guess)
{
  double rho_orig = matl->getInitialDensity();
  double bulk = d_initialData.d_Bulk;

  double p_gauge = pressure - p_ref;
  double rho_cur;

  double p_g_over_bulk = p_gauge / bulk;
  rho_cur =
    rho_orig * (p_g_over_bulk + sqrt(p_g_over_bulk * p_g_over_bulk + 1.));

  return rho_cur;
}

void
Water::computePressEOSCM(const double rho_cur, double& pressure,
                         const double p_ref, double& dp_drho, double& tmp,
                         const MPMMaterial* matl, double temperature)
{
  double bulk = d_initialData.d_Bulk;
  double rho_orig = matl->getInitialDensity();

  double p_g = .5 * bulk * (rho_cur / rho_orig - rho_orig / rho_cur);
  pressure = p_ref + p_g;
  dp_drho = .5 * bulk * (rho_orig / (rho_cur * rho_cur) + 1. / rho_orig);
  tmp = bulk / rho_cur; // speed of sound squared
}

double
Water::getCompressibility()
{
  return 1.0 / d_initialData.d_Bulk;
}

