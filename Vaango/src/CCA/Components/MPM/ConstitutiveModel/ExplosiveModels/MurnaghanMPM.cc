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

#include <CCA/Components/MPM/ConstitutiveModel/ExplosiveModels/MurnaghanMPM.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/Constants.h>
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
#include <Core/ProblemSpec/ProblemSpec.h>
#include <fstream>
#include <iostream>

using namespace Uintah;

MurnaghanMPM::MurnaghanMPM(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  d_useModifiedEOS = false;
  ps->require("K", d_modelParam.d_K);
  ps->require("Kprime", d_modelParam.d_Kprime);
  ps->require("pressure", d_modelParam.d_P0);
  ps->require("rho0", d_modelParam.d_rho0);
  ps->require("viscosity", d_modelParam.d_viscosity);
}

MurnaghanMPM::MurnaghanMPM(const MurnaghanMPM* cm)
  : ConstitutiveModel(cm)
{
  d_useModifiedEOS         = cm->d_useModifiedEOS;
  d_modelParam.d_K         = cm->d_modelParam.d_K;
  d_modelParam.d_Kprime    = cm->d_modelParam.d_Kprime;
  d_modelParam.d_P0        = cm->d_modelParam.d_P0;
  d_modelParam.d_rho0      = cm->d_modelParam.d_rho0;
  d_modelParam.d_viscosity = cm->d_modelParam.d_viscosity;
}

void
MurnaghanMPM::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "murnaghanMPM");
  }

  cm_ps->appendElement("K", d_modelParam.d_K);
  cm_ps->appendElement("Kprime", d_modelParam.d_Kprime);
  cm_ps->appendElement("pressure", d_modelParam.d_P0);
  cm_ps->appendElement("rho0", d_modelParam.d_rho0);
  cm_ps->appendElement("viscosity", d_modelParam.d_viscosity);
}

MurnaghanMPM*
MurnaghanMPM::clone()
{
  return scinew MurnaghanMPM(*this);
}

void
MurnaghanMPM::addParticleState(std::vector<const VarLabel*>&,
                               std::vector<const VarLabel*>&)
{
  // Add the local particle state data for this constitutive model.
}

void
MurnaghanMPM::initializeCMData(const Patch* patch,
                               const MPMMaterial* matl,
                               DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);
  computeStableTimestep(patch, matl, new_dw);
}

void
MurnaghanMPM::computeStableTimestep(const Patch* patch,
                                    const MPMMaterial* matl,
                                    DataWarehouse* new_dw)
{
  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of computeStressTensor
  int matID            = matl->getDWIndex();
  ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);

  constParticleVariable<double> pMass, pVolume;
  constParticleVariable<Vector> pVelocity;
  new_dw->get(pMass, lb->pMassLabel, pset);
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pVelocity, lb->pVelocityLabel, pset);

  // Compute wave speed at each particle, store the maximum
  Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);
  double K = d_modelParam.d_K;
  for (int idx : *pset) {
    double c_dil  = std::sqrt((K)*pVolume[idx] / pMass[idx]);
    Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
    waveSpeed     = Max(velMax, waveSpeed);
  }
  Vector dx       = patch->dCell();
  waveSpeed       = dx / waveSpeed;
  double delT_new = waveSpeed.minComponent();

  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

void
MurnaghanMPM::addComputesAndRequires(Task* task,
                                     const MPMMaterial* matl,
                                     const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForExplicit(task, matlset, patches);
}

void
MurnaghanMPM::computeStressTensor(const PatchSubset* patches,
                                  const MPMMaterial* matl,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw)
{
  double K         = d_modelParam.d_K;
  double Kprime    = d_modelParam.d_Kprime;
  double rho_orig  = d_modelParam.d_rho0; // matl->getInitialDensity();
  double viscosity = d_modelParam.d_viscosity;

  int matID = matl->getDWIndex();

  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  for (int pp = 0; pp < patches->size(); pp++) {
    const Patch* patch   = patches->get(pp);
    Vector dx            = patch->dCell();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    //constParticleVariable<double> pVolume_new;
    constParticleVariable<Vector> pVelocity;
    constParticleVariable<Matrix3> pVelGrad_mid, pDefGrad_new;
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);
    //new_dw->get(pVolume_new, lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pVelGrad_mid, lb->pVelGradLabel_preReloc, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    ParticleVariable<double> pdTdt, p_q;
    ParticleVariable<Matrix3> pStress_new;
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);

    if (!flag->d_doGridReset) {
      std::cerr << "The Murnaghan model doesn't work without resetting the grid"
                << "\n";
    }

    Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);
    double strainEnergy = 0.0;
    for (int idx : *pset) {
      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;

      // Viscous part of the stress
      Matrix3 D = (pVelGrad_mid[idx] + pVelGrad_mid[idx].Transpose()) * 0.5;
      double trD = D.Trace();
      Matrix3 DPrime =
        D - Vaango::Util::Identity * Vaango::Util::one_third * trD;
      Matrix3 shearStress = DPrime * (2. * viscosity);

      // get the hydrostatic part of the stress
      double J        = pDefGrad_new[idx].Determinant();
      double pressure = (K/Kprime) * (std::pow(J, -Kprime) - 1.0);

      // compute the total stress (volumetric + deviatoric)
      pStress_new[idx] = Vaango::Util::Identity * (-pressure) + shearStress;

      // compute wave speed
      //double V0      = pVolume_new[idx] / J;
      double rho_cur = rho_orig / J;
      double bulk    = K * std::pow(J, -(Kprime+1));
      double c_bulk  = std::sqrt(bulk / rho_cur);
      Vector velMax  = pVelocity[idx].cwiseAbs() + c_bulk;
      waveSpeed      = Max(velMax, waveSpeed);

      //std::cout << "J = " << J << " bulk = " << bulk << " c = " << c_bulk << "\n";

      // Compute artificial viscosity term
      if (flag->d_artificialViscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;
        p_q[idx] = artificialBulkViscosity(trD, c_bulk, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.;
      }
    } // end loop over particles

    waveSpeed       = dx / waveSpeed;
    double delT_new = waveSpeed.minComponent();
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());

    // No strain energy is calculated for this model
    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(strainEnergy), lb->StrainEnergyLabel);
    }
  }
}

void
MurnaghanMPM::addComputesAndRequires(Task*,
                                     const MPMMaterial*,
                                     const PatchSet*,
                                     const bool,
                                     const bool) const
{
}

void
MurnaghanMPM::allocateCMDataAddRequires(Task* task,
                                        const MPMMaterial* matl,
                                        const PatchSet* patches,
                                        MPMLabel*) const
{
  const MaterialSubset* matlset = matl->thisMaterial();

  // Allocate the variables shared by all constitutive models
  // for the particle convert operation
  // This method is defined in the ConstitutiveModel base class.
  addSharedRForConvertExplicit(task, matlset, patches);
}

void
MurnaghanMPM::allocateCMDataAdd(DataWarehouse* new_dw,
                                ParticleSubset* addset,
                                ParticleLabelVariableMap* newState,
                                ParticleSubset* delset,
                                DataWarehouse*)
{
  // Copy the data common to all constitutive models from the particle to be
  // deleted to the particle to be added.
  // This method is defined in the ConstitutiveModel base class.
  copyDelToAddSetForConvertExplicit(new_dw, delset, addset, newState);

  // Copy the data local to this constitutive model from the particles to
  // be deleted to the particles to be added
}

void
MurnaghanMPM::carryForward(const PatchSubset* patches,
                           const MPMMaterial* matl,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch   = patches->get(p);
    int matID            = matl->getDWIndex();
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
MurnaghanMPM::computeRhoMicroCM(double pressure,
                                const double p_ref,
                                const MPMMaterial* matl,
                                double temperature,
                                double rho_guess)
{
  double rhoM;
  double rho_orig = d_modelParam.d_rho0; // matl->getInitialDensity();
  double K        = d_modelParam.d_K;
  double Kprime   = d_modelParam.d_Kprime;
  double P0       = d_modelParam.d_P0;

  if (pressure >= P0) {
    rhoM = rho_orig * std::pow((1.0 + Kprime/K * (pressure - P0)), 1. / Kprime);
  } else {
    rhoM = rho_orig * std::pow((pressure / P0), (K / P0));
  }

  return rhoM;
}

void
MurnaghanMPM::computePressEOSCM(const double rhoM,
                                double& pressure,
                                const double p_ref,
                                double& dp_drho,
                                double& tmp,
                                const MPMMaterial* matl,
                                double temperature)
{
  double rho_orig = matl->getInitialDensity();
  double K        = d_modelParam.d_K;
  double Kprime   = d_modelParam.d_Kprime;
  double P0       = d_modelParam.d_P0;

  // Pointwise computation of thermodynamic quantities
  if (rhoM >= rho_orig) {
    pressure = P0 + (K / Kprime) * (std::pow(rhoM / rho_orig, Kprime) - 1.);
    dp_drho  = (K / rho_orig) * std::pow((rhoM / rho_orig), Kprime - 1.);
  } else {
    pressure = P0 * std::pow(rhoM / rho_orig, (K / P0));
    dp_drho = ((K / rho_orig)) * std::pow(rhoM / rho_orig, (K / P0 - 1));
  }
}

double
MurnaghanMPM::getCompressibility()
{
  return 1.0 / d_modelParam.d_K;
}
