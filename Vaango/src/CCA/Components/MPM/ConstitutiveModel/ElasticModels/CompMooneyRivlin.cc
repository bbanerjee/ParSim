/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
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

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/CompMooneyRivlin.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Grid/Box.h>
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
#include <sci_values.h>


using namespace Uintah;

// Material Constants are C1, C2 and PR (poisson's ratio).
// The shear modulus = 2(C1 + C2).

CompMooneyRivlin::CompMooneyRivlin(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  ps->require("he_constant_1", d_initialData.C1);
  ps->require("he_constant_2", d_initialData.C2);
  ps->require("he_PR", d_initialData.PR);
}

CompMooneyRivlin::CompMooneyRivlin(const CompMooneyRivlin* cm)
  : ConstitutiveModel(cm)
{
  d_initialData.C1 = cm->d_initialData.C1;
  d_initialData.C2 = cm->d_initialData.C2;
  d_initialData.PR = cm->d_initialData.PR;
}

void
CompMooneyRivlin::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "comp_mooney_rivlin");
  }

  cm_ps->appendElement("he_constant_1", d_initialData.C1);
  cm_ps->appendElement("he_constant_2", d_initialData.C2);
  cm_ps->appendElement("he_PR", d_initialData.PR);
}

std::unique_ptr<ConstitutiveModel>
CompMooneyRivlin::clone()
{
  return std::make_unique<CompMooneyRivlin>(*this);
}

void
CompMooneyRivlin::addParticleState(std::vector<const VarLabel*>& from,
                                   std::vector<const VarLabel*>& to)
{
}

void
CompMooneyRivlin::initializeCMData(const Patch* patch, const MPMMaterial* matl,
                                   DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);

  computeStableTimestep(patch, matl, new_dw);
}

void
CompMooneyRivlin::computeStableTimestep(const Patch* patch,
                                        const MPMMaterial* matl,
                                        DataWarehouse* new_dw)
{
  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of computeStressTensor
  Vector dx = patch->dCell();
  int dwi = matl->getDWIndex();
  // Retrieve the array of constitutive parameters
  ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
  constParticleVariable<double> pMass, pVolume;
  constParticleVariable<Vector> pVelocity;

  new_dw->get(pMass, lb->pMassLabel, pset);
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pVelocity, lb->pVelocityLabel, pset);

  double c_dil = 0.0;
  Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);
  double C1 = d_initialData.C1;
  double C2 = d_initialData.C2;
  double PR = d_initialData.PR;

  for (int idx : *pset) {
    // Compute wave speed + particle velocity at each particle,
    // store the maximum
    double mu = 2. * (C1 + C2);
    // double C4 = .5*(C1*(5.*PR-2) + C2*(11.*PR-5)) / (1. - 2.*PR);
    c_dil =
      sqrt(2. * mu * (1. - PR) * pVolume[idx] / ((1. - 2. * PR) * pMass[idx]));
    Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
    waveSpeed = Max(velMax, waveSpeed);
  }
  waveSpeed = dx / waveSpeed;
  double delT_new = waveSpeed.minComponent();
  if (delT_new < 1.e-12)
    new_dw->put(delt_vartype(DBL_MAX), lb->delTLabel, patch->getLevel());
  else
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

void
CompMooneyRivlin::addComputesAndRequires(Task* task, const MPMMaterial* matl,
                                         const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForExplicit(task, matlset, patches);
}

void
CompMooneyRivlin::computeStressTensor(const PatchSubset* patches,
                                      const MPMMaterial* matl,
                                      DataWarehouse* old_dw,
                                      DataWarehouse* new_dw)
{
  Matrix3 Identity;
  Identity.Identity();

  double rho_orig = matl->getInitialDensity();
  double C1 = d_initialData.C1;
  double C2 = d_initialData.C2;
  double C3 = .5 * C1 + C2;
  double PR = d_initialData.PR;
  double C4 =
    .5 * (C1 * (5. * PR - 2) + C2 * (11. * PR - 5)) / (1. - 2. * PR);

  int dwi = matl->getDWIndex();

  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    Vector dx = patch->dCell();
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

    double strainEnergy = 0.0;
    Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);

    // Create array for the particle position
    constParticleVariable<double> pMass, pVolume;
    constParticleVariable<Vector> pVelocity;
    constParticleVariable<Matrix3> pVelGrad, pDefGrad_new;

    old_dw->get(pMass, lb->pMassLabel, pset);
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);
    new_dw->get(pVolume, lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pVelGrad, lb->pVelGradLabel_preReloc, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    ParticleVariable<double> pdTdt, p_q;
    ParticleVariable<Matrix3> pStress;
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress, lb->pStressLabel_preReloc, pset);

    for (int idx : *pset) {
      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;

      // Compute the left Cauchy-Green deformation tensor
      Matrix3 B = pDefGrad_new[idx] * pDefGrad_new[idx].Transpose();

      // Compute the invariants
      double I1_b = B.Trace();
      double I2_b = 0.5 * ((I1_b * I1_b) - (B * B).Trace());
      double J = pDefGrad_new[idx].Determinant();
      double I3_b = J * J;

      double w3 = -2.0 * C3 / (I3_b * I3_b * I3_b) + 2.0 * C4 * (I3_b - 1.0);

      // Compute T = 2/sqrt(I3)*(I3*W3*Identity + (W1+I1*W2)*B - W2*B^2)
      // W1 = C1, W2 = C2
      double C1pi1C2 = C1 + I1_b * C2;
      double i3w3 = I3_b * w3;

      pStress[idx] = (B * C1pi1C2 - (B * B) * C2 + Identity * i3w3) * 2.0 / J;

      // Compute wave speed + particle velocity at each particle,
      // store the maximum
      double c_dil = std::sqrt(
        (4. * (C1 + C2 * I2_b) / J +
         8. * (2. * C3 / (I3_b * I3_b * I3_b) + C4 * (2. * I3_b - 1.)) -
         Min((pStress[idx])(0, 0), (pStress[idx])(1, 1), (pStress[idx])(2, 2)) /
           J) *
        pVolume[idx] / pMass[idx]);
      Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
      waveSpeed = Max(velMax, waveSpeed);

      // Compute artificial viscosity term
      if (flag->d_artificialViscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;
        double bulk =
          (4. * (C1 + C2 * I2_b) / J); // I'm a little fuzzy here - JG
        double rho_cur = rho_orig / J;
        double c_bulk = sqrt(bulk / rho_cur);
        Matrix3 D = (pVelGrad[idx] + pVelGrad[idx].Transpose()) * 0.5;
        p_q[idx] = artificialBulkViscosity(D.Trace(), c_bulk, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.;
      }

      // Compute the strain energy for all the particles
      double pStrainEnergy = (C1 * (I1_b - 3.0) + C2 * (I2_b - 3.0) +
                              C3 * (1.0 / (I3_b * I3_b) - 1.0) +
                              C4 * (I3_b - 1.0) * (I3_b - 1.0)) *
                             pVolume[idx] / J;

      strainEnergy += pStrainEnergy;
    } // end loop over particles

    waveSpeed = dx / waveSpeed;
    double delT_new = waveSpeed.minComponent();

    if (delT_new < 1.e-12)
      new_dw->put(delt_vartype(DBL_MAX), lb->delTLabel);
    else
      new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(strainEnergy), lb->StrainEnergyLabel);
    }
  }
}

void
CompMooneyRivlin::addComputesAndRequires(Task*, const MPMMaterial*,
                                         const PatchSet*, const bool,
                                         const bool) const
{
}

void
CompMooneyRivlin::carryForward(const PatchSubset* patches,
                               const MPMMaterial* matl, DataWarehouse* old_dw,
                               DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    int dwi = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

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

///////////////////////////////////////////////////////////////////////////
/*! Allocate data required during the conversion of failed particles
    from one material to another */
///////////////////////////////////////////////////////////////////////////
void
CompMooneyRivlin::allocateCMDataAddRequires(Task* task, const MPMMaterial* matl,
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
CompMooneyRivlin::allocateCMDataAdd(
  DataWarehouse* new_dw, ParticleSubset* addset,
  ParticleLabelVariableMap* newState, ParticleSubset* delset,
  DataWarehouse*)
{
  // Copy the data common to all constitutive models from the particle to be
  // deleted to the particle to be added.
  // This method is defined in the ConstitutiveModel base class.
  copyDelToAddSetForConvertExplicit(new_dw, delset, addset, newState);
}

double
CompMooneyRivlin::computeRhoMicroCM(double /*pressure*/, const double /*p_ref*/,
                                    const MPMMaterial* /*matl*/,
                                    double temperature, double rho_guess)
{
#if 0
  double rho_orig = matl->getInitialDensity();
  double bulk = d_initialData.Bulk;

  double p_gauge = pressure - p_ref;
  double rho_cur;

  rho_cur = rho_orig*(p_gauge/bulk + sqrt((p_gauge/bulk)*(p_gauge/bulk) +1));
#endif

  std::cout << "NO VERSION OF computeRhoMicroCM EXISTS YET FOR CompMooneyRivlin"
       << endl;

  double rho_cur = 0.;

  return rho_cur;
}

void
CompMooneyRivlin::computePressEOSCM(double /*rho_cur*/, double& /*pressure*/,
                                    double /*p_ref*/, double& /*dp_drho*/,
                                    double& /*tmp*/,
                                    const MPMMaterial* /*matl*/,
                                    double temperature)
{
#if 0
  double bulk = d_initialData.Bulk;
  double shear = d_initialData.Shear;
  double rho_orig = matl->getInitialDensity();

  double p_g = .5*bulk*(rho_cur/rho_orig - rho_orig/rho_cur);
  pressure = p_ref + p_g;
  dp_drho  = .5*bulk*(rho_orig/(rho_cur*rho_cur) + 1./rho_orig);
  tmp = (bulk + 4.*shear/3.)/rho_cur;  // speed of sound squared
#endif

  std::cout << "NO VERSION OF computePressEOSCM EXISTS YET FOR CompMooneyRivlin"
       << endl;
}

double
CompMooneyRivlin::getCompressibility()
{
  std::cout << "NO VERSION OF computePressEOSCM EXISTS YET FOR CompMooneyRivlin"
       << endl;
  return 1.0;
}
