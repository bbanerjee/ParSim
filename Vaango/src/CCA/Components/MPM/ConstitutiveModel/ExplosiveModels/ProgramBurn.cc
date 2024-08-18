/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#include <CCA/Components/MPM/ConstitutiveModel/Utilities/Constants.h>
#include <CCA/Components/MPM/ConstitutiveModel/ExplosiveModels/ProgramBurn.h>
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
#include<CCA/Components/MPM/Core/MPMLabel.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <fstream>
#include <iostream>

using namespace Uintah;

ProgramBurn::ProgramBurn(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  d_useModifiedEOS = false;

  // These two parameters are used for the unburned Murnahan EOS
  ps->require("K", d_cm.d_K);
  ps->require("n", d_cm.d_n);

  // These parameters are used for the product JWL EOS
  ps->require("A", d_cm.d_A);
  ps->require("B", d_cm.d_B);
  ps->require("C", d_cm.d_C);
  ps->require("R1", d_cm.d_R1);
  ps->require("R2", d_cm.d_R2);
  ps->require("om", d_cm.d_om);
  ps->require("rho0", d_cm.d_rho0);

  // These parameters are needed for the reaction model
  ps->require("starting_location", d_cm.d_start_place);
  ps->require("D", d_cm.d_D); // Detonation velocity
  ps->getWithDefault(
    "direction_if_plane", d_cm.d_direction, Vector(0., 0., 0.));
  ps->getWithDefault("T0", d_cm.d_T0, 0.0);

  pProgressFLabel = VarLabel::create(
    "p.progressF", ParticleVariable<double>::getTypeDescription());
  pProgressFLabel_preReloc = VarLabel::create(
    "p.progressF+", ParticleVariable<double>::getTypeDescription());
  pLocalizedLabel = VarLabel::create(
    "p.localized", ParticleVariable<int>::getTypeDescription());
  pLocalizedLabel_preReloc = VarLabel::create(
    "p.localized+", ParticleVariable<int>::getTypeDescription());
}

ProgramBurn::ProgramBurn(const ProgramBurn* cm)
  : ConstitutiveModel(cm)
{
  d_useModifiedEOS = cm->d_useModifiedEOS;

  d_cm.d_K = cm->d_cm.d_K;
  d_cm.d_n = cm->d_cm.d_n;

  d_cm.d_A    = cm->d_cm.d_A;
  d_cm.d_B    = cm->d_cm.d_B;
  d_cm.d_C    = cm->d_cm.d_C;
  d_cm.d_R1   = cm->d_cm.d_R1;
  d_cm.d_R2   = cm->d_cm.d_R2;
  d_cm.d_om   = cm->d_cm.d_om;
  d_cm.d_rho0 = cm->d_cm.d_rho0;

  d_cm.d_start_place = cm->d_cm.d_start_place;
  d_cm.d_direction   = cm->d_cm.d_direction;
  d_cm.d_D           = cm->d_cm.d_D;

  pProgressFLabel = VarLabel::create(
    "p.progressF", ParticleVariable<double>::getTypeDescription());
  pProgressFLabel_preReloc = VarLabel::create(
    "p.progressF+", ParticleVariable<double>::getTypeDescription());
  pLocalizedLabel = VarLabel::create(
    "p.localized", ParticleVariable<int>::getTypeDescription());
  pLocalizedLabel_preReloc = VarLabel::create(
    "p.localized+", ParticleVariable<int>::getTypeDescription());
}

ProgramBurn::~ProgramBurn()
{
  VarLabel::destroy(pProgressFLabel);
  VarLabel::destroy(pProgressFLabel_preReloc);
  VarLabel::destroy(pLocalizedLabel);
  VarLabel::destroy(pLocalizedLabel_preReloc);
}

std::unique_ptr<ConstitutiveModel>
ProgramBurn::clone()
{
  return std::make_unique<ProgramBurn>(*this);
}

void
ProgramBurn::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "program_burn");
  }

  cm_ps->appendElement("K", d_cm.d_K);
  cm_ps->appendElement("n", d_cm.d_n);

  cm_ps->appendElement("A", d_cm.d_A);
  cm_ps->appendElement("B", d_cm.d_B);
  cm_ps->appendElement("C", d_cm.d_C);
  cm_ps->appendElement("R1", d_cm.d_R1);
  cm_ps->appendElement("R2", d_cm.d_R2);
  cm_ps->appendElement("om", d_cm.d_om);
  cm_ps->appendElement("rho0", d_cm.d_rho0);

  cm_ps->appendElement("starting_location", d_cm.d_start_place);
  cm_ps->appendElement("direction_if_plane", d_cm.d_direction);
  cm_ps->appendElement("D", d_cm.d_D);
  cm_ps->appendElement("T0", d_cm.d_T0);
}

void
ProgramBurn::addParticleState(std::vector<const VarLabel*>& from,
                              std::vector<const VarLabel*>& to)
{
  // Add the local particle state data for this constitutive model.
  from.push_back(pProgressFLabel);
  to.push_back(pProgressFLabel_preReloc);

  from.push_back(pLocalizedLabel);
  to.push_back(pLocalizedLabel_preReloc);
}

void
ProgramBurn::addInitialComputesAndRequires(Task* task,
                                           const MPMMaterial* matl,
                                           const PatchSet*) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pProgressFLabel, matlset);
  task->computes(pLocalizedLabel, matlset);
}

void
ProgramBurn::initializeCMData(const Patch* patch,
                              const MPMMaterial* matl,
                              DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  ParticleVariable<double> pProgress;
  ParticleVariable<int> pLocalized;
  new_dw->allocateAndPut(pProgress, pProgressFLabel, pset);
  new_dw->allocateAndPut(pLocalized, pLocalizedLabel, pset);

  for (int& iter : *pset) {
    pProgress[iter]  = 0.;
    pLocalized[iter] = 0;
  }

  computeStableTimestep(patch, matl, new_dw);
}

void
ProgramBurn::computeStableTimestep(const Patch* patch,
                                   const MPMMaterial* matl,
                                   DataWarehouse* new_dw)
{
  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of computeStressTensor
  Vector dx = patch->dCell();
  int matID = matl->getDWIndex();

  ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);
  constParticleVariable<double> pMass, pVolume, ptemperature;
  constParticleVariable<Vector> pVelocity;

  new_dw->get(pMass, lb->pMassLabel, pset);
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pVelocity, lb->pVelocityLabel, pset);
  new_dw->get(ptemperature, lb->pTemperatureLabel, pset);

  double K    = d_cm.d_K;
  double n    = d_cm.d_n;
  double rho0 = d_cm.d_rho0;

  Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);
  for (int idx : *pset) {
    // Compute wave speed at each particle, store the maximum
    double rhoM    = pMass[idx] / pVolume[idx];
    double dp_drho = (1. / (K * rho0)) * pow((rhoM / rho0), n - 1.);
    double c_dil   = std::sqrt(dp_drho);
    Vector velMax  = pVelocity[idx].cwiseAbs() + c_dil;
    waveSpeed      = Max(velMax, waveSpeed);
  }
  waveSpeed       = dx / waveSpeed;
  double delT_new = waveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

void
ProgramBurn::addComputesAndRequires(Task* task,
                                    const MPMMaterial* matl,
                                    const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForExplicit(task, matlset, patches);

  task->requires(Task::OldDW, lb->simulationTimeLabel);
  task->requires(Task::OldDW, lb->pParticleIDLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pProgressFLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pLocalizedLabel, matlset, Ghost::None);
  task->computes(pProgressFLabel_preReloc, matlset);
  task->computes(pLocalizedLabel_preReloc, matlset);
}

void
ProgramBurn::computeStressTensor(const PatchSubset* patches,
                                 const MPMMaterial* matl,
                                 DataWarehouse* old_dw,
                                 DataWarehouse* new_dw)
{
  int matID = matl->getDWIndex();

  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  simTime_vartype simTimeVar;
  old_dw->get(simTimeVar, lb->simulationTimeLabel);

  double time = simTimeVar - d_cm.d_T0;

  double K    = d_cm.d_K;
  double n    = d_cm.d_n;
  double A    = d_cm.d_A;
  double B    = d_cm.d_B;
  double C    = d_cm.d_C;
  double R1   = d_cm.d_R1;
  double R2   = d_cm.d_R2;
  double om   = d_cm.d_om;
  double rho0 = d_cm.d_rho0; // matl->getInitialDensity();

  double A_d = d_cm.d_direction.x();
  double B_d = d_cm.d_direction.y();
  double C_d = d_cm.d_direction.z();

  double x0 = d_cm.d_start_place.x();
  double y0 = d_cm.d_start_place.y();
  double z0 = d_cm.d_start_place.z();

  double D_d = -A_d * x0 - B_d * y0 - C_d * z0;

  double denom = 1.0;
  double plane = 0.;

  if (d_cm.d_direction.length() > 0.0) {
    plane = 1.0;
    denom = std::sqrt(A_d * A_d + B_d * B_d + C_d * C_d);
  }

  for (int pp = 0; pp < patches->size(); pp++) {
    const Patch* patch   = patches->get(pp);
    Vector dx            = patch->dCell();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    constParticleVariable<int> pLocalized;
    constParticleVariable<long64> pParticleID;
    constParticleVariable<double> pMass;
    constParticleVariable<Point> pX;
    constParticleVariable<Vector> pVelocity;
    constParticleVariable<Matrix3> pVelGrad, pDefGrad_new;
    old_dw->get(pLocalized, pLocalizedLabel, pset);
    old_dw->get(pParticleID, lb->pParticleIDLabel, pset);
    old_dw->get(pMass, lb->pMassLabel, pset);
    old_dw->get(pX, lb->pXLabel, pset);
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);
    new_dw->get(pVelGrad, lb->pVelGradLabel_preReloc, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    ParticleVariable<int> pLocalized_new;
    ParticleVariable<double> pdTdt, p_q, pProgressF_new;
    ParticleVariable<Matrix3> pStress;
    new_dw->allocateAndPut(pLocalized_new, pLocalizedLabel_preReloc, pset);
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);
    new_dw->allocateAndPut(pProgressF_new, pProgressFLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress, lb->pStressLabel_preReloc, pset);

    if (!flag->d_doGridReset) {
      std::cerr
        << "The program_burn model doesn't work without resetting the grid"
        << "\n";
    }

    double strainEnergy = 0.;
    Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);
    for (int idx : *pset) {
      pLocalized_new[idx] = pLocalized[idx];

      Point pt = pX[idx];

      double dist_plane =
        fabs(A_d * pt.x() + B_d * pt.y() + C_d * pt.z() + D_d) / denom;

      double dist_straight = (pt - d_cm.d_start_place).length();

      double dist = dist_plane * plane + dist_straight * (1. - plane);

      double t_b = dist / d_cm.d_D;

      double delta_L = 1.5 * std::pow(pMass[idx] / rho0, 1. / 3.) / d_cm.d_D;

      if (time >= t_b) {
        pProgressF_new[idx] = (time - t_b) / delta_L;
        if (pProgressF_new[idx] > 0.96) {
          pProgressF_new[idx] = 1.0;
        }
      } else {
        pProgressF_new[idx] = 0.0;
      }

      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;

      double J = pDefGrad_new[idx].Determinant();

      if (J <= 0.0) {
        std::cout << "negative J in ProgramBurn, J =" << J << "\n";
        std::cout << "pos = " << pX[idx] << "\n";
        pLocalized_new[idx] = -999;
        std::cout << "localizing (deleting) particle " << pParticleID[idx]
                  << "\n";
        std::cout << "material = " << matID << "\n"
                  << "Momentum deleted = " << pVelocity[idx] * pMass[idx]
                  << "\n";
        J = 1;
      }

      //  The following computes a pressure for partially burned particles
      //  as a mixture of Murnahan and JWL pressures, based on pProgressF
      //  This is as described in Eq. 5 of "JWL++: ..." by Souers, et al.
      double pM   = (1. / (n * K)) * (std::pow(J, -n) - 1.);
      double pJWL = pM;

      // For computing speed of sound if not yet detonating
      double rho_cur = rho0 / J;
      double dp_drho = (1. / (K * rho0)) * std::pow((rho_cur / rho0), n - 1.);
      if (pProgressF_new[idx] > 0.0) {
        double one_plus_omega               = 1. + om;
        double inv_rho_rat                  = J;      // rho0/rhoM;
        double rho_rat                      = 1. / J; // rhoM/rho0;
        double A_e_to_the_R1_rho0_over_rhoM = A * exp(-R1 * inv_rho_rat);
        double B_e_to_the_R2_rho0_over_rhoM = B * exp(-R2 * inv_rho_rat);
        double C_rho_rat_tothe_one_plus_omega =
          C * pow(rho_rat, one_plus_omega);

        pJWL = A_e_to_the_R1_rho0_over_rhoM + B_e_to_the_R2_rho0_over_rhoM +
               C_rho_rat_tothe_one_plus_omega;

        // For computing speed of sound if detonat(ing/ed)
        double rho0_rhoMsqrd = rho0 / (rho_cur * rho_cur);
        dp_drho = R1 * rho0_rhoMsqrd * A_e_to_the_R1_rho0_over_rhoM +
                  R2 * rho0_rhoMsqrd * B_e_to_the_R2_rho0_over_rhoM +
                  (one_plus_omega / rho_cur) * C_rho_rat_tothe_one_plus_omega;
      }

      double p = pM * (1.0 - pProgressF_new[idx]) + pJWL * pProgressF_new[idx];

      // compute the total stress
      pStress[idx] = Vaango::Util::Identity * (-p);

      // Compute wave speed at each particle, store the maximum
      double c_dil  = std::sqrt(dp_drho);
      Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
      waveSpeed     = Max(velMax, waveSpeed);

      // Compute artificial viscosity term
      if (flag->d_artificialViscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;
        double c_bulk = sqrt(1. / (K * rho_cur));
        Matrix3 D     = (pVelGrad[idx] + pVelGrad[idx].Transpose()) * 0.5;
        p_q[idx] = artificialBulkViscosity(D.Trace(), c_bulk, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.;
      }
    } // end loop over particles

    waveSpeed       = dx / waveSpeed;
    double delT_new = waveSpeed.minComponent();
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(strainEnergy), lb->StrainEnergyLabel);
    }
  }
}

void
ProgramBurn::addComputesAndRequires(Task*,
                                    const MPMMaterial*,
                                    const PatchSet*,
                                    const bool,
                                    const bool) const
{
}

void
ProgramBurn::addRequiresDamageParameter(Task* task,
                                        const MPMMaterial* matl,
                                        const PatchSet*) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, pLocalizedLabel_preReloc, matlset, Ghost::None);
}

void
ProgramBurn::getDamageParameter(const Patch* patch,
                                ParticleVariable<int>& damage,
                                int matID,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw)
{
  ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);
  constParticleVariable<int> pLocalized;
  new_dw->get(pLocalized, pLocalizedLabel_preReloc, pset);

  // Only update the damage variable if it hasn't been modified by a damage
  // model earlier
  for (auto particle : *pset) {
    if (damage[particle] == 0) {
      damage[particle] = pLocalized[particle];
    }
  }
}

void
ProgramBurn::allocateCMDataAddRequires(Task* task,
                                       const MPMMaterial* matl,
                                       const PatchSet* patches,
                                       MPMLabel*) const
{
  const MaterialSubset* matlset = matl->thisMaterial();

  // Allocate the variables shared by all constitutive models
  // for the particle convert operation
  // This method is defined in the ConstitutiveModel base class.
  addSharedRForConvertExplicit(task, matlset, patches);
  task->requires(Task::NewDW, pProgressFLabel_preReloc, matlset, Ghost::None);
  task->requires(Task::NewDW, pLocalizedLabel_preReloc, matlset, Ghost::None);
}

void
ProgramBurn::allocateCMDataAdd(DataWarehouse* new_dw,
                               ParticleSubset* addset,
                               ParticleLabelVariableMap* newState,
                               ParticleSubset* delset,
                               DataWarehouse*)
{
  // Copy the data common to all constitutive models from the particle to be
  // deleted to the particle to be added.
  // This method is defined in the ConstitutiveModel base class.
  copyDelToAddSetForConvertExplicit(new_dw, delset, addset, newState);

  ParticleVariable<int> pLocalized;
  constParticleVariable<int> o_Localized;
  new_dw->allocateTemporary(pLocalized, addset);
  new_dw->get(o_Localized, pLocalizedLabel_preReloc, delset);

  ParticleVariable<int> pProgressF;
  constParticleVariable<int> o_ProgressF;
  new_dw->allocateTemporary(pProgressF, addset);
  new_dw->get(o_ProgressF, pProgressFLabel_preReloc, delset);

  // Copy the data local to this constitutive model from the particles to
  // be deleted to the particles to be added
  ParticleSubset::iterator o, n = addset->begin();
  for (o = delset->begin(); o != delset->end(); o++, n++) {
    pLocalized[*n] = o_Localized[*o];
    pProgressF[*n] = o_ProgressF[*o];
  }
  (*newState)[pLocalizedLabel] = pLocalized.clone();
  (*newState)[pProgressFLabel] = pProgressF.clone();
}

void
ProgramBurn::carryForward(const PatchSubset* patches,
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

// This is not yet implemented - JG- 7/26/10
double
ProgramBurn::computeRhoMicroCM(double pressure,
                               const double p_ref,
                               const MPMMaterial* matl,
                               double temperature,
                               double rho_guess)
{
  std::cout << "NO VERSION OF computeRhoMicroCM EXISTS YET FOR ProgramBurn"
            << "\n";
  double rho_orig = d_cm.d_rho0; // matl->getInitialDensity();

  return rho_orig;
}

void
ProgramBurn::computePressEOSCM(const double rhoM,
                               double& pressure,
                               const double p_ref,
                               double& dp_drho,
                               double& tmp,
                               const MPMMaterial* matl,
                               double temperature)
{
  double A    = d_cm.d_A;
  double B    = d_cm.d_B;
  double R1   = d_cm.d_R1;
  double R2   = d_cm.d_R2;
  double om   = d_cm.d_om;
  double rho0 = d_cm.d_rho0;
  double cv   = matl->getSpecificHeat();
  double V    = rho0 / rhoM;
  double P1   = A * exp(-R1 * V);
  double P2   = B * exp(-R2 * V);
  double P3   = om * cv * tmp * rhoM;

  pressure = P1 + P2 + P3;

  dp_drho = (R1 * rho0 * P1 + R2 * rho0 * P2) / (rhoM * rhoM) + om * cv * tmp;
}

// This is not yet implemented - JG- 7/26/10
double
ProgramBurn::getCompressibility()
{
  std::cout << "NO VERSION OF getCompressibility EXISTS YET FOR ProgramBurn"
            << "\n";
  return 1.0;
}
