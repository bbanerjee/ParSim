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

#include <CCA/Components/MPM/ConstitutiveModel/ExplosiveModels/JWLppMPM.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/Constants.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Grid/Variables/VarTypes.h>
#include<CCA/Components/MPM/Core/MPMLabel.h>

#define CHECK_ISFINITE

using namespace Uintah;

JWLppMPM::JWLppMPM(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  d_useModifiedEOS = false;

  // Read the ignition pressure
  ps->require("ignition_pressure", d_cm.ignition_pressure);

  // Shear viscosity
  ps->require("viscosity", d_cm.mu);

  // These two parameters are used for the unburned Murnahan EOS
  ps->require("murnaghan_K", d_cm.K);
  ps->require("murnaghan_n", d_cm.n);

  // These parameters are used for the product JWL EOS
  ps->require("jwl_A", d_cm.A);
  ps->require("jwl_B", d_cm.B);
  ps->require("jwl_C", d_cm.C);
  ps->require("jwl_R1", d_cm.R1);
  ps->require("jwl_R2", d_cm.R2);
  ps->require("jwl_om", d_cm.omega);
  // ps->require("jwl_rho0", d_cm.rho0);  // Get from matl->getInitialDensity()

  // These parameters are needed for the reaction model
  ps->require("reaction_G", d_cm.G); // Rate coefficient
  ps->require("reaction_b", d_cm.b); // Pressure exponent
  // Maximum time increment for burn model subcycling
  ps->getWithDefault("max_burn_timestep_size", d_cm.max_burn_timestep, 1.0e-12);
  // Limit on the fraction that remains unburned
  ps->getWithDefault("max_burned_fraction", d_cm.max_burned_frac, 1.0);

  // Initial stress
  // Fix: Need to make it more general.  Add gravity turn-on option and
  //      read from file option etc.
  ps->getWithDefault("useInitialStress", d_useInitialStress, false);
  d_init_pressure = 0.0;
  if (d_useInitialStress) {
    ps->getWithDefault("initial_pressure", d_init_pressure, 0.0);
  }

  // Use Newton iterations for stress update by default.  Else use two step
  // algorithm.
  ps->getWithDefault("doFastStressCompute", d_fastCompute, false);
  ps->getWithDefault(
    "tolerance_for_Newton_iterations", d_newtonIterTol, 1.0e-3);
  ps->getWithDefault("max_number_of_Newton_iterations", d_newtonIterMax, 20);

  pProgressFLabel = VarLabel::create(
    "p.progressF", ParticleVariable<double>::getTypeDescription());
  pProgressFLabel_preReloc = VarLabel::create(
    "p.progressF+", ParticleVariable<double>::getTypeDescription());
  pProgressdelFLabel = VarLabel::create(
    "p.progressdelF", ParticleVariable<double>::getTypeDescription());
  pProgressdelFLabel_preReloc = VarLabel::create(
    "p.progressdelF+", ParticleVariable<double>::getTypeDescription());
  pLocalizedLabel = VarLabel::create(
    "p.localized", ParticleVariable<int>::getTypeDescription());
  pLocalizedLabel_preReloc = VarLabel::create(
    "p.localized+", ParticleVariable<int>::getTypeDescription());
}

JWLppMPM::JWLppMPM(const JWLppMPM* cm)
  : ConstitutiveModel(cm)
{
  d_useModifiedEOS = cm->d_useModifiedEOS;

  d_cm.ignition_pressure = cm->d_cm.ignition_pressure;

  d_cm.mu = cm->d_cm.mu;

  d_cm.K = cm->d_cm.K;
  d_cm.n = cm->d_cm.n;

  d_cm.A     = cm->d_cm.A;
  d_cm.B     = cm->d_cm.B;
  d_cm.C     = cm->d_cm.C;
  d_cm.R1    = cm->d_cm.R1;
  d_cm.R2    = cm->d_cm.R2;
  d_cm.omega = cm->d_cm.omega;
  // d_cm.rho0 = cm->d_cm.rho0;

  d_cm.G                 = cm->d_cm.G;
  d_cm.b                 = cm->d_cm.b;
  d_cm.max_burn_timestep = cm->d_cm.max_burn_timestep;
  d_cm.max_burned_frac   = cm->d_cm.max_burned_frac;

  // Initial stress
  d_useInitialStress = cm->d_useInitialStress;
  d_init_pressure    = cm->d_init_pressure;

  // Stress compute algorithms
  d_fastCompute   = cm->d_fastCompute;
  d_newtonIterTol = cm->d_newtonIterTol;
  d_newtonIterMax = cm->d_newtonIterMax;

  pProgressFLabel = VarLabel::create(
    "p.progressF", ParticleVariable<double>::getTypeDescription());
  pProgressFLabel_preReloc = VarLabel::create(
    "p.progressF+", ParticleVariable<double>::getTypeDescription());
  pProgressdelFLabel = VarLabel::create(
    "p.progressdelF", ParticleVariable<double>::getTypeDescription());
  pProgressdelFLabel_preReloc = VarLabel::create(
    "p.progressdelF+", ParticleVariable<double>::getTypeDescription());
  pLocalizedLabel = VarLabel::create(
    "p.localized", ParticleVariable<int>::getTypeDescription());
  pLocalizedLabel_preReloc = VarLabel::create(
    "p.localized+", ParticleVariable<int>::getTypeDescription());
}

JWLppMPM::~JWLppMPM()
{
  VarLabel::destroy(pProgressFLabel);
  VarLabel::destroy(pProgressFLabel_preReloc);
  VarLabel::destroy(pProgressdelFLabel);
  VarLabel::destroy(pProgressdelFLabel_preReloc);
  VarLabel::destroy(pLocalizedLabel);
  VarLabel::destroy(pLocalizedLabel_preReloc);
}

std::unique_ptr<ConstitutiveModel>
JWLppMPM::clone()
{
  return std::make_unique<JWLppMPM>(*this);
}

void
JWLppMPM::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "jwlpp_mpm");
  }

  cm_ps->appendElement("ignition_pressure", d_cm.ignition_pressure);

  cm_ps->appendElement("viscosity", d_cm.mu);

  cm_ps->appendElement("murnaghan_K", d_cm.K);
  cm_ps->appendElement("murnaghan_n", d_cm.n);

  cm_ps->appendElement("jwl_A", d_cm.A);
  cm_ps->appendElement("jwl_B", d_cm.B);
  cm_ps->appendElement("jwl_C", d_cm.C);
  cm_ps->appendElement("jwl_R1", d_cm.R1);
  cm_ps->appendElement("jwl_R2", d_cm.R2);
  cm_ps->appendElement("jwl_om", d_cm.omega);
  // cm_ps->appendElement("jwl_rho0", d_cm.rho0);

  cm_ps->appendElement("reaction_b", d_cm.b);
  cm_ps->appendElement("reaction_G", d_cm.G);
  cm_ps->appendElement("max_burn_timestep_size", d_cm.max_burn_timestep);
  cm_ps->appendElement("max_burned_fraction", d_cm.max_burned_frac);

  cm_ps->appendElement("useInitialStress", d_useInitialStress);
  if (d_useInitialStress) {
    cm_ps->appendElement("initial_pressure", d_init_pressure);
  }

  cm_ps->appendElement("doFastStressCompute", d_fastCompute);
  cm_ps->appendElement("tolerance_for_Newton_iterations", d_newtonIterTol);
  cm_ps->appendElement("max_number_of_Newton_iterations", d_newtonIterMax);
}

void
JWLppMPM::addParticleState(std::vector<const VarLabel*>& from,
                           std::vector<const VarLabel*>& to)
{
  // Add the local particle state data for this constitutive model.
  from.push_back(pProgressFLabel);
  to.push_back(pProgressFLabel_preReloc);
  from.push_back(pProgressdelFLabel);
  to.push_back(pProgressdelFLabel_preReloc);
  from.push_back(pLocalizedLabel);
  to.push_back(pLocalizedLabel_preReloc);
}

void
JWLppMPM::addInitialComputesAndRequires(Task* task,
                                        const MPMMaterial* matl,
                                        const PatchSet*) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pProgressFLabel, matlset);
  task->computes(pProgressdelFLabel, matlset);
  task->computes(pLocalizedLabel, matlset);
}

void
JWLppMPM::initializeCMData(const Patch* patch,
                           const MPMMaterial* matl,
                           DataWarehouse* new_dw)
{
  // Initialize local variables
  Matrix3 zero(0.0);
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  ParticleVariable<int> pLocalized;
  ParticleVariable<double> pProgress, pProgressdelF;
  new_dw->allocateAndPut(pProgress, pProgressFLabel, pset);
  new_dw->allocateAndPut(pProgressdelF, pProgressdelFLabel, pset);
  new_dw->allocateAndPut(pLocalized, pLocalizedLabel, pset);
  for (auto pidx : *pset) {
    pProgress[pidx]     = 0.0;
    pProgressdelF[pidx] = 0.0;
    pLocalized[pidx]    = 0;
  }

  // Initialize the variables shared by all constitutive models
  if (!d_useInitialStress) {
    // This method is defined in the ConstitutiveModel base class.
    initSharedDataForExplicit(patch, matl, new_dw);

  } else {
    ParticleVariable<double> pdTdt;
    ParticleVariable<Matrix3> pDefGrad;
    ParticleVariable<Matrix3> pStress;

    new_dw->getModifiable(pDefGrad, lb->pDefGradLabel, pset);

    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel, pset);
    new_dw->allocateAndPut(pStress, lb->pStressLabel, pset);

    // Set the initial pressure
    double p        = d_init_pressure;
    Matrix3 sigInit = Vaango::Util::Identity * (-p);

    // Compute deformation gradient
    //  using the Murnaghan eos
    //     p = (1/nK) [J^(-n) - 1]
    //     =>
    //     det(F) = (1 + nKp)^(-1/n)
    //     =>
    //     F_{11} = F_{22} = F_{33} = (1 + nKp)^(-1/3n)
    double F11 = std::pow((1.0 + d_cm.K * d_cm.n * p), (-1.0 / (3.0 * d_cm.n)));
    Matrix3 defGrad = Vaango::Util::Identity * F11;

    for (auto pidx : *pset) {
      pdTdt[pidx]    = 0.0;
      pStress[pidx]  = sigInit;
      pDefGrad[pidx] = defGrad;
    }
  }

  computeStableTimestep(patch, matl, new_dw);
}

void
JWLppMPM::computeStableTimestep(const Patch* patch,
                                const MPMMaterial* matl,
                                DataWarehouse* new_dw)
{
  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of computeStressTensor
  Vector dx = patch->dCell();
  int matID = matl->getDWIndex();
  // Retrieve the array of constitutive parameters
  ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);
  constParticleVariable<double> pMass, pVolume, ptemperature;
  constParticleVariable<Vector> pVelocity;

  new_dw->get(pMass, lb->pMassLabel, pset);
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pVelocity, lb->pVelocityLabel, pset);
  new_dw->get(ptemperature, lb->pTemperatureLabel, pset);

  double K = d_cm.K;
  double n = d_cm.n;
  // double rho0 = d_cm.rho0;
  double rho0 = matl->getInitialDensity();

  Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);
  for (int idx : *pset) {
    // Compute wave speed at each particle, store the maximum
    double rhoM    = pMass[idx] / pVolume[idx];
    double dp_drho = (1. / (K * rho0)) * std::pow((rhoM / rho0), n - 1.);
    double c_dil   = std::sqrt(dp_drho);
    Vector velMax  = pVelocity[idx].cwiseAbs() + c_dil;
    waveSpeed      = Max(velMax, waveSpeed);
  }
  waveSpeed       = dx / waveSpeed;
  double delT_new = waveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

void
JWLppMPM::addComputesAndRequires(Task* task,
                                 const MPMMaterial* matl,
                                 const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForHypoExplicit(task, matlset, patches);

  task->needs(Task::OldDW, pProgressFLabel, matlset, Ghost::None);
  task->needs(Task::OldDW, pProgressdelFLabel, matlset, Ghost::None);
  task->needs(Task::OldDW, pLocalizedLabel, matlset, Ghost::None);

  task->computes(pProgressFLabel_preReloc, matlset);
  task->computes(pProgressdelFLabel_preReloc, matlset);
  task->computes(pLocalizedLabel_preReloc, matlset);
}

void
JWLppMPM::computeStressTensor(const PatchSubset* patches,
                              const MPMMaterial* matl,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw)
{

  // Material parameters
  // double d_rho0 = d_cm.rho0;
  double d_rho0 = matl->getInitialDensity();

  int matID = matl->getDWIndex();

  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  // Loop through patches
  for (int pp = 0; pp < patches->size(); pp++) {
    const Patch* patch   = patches->get(pp);
    Vector dx            = patch->dCell();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    // variables to hold this timestep's values
    constParticleVariable<int> pLocalized_old;
    constParticleVariable<double> pMass, pProgressF_old, pProgressdelF_old;
    constParticleVariable<Vector> pVelocity;
    constParticleVariable<Matrix3> pVelGrad_new, pDefGrad_old, pStress_old;
    old_dw->get(pLocalized_old, pLocalizedLabel, pset);
    old_dw->get(pMass, lb->pMassLabel, pset);
    old_dw->get(pProgressF_old, pProgressFLabel, pset);
    old_dw->get(pProgressdelF_old, pProgressdelFLabel, pset);
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);
    old_dw->get(pDefGrad_old, lb->pDefGradLabel, pset);
    old_dw->get(pStress_old, lb->pStressLabel, pset);
    new_dw->get(pVelGrad_new, lb->pVelGradLabel_preReloc, pset);

    ParticleVariable<double> pVolume_new;
    ParticleVariable<Matrix3> pDefGrad_new;
    new_dw->getModifiable(pVolume_new, lb->pVolumeLabel_preReloc, pset);
    new_dw->getModifiable(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    ParticleVariable<int> pLocalized_new;
    ParticleVariable<double> pdTdt, p_q, pProgressF_new, pProgressdelF_new;
    ParticleVariable<Matrix3> pStress_new;
    new_dw->allocateAndPut(pLocalized_new, pLocalizedLabel_preReloc, pset);
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);
    new_dw->allocateAndPut(pProgressF_new, pProgressFLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pProgressdelF_new, pProgressdelFLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);

    if (!flag->d_doGridReset) {
      std::cerr << "The jwlpp_mpm model doesn't work without resetting the grid"
                << "\n";
    }

    // Compute deformation gradient and velocity gradient at each
    // particle before pressure stabilization
    double strainEnergy = 0.0;
    Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);
    for (int idx : *pset) {
      // If the particle has already failed just ignore
      pLocalized_new[idx] = 0;
      if (pLocalized_old[idx]) {
        pStress_new[idx] = pStress_old[idx];
        pdTdt[idx]       = 0.0;
        p_q[idx]         = 0.0;

        pProgressF_new[idx]    = pProgressF_old[idx];
        pProgressdelF_new[idx] = pProgressdelF_old[idx];
        pLocalized_new[idx]    = pLocalized_old[idx];
        continue;
      }

      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;

      double J = pDefGrad_new[idx].Determinant();

      // If the particle has already failed just ignore
      if (pLocalized_old[idx])
        continue;

      if (!(J > 0.0)) {
        std::cerr
          << "**ERROR in JWL++MPM** Negative Jacobian of deformation gradient"
          << "\n";
        std::cerr << "idx = " << idx << " J = " << J << " matl = " << matl
                  << "\n";
        std::cerr << "F_new = " << pDefGrad_new[idx] << "\n";
        std::cerr << "VelGrad = " << pVelGrad_new[idx] << "\n";
        std::cerr << "**Particle is being removed from the computation**"
                  << "\n";
        // throw InvalidValue("**ERROR**: Error in deformation gradient",
        // __FILE__, __LINE__);

        pStress_new[idx]  = Vaango::Util::Zero;
        pVolume_new[idx]  = d_rho0 / pMass[idx];
        pdTdt[idx]        = 0.0;
        p_q[idx]          = 0.0;
        pDefGrad_new[idx] = Vaango::Util::Identity;

        pProgressF_new[idx]    = pProgressF_old[idx];
        pProgressdelF_new[idx] = pProgressdelF_old[idx];
        pLocalized_new[idx]    = 1;
        continue;
      }

      // Compute new mass density and update the deformed volume
      double rho_cur = d_rho0 / J;

      // Update the burn fraction and pressure
      double J_old = pDefGrad_old[idx].Determinant();
      double p_old = -(1.0 / 3.0) * pStress_old[idx].Trace();
      double f_old = pProgressF_old[idx];
      double f_new = f_old;
      double p_new = p_old;
      computeUpdatedFractionAndPressure(
        J_old, J, f_old, p_old, delT, f_new, p_new);

      // Compute a viscous stress
      Matrix3 D = (pVelGrad_new[idx] + pVelGrad_new[idx].Transpose()) * 0.5;
      double trD = D.Trace();
      Matrix3 devD = D - Vaango::Util::Identity * (trD * Vaango::Util::one_third);

      // Update the volume fraction and the stress in the data warehouse
      pProgressdelF_new[idx] = f_new - f_old;
      pProgressF_new[idx]    = f_new;
      pStress_new[idx]       = Vaango::Util::Identity * (-p_new) + devD * (2.0 * d_cm.mu) ;

      // if (std::isnan(pStress_new[idx].Norm()) || pStress_new[idx].Norm() >
      // 1.0e20) {
      //  std::cerr << "particle = " << idx << " velGrad = " <<
      //  pVelGrad_new[idx] <<
      //  " stress_old = " << pStress_old[idx] << "\n";
      //  std::cerr << " stress = " << pStress_new[idx]
      //       << "  pProgressdelF_new = " << pProgressdelF_new[idx]
      //       << "  pProgressF_new = " << pProgressF_new[idx]
      //       << " pm = " << pM << " pJWL = " << pJWL <<  " rho_cur = " <<
      //       rho_cur << "\n";
      //  std::cerr << " pMass = " << pMass[idx] << " pvol = " <<
      //  pVolume_new[idx] <<
      //  "\n";
      //  throw InvalidValue("**JWLppMPM ERROR**: Nan in stress value",
      //  __FILE__, __LINE__);
      //}

      // Compute wave speed at each particle, store the maximum
      double dp_drho =
        (1. / (d_cm.K * d_rho0)) * std::pow((rho_cur / d_rho0), d_cm.n - 1.);
      double c_dil  = std::sqrt(dp_drho);
      Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
      waveSpeed     = Max(velMax, waveSpeed);

      // Compute artificial viscosity term
      if (flag->d_artificialViscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;
        double c_bulk = std::sqrt(1.0 / (d_cm.K * rho_cur));
        p_q[idx]  = artificialBulkViscosity(trD, c_bulk, rho_cur, dx_ave);

#ifdef CHECK_ISFINITE
        if (!std::isfinite(p_q[idx])) {
          std::cout << "K = " << d_cm.K << " rho_cur = " << rho_cur
                    << " c_bulk = " << c_bulk << " D = " << D << "\n";
        }
#endif
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
JWLppMPM::addComputesAndRequires(Task*,
                                 const MPMMaterial*,
                                 const PatchSet*,
                                 const bool,
                                 const bool) const
{
}

void
JWLppMPM::allocateCMDataAddRequires(Task* task,
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
JWLppMPM::allocateCMDataAdd(DataWarehouse* new_dw,
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
JWLppMPM::carryForward(const PatchSubset* patches,
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
JWLppMPM::computeRhoMicroCM([[maybe_unused]] double pressure,
                            [[maybe_unused]] const double p_ref,
                            const MPMMaterial* matl,
                            [[maybe_unused]] double temperature,
                            [[maybe_unused]] double rho_guess)
{
  std::cout << "NO VERSION OF computeRhoMicroCM EXISTS YET FOR JWLppMPM"
            << "\n";
  // double rho_orig = d_cm.rho0; //matl->getInitialDensity();
  double rho_orig = matl->getInitialDensity();

  return rho_orig;
}

void
JWLppMPM::computePressEOSCM(const double rhoM,
                            double& pressure,
                            [[maybe_unused]] const double p_ref,
                            double& dp_drho,
                            double& tmp,
                            const MPMMaterial* matl,
                            [[maybe_unused]] double temperature)
{
  double A     = d_cm.A;
  double B     = d_cm.B;
  double R1    = d_cm.R1;
  double R2    = d_cm.R2;
  double omega = d_cm.omega;
  // double rho0 = d_cm.rho0;
  double rho0 = matl->getInitialDensity();
  double cv   = matl->getSpecificHeat();
  double V    = rho0 / rhoM;
  double P1   = A * std::exp(-R1 * V);
  double P2   = B * std::exp(-R2 * V);
  double P3   = omega * cv * tmp * rhoM;

  pressure = P1 + P2 + P3;

  dp_drho =
    (R1 * rho0 * P1 + R2 * rho0 * P2) / (rhoM * rhoM) + omega * cv * tmp;
}

// This is not yet implemented - JG- 7/26/10
double
JWLppMPM::getCompressibility()
{
  std::cout << "NO VERSION OF getCompressibility EXISTS YET FOR JWLppMPM"
            << "\n";
  return 1.0;
}

// This is the burn logic used in the reaction model  (more complex versions
//   are available -- see LS-DYNA manual)
//       df/dt = G (1-f) p^b
//       Forward Euler: f_{n+1} = f_n + G*(1-f_n)*p_n^b*delT
//       Backward Euler: f_{n+1} = f_n + G*(1-f_{n+1})*p_n^b*delT
//                       or, f_{n+1} = (f_n + G*p_n^b*delT)/(1 + G*p_n^b*delT)
//       Fourth-order R-K: f_{n+1} = f_n + 1/6(k1 + 2k2 + 2k3 + k4)
//         k1 = G*(1-f_n)*p_n^b*delT
//         k2 = G*(1-f_n-k1/2)*p_n^b*delT
//         k3 = G*(1-f_n-k2/2)*p_n^b*delT
//         k4 = G*(1-f_n-k3)*p_n^b*delT
// (ignition_pressure in previous versions hardcoded to 2.0e8 Pa)
void
JWLppMPM::computeUpdatedFractionAndPressure(const double& J_old,
                                            const double& J,
                                            const double& f_old_orig,
                                            const double& p_old_orig,
                                            const double& delT,
                                            double& f_new,
                                            double& p_new) const
{
  if ((p_old_orig > d_cm.ignition_pressure) &&
      (f_old_orig < d_cm.max_burned_frac)) {

    // std::cerr << setprecision(10) << scientific
    //     << " p_old = " << p_old_orig << " ignition = " <<
    //     d_cm.ignition_pressure
    //     << " f_old = " << f_old_orig << " max_burn = " <<
    //     d_cm.max_burned_frac
    //     << " f_old - max_f = " << f_old_orig - d_cm.max_burned_frac << "\n";
    int numCycles  = std::max(1, (int)std::ceil(delT / d_cm.max_burn_timestep));
    double delTinc = delT / ((double)numCycles);
    double delJ    = J / J_old;
    double delJinc = std::pow(delJ, 1.0 / ((double)numCycles));
    double p_old   = p_old_orig;
    double f_old   = f_old_orig;
    double J_new   = J_old;
    f_new          = f_old_orig;
    p_new          = p_old_orig;

    if (d_fastCompute) {
      // std::cerr << "Using Fast" << "\n";
      for (int ii = 0; ii < numCycles; ++ii) {

        // Compute Murnaghan and JWL pressures
        J_new *= delJinc;
        double pM   = computePressureMurnaghan(J_new);
        double pJWL = computePressureJWL(J_new);

        computeWithTwoStageBackwardEuler(
          J_new, f_old, p_old, delTinc, pM, pJWL, f_new, p_new);
        f_old = f_new;
        p_old = p_new;
      }
    } else {
      // std::cerr << "Using Newton" << "\n";
      for (int ii = 0; ii < numCycles; ++ii) {

        // Compute Murnaghan and JWL pressures
        J_new *= delJinc;
        double pM   = computePressureMurnaghan(J_new);
        double pJWL = computePressureJWL(J_new);

        computeWithNewtonIterations(
          J_new, f_old, p_old, delTinc, pM, pJWL, f_new, p_new);
        f_old = f_new;
        p_old = p_new;
      }
    }
    if (f_new > d_cm.max_burned_frac) {
      f_new       = d_cm.max_burned_frac;
      double pM   = computePressureMurnaghan(J);
      double pJWL = computePressureJWL(J);
      p_new       = pM * (1.0 - f_new) + pJWL * f_new;
    }
  } else {
    // Compute Murnaghan and JWL pressures
    double pM   = computePressureMurnaghan(J);
    double pJWL = computePressureJWL(J);

    //  The following computes a pressure for partially burned particles
    //  as a mixture of Murnaghan and JWL pressures, based on pProgressF
    //  This is as described in Eq. 5 of "JWL++: ..." by Souers, et al.
    f_new = f_old_orig;
    p_new = pM * (1.0 - f_new) + pJWL * f_new;
  }

  return;
}

//  This is the original two stage Backward Euler
void
JWLppMPM::computeWithTwoStageBackwardEuler([[maybe_unused]] const double& J,
                                           const double& f_old,
                                           const double& p_old,
                                           const double& delT,
                                           const double& pM,
                                           const double& pJWL,
                                           double& f_new,
                                           double& p_new) const
{
  double fac = (delT * d_cm.G) * std::pow(p_old, d_cm.b);

  // Backward Euler
  f_new = (f_old + fac) / (1.0 + fac);

  // Forward Euler
  // f_new = f_old + (1.0 - f_old)*fac;

  // Fourth-order R-K
  // double k1 = (1.0 - f_old)*fac;
  // double k2 = (1.0 - f_old - 0.5*k1)*fac;
  // double k3 = (1.0 - f_old - 0.5*k2)*fac;
  // double k4 = (1.0 - f_old - k3)*fac;
  // f_new = f_old + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);

  // if (f_new < 0.0) f_new = 0.0;
  if (f_new > d_cm.max_burned_frac)
    f_new = d_cm.max_burned_frac; // Max burned volume fraction

  //  The following computes a pressure for partially burned particles
  //  as a mixture of Murnaghan and JWL pressures, based on pProgressF
  //  This is as described in Eq. 5 of "JWL++: ..." by Souers, et al.
  p_new = pM * (1.0 - f_new) + pJWL * f_new;

  return;
}

//  This is the Newton iteration with Backward Euler
void
JWLppMPM::computeWithNewtonIterations(const double& J,
                                      const double& f_old,
                                      const double& p_old,
                                      const double& delT,
                                      const double& pM,
                                      const double& pJWL,
                                      double& f_new,
                                      double& p_new) const
{
  // Initialize matrices
  std::vector<double> G(2); // The vector [F_n+1 P_n+1]^T = 0
  FastMatrix JacobianG(2, 2);

  // Initial values of f and p
  f_new = f_old;
  p_new = p_old;

  // Set iteration controls
  int iter    = 0;
  double norm = 0.0;

  // Compute G
  computeG(J, f_old, f_new, p_new, pM, pJWL, delT, G);

  // Do Newton iterations
  FastMatrix Jinv(2, 2);
  std::vector<double> Finc(2);
  do {

    // Compute Jacobian of G
    computeJacobianG(J, f_new, p_new, pM, pJWL, delT, JacobianG);

    // Invert Jacobian of G
    Jinv.destructiveInvert(JacobianG);

    // Compute increment
    Jinv.multiply(G, Finc);

    // Update the variables
    f_new -= Finc[0];
    p_new -= Finc[1];

    // Compute G
    computeG(J, f_old, f_new, p_new, pM, pJWL, delT, G);

    // Compute L2 norm and increment iter
    norm = std::sqrt(G[0] * G[0] + G[1] * G[1]);
    iter++;

  } while ((norm > d_newtonIterTol) && (iter < d_newtonIterMax));

  if (iter > d_newtonIterMax) {
    std::cerr << "**JWLppMPM** Newton iterations failed to converge."
              << "\n";
    std::cerr << "iter = " << iter << " norm = " << norm
              << " tol = " << d_newtonIterTol << " p_new = " << p_new
              << " f_new = " << f_new << " p_old = " << p_old
              << " f_old = " << f_old << " J = " << J << "\n";
    std::cerr << " pM = " << pM << " pJWL = " << pJWL << " G = [" << G[0] << ","
              << G[1] << "]"
              << " JacobianG = [[" << JacobianG(0, 0) << "," << JacobianG(0, 1)
              << "],[" << JacobianG(1, 0) << "," << JacobianG(1, 1) << "]]"
              << "\n";
    std::cerr << " Jinv = [[" << Jinv(0, 0) << "," << Jinv(0, 1) << "],["
              << Jinv(1, 0) << "," << Jinv(1, 1) << "]]"
              << " Finc = [" << Finc[0] << "," << Finc[1] << "]"
              << "\n";
  }
  if (std::isnan(p_new) || std::isnan(f_new)) {
    std::cerr << "iter = " << iter << " norm = " << norm
              << " tol = " << d_newtonIterTol << " p_new = " << p_new
              << " f_new = " << f_new << " p_old = " << p_old
              << " f_old = " << f_old << " J = " << J << "\n";
    std::cerr << " pM = " << pM << " pJWL = " << pJWL << " G = [" << G[0] << ","
              << G[1] << "]"
              << " JacobianG = [[" << JacobianG(0, 0) << "," << JacobianG(0, 1)
              << "],[" << JacobianG(1, 0) << "," << JacobianG(1, 1) << "]]"
              << "\n";
    std::cerr << " Jinv = [[" << Jinv(0, 0) << "," << Jinv(0, 1) << "],["
              << Jinv(1, 0) << "," << Jinv(1, 1) << "]]"
              << " Finc = [" << Finc[0] << "," << Finc[1] << "]"
              << "\n";
    throw InvalidValue(
      "**JWLppMPM ERROR**: Nan in p_new/f_new value or no convergence",
      __FILE__,
      __LINE__);
  }

  return;
}

//------------------------------------------------------------------
// Compute G
//  G = [F_n+1 P_n+1]^T
//   F_n+1 = 0 = f_n+1 - f_n - G*(1 - f_n+1)*(p_n+1)^b*Delta t
//   P_n+1 = 0 = p_n+1 - (1 - f_n+1) p_m - f_n+1 p_jwl
//------------------------------------------------------------------
void
JWLppMPM::computeG([[maybe_unused]] const double& J,
                   const double& f_old,
                   const double& f_new,
                   const double& p_new,
                   const double& pM,
                   const double& pJWL,
                   const double& delT,
                   std::vector<double>& G) const
{
  double dfdt_new = computeBurnRate(f_new, p_new);
  double f_func   = f_new - f_old - dfdt_new * delT;
  double p_func   = p_new - (1.0 - f_new) * pM - f_new * pJWL;
  G[0]            = f_func;
  G[1]            = p_func;
  return;
}

//------------------------------------------------------------------
// Compute the Jacobian of G
//  J_G = [[dF_n+1/df_n+1 dF_n+1/dp_n+1];[dP_n+1/df_n+1 dP_n+1/dp_n+1]]
//   F_n+1 = 0 = f_n+1 - f_n - G*(1 - f_n+1)*(p_n+1)^b*Delta t
//   P_n+1 = 0 = p_n+1 - (1 - f_n+1) p_m - f_n+1 p_jwl
//   dF_n+1/df_n+1 = 1 + G*(p_n+1)^b*Delta t
//   dF_n+1/dp_n+1 =  b*G*(1 - f_n+1)*(p_n+1)^(b-1)*Delta t
//   dP_n+1/df_n+1 =  p_m - p_jwl
//   dP_n+1/dp_n+1 = 1
//------------------------------------------------------------------
void
JWLppMPM::computeJacobianG([[maybe_unused]] const double& J,
                           const double& f_new,
                           const double& p_new,
                           const double& pM,
                           const double& pJWL,
                           const double& delT,
                           FastMatrix& JacobianG) const
{
  double fac   = d_cm.G * std::pow(p_new, d_cm.b) * delT;
  double dF_df = 1.0 + fac;
  double dF_dp = d_cm.b * (1.0 - f_new) * (fac / p_new);
  double dP_df = pM - pJWL;
  JacobianG(0, 0) = dF_df;
  JacobianG(0, 1) = dF_dp;
  JacobianG(1, 0) = dP_df;
  JacobianG(1, 1) = 1.0;
  return;
}

//------------------------------------------------------------------
//  df/dt = G (1-f) p^b
//------------------------------------------------------------------
double
JWLppMPM::computeBurnRate(const double& f, const double& p) const
{
  double dfdt = d_cm.G * (1.0 - f) * std::pow(p, d_cm.b);
  return dfdt;
}

//------------------------------------------------------------------
//  p_m = (1/nK) [J^(-n) - 1]
//------------------------------------------------------------------
double
JWLppMPM::computePressureMurnaghan(const double& J) const
{
  double pM = (1.0 / (d_cm.n * d_cm.K)) * (std::pow(J, -d_cm.n) - 1.0);
  return pM;
}

//------------------------------------------------------------------
// p_jwl = A exp(-R1 J) + B exp(-R2 J) + C J^[-(1+omega)]
//------------------------------------------------------------------
double
JWLppMPM::computePressureJWL(const double& J) const
{
  double one_plus_omega                 = 1.0 + d_cm.omega;
  double A_e_to_the_R1_rho0_over_rhoM   = d_cm.A * std::exp(-d_cm.R1 * J);
  double B_e_to_the_R2_rho0_over_rhoM   = d_cm.B * std::exp(-d_cm.R2 * J);
  double C_rho_rat_tothe_one_plus_omega = d_cm.C * std::pow(J, -one_plus_omega);

  double pJWL = A_e_to_the_R1_rho0_over_rhoM + B_e_to_the_R2_rho0_over_rhoM +
                C_rho_rat_tothe_one_plus_omega;
  return pJWL;
}
