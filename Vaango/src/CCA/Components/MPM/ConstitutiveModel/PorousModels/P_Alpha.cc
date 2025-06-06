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
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/PorousModels/P_Alpha.h>
#include <Core/Grid/Variables/VarTypes.h>
#include<CCA/Components/MPM/Core/MPMLabel.h>

using namespace Uintah;

P_Alpha::P_Alpha(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  //  For P-alpha part of response
  ps->require("Ps", d_modelParam.Ps);
  ps->require("Pe", d_modelParam.Pe);
  ps->require("rhoS", d_modelParam.rhoS);
  // Compute alpha0 from material density and rhoS - Jim 9/8/2011
  // ps->require("alpha0",    d_modelParam.alpha0);
  ps->require("K0", d_modelParam.K0);
  ps->require("Ks", d_modelParam.Ks);
  //  For M-G part of response
  ps->require("T_0", d_modelParam.T_0);
  ps->require("C_0", d_modelParam.C_0);
  ps->require("Gamma_0", d_modelParam.Gamma_0);
  ps->require("S_alpha", d_modelParam.S_alpha);
  // For the unloading response
  ps->getWithDefault("Ku", d_modelParam.Ku, .1 * d_modelParam.K0);

  pAlphaLabel =
    VarLabel::create("p.alpha", ParticleVariable<double>::getTypeDescription());
  pAlphaMinLabel = VarLabel::create(
    "p.pAlphaMin", ParticleVariable<double>::getTypeDescription());
  pAlphaMinLabel_preReloc = VarLabel::create(
    "p.pAlphaMin+", ParticleVariable<double>::getTypeDescription());
  pTempAlpha1Label = VarLabel::create(
    "p.pTempAlpha1", ParticleVariable<double>::getTypeDescription());
  pTempAlpha1Label_preReloc = VarLabel::create(
    "p.pTempAlpha1+", ParticleVariable<double>::getTypeDescription());
}

P_Alpha::P_Alpha(const P_Alpha* cm)
  : ConstitutiveModel(cm)
{
  d_modelParam.Ps   = cm->d_modelParam.Ps;
  d_modelParam.Pe   = cm->d_modelParam.Pe;
  d_modelParam.rhoS = cm->d_modelParam.rhoS;
  // Compute alpha0 from material density and rhoS - Jim 9/8/2011
  // d_modelParam.alpha0 = cm->d_modelParam.alpha0;
  d_modelParam.K0 = cm->d_modelParam.K0;
  d_modelParam.Ks = cm->d_modelParam.Ks;
  d_modelParam.Ku = cm->d_modelParam.Ku;

  d_modelParam.T_0     = cm->d_modelParam.T_0;
  d_modelParam.C_0     = cm->d_modelParam.C_0;
  d_modelParam.Gamma_0 = cm->d_modelParam.Gamma_0;
  d_modelParam.S_alpha = cm->d_modelParam.S_alpha;

  pAlphaLabel =
    VarLabel::create("p.alpha", ParticleVariable<double>::getTypeDescription());
  pAlphaMinLabel = VarLabel::create(
    "p.pAlphaMin", ParticleVariable<double>::getTypeDescription());
  pAlphaMinLabel_preReloc = VarLabel::create(
    "p.pAlphaMin+", ParticleVariable<double>::getTypeDescription());
  pTempAlpha1Label = VarLabel::create(
    "p.pTempAlpha1", ParticleVariable<double>::getTypeDescription());
  pTempAlpha1Label_preReloc = VarLabel::create(
    "p.pTempAlpha1+", ParticleVariable<double>::getTypeDescription());
}

P_Alpha::~P_Alpha()
{
  VarLabel::destroy(pAlphaLabel);
  VarLabel::destroy(pAlphaMinLabel);
  VarLabel::destroy(pAlphaMinLabel_preReloc);
  VarLabel::destroy(pTempAlpha1Label);
  VarLabel::destroy(pTempAlpha1Label_preReloc);
}

std::unique_ptr<ConstitutiveModel>
P_Alpha::clone()
{
  return std::make_unique<P_Alpha>(*this);
}

void
P_Alpha::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "p_alpha");
  }

  cm_ps->appendElement("Ps", d_modelParam.Ps);
  cm_ps->appendElement("Pe", d_modelParam.Pe);
  cm_ps->appendElement("rhoS", d_modelParam.rhoS);
  // Compute alpha0 from material density and rhoS - Jim 9/8/2011
  // cm_ps->appendElement("alpha0",  d_modelParam.alpha0);
  cm_ps->appendElement("K0", d_modelParam.K0);
  cm_ps->appendElement("Ks", d_modelParam.Ks);
  cm_ps->appendElement("Ku", d_modelParam.Ku);
  cm_ps->appendElement("T_0", d_modelParam.T_0);
  cm_ps->appendElement("C_0", d_modelParam.C_0);
  cm_ps->appendElement("Gamma_0", d_modelParam.Gamma_0);
  cm_ps->appendElement("S_alpha", d_modelParam.S_alpha);
}

void
P_Alpha::addParticleState(std::vector<const VarLabel*>& from,
                          std::vector<const VarLabel*>& to)
{
  // Add the local particle state data for this constitutive model.
  from.push_back(pAlphaMinLabel);
  from.push_back(pTempAlpha1Label);
  to.push_back(pAlphaMinLabel_preReloc);
  to.push_back(pTempAlpha1Label_preReloc);
}

void
P_Alpha::addInitialComputesAndRequires(Task* task,
                                       const MPMMaterial* matl,
                                       const PatchSet*) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pAlphaMinLabel, matlset);
  task->computes(pTempAlpha1Label, matlset);
}

void
P_Alpha::initializeCMData(const Patch* patch,
                          const MPMMaterial* matl,
                          DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);

  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  ParticleVariable<double> pAlphaMin, pTempAlpha1;
  new_dw->allocateAndPut(pAlphaMin, pAlphaMinLabel, pset);
  new_dw->allocateAndPut(pTempAlpha1, pTempAlpha1Label, pset);

  double rhoS       = d_modelParam.rhoS;
  double rho_orig   = matl->getInitialDensity();
  double pAlphaMin0 = rhoS / rho_orig;

  for (auto pidx : *pset) {
    // Compute alpha0 from material density and rhoS - Jim 9/8/2011
    // pAlphaMin[pidx]    = d_modelParam.alpha0;
    pAlphaMin[pidx]   = pAlphaMin0;
    pTempAlpha1[pidx] = d_modelParam.T_0;
  }

  computeStableTimestep(patch, matl, new_dw);
}

void
P_Alpha::computeStableTimestep(const Patch* patch,
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

  double K0 = d_modelParam.K0;
  Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);

  for (int idx : *pset) {
    double rhoM = pMass[idx] / pVolume[idx];
    double tmp  = K0 / rhoM;

    // Compute wave speed at each particle, store the maximum
    double c_dil  = std::sqrt(tmp);
    Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
    waveSpeed     = Max(velMax, waveSpeed);
  }
  waveSpeed       = dx / waveSpeed;
  double delT_new = waveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

void
P_Alpha::addComputesAndRequires(Task* task,
                                const MPMMaterial* matl,
                                const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForExplicit(task, matlset, patches);

  task->needs(Task::OldDW, pAlphaMinLabel, matlset, Ghost::None);
  task->needs(Task::OldDW, pTempAlpha1Label, matlset, Ghost::None);

  task->computes(pAlphaLabel, matlset);
  task->computes(pAlphaMinLabel_preReloc, matlset);
  task->computes(pTempAlpha1Label_preReloc, matlset);
}

void
P_Alpha::computeStressTensor(const PatchSubset* patches,
                             const MPMMaterial* matl,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw)
{
  double rho_orig = matl->getInitialDensity();
  double Ps       = d_modelParam.Ps;
  double Pe       = d_modelParam.Pe;
  // Compute alpha0 from material density and rhoS - Jim 9/8/2011
  // double alpha0 = d_modelParam.alpha0;
  double alpha0 = d_modelParam.rhoS / rho_orig;
  double K0     = d_modelParam.K0;
  double Ks     = d_modelParam.Ks;
  double Ku     = d_modelParam.Ku;
  double rhoS   = d_modelParam.rhoS;

  double cv     = matl->getSpecificHeat();
  double rhoP   = rho_orig / (1. - Pe / K0);
  double alphaP = rhoS / rhoP;

  double cs = std::sqrt(Ks / rhoS);
  double ce = std::sqrt(K0 / rho_orig);

  int matID = matl->getDWIndex();

  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  for (int pp = 0; pp < patches->size(); pp++) {

    const Patch* patch = patches->get(pp);
    Vector dx          = patch->dCell();

    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    constParticleVariable<double> pAlphaMin_old, pTemperature, pTempAlpha1_old;
    constParticleVariable<Vector> pVelocity;
    constParticleVariable<Matrix3> pDefGrad_old, pDefGrad_new, pVelGrad;
    old_dw->get(pAlphaMin_old, pAlphaMinLabel, pset);
    old_dw->get(pTemperature, lb->pTemperatureLabel, pset);
    old_dw->get(pTempAlpha1_old, pTempAlpha1Label, pset);
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);
    old_dw->get(pDefGrad_old, lb->pDefGradLabel, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);
    new_dw->get(pVelGrad, lb->pVelGradLabel_preReloc, pset);

    ParticleVariable<double> pAlpha_new, pAlphaMin_new, pTempAlpha1_new;
    ParticleVariable<double> pdTdt, p_q;
    ParticleVariable<Matrix3> pStress_new;
    new_dw->allocateAndPut(pAlpha_new, pAlphaLabel, pset);
    new_dw->allocateAndPut(pAlphaMin_new, pAlphaMinLabel_preReloc, pset);
    new_dw->allocateAndPut(pTempAlpha1_new, pTempAlpha1Label_preReloc, pset);
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);

    double strainEnergy = 0.;
    Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);

    for (int idx : *pset) {

      double Jold = pDefGrad_old[idx].Determinant();
      double Jnew = pDefGrad_new[idx].Determinant();
      double Jinc = Jnew / Jold;
      double rhoM = rho_orig / Jnew;

      double alpha       = rhoS / rhoM;
      pAlphaMin_new[idx] = std::min(alpha, pAlphaMin_old[idx]);
      pAlpha_new[idx]    = alpha;

      double p       = 0.;
      double dAel_dp = 0.;
      double c = cs; // default to unstressed solid material speed of sound
      double dTdt_plas = 0., dTdt_MG = 0.;

      if (alpha < alpha0 && alpha >= 1.0) {
        if (alpha <= pAlphaMin_old[idx]) { // loading
          if (alpha <= alpha0 && alpha > alphaP) {
            // elastic response
            p = K0 * (1. - rho_orig / rhoM);
            c = std::sqrt(K0 / rhoM);
          } else if (alpha <= alphaP && alpha > 1.0) {
            // crushing out the voids
            p = Ps - (Ps - Pe) * std::sqrt((alpha - 1.) / (alphaP - 1.0));
            c = cs + (ce - cs) * ((alpha - 1.) / (alpha0 - 1.));
            dTdt_plas = (-p) * (Jinc - 1.) * (1. / (rhoM * cv)) / delT;
          }
        } else { // alpha < pAlphaMin, unloading
          if (alpha < alpha0 && alpha >= alphaP &&
              pAlphaMin_old[idx] >= alphaP) {
            // still in initial elastic response
            p = K0 * (1. - rho_orig / rhoM);
            c = std::sqrt(K0 / rhoM);
          } else if ((alpha < alphaP && alpha > 1.0) ||
                     pAlphaMin_old[idx] < alphaP) {
            // First, get plastic pressure
            p = Ps - (Ps - Pe) * std::sqrt((alpha - 1.) / (alphaP - 1.0));
            double h = 1. + (ce - cs) * (alpha - 1.0) / (cs * (alpha0 - 1.));
            dAel_dp  = ((alpha * alpha) / Ks) * (1. - 1. / (h * h));
            // Limit the unloading modulus, mostly to avoid numerical issues
            dAel_dp     = std::min(dAel_dp, -1. / Ks);
            double dPel = (alpha - pAlphaMin_old[idx]) / dAel_dp;
            p += dPel;
            c = cs + (ce - cs) * ((alpha - 1.) / (alpha0 - 1.));
          }
        }
        pTempAlpha1_new[idx] = pTemperature[idx];
      }

      // Mie-Gruneisen response for fully densified solid
      if (alpha < 1.0 || pAlphaMin_new[idx] < 1.0) {
        // Get the state data
        double Gamma_0 = d_modelParam.Gamma_0; // 1.54
        double C_0     = d_modelParam.C_0;     // 4029.
        double S_alpha = d_modelParam.S_alpha; // 1.237;

        // Calc. zeta
        double zeta = (rhoM / rhoS - 1.0);

        // Calculate internal energy E
        double E = (cv) * (pTemperature[idx] - pTempAlpha1_old[idx]) * rhoS;

        // Calculate the pressure
        p = Gamma_0 * E;
        if (rhoM != rhoS) {
          double numer =
            rhoS * (C_0 * C_0) * (1.0 / zeta + (1.0 - 0.5 * Gamma_0));
          double denom = 1.0 / zeta - (S_alpha - 1.0);
          if (denom == 0.0) {
            std::cout << "rho_0 = " << rhoS << " zeta = " << zeta
                      << " numer = " << numer << std::endl;
            denom = 1.0e-5;
          }
          p += numer / (denom * denom);
        }
        p = Ps + p;

        double DTrace = pVelGrad[idx].Trace();
        dTdt_MG       = -pTemperature[idx] * Gamma_0 * rhoS * DTrace / rhoM;
        pTempAlpha1_new[idx] = pTempAlpha1_old[idx];
      }

      // Unloading cases that get into either alpha > alpha0, or negative P
      if (alpha > alpha0 || p < 0) {
        // This still may need some work - Jim (2/11/2011)
        // I think I've improved things, but
        // the plastic work (dTdt_plas) still needs work - Jim (9/8/2011)
        double rho_max = std::min(rhoS, rhoS / pAlphaMin_new[idx]);
        p              = Ku * (1. - rho_max / rhoM);
      }

      //      double etime = d_mat_manager->getElapsedTime();
      //      std::cout << "12345 " << " " << etime << " " << alpha << " " <<
      //      pTemperature[idx] << " " << p << std::endl;

      // Compute artificial viscosity term
      if (flag->d_artificialViscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;
        double c_bulk = c;
        double DTrace = pVelGrad[idx].Trace();
        p_q[idx]      = artificialBulkViscosity(DTrace, c_bulk, rhoM, dx_ave);
      } else {
        p_q[idx] = 0.;
      }

      pStress_new[idx] = Vaango::Util::Identity * (-p);

      // Temp increase
      pdTdt[idx] = dTdt_plas + dTdt_MG;

      double c_dil  = c;
      Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
      waveSpeed     = Max(velMax, waveSpeed);
    }

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
P_Alpha::allocateCMDataAddRequires(Task* task,
                                   const MPMMaterial* matl,
                                   const PatchSet* patches,
                                   MPMLabel*) const
{
  std::ostringstream err;
  err << "**INTERNAL ERROR** No version of allocateCMDataAddRequires "
      << "exists yet for P_Alpha\n";
  throw InternalError(err.str(), __FILE__, __LINE__);
}

void
P_Alpha::allocateCMDataAdd(DataWarehouse* new_dw,
                           ParticleSubset* addset,
                           ParticleLabelVariableMap* newState,
                           ParticleSubset* delset,
                           DataWarehouse*)
{
  std::ostringstream err;
  err << "**INTERNAL ERROR** No version of allocateCMDataAdd "
      << "exists yet for P_Alpha\n";
  throw InternalError(err.str(), __FILE__, __LINE__);
}

double
P_Alpha::computeRhoMicroCM(double press,
                           const double Temp,
                           const MPMMaterial* matl,
                           double temperature,
                           double rho_guess)
{
  std::ostringstream err;
  err << "**INTERNAL ERROR** No version of computeRhoMicroCM exists yet for "
         "P_Alpha\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return 0.0;
}

void
P_Alpha::computePressEOSCM(double rhoM,
                           double& pressure,
                           double Temp,
                           double& dp_drho,
                           double& tmp,
                           const MPMMaterial*,
                           double temperature)
{
  std::ostringstream err;
  err << "**INTERNAL ERROR** No version of computePressEOSCM exists yet for "
         "P_Alpha\n";
  throw InternalError(err.str(), __FILE__, __LINE__);
}

double
P_Alpha::getCompressibility()
{
  std::ostringstream err;
  err << "**INTERNAL ERROR** No version of getCompressibility exists yet for "
         "P_Alpha\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return 0.0;
}
