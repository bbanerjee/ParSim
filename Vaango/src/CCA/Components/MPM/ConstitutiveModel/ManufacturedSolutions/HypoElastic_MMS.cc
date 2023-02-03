/*
 * The MIT License
 *
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

#include <CCA/Components/MPM/ConstitutiveModel/ManufacturedSolutions/HypoElastic_MMS.h>
#include <CCA/Components/MPM/Core/MPMMaterial.h>

#include <CCA/Ports/DataWarehouse.h>
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
#include <Core/ProblemSpec/ProblemSpec.h>

#include <Core/Malloc/Allocator.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/Math/Short27.h>
#include <fstream>
#include <iostream>

using namespace Uintah;

HypoElastic_MMS::HypoElastic_MMS(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  // Hardcoded density (1.7 gm/cc) and elastic properties (K = 60 MPa, G = 100
  // MPa)
  d_cm.rho0 = 1700.0;
  d_cm.kappa = 6.0e7;
  d_cm.mu = 1.0e8;
  d_cm.lambda = d_cm.kappa - (d_cm.mu * 2.0) / 3.0;
  d_cm.cp = std::sqrt((d_cm.lambda + 2.0 * d_cm.mu) / d_cm.rho0);

  // Hardcoded amplitude (0.01 m) and frequency (10000 rad/s)
  d_cm.alpha = 0.01;
  d_cm.omega = 1000;
}

HypoElastic_MMS::HypoElastic_MMS(const HypoElastic_MMS* cm)
  : ConstitutiveModel(cm)
{
  d_cm.rho0 = cm->d_cm.rho0;
  d_cm.kappa = cm->d_cm.kappa;
  d_cm.mu = cm->d_cm.mu;
  d_cm.lambda = cm->d_cm.lambda;
  d_cm.cp = cm->d_cm.cp;

  // Hardcoded amplitude (150 m/s) and frequency (1000 rad/s)
  d_cm.alpha = cm->d_cm.alpha;
  d_cm.omega = cm->d_cm.omega;
}

HypoElastic_MMS::~HypoElastic_MMS() = default;

void
HypoElastic_MMS::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "hypo_elastic_mms");
  }
}

std::unique_ptr<ConstitutiveModel>
HypoElastic_MMS::clone()
{
  return std::make_unique<HypoElastic_MMS>(*this);
}

void
HypoElastic_MMS::initializeCMData(const Patch* patch, const MPMMaterial* matl,
                                  DataWarehouse* new_dw)
{
  std::string mms_type = flag->d_mmsType;
  if (mms_type != "none") {

    if (mms_type == "UniaxialStrainHarmonic") {

      initStressAndDefGradUniaxialStrainHarmonic(patch, matl, new_dw);

    } else if (mms_type == "UniaxialStrainHomogeneousLinear") {

      initStressAndDefGradUniaxialStrainHomogeneous(patch, matl, new_dw);

    } else if (mms_type == "UniaxialStrainHomogeneousQuadratic") {

      initStressAndDefGradUniaxialStrainHomogeneous(patch, matl, new_dw);

    } else {

      std::ostringstream out;
      out << "**ERROR** Hypoelastic version of MMS:" << mms_type
          << " does not exist";
      throw InvalidValue(out.str(), __FILE__, __LINE__);
    }
  }
  computeStableTimestep(patch, matl, new_dw);
}

void
HypoElastic_MMS::initStressAndDefGradUniaxialStrainHarmonic(
  const Patch* patch, const MPMMaterial* matl, DataWarehouse* new_dw)
{
  Matrix3 Zero(0.0), One;
  One.Identity();

  double omega_cp = d_cm.omega / d_cm.cp;
  double M = (d_cm.lambda + 2.0 * d_cm.mu);

  int matIndex = matl->getDWIndex();
  ParticleSubset* pset = new_dw->getParticleSubset(matIndex, patch);

  constParticleVariable<Point> pX;
  constParticleVariable<Vector> pDisp, pVelocity;
  new_dw->get(pX, lb->pXLabel, pset);
  new_dw->get(pDisp, lb->pDispLabel, pset);
  new_dw->get(pVelocity, lb->pVelocityLabel, pset);

  ParticleVariable<double> pVolume, pdTdt;
  ParticleVariable<Matrix3> pDefGrad, pStress;
  new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel, pset);
  new_dw->getModifiable(pDefGrad, lb->pDefGradLabel, pset);
  new_dw->getModifiable(pVolume, lb->pVolumeLabel, pset);
  new_dw->allocateAndPut(pStress, lb->pStressLabel, pset);

  for (int idx : *pset) {
    pdTdt[idx] = 0.0;

    // Get the current position
    double X = pX[idx].x() - pDisp[idx].x();

    // Deformation gradient
    double omega_X_cp = omega_cp * X;
    double U0i = d_cm.alpha * std::sin(omega_X_cp);
    double F11 = 1.0 - omega_cp * U0i;
    // std::cout << " F11 = " << F11 << std::endl;
    pDefGrad[idx] = Matrix3(F11, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

    // Volume
    double J = pDefGrad[idx].Determinant();
    pVolume[idx] *= J;

    // Cauchy stress
    double sigma11 = M * std::log(F11);
    pStress[idx] = Matrix3(sigma11, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  }
}

void
HypoElastic_MMS::initStressAndDefGradUniaxialStrainHomogeneous(
  const Patch* patch, const MPMMaterial* matl, DataWarehouse* new_dw)
{
  Matrix3 Zero(0.0), One;
  One.Identity();

  double M = (d_cm.lambda + 2.0 * d_cm.mu);

  int matIndex = matl->getDWIndex();
  ParticleSubset* pset = new_dw->getParticleSubset(matIndex, patch);

  ParticleVariable<double> pVolume, pdTdt;
  ParticleVariable<Matrix3> pDefGrad, pStress;
  new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel, pset);
  new_dw->getModifiable(pDefGrad, lb->pDefGradLabel, pset);
  new_dw->getModifiable(pVolume, lb->pVolumeLabel, pset);
  new_dw->allocateAndPut(pStress, lb->pStressLabel, pset);

  for (int idx : *pset) {
    pdTdt[idx] = 0.0;

    // Deformation gradient
    double F11 = 1.0;
    pDefGrad[idx] = Matrix3(F11, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

    // Volume
    double J = pDefGrad[idx].Determinant();
    pVolume[idx] *= J;

    // Cauchy stress
    double sigma11 = M * std::log(F11);
    pStress[idx] = Matrix3(sigma11, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  }
}

void
HypoElastic_MMS::allocateCMDataAddRequires(Task* task, const MPMMaterial* matl,
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
HypoElastic_MMS::allocateCMDataAdd(DataWarehouse* new_dw,
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
HypoElastic_MMS::addParticleState(std::vector<const VarLabel*>& from,
                                  std::vector<const VarLabel*>& to)
{
}

void
HypoElastic_MMS::computeStableTimestep(const Patch* patch,
                                       const MPMMaterial* matl,
                                       DataWarehouse* new_dw)
{
  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of computeStressTensor
  Vector dx = patch->dCell();
  int matIndex = matl->getDWIndex();
  ParticleSubset* pset = new_dw->getParticleSubset(matIndex, patch);
  constParticleVariable<double> pMass, pVolume;
  constParticleVariable<Vector> pVelocity;

  new_dw->get(pMass, lb->pMassLabel, pset);
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pVelocity, lb->pVelocityLabel, pset);

  double c_dil = 0.0;
  Vector WaveSpeed(1.e-12, 1.e-12, 1.e-12);

  double G = d_cm.mu;
  double bulk = d_cm.kappa;
  for (int idx : *pset) {
    // Compute wave speed at each particle, store the maximum
    c_dil = sqrt((bulk + 4. * G / 3.) * pVolume[idx] / pMass[idx]);
    WaveSpeed = Vector(Max(c_dil + fabs(pVelocity[idx].x()), WaveSpeed.x()),
                       Max(c_dil + fabs(pVelocity[idx].y()), WaveSpeed.y()),
                       Max(c_dil + fabs(pVelocity[idx].z()), WaveSpeed.z()));
  }
  WaveSpeed = dx / WaveSpeed;
  double delT_new = WaveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

void
HypoElastic_MMS::computeStressTensor(const PatchSubset* patches,
                                     const MPMMaterial* matl,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw)
{
  double rho_orig = d_cm.rho0;
  double G = d_cm.mu;
  double bulk = d_cm.kappa;
  double M = bulk + 4.0 / 3.0 * G;

  double onethird = (1.0 / 3.0);
  Matrix3 Identity;
  Identity.Identity();

  int matIndex = matl->getDWIndex();

  for (int p = 0; p < patches->size(); p++) {
    double se = 0.0;
    const Patch* patch = patches->get(p);
    Vector dx = patch->dCell();

    double c_dil = 0.0;
    Vector WaveSpeed(1.e-12, 1.e-12, 1.e-12);

    // Get particles in this patch
    ParticleSubset* pset = old_dw->getParticleSubset(matIndex, patch);

    // Get particle data
    delt_vartype delT;
    old_dw->get(delT, lb->delTLabel, getLevel(patches));
    constParticleVariable<Point> px;
    constParticleVariable<double> pMass, pVolume, pVolume_new;
    constParticleVariable<Vector> pVelocity;
    constParticleVariable<Matrix3> pVelGrad, pDefGrad, pDefGrad_new, pStress;
    constParticleVariable<Matrix3> pSize;

    old_dw->get(px, lb->pXLabel, pset);
    old_dw->get(pStress, lb->pStressLabel, pset);
    old_dw->get(pSize, lb->pSizeLabel, pset);
    old_dw->get(pMass, lb->pMassLabel, pset);
    old_dw->get(pVolume, lb->pVolumeLabel, pset);
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);
    old_dw->get(pDefGrad, lb->pDefGradLabel, pset);

    new_dw->get(pVolume_new, lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);
    new_dw->get(pVelGrad, lb->pVelGradLabel_preReloc, pset);

    // Create space for updated particle data
    ParticleVariable<double> pdTdt, p_q;
    ParticleVariable<Matrix3> pStress_new;
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);

    // Loop thru particles
    for (int idx : *pset) {
      pdTdt[idx] = 0.0;

      // Calculate rate of deformation D, and deviatoric rate DPrime,
      Matrix3 D = (pVelGrad[idx] + pVelGrad[idx].Transpose()) * .5;
      double DTrace = D.Trace();
      Matrix3 DPrime = D - Identity * onethird * DTrace;

      // get the volumetric part of the deformation
      double J = pDefGrad[idx].Determinant();

      // Compute the local sound speed
      double rho_cur = rho_orig / J;
      c_dil = sqrt(M / rho_cur);

      // This is the (updated) Cauchy stress
      pStress_new[idx] =
        pStress[idx] + (DPrime * 2. * G + Identity * bulk * DTrace) * delT;

      // Compute the strain energy for all the particles
      Matrix3 AvgStress = (pStress_new[idx] + pStress[idx]) * .5;

      double e = (D(0, 0) * AvgStress(0, 0) + D(1, 1) * AvgStress(1, 1) +
                  D(2, 2) * AvgStress(2, 2) +
                  2. * (D(0, 1) * AvgStress(0, 1) + D(0, 2) * AvgStress(0, 2) +
                        D(1, 2) * AvgStress(1, 2))) *
                 pVolume_new[idx] * delT;

      se += e;

      // Compute wave speed at each particle, store the maximum
      Vector pVelocity_idx = pVelocity[idx];
      WaveSpeed = Vector(Max(c_dil + fabs(pVelocity_idx.x()), WaveSpeed.x()),
                         Max(c_dil + fabs(pVelocity_idx.y()), WaveSpeed.y()),
                         Max(c_dil + fabs(pVelocity_idx.z()), WaveSpeed.z()));

      // Compute artificial viscosity term
      if (flag->d_artificialViscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;
        double c_bulk = sqrt(bulk / rho_cur);
        p_q[idx] = artificialBulkViscosity(DTrace, c_bulk, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.;
      }
    } // end loop over particles

    WaveSpeed = dx / WaveSpeed;
    double delT_new = WaveSpeed.minComponent();
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(se), lb->StrainEnergyLabel);
    }
  }
}

void
HypoElastic_MMS::carryForward(const PatchSubset* patches,
                              const MPMMaterial* matl, DataWarehouse* old_dw,
                              DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    int matIndex = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(matIndex, patch);

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
HypoElastic_MMS::addComputesAndRequires(Task* task, const MPMMaterial* matl,
                                        const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForHypoExplicit(task, matlset, patches);
}

void
HypoElastic_MMS::addComputesAndRequires(Task*, const MPMMaterial*,
                                        const PatchSet*, const bool,
                                        const bool) const
{
}

double
HypoElastic_MMS::computeRhoMicroCM(double pressure, const double p_ref,
                                   const MPMMaterial* matl, double temperature,
                                   double rho_guess)
{
  double rho_orig = d_cm.rho0;
  double p_gauge = pressure - p_ref;
  double bulk = d_cm.kappa;
  double rho_cur = rho_orig / (1 - p_gauge / bulk);

  return rho_cur;
}

void
HypoElastic_MMS::computePressEOSCM(double rho_cur, double& pressure,
                                   double p_ref, double& dp_drho, double& tmp,
                                   const MPMMaterial* matl, double temperature)
{
  double bulk = d_cm.kappa;
  double rho_orig = d_cm.rho0;
  double p_g = bulk * (1.0 - rho_orig / rho_cur);
  pressure = p_ref + p_g;
  dp_drho = bulk * rho_orig / (rho_cur * rho_cur);
  tmp = bulk / rho_cur; // speed of sound squared
}

double
HypoElastic_MMS::getCompressibility()
{
  return 1.0 / d_cm.kappa;
}
