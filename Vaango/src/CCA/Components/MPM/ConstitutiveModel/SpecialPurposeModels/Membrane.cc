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

#include <CCA/Components/MPM/ConstitutiveModel/Utilities/Constants.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/SpecialPurposeModels/Membrane.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Labels/MPMLabel.h>
#include <Core/Math/Matrix3.h>
#include <fstream>
#include <iostream>

using namespace Uintah;

Membrane::Membrane(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  ps->require("bulk_modulus", d_modelParam.bulk);
  ps->require("shear_modulus", d_modelParam.shear);

  pDefGradInPlaneLabel = VarLabel::create(
    "p.defgrad_in_plane", ParticleVariable<Matrix3>::getTypeDescription());

  pDefGradInPlaneLabel_preReloc = VarLabel::create(
    "p.defgrad_in_plane+", ParticleVariable<Matrix3>::getTypeDescription());

  pTang1Label =
    VarLabel::create("p.tang1", ParticleVariable<Vector>::getTypeDescription());

  pTang2Label =
    VarLabel::create("p.tang2", ParticleVariable<Vector>::getTypeDescription());

  pNormLabel =
    VarLabel::create("p.norm", ParticleVariable<Vector>::getTypeDescription());

  pTang1Label_preReloc = VarLabel::create(
    "p.tang1+", ParticleVariable<Vector>::getTypeDescription());

  pTang2Label_preReloc = VarLabel::create(
    "p.tang2+", ParticleVariable<Vector>::getTypeDescription());

  pNormLabel_preReloc =
    VarLabel::create("p.norm+", ParticleVariable<Vector>::getTypeDescription());
}

Membrane::Membrane(const Membrane* cm)
  : ConstitutiveModel(cm)
{
  d_modelParam.bulk  = cm->d_modelParam.bulk;
  d_modelParam.shear = cm->d_modelParam.shear;

  pDefGradInPlaneLabel = VarLabel::create(
    "p.defgrad_in_plane", ParticleVariable<Matrix3>::getTypeDescription());

  pDefGradInPlaneLabel_preReloc = VarLabel::create(
    "p.defgrad_in_plane+", ParticleVariable<Matrix3>::getTypeDescription());

  pTang1Label =
    VarLabel::create("p.tang1", ParticleVariable<Vector>::getTypeDescription());

  pTang2Label =
    VarLabel::create("p.tang2", ParticleVariable<Vector>::getTypeDescription());

  pNormLabel =
    VarLabel::create("p.norm", ParticleVariable<Vector>::getTypeDescription());

  pTang1Label_preReloc = VarLabel::create(
    "p.tang1+", ParticleVariable<Vector>::getTypeDescription());

  pTang2Label_preReloc = VarLabel::create(
    "p.tang2+", ParticleVariable<Vector>::getTypeDescription());

  pNormLabel_preReloc =
    VarLabel::create("p.norm+", ParticleVariable<Vector>::getTypeDescription());
}

Membrane::~Membrane()
{
  // Destructor
  VarLabel::destroy(pDefGradInPlaneLabel);
  VarLabel::destroy(pDefGradInPlaneLabel_preReloc);
  VarLabel::destroy(pTang1Label);
  VarLabel::destroy(pTang1Label_preReloc);
  VarLabel::destroy(pTang2Label);
  VarLabel::destroy(pTang2Label_preReloc);
  VarLabel::destroy(pNormLabel);
  VarLabel::destroy(pNormLabel_preReloc);
}

void
Membrane::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "membrane");
  }

  cm_ps->appendElement("bulk_modulus", d_modelParam.bulk);
  cm_ps->appendElement("shear_modulus", d_modelParam.shear);
}

Membrane*
Membrane::clone()
{
  return scinew Membrane(*this);
}

void
Membrane::addParticleState(std::vector<const VarLabel*>& from,
                           std::vector<const VarLabel*>& to)
{
  // Add the local particle state data for this constitutive model.
  from.push_back(pDefGradInPlaneLabel);
  to.push_back(pDefGradInPlaneLabel_preReloc);
}

void
Membrane::addInitialComputesAndRequires(Task* task,
                                        const MPMMaterial* matl,
                                        const PatchSet*) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pDefGradInPlaneLabel, matlset);
}

void
Membrane::initializeCMData(const Patch* patch,
                           const MPMMaterial* matl,
                           DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);

  // Put stuff in here to initialize each particle's
  // constitutive model parameters and deformationMeasure
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);
  ParticleVariable<Matrix3> pDefGradIP;
  new_dw->allocateAndPut(pDefGradIP, pDefGradInPlaneLabel, pset);

  for (auto pidx : *pset) {
    pDefGradIP[pidx] = Vaango::Util::Identity;
  }
  computeStableTimestep(patch, matl, new_dw);
}

void
Membrane::computeStableTimestep(const Patch* patch,
                                const MPMMaterial* matl,
                                DataWarehouse* new_dw)
{
  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of cSTensor
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

  double mu   = d_modelParam.shear;
  double bulk = d_modelParam.bulk;
  for (int idx : *pset) {
    // Compute wave speed at each particle, store the maximum
    c_dil         = sqrt((bulk + 4. * mu / 3.) * pVolume[idx] / pMass[idx]);
    Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
    waveSpeed     = Max(velMax, waveSpeed);
  }
  waveSpeed       = dx / waveSpeed;
  double delT_new = waveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

void
Membrane::addComputesAndRequires(Task* task,
                                 const MPMMaterial* matl,
                                 const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForHypoExplicit(task, matlset, patches);

  // Other constitutive model and input dependent computes and requires
  Ghost::GhostType gnone = Ghost::None;

  task->requires(Task::OldDW, pDefGradInPlaneLabel, matlset, gnone);
  task->requires(Task::OldDW, pTang1Label, matlset, gnone);
  task->requires(Task::OldDW, pTang2Label, matlset, gnone);
  task->requires(Task::OldDW, pNormLabel, matlset, gnone);

  task->computes(pDefGradInPlaneLabel_preReloc, matlset);
  task->computes(pTang1Label_preReloc, matlset);
  task->computes(pTang2Label_preReloc, matlset);
  task->computes(pNormLabel_preReloc, matlset);
}

void
Membrane::computeStressTensor(const PatchSubset* patches,
                              const MPMMaterial* matl,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw)
{
  double shear = d_modelParam.shear;
  double bulk  = d_modelParam.bulk;
  int matID    = matl->getDWIndex();

  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  for (int pp = 0; pp < patches->size(); pp++) {
    const Patch* patch = patches->get(pp);
    Vector dx          = patch->dCell();

    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    constParticleVariable<double> pMass, pVolume;
    constParticleVariable<Vector> pVelocity, pTang1, pTang2, pNorm;
    constParticleVariable<Matrix3> pDefGrad_old, pDefGrad_new, pVelGrad;
    constParticleVariable<Matrix3> pDefGradIPOld;
    old_dw->get(pMass, lb->pMassLabel, pset);
    old_dw->get(pVolume, lb->pVolumeLabel, pset);
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);
    old_dw->get(pTang1, pTang1Label, pset);
    old_dw->get(pTang2, pTang2Label, pset);
    old_dw->get(pNorm, pNormLabel, pset);
    old_dw->get(pDefGrad_old, lb->pDefGradLabel, pset);
    old_dw->get(pDefGradIPOld, pDefGradInPlaneLabel, pset);
    new_dw->get(pVelGrad, lb->pVelGradLabel_preReloc, pset);

    ParticleVariable<double> pdTdt, p_q;
    ParticleVariable<Vector> T1, T2, T3;
    ParticleVariable<Matrix3> pStress_new, pDefGradIP;
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);
    new_dw->allocateAndPut(T1, pTang1Label_preReloc, pset);
    new_dw->allocateAndPut(T2, pTang2Label_preReloc, pset);
    new_dw->allocateAndPut(T3, pNormLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);
    new_dw->allocateAndPut(pDefGradIP, pDefGradInPlaneLabel_preReloc, pset);

    double strainEnergy = 0.0;
    Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);
    for (int idx : *pset) {

      pdTdt[idx] = 0.0;
      p_q[idx]   = 0.0;

      T1[idx] = pTang1[idx];
      T2[idx] = pTang2[idx];

      Matrix3 UU = Vaango::Util::Identity;
      Matrix3 RR = Vaango::Util::Identity;
      pDefGrad_old[idx].polarDecompositionRMB(UU, RR);
      // std::cout << RR << endl << endl;

      T1[idx] = RR * pTang1[idx];
      T2[idx] = RR * pTang2[idx];

      // T3 = T1 X T2
      T3[idx] = Vector(T1[idx].y() * T2[idx].z() - T1[idx].z() * T2[idx].y(),
                       -(T1[idx].x() * T2[idx].z() - T1[idx].z() * T2[idx].x()),
                       T1[idx].x() * T2[idx].y() - T1[idx].y() * T2[idx].x());

      // The following code is carrying out:
      // Matrix3 Q(Dot(T1[idx],I), Dot(T1[idx],J), Dot(T1[idx],K),
      //          Dot(T2[idx],I), Dot(T2[idx],J), Dot(T2[idx],K),
      //          Dot(T3[idx],I), Dot(T3[idx],J), Dot(T3[idx],K));
      // assuming that I, J and K are the (1,0,0), (0,1,0) and (0,0,1)
      Matrix3 Q(T1[idx].x(), T1[idx].y(), T1[idx].z(),
                T2[idx].x(), T2[idx].y(), T2[idx].z(),
                T3[idx].x(), T3[idx].y(), T3[idx].z());

      Vector vGT1 = pVelGrad[idx] * T1[idx];
      Vector vGT2 = pVelGrad[idx] * T2[idx];
      Vector vGT3 = pVelGrad[idx] * T3[idx];

      Matrix3 L_ij_ip(0.0);
      L_ij_ip(0, 0) = Dot(T1[idx], vGT1);
      L_ij_ip(0, 1) = Dot(T1[idx], vGT2);
      L_ij_ip(0, 2) = Dot(T1[idx], vGT3);
      L_ij_ip(1, 0) = Dot(T2[idx], vGT1);
      L_ij_ip(1, 1) = Dot(T2[idx], vGT2);
      L_ij_ip(1, 2) = Dot(T2[idx], vGT3);
      L_ij_ip(2, 0) = Dot(T3[idx], vGT1);
      L_ij_ip(2, 1) = Dot(T3[idx], vGT2);
      L_ij_ip(2, 2) = Dot(T3[idx], vGT3);

      Matrix3 T1T1, T1T2, T2T1, T2T2;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          T1T1(i, j) = T1[idx][i] * T1[idx][j];
          T1T2(i, j) = T1[idx][i] * T2[idx][j];
          T2T1(i, j) = T2[idx][i] * T1[idx][j];
          T2T2(i, j) = T2[idx][i] * T2[idx][j];
        }
      }

      Matrix3 L_ip = T1T1 * L_ij_ip(0, 0) + T1T2 * L_ij_ip(0, 1) +
                     T2T1 * L_ij_ip(1, 0) + T2T2 * L_ij_ip(1, 1);

      Matrix3 L_local = Q * L_ip * Q.Transpose();

      // BE SURE TO FIX THIS
      // I'm currently setting the tangent and normals back to
      // their original position, and then each timestep rotating
      // them by the total R.  It should be possible to do this
      // incrementally, but the error might be greater.
      T1[idx] = pTang1[idx];
      T2[idx] = pTang2[idx];
      T3[idx] = pNorm[idx];

      Matrix3 pDefGradIPInc = L_local * delT + Vaango::Util::Identity;

      // Update the deformation gradient tensor to its time n+1 value.
      pDefGradIP[idx] = pDefGradIPInc * pDefGradIPOld[idx];

      // Use Newton's method to determine pDefGradIP(3,3)
      Matrix3 F = pDefGradIP[idx];

      double epsilon = 1.e-14;
      double delta   = 1.;
      double f33, f33p, f33m, jv, jvp, jvm, sig33, sig33p, sig33m;
      f33 = 1. / (F(0, 0) * F(1, 1));
      while (std::abs(delta) > epsilon) {
        double detF2 = (F(0, 0) * F(1, 1) - F(1, 0) * F(0, 1));
        jv           = f33 * detF2;
        double FinF  = F(0, 0) * F(0, 0) + F(0, 1) * F(0, 1) +
                      F(1, 0) * F(1, 0) + F(1, 1) * F(1, 1);
        sig33 =
          (shear / (3. * std::pow(jv, 2. / 3.))) * (2. * f33 * f33 - FinF) +
          (.5 * bulk) * (jv - 1. / jv);

        f33p = 1.01 * f33;
        f33m = 0.99 * f33;
        jvp  = f33p * detF2;
        jvm  = f33m * detF2;

        sig33p =
          (shear / (3. * std::pow(jvp, 2. / 3.))) * (2. * f33p * f33p - FinF) +
          (.5 * bulk) * (jvp - 1. / jvp);

        sig33m =
          (shear / (3. * std::pow(jvm, 2. / 3.))) * (2. * f33m * f33m - FinF) +
          (.5 * bulk) * (jvm - 1. / jvm);

        delta = -sig33 / ((sig33p - sig33m) / (f33p - f33m));

        f33 = f33 + delta;
      }

      // get the volumetric part of the deformation
      jv = f33 * (F(0, 0) * F(1, 1) - F(1, 0) * F(0, 1));
      pDefGradIP[idx](2, 2) = f33;

      Matrix3 B_bar_new = pDefGradIP[idx] * pDefGradIP[idx].Transpose() *
                          std::pow(jv, -(2. / 3.));

      double IEl = Vaango::Util::one_third * B_bar_new.Trace();

      // shear is equal to the shear modulus times dev(B_bar)
      Matrix3 shearStress = (B_bar_new - Vaango::Util::Identity * IEl) * shear;

      // get the hydrostatic part of the stress
      double p = 0.5 * bulk * (jv - 1.0 / jv);

      // compute the total stress (volumetric + deviatoric)
      pStress_new[idx] = Vaango::Util::Identity * p + shearStress / jv;
      pStress_new[idx](2, 2) = 0.;

      pStress_new[idx] = Q.Transpose() * pStress_new[idx] * Q;

      // Compute the strain energy for all the particles
      double U = .5 * bulk * (.5 * (jv * jv - 1.0) - log(jv));
      double W = .5 * shear * (B_bar_new.Trace() - 3.0);

      double e = (U + W) * pVolume[idx] / jv;

      strainEnergy += e;

      double c_dil = sqrt((bulk + 4. * shear / 3.) * pVolume[idx] / pMass[idx]);
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
Membrane::addComputesAndRequires(Task*,
                                 const MPMMaterial*,
                                 const PatchSet*,
                                 const bool,
                                 const bool) const
{
}

void
Membrane::allocateCMDataAddRequires(Task* task,
                                    const MPMMaterial* matl,
                                    const PatchSet* patches,
                                    MPMLabel*) const
{
  const MaterialSubset* matlset = matl->thisMaterial();

  // Allocate the variables shared by all constitutive models
  // for the particle convert operation
  // This method is defined in the ConstitutiveModel base class.
  addSharedRForConvertExplicit(task, matlset, patches);

  // Add requires local to this model
  Ghost::GhostType gnone = Ghost::None;
  task->requires(Task::NewDW, pDefGradInPlaneLabel_preReloc, matlset, gnone);
}

void
Membrane::allocateCMDataAdd(DataWarehouse* new_dw,
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
  ParticleVariable<Matrix3> pDefGradIP;
  constParticleVariable<Matrix3> o_pDefGradIP;

  new_dw->allocateTemporary(pDefGradIP, addset);
  new_dw->get(o_pDefGradIP, pDefGradInPlaneLabel_preReloc, delset);

  ParticleSubset::iterator o, n = addset->begin();
  for (o = delset->begin(); o != delset->end(); o++, n++) {
    pDefGradIP[*n] = o_pDefGradIP[*o];
  }

  (*newState)[pDefGradInPlaneLabel] = pDefGradIP.clone();
}

// The "CM" versions use the pressure-volume relationship of the CNH model
double
Membrane::computeRhoMicroCM(double pressure,
                            const double p_ref,
                            const MPMMaterial* matl,
                            double temperature,
                            double rho_guess)
{
  double rho_orig = matl->getInitialDensity();
  // double p_ref=101325.0;
  double bulk = d_modelParam.bulk;

  double p_gauge = pressure - p_ref;
  double rho_cur;

  rho_cur =
    rho_orig * (p_gauge / bulk + sqrt((p_gauge / bulk) * (p_gauge / bulk) + 1));

  return rho_cur;
}

void
Membrane::computePressEOSCM(double rho_cur,
                            double& pressure,
                            double p_ref,
                            double& dp_drho,
                            double& tmp,
                            const MPMMaterial* matl,
                            double temperature)
{
  // double p_ref=101325.0;
  double bulk = d_modelParam.bulk;
  // double shear = d_modelParam.shear;
  double rho_orig = matl->getInitialDensity();

  double p_g = .5 * bulk * (rho_cur / rho_orig - rho_orig / rho_cur);
  pressure   = p_ref + p_g;
  dp_drho    = .5 * bulk * (rho_orig / (rho_cur * rho_cur) + 1. / rho_orig);
  tmp        = bulk / rho_cur; // speed of sound squared
}

double
Membrane::getCompressibility()
{
  return 1.0 / d_modelParam.bulk;
}
