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

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/TransIsoHyper.h>
#include <CCA/Ports/DataWarehouse.h>
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
#include <Core/Math/Short27.h>              //for Fracture
#include <Core/Math/TangentModulusTensor.h> //added this for stiffness
#include <Core/ProblemSpec/ProblemSpec.h>
#include <fstream>
#include <iostream>

using namespace Uintah;

// _________________transversely isotropic hyperelastic material [Jeff Weiss's]

TransIsoHyper::TransIsoHyper(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  d_useModifiedEOS = false;

  ps->require("bulk_modulus", d_modelParam.bulkModulus);
  ps->require("c1", d_modelParam.c1); // Mooney Rivlin constant 1
  ps->require("c2", d_modelParam.c2); // Mooney Rivlin constant 2
  ps->require("c3", d_modelParam.c3); // scales exponential stresses
  ps->require("c4", d_modelParam.c4); // controls uncrimping of fibers
  ps->require("c5", d_modelParam.c5); // straightened fibers modulus
  ps->require("fiber_pStretch", d_modelParam.lambdaStar);
  ps->require("direction_of_symm", d_modelParam.a0);
  ps->require("failure_option",
              d_modelParam.failure); // failure flag True/False
  ps->require("max_fiber_strain", d_modelParam.critStretch);
  ps->require("max_matrix_strain", d_modelParam.critShear);
  ps->get("useModifiedEOS", d_useModifiedEOS); // no negative pressure for
                                               // solids

  pStretchLabel = VarLabel::create(
    "p.pStretch", ParticleVariable<double>::getTypeDescription());
  pStretchLabel_preReloc = VarLabel::create(
    "p.pStretch+", ParticleVariable<double>::getTypeDescription());

  pFailureLabel =
    VarLabel::create("p.fail", ParticleVariable<double>::getTypeDescription());
  pFailureLabel_preReloc =
    VarLabel::create("p.fail+", ParticleVariable<double>::getTypeDescription());
}

TransIsoHyper::TransIsoHyper(const TransIsoHyper* cm)
  : ConstitutiveModel(cm)
{
  d_useModifiedEOS = cm->d_useModifiedEOS;

  d_modelParam.bulkModulus = cm->d_modelParam.bulkModulus;
  d_modelParam.c1 = cm->d_modelParam.c1;
  d_modelParam.c2 = cm->d_modelParam.c2;
  d_modelParam.c3 = cm->d_modelParam.c3;
  d_modelParam.c4 = cm->d_modelParam.c4;
  d_modelParam.c5 = cm->d_modelParam.c5;
  d_modelParam.lambdaStar = cm->d_modelParam.lambdaStar;
  d_modelParam.a0 = cm->d_modelParam.a0;
  d_modelParam.failure = cm->d_modelParam.failure;
  d_modelParam.critStretch = cm->d_modelParam.critStretch;
  d_modelParam.critShear = cm->d_modelParam.critShear;
}

TransIsoHyper::~TransIsoHyper()
{
  VarLabel::destroy(pStretchLabel);
  VarLabel::destroy(pStretchLabel_preReloc);
  VarLabel::destroy(pFailureLabel);
  VarLabel::destroy(pFailureLabel_preReloc);
}

void
TransIsoHyper::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "trans_iso_hyper");
  }

  cm_ps->appendElement("bulk_modulus", d_modelParam.bulkModulus);
  cm_ps->appendElement("c1", d_modelParam.c1);
  cm_ps->appendElement("c2", d_modelParam.c2);
  cm_ps->appendElement("c3", d_modelParam.c3);
  cm_ps->appendElement("c4", d_modelParam.c4);
  cm_ps->appendElement("c5", d_modelParam.c5);
  cm_ps->appendElement("fiber_pStretch", d_modelParam.lambdaStar);
  cm_ps->appendElement("direction_of_symm", d_modelParam.a0);
  cm_ps->appendElement("failure_option", d_modelParam.failure);
  cm_ps->appendElement("max_fiber_strain", d_modelParam.critStretch);
  cm_ps->appendElement("max_matrix_strain", d_modelParam.critShear);
  cm_ps->appendElement("useModifiedEOS", d_useModifiedEOS);
}

TransIsoHyper*
TransIsoHyper::clone()
{
  return scinew TransIsoHyper(*this);
}

Vector
TransIsoHyper::getInitialFiberDir()
{
  return d_modelParam.a0;
}

void
TransIsoHyper::addParticleState(std::vector<const VarLabel*>& from,
                                std::vector<const VarLabel*>& to)
{
  // Add the local particle state data for this constitutive model.
  from.push_back(lb->pFiberDirLabel);
  from.push_back(pStretchLabel);
  from.push_back(pFailureLabel);

  to.push_back(lb->pFiberDirLabel_preReloc);
  to.push_back(pStretchLabel_preReloc);
  to.push_back(pFailureLabel_preReloc);
}

/**
 * Initialiazation:
 * Assumption: STRESS FREE REFERENCE CONFIG
 */
void
TransIsoHyper::addInitialComputesAndRequires(Task* task,
                                             const MPMMaterial* matl,
                                             const PatchSet*) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pFailureLabel, matlset);
  task->computes(pStretchLabel, matlset);
  task->computes(lb->pStressLabel_preReloc, matlset);
}

void
TransIsoHyper::initializeCMData(const Patch* patch, const MPMMaterial* matl,
                                DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);

  Matrix3 Identity, zero(0.);
  Identity.Identity();

  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);
  ParticleVariable<double> pStretch, pFail; // pFail_label
  new_dw->allocateAndPut(pFail, pFailureLabel, pset);
  new_dw->allocateAndPut(pStretch, pStretchLabel, pset);

  for (auto particle : *pset) {
    pFail[particle] = 0.0;
    pStretch[particle] = 1.0;
  }
  computeStableTimestep(patch, matl, new_dw);
}

void
TransIsoHyper::computeStableTimestep(const Patch* patch,
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

  double c_dil = 0.0;
  Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);

  double bulkModulus = d_modelParam.bulkModulus;
  double c1 = d_modelParam.c1;

  for (int idx : *pset) {
    // this is valid only for F=Identity
    c_dil = sqrt((bulkModulus + 2. / 3. * c1) * pVolume[idx] / pMass[idx]);
    Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
    waveSpeed = Max(velMax, waveSpeed);
  }
  waveSpeed = dx / waveSpeed;
  double delT_new = waveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

void
TransIsoHyper::addComputesAndRequires(Task* task, const MPMMaterial* matl,
                                      const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForExplicit(task, matlset, patches);

  Ghost::GhostType gnone = Ghost::None;
  task->requires(Task::OldDW, lb->pFiberDirLabel, matlset, gnone);
  task->requires(Task::OldDW, pFailureLabel, matlset, gnone);

  task->computes(lb->pFiberDirLabel_preReloc, matlset);
  task->computes(pStretchLabel_preReloc, matlset);
  task->computes(pFailureLabel_preReloc, matlset);
}

void
TransIsoHyper::computeStressTensor(const PatchSubset* patches,
                                   const MPMMaterial* matl,
                                   DataWarehouse* old_dw, DataWarehouse* new_dw)
{
  Matrix3 Identity; Identity.Identity();
  Matrix3 Zero(0.0);

  double bulkModulus = d_modelParam.bulkModulus;
  double c1 = d_modelParam.c1;
  double c2 = d_modelParam.c2;
  double c3 = d_modelParam.c3;
  double c4 = d_modelParam.c4;
  double c5 = d_modelParam.c5;
  double lambdaStar = d_modelParam.lambdaStar;
  double c6 = c3 * (std::exp(c4 * (lambdaStar - 1.0)) - 1.0) -
              c5 * lambdaStar; // c6 = y-intercept
  double rho_orig = matl->getInitialDensity();
  double failure = d_modelParam.failure;
  double critShear = d_modelParam.critShear;
  double critStretch = d_modelParam.critStretch;

  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  int matID = matl->getDWIndex();

  for (int pp = 0; pp < patches->size(); pp++) {
    const Patch* patch = patches->get(pp);
    Vector dx = patch->dCell();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    double strainEnergy = 0.;
    Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);

    constParticleVariable<double> pVolume_new, pFail_old;
    constParticleVariable<Vector> pVelocity, pFiberDir;
    constParticleVariable<Matrix3> pDefGrad, pVelGrad, pDefGrad_new;

    old_dw->get(pFail_old, pFailureLabel, pset);
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);
    old_dw->get(pFiberDir, lb->pFiberDirLabel, pset);
    old_dw->get(pDefGrad, lb->pDefGradLabel, pset);

    new_dw->get(pVolume_new, lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pVelGrad, lb->pVelGradLabel_preReloc, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    ParticleVariable<double> pFail, pdTdt, pStretch, p_q;
    ParticleVariable<Vector> pFiberDir_carry;
    ParticleVariable<Matrix3> pStress;

    new_dw->allocateAndPut(pFail, pFailureLabel_preReloc, pset);
    new_dw->allocateAndPut(pStretch, pStretchLabel_preReloc, pset);
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);
    new_dw->allocateAndPut(pFiberDir_carry, lb->pFiberDirLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress, lb->pStressLabel_preReloc, pset);

    for (int idx : *pset) {
      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;

      // get the volumetric part of the deformation
      double J = pDefGrad_new[idx].Determinant();
      double J_onethird = std::pow(J, (1.0/3.0));
      double J_twothird = J_onethird * J_onethird;
      double inv_J_twothird = 1.0 / J_twothird;

      // carry forward fiber direction
      pFiberDir_carry[idx] = pFiberDir[idx];

      //_______________________UNCOUPLE DEVIATORIC AND DILATIONAL PARTS
      //_______________________Fbar=J^(-1/3)*F
      //_______________________Fvol=J^1/3*Identity
      //_______________________right Cauchy Green (C) bar and invariants
      Matrix3 C = pDefGrad_new[idx].Transpose() * pDefGrad_new[idx];
      Matrix3 B = pDefGrad_new[idx] * pDefGrad_new[idx].Transpose();
      Matrix3 C_bar = C * inv_J_twothird;
      Matrix3 B_bar = B * inv_J_twothird;

      double I1_bar = C_bar.Trace();
      double I2_bar = .5 * (I1_bar * I1_bar - (C_bar * C_bar).Trace());
      double I4_bar = Dot(pFiberDir[idx], (C_bar * pFiberDir[idx]));
      double lambda_bar = std::sqrt(I4_bar);

      double I4 = I4_bar * J_twothird; // For diagnostics only
      pStretch[idx] = std::sqrt(I4);

      Vector fiberDir_new = pDefGrad_new[idx] * pFiberDir[idx] *
                              (1. / (lambda_bar * J_onethird));
      Matrix3 ff_dyad(fiberDir_new, fiberDir_new);

      //________________________________strain energy derivatives
      double dWdI4_bar = 0.;
      double d2WdI4_bar2 = 0.;
      double shearModulus = 2.0 * c1 + c2;
      if (lambda_bar >= 1.) {
        double lam_bar_sq = I4_bar;
        double lam_bar_cub = I4_bar * lambda_bar;
        if (lambda_bar < lambdaStar) {
          double dW_dlambda_bar = c3 * c4 * std::exp(c4 * (lambda_bar - 1));
          double dlambda_bar_dI4_bar = 1.0 /(2.0 * lambda_bar);
          double d2lambda_bar_dI4_bar2 = - 1.0 /(4.0 * lam_bar_cub);
          dWdI4_bar = dW_dlambda_bar * dlambda_bar_dI4_bar;
          d2WdI4_bar2 = dW_dlambda_bar * d2lambda_bar_dI4_bar2;
          
          //dWdI4_bar =
          //  0.5 * c3 * (std::exp(c4 * (lambda_bar - 1.)) - 1.) / lam_bar_sq;
          //d2WdI4_bar2 = 0.25 * c3 * (c4 * std::exp(c4 * (lambda_bar - 1.)) -
          //                            1. / lambda_bar *
          //                              (std::exp(c4 * (lambda_bar - 1.)) - 1.)) /
          //               lam_bar_cub;

        } else {
          double lam_bar_4th = I4_bar * I4_bar;
          dWdI4_bar = 0.5 * (c5 + c6 / lambda_bar) / lambda_bar;
          d2WdI4_bar2 = -0.25 * c6 / lam_bar_4th;
        }
        shearModulus += I4_bar * (4. * d2WdI4_bar2 * lam_bar_sq -
                           2. * dWdI4_bar * lambda_bar);
      }

      // Compute deformed volume and local wave speed
      double rho_cur = rho_orig / J;
      double c_dil = sqrt((bulkModulus + 1. / 3. * shearModulus) / rho_cur);
      double p = bulkModulus * log(J) / J; // p -= qVisco;
      if (p >= -1.e-5 && p <= 1.e-5) {
        p = 0.;
      }

      //________________________________Failure and stress terms
      Matrix3 deviatoric_stress(0.0), fiber_stress(0.0);
      pFail[idx] = 0.;
      if (failure != 1) {
        deviatoric_stress =
          (B_bar * (c1 + c2 * I1_bar) - B_bar * B_bar * c2 -
           Identity * (1. / 3.) * (c1 * I1_bar + 2. * c2 * I2_bar)) *
          2. / J;
        fiber_stress =
          (ff_dyad - Identity * (1. / 3.)) * (2. / J) * dWdI4_bar * I4_bar;
      } else {
        double matrix_failed = 0.;
        double fiber_failed = 0.;
        Vector eigVal = computeEigenvalues(C);
        double e1 = eigVal[0];
        double e3 = eigVal[2];
        double max_shear_strain = (e1 - e3) / 2.;
        if (max_shear_strain > critShear || pFail_old[idx] == 1.0 ||
            pFail_old[idx] == 3.0) {
          deviatoric_stress = Zero;
          pFail[idx] = 1.;
          matrix_failed = 1.;
        } else {
          deviatoric_stress =
            (B_bar * (c1 + c2 * I1_bar) - B_bar * B_bar * c2 -
             Identity * (1. / 3.) * (c1 * I1_bar + 2. * c2 * I2_bar)) *
            2. / J;
        }
        //________________________________fiber stress term + failure of fibers
        if (pStretch[idx] > critStretch || pFail_old[idx] == 2. ||
            pFail_old[idx] == 3.) {
          fiber_stress = Zero;
          pFail[idx] = 2.;
          fiber_failed = 1.;
        } else {
          fiber_stress = (ff_dyad * dWdI4_bar * I4_bar -
                          Identity * (1. / 3.) * dWdI4_bar * I4_bar) *
                         2. / J;
        }
        if ((matrix_failed + fiber_failed) == 2. || pFail_old[idx] == 3.) {
          pFail[idx] = 3.;
        }
        //________________________________hydrostatic pressure term
        if (pFail[idx] == 1.0 || pFail[idx] == 3.0) {
          p = 0.;
        }
        //_______________________________Cauchy stress
      }

      Matrix3 pressure = Identity * p;
      // Cauchy stress
      pStress[idx] = pressure + deviatoric_stress + fiber_stress;
      //________________________________end stress

      // Compute the strain energy for all the particles
      double U = .5 * std::log(J) * std::log(J) * bulkModulus;
      double W = 0.0;
      if (lambda_bar < lambdaStar) {
        W = c1 * (I1_bar - 3.) + c2 * (I2_bar - 3.) +
            (std::exp(c4 * (lambda_bar - 1.)) - 1.) * c3;
      } else {
        W = c1 * (I1_bar - 3.) + c2 * (I2_bar - 3.) + c5 * lambda_bar +
            c6 * std::log(lambda_bar);
      }

      double e = (U + W) * pVolume_new[idx] / J;
      strainEnergy += e;

      // Compute artificial viscosity term
      if (flag->d_artificialViscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;
        double c_bulk = sqrt(bulkModulus / rho_cur);
        Matrix3 D = (pVelGrad[idx] + pVelGrad[idx].Transpose()) * 0.5;
        p_q[idx] = artificialBulkViscosity(D.Trace(), c_bulk, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.;
      }

      Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
      waveSpeed = Max(velMax, waveSpeed);
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

Vector 
TransIsoHyper::computeEigenvalues(const Matrix3& C) const
{
  double e1, e2, e3; // eigenvalues of C=symm.+pos.def.->Dis<=0
  double pi = 3.1415926535897932384;
  double I1 = C.Trace();
  double I2 = .5 * (I1 * I1 - (C * C).Trace());
  double I3 = C.Determinant();
  double Q = (1. / 9.) * (3. * I2 - I1 * I1);
  double R =
    (1. / 54.) * (-9. * I1 * I2 + 27. * I3 + 2. * (I1 * I1 * I1));
  double Dis = Q * Q * Q + R * R;
  if (Dis <= 1.e-5 && Dis >= 0.) {
    if (R >= -1.e-5 && R <= 1.e-5) {
      e1 = e2 = e3 = I1 / 3.;
    } else {
      e1 = 2. * std::pow(R, 1. / 3.) + I1 / 3.;
      e3 = -std::pow(R, 1. / 3.) + I1 / 3.;
      if (e1 < e3)
        std::swap(e1, e3);
      e2 = e3;
    }
  } else {
    double theta = std::acos(R / std::pow(-Q, 3. / 2.));
    double sqrt_negQ = std::sqrt(-Q);
    e1 = 2. * sqrt_negQ * std::cos(theta / 3.) + I1 / 3.;
    e2 = 2. * sqrt_negQ * std::cos(theta / 3. + 2. * pi / 3.) + I1 / 3.;
    e3 = 2. * sqrt_negQ * std::cos(theta / 3. + 4. * pi / 3.) + I1 / 3.;
    if (e1 < e2)
      std::swap(e1, e2);
    if (e1 < e3)
      std::swap(e1, e3);
    if (e2 < e3)
      std::swap(e2, e3);
  }
  return Vector(e1, e2, e3);
}

void
TransIsoHyper::addComputesAndRequires(Task*, const MPMMaterial*,
                                      const PatchSet*, const bool,
                                      const bool) const
{
}

/**
  * Used with RigidMPM
  */
void
TransIsoHyper::carryForward(const PatchSubset* patches, const MPMMaterial* matl,
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
    constParticleVariable<Vector> pFiberDir;
    ParticleVariable<double> pStretch;
    constParticleVariable<Vector> pFail_old;
    ParticleVariable<double> pFail;

    ParticleVariable<Vector> pFiberDir_new;
    old_dw->get(pFiberDir, lb->pFiberDirLabel, pset);
    old_dw->get(pFail_old, pFailureLabel, pset);

    new_dw->allocateAndPut(pFiberDir_new, lb->pFiberDirLabel_preReloc, pset);
    new_dw->allocateAndPut(pStretch, pStretchLabel_preReloc, pset);
    new_dw->allocateAndPut(pFail, pFailureLabel_preReloc, pset);

    for (int idx : *pset) {
      pFiberDir_new[idx] = pFiberDir[idx];
      pStretch[idx] = 1.0;
      pFail[idx] = 0.0;
    }
    new_dw->put(delt_vartype(1.e10), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(0.), lb->StrainEnergyLabel);
    }
  }
}

void
TransIsoHyper::allocateCMDataAddRequires(Task* task, const MPMMaterial* matl,
                                         const PatchSet* patches,
                                         MPMLabel*) const
{
  const MaterialSubset* matlset = matl->thisMaterial();

  // Allocate the variables shared by all constitutive models
  // for the particle convert operation
  // This method is defined in the ConstitutiveModel base class.
  addSharedRForConvertExplicit(task, matlset, patches);
  // Add requires local to this model
  task->requires(Task::NewDW, pFailureLabel_preReloc, matlset, Ghost::None);
  task->requires(Task::NewDW, pStretchLabel_preReloc, matlset, Ghost::None);
}

void
TransIsoHyper::allocateCMDataAdd(DataWarehouse* new_dw, ParticleSubset* addset,
                                 ParticleLabelVariableMap* newState,
                                 ParticleSubset* delset, DataWarehouse*)
{
  // Copy the data common to all constitutive models from the particle to be
  // deleted to the particle to be added.
  // This method is defined in the ConstitutiveModel base class.
  copyDelToAddSetForConvertExplicit(new_dw, delset, addset, newState);

  // Copy the data local to this constitutive model from the particles to
  // be deleted to the particles to be added
  // ParticleVariable<Matrix3> pDefGrad, pStress;
  // constParticleVariable<Matrix3> o_defGrad, o_stress;
  ParticleVariable<double> pStretch, pFail;
  constParticleVariable<double> o_pStretch, o_pFail;

  // new_dw->allocateTemporary(pDefGrad,addset);
  // new_dw->allocateTemporary(pStress,            addset);
  new_dw->allocateTemporary(pStretch, addset);
  new_dw->allocateTemporary(pFail, addset);

  // new_dw->get(o_defGrad,lb->pDeformationMeasureLabel_preReloc,   delset);
  // new_dw->get(o_stress, lb->pStressLabel_preReloc,               delset);
  new_dw->get(o_pStretch, pStretchLabel_preReloc, delset);
  new_dw->get(o_pFail, pFailureLabel_preReloc, delset);

  ParticleSubset::iterator o, n = addset->begin();
  for (o = delset->begin(); o != delset->end(); o++, n++) {
    // pDefGrad[*n] = o_defGrad[*o];
    // pStress[*n] = o_stress[*o];
    pStretch[*n] = o_pStretch[*o];
    pFail[*n] = o_pFail[*o];
  }
  //(*newState)[lb->pDeformationMeasureLabel]=pDefGrad.clone();
  //(*newState)[lb->pStressLabel]=pStress.clone();
  (*newState)[pStretchLabel] = pStretch.clone();
  (*newState)[pFailureLabel] = pFail.clone();
}

// The "CM" versions use the pressure-volume relationship of the CNH model
double
TransIsoHyper::computeRhoMicroCM(double pressure, const double p_ref,
                                 const MPMMaterial* matl, double temperature,
                                 double rho_guess)
{
  double rho_orig = matl->getInitialDensity();
  double bulkModulus = d_modelParam.bulkModulus;

  double p_gauge = pressure - p_ref;
  double rho_cur;

  if (d_useModifiedEOS && p_gauge < 0.0) {
    double A = p_ref; // MODIFIED EOS
    double n = p_ref / bulkModulus;
    rho_cur = rho_orig * pow(pressure / A, n);
  } else { // STANDARD EOS
    rho_cur = rho_orig *
              (p_gauge / bulkModulus + sqrt((p_gauge / bulkModulus) * (p_gauge / bulkModulus) + 1));
  }
  return rho_cur;
}

void
TransIsoHyper::computePressEOSCM(const double rho_cur, double& pressure,
                                 const double p_ref, double& dp_drho,
                                 double& tmp, const MPMMaterial* matl,
                                 double temperature)
{
  double bulkModulus = d_modelParam.bulkModulus;
  double rho_orig = matl->getInitialDensity();

  if (d_useModifiedEOS && rho_cur < rho_orig) {
    double A = p_ref; // MODIFIED EOS
    double n = bulkModulus / p_ref;
    pressure = A * pow(rho_cur / rho_orig, n);
    dp_drho = (bulkModulus / rho_orig) * pow(rho_cur / rho_orig, n - 1);
    tmp = dp_drho; // speed of sound squared
  } else {         // STANDARD EOS
    double p_g = .5 * bulkModulus * (rho_cur / rho_orig - rho_orig / rho_cur);
    pressure = p_ref + p_g;
    dp_drho = .5 * bulkModulus * (rho_orig / (rho_cur * rho_cur) + 1. / rho_orig);
    tmp = bulkModulus / rho_cur; // speed of sound squared
  }
}

double
TransIsoHyper::getCompressibility()
{
  return 1.0 / d_modelParam.bulkModulus;
}

