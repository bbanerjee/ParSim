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

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/MPMFlags.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Labels/MPMLabel.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/FastMatrix.h>
#include <Core/Math/Matrix3.h>
#include <cmath>
#include <iostream>

using namespace Uintah;
using std::map;
using std::ostringstream;

#ifndef M_PI
#define M_PI 3.14159265358979323846 /* pi */
#endif

std::ostream&
Uintah::operator<<(std::ostream& out, const ConstitutiveModel::ModelType& mt)
{
  static std::array<const char*, 3> names = {
    { "TOTAL_FORM", "RATE_FORM", "INCREMENTAL" }
  };
  out << names[static_cast<int>(mt)];
  return out;
}

ConstitutiveModel::ConstitutiveModel(MPMFlags* Mflag)
{
  lb = scinew MPMLabel();
  flag = Mflag;
  if (flag->d_8or27 == 8) {
    NGN = 1;
  } else {
    NGN = 2;
  }
}

ConstitutiveModel::ConstitutiveModel(const ConstitutiveModel* cm)
{
  lb = scinew MPMLabel();
  flag = cm->flag;
  NGN = cm->NGN;
  NGP = cm->NGP;
  d_sharedState = cm->d_sharedState;
}

ConstitutiveModel::~ConstitutiveModel()
{
  if (lb) {
    delete lb;
  }
}

void
ConstitutiveModel::addInitialComputesAndRequires(Task*, const MPMMaterial*,
                                                 const PatchSet*) const
{
  // Do nothing unless this function has been implemented by the
  // constitutive model
}

///////////////////////////////////////////////////////////////////////
/*! Initialize the common quantities that all the explicit constituive
 *  models compute */
///////////////////////////////////////////////////////////////////////
void
ConstitutiveModel::initSharedDataForExplicit(const Patch* patch,
                                             const MPMMaterial* matl,
                                             DataWarehouse* new_dw)
{
  Matrix3 I;
  I.Identity();
  Matrix3 zero(0.);
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  ParticleVariable<double> pdTdt;
  ParticleVariable<Matrix3> pStress;

  ParticleVariable<Matrix3> pDefGrad;
  new_dw->getModifiable(pDefGrad, lb->pDefGradLabel, pset);

  new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel, pset);
  new_dw->allocateAndPut(pStress, lb->pStressLabel, pset);

  // To fix : For a material that is initially stressed we need to
  // modify the stress tensors to comply with the initial stress state
  ParticleSubset::iterator iter = pset->begin();
  for (; iter != pset->end(); iter++) {
    particleIndex idx = *iter;
    pdTdt[idx] = 0.0;
    pDefGrad[idx] = I;
    pStress[idx] = zero;
  }
}

///////////////////////////////////////////////////////////////////////
/*!
* Actually initialize the stress and deformation gradient assuming linear
* elastic behavior after computing the body force acceleration
*
* **WARNING** 1) Assumes zero shear stresses and that body forces are aligned
*                with coordinate directions
*             2) Needs the model to have a "initializeWithBodyForce" flag
*                set as true.  A more general implementation is not worth
*                the significant extra effort.
*/
///////////////////////////////////////////////////////////////////////
void
ConstitutiveModel::initializeStressAndDefGradFromBodyForce(const Patch*,
                                                           const MPMMaterial*,
                                                           DataWarehouse*) const
{
  // Do nothing unless this function has been implemented by the
  // constitutive model
}

void
ConstitutiveModel::addComputesAndRequires(Task*, const MPMMaterial*,
                                          const PatchSet*) const
{
  throw InternalError("Stub Task: ConstitutiveModel::addComputesAndRequires ",
                      __FILE__, __LINE__);
}

void
ConstitutiveModel::addComputesAndRequires(Task*, const MPMMaterial*,
                                          const PatchSet*, const bool,
                                          const bool) const
{
  throw InternalError("Stub Task: ConstitutiveModel::addComputesAndRequires ",
                      __FILE__, __LINE__);
}

void
ConstitutiveModel::scheduleCheckNeedAddMPMMaterial(Task* task,
                                                   const MPMMaterial*,
                                                   const PatchSet*) const
{
  task->computes(lb->NeedAddMPMMaterialLabel);
}

void
ConstitutiveModel::addSharedCRForHypoExplicit(Task* task,
                                              const MaterialSubset* matlset,
                                              const PatchSet* p) const
{
  Ghost::GhostType gnone = Ghost::None;
  addSharedCRForExplicit(task, matlset, p);
  task->requires(Task::OldDW, lb->pStressLabel, matlset, gnone);
}

void
ConstitutiveModel::addSharedCRForExplicit(Task* task,
                                          const MaterialSubset* matlset,
                                          const PatchSet*) const
{
  Ghost::GhostType gnone = Ghost::None;
  Ghost::GhostType gac = Ghost::AroundCells;

  task->requires(Task::OldDW, lb->delTLabel);
  task->requires(Task::OldDW, lb->pXLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pMassLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pVolumeLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pTemperatureLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pVelocityLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pAccelerationLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pDefGradLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pVelGradLabel, matlset, gnone);
  task->requires(Task::NewDW, lb->gVelocityStarLabel, matlset, gac, NGN);
  if (!flag->d_doGridReset) {
    task->requires(Task::NewDW, lb->gDisplacementLabel, matlset, gac, NGN);
  }
  task->requires(Task::OldDW, lb->pSizeLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pTempPreviousLabel, matlset, gnone);
  if (flag->d_fracture) {
    task->requires(Task::NewDW, lb->pgCodeLabel, matlset, gnone);
    task->requires(Task::NewDW, lb->GVelocityStarLabel, matlset, gac, NGN);
  }

  task->modifies(lb->pDefGradLabel_preReloc, matlset);
  task->modifies(lb->pVelGradLabel_preReloc, matlset);
  task->computes(lb->pStressLabel_preReloc, matlset);
  task->modifies(lb->pVolumeLabel_preReloc, matlset);
  task->computes(lb->pdTdtLabel_preReloc, matlset);
}

void
ConstitutiveModel::addComputesAndRequiresForRotatedExplicit(Task* task,
                                          const MaterialSubset* matlset,
                                          const PatchSet* patches) const
{
  Ghost::GhostType gnone = Ghost::None;
  task->requires(Task::OldDW, lb->delTLabel);
  task->requires(Task::OldDW, lb->pParticleIDLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pXLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pMassLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pVolumeLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pTempPreviousLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pTemperatureLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pVelocityLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pAccelerationLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pSizeLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pDefGradLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pVelGradLabel, matlset, gnone);

  task->requires(Task::NewDW, lb->pVolumeLabel_preReloc, matlset, Ghost::None);
  task->requires(Task::NewDW, lb->pDeformRateMidLabel, matlset, Ghost::None);
  task->requires(Task::NewDW, lb->pStressUnrotatedLabel, matlset, Ghost::None);
  task->requires(Task::NewDW, lb->pDefGradLabel_preReloc, matlset, Ghost::None);

  task->computes(lb->pStressLabel_preReloc, matlset);
  task->computes(lb->pdTdtLabel_preReloc, matlset);
}

void
ConstitutiveModel::computeStressTensor(const PatchSubset*, const MPMMaterial*,
                                       DataWarehouse*, DataWarehouse*)
{
  throw InternalError("Stub Task: ConstitutiveModel::computeStressTensor ",
                      __FILE__, __LINE__);
}

void
ConstitutiveModel::computeStressTensorImplicit(const PatchSubset*,
                                               const MPMMaterial*,
                                               DataWarehouse*, DataWarehouse*)
{
  throw InternalError(
    "Stub Task: ConstitutiveModel::computeStressTensorImplicit ", __FILE__,
    __LINE__);
}

void
ConstitutiveModel::checkNeedAddMPMMaterial(const PatchSubset*,
                                           const MPMMaterial*,
                                           DataWarehouse* new_dw,
                                           DataWarehouse*)
{
  double need_add = 0.;

  new_dw->put(sum_vartype(need_add), lb->NeedAddMPMMaterialLabel);
}

void
ConstitutiveModel::carryForward(const PatchSubset*, const MPMMaterial*,
                                DataWarehouse*, DataWarehouse*)
{
  throw InternalError("Stub Task: ConstitutiveModel::carryForward ", __FILE__,
                      __LINE__);
}

void
ConstitutiveModel::carryForwardSharedData(ParticleSubset* pset,
                                          DataWarehouse* old_dw,
                                          DataWarehouse* new_dw,
                                          const MPMMaterial* matl)
{
  double rho_orig = matl->getInitialDensity();
  Matrix3 Id, Zero(0);
  Id.Identity();

  constParticleVariable<double> pMass;
  constParticleVariable<Matrix3> pDefGrad_old;
  old_dw->get(pMass, lb->pMassLabel, pset);
  old_dw->get(pDefGrad_old, lb->pDefGradLabel, pset);

  ParticleVariable<double> pVol_new;
  ParticleVariable<Matrix3> pDefGrad_new;
  new_dw->getModifiable(pVol_new, lb->pVolumeLabel_preReloc, pset);
  new_dw->getModifiable(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

  ParticleVariable<double> pIntHeatRate_new, p_q;
  ParticleVariable<Matrix3> pStress_new;
  new_dw->allocateAndPut(pIntHeatRate_new, lb->pdTdtLabel_preReloc, pset);
  new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);
  new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);

  ParticleSubset::iterator iter = pset->begin();
  for (; iter != pset->end(); iter++) {
    particleIndex idx = *iter;
    pVol_new[idx] = (pMass[idx] / rho_orig);
    pIntHeatRate_new[idx] = 0.0;
    pDefGrad_new[idx] = pDefGrad_old[idx];
    pStress_new[idx] = Zero;
    p_q[idx] = 0.;
  }
}

void
ConstitutiveModel::allocateCMDataAddRequires(Task*, const MPMMaterial*,
                                             const PatchSet*, MPMLabel*) const
{
  throw InternalError(
    "Stub Task: ConstitutiveModel::allocateCMDataAddRequires ", __FILE__,
    __LINE__);
}

void
ConstitutiveModel::addSharedRForConvertExplicit(Task* task,
                                                const MaterialSubset* mset,
                                                const PatchSet*) const
{
  Ghost::GhostType gnone = Ghost::None;
  task->requires(Task::NewDW, lb->pdTdtLabel_preReloc, mset, gnone);
  // task->requires(Task::NewDW,lb->pDeformationMeasureLabel_preReloc,mset,gnone);
  task->requires(Task::NewDW, lb->pStressLabel_preReloc, mset, gnone);
}

void
ConstitutiveModel::copyDelToAddSetForConvertExplicit(
  DataWarehouse* new_dw, ParticleSubset* delset, ParticleSubset* addset,
  ParticleLabelVariableMap* newState)
{
  constParticleVariable<double> pIntHeatRate_del;
  // constParticleVariable<Matrix3> pDefGrad_del;
  constParticleVariable<Matrix3> pStress_del;

  new_dw->get(pIntHeatRate_del, lb->pdTdtLabel_preReloc, delset);
  // new_dw->get(pDefGrad_del,     lb->pDeformationMeasureLabel_preReloc,
  // delset);
  new_dw->get(pStress_del, lb->pStressLabel_preReloc, delset);

  ParticleVariable<double> pIntHeatRate_add;
  // ParticleVariable<Matrix3> pDefGrad_add;
  ParticleVariable<Matrix3> pStress_add;

  new_dw->allocateTemporary(pIntHeatRate_add, addset);
  // new_dw->allocateTemporary(pDefGrad_add,     addset);
  new_dw->allocateTemporary(pStress_add, addset);

  ParticleSubset::iterator del = delset->begin();
  ParticleSubset::iterator add = addset->begin();
  for (; del != delset->end(); del++, add++) {
    pIntHeatRate_add[*add] = pIntHeatRate_del[*del];
    // pDefGrad_add[*add] = pDefGrad_del[*del];
    pStress_add[*add] = pStress_del[*del];
  }

  (*newState)[lb->pdTdtLabel] = pIntHeatRate_add.clone();
  //(*newState)[lb->pDeformationMeasureLabel] = pDefGrad_add.clone();
  (*newState)[lb->pStressLabel] = pStress_add.clone();
}

void
ConstitutiveModel::addRequiresDamageParameter(Task*, const MPMMaterial*,
                                              const PatchSet*) const
{
}

void
ConstitutiveModel::getDamageParameter(const Patch*, ParticleVariable<int>&, int,
                                      DataWarehouse*, DataWarehouse*)
{
}

double 
ConstitutiveModel::computeRateOfWork(const Matrix3& stress, const Matrix3& rateOfDeformation) const
{
  double rate = rateOfDeformation(0, 0) * stress(0, 0) + 
                rateOfDeformation(1, 1) * stress(1, 1) +
                rateOfDeformation(2, 2) * stress(2, 2) +
                2. * (rateOfDeformation(0, 1) * stress(0, 1) + 
                      rateOfDeformation(0, 2) * stress(0, 2) +
                      rateOfDeformation(1, 2) * stress(1, 2));
  return rate;
}

Vector
ConstitutiveModel::getInitialFiberDir()
{
  return Vector(0., 0., 1);
}

//______________________________________________________________________
//______________________________________________________________________
//          HARDWIRE FOR AN IDEAL GAS -Todd
double
ConstitutiveModel::computeRhoMicro(double press, double gamma, double cv,
                                   double Temp, double rho_guess)
{
  // Pointwise computation of microscopic density
  return press / ((gamma - 1.0) * cv * Temp);
}

void
ConstitutiveModel::computePressEOS(double rhoM, double gamma, double cv,
                                   double Temp, double& press, double& dp_drho,
                                   double& dp_de)
{
  // Pointwise computation of thermodynamic quantities
  press = (gamma - 1.0) * rhoM * cv * Temp;
  dp_drho = (gamma - 1.0) * cv * Temp;
  dp_de = (gamma - 1.0) * rhoM;
}
//______________________________________________________________________

// Convert J-integral into stress intensity (for FRACTURE)
void
ConstitutiveModel::convertJToK(const MPMMaterial*, const string&, const Vector&,
                               const double&, const Vector&, Vector& SIF)
{
  SIF = Vector(-9999., -9999., -9999.);
}

// Detect if crack propagtes and the propagation direction (for FRACTURE)
short
ConstitutiveModel::crackPropagates(const double&, const double&, const double&,
                                   double& theta)
{
  enum
  {
    NO = 0,
    YES
  };
  theta = 0.0;
  return NO;
}

double
ConstitutiveModel::artificialBulkViscosity(double Dkk, double c_bulk,
                                           double rho, double dx) const
{
  double q = 0.0;
  if (Dkk < 0.0) {
    double A1 = flag->d_artificialViscCoeff1;
    double A2 = flag->d_artificialViscCoeff2;
    // double c_bulk = sqrt(K/rho);
    q = (A1 * fabs(c_bulk * Dkk * dx) + A2 * (Dkk * Dkk * dx * dx)) * rho;
  }
  return q;
}

void 
ConstitutiveModel::computeGradAndBmats(Matrix3& grad, vector<IntVector>& ni,
                                std::vector<Vector>& d_S, const double* oodx,
                                constNCVariable<Vector>& gVec,
                                const Array3<int>& l2g, double B[6][24],
                                double Bnl[3][24], int* dof)
{
  int l2g_node_num = -1;

  // Compute gradient matrix
  grad.set(0.0);
  for (int k = 0; k < flag->d_8or27; k++) {
    const Vector& vec = gVec[ni[k]];
    for (int j = 0; j < 3; j++) {
      double fac = d_S[k][j] * oodx[j];
      for (int i = 0; i < 3; i++) {
        grad(i, j) += vec[i] * fac;
      }
    }
  }

  for (int k = 0; k < 8; k++) {
    B[0][3 * k] = d_S[k][0] * oodx[0];
    B[3][3 * k] = d_S[k][1] * oodx[1];
    B[5][3 * k] = d_S[k][2] * oodx[2];
    B[1][3 * k] = 0.;
    B[2][3 * k] = 0.;
    B[4][3 * k] = 0.;

    B[1][3 * k + 1] = d_S[k][1] * oodx[1];
    B[3][3 * k + 1] = d_S[k][0] * oodx[0];
    B[4][3 * k + 1] = d_S[k][2] * oodx[2];
    B[0][3 * k + 1] = 0.;
    B[2][3 * k + 1] = 0.;
    B[5][3 * k + 1] = 0.;

    B[2][3 * k + 2] = d_S[k][2] * oodx[2];
    B[4][3 * k + 2] = d_S[k][1] * oodx[1];
    B[5][3 * k + 2] = d_S[k][0] * oodx[0];
    B[0][3 * k + 2] = 0.;
    B[1][3 * k + 2] = 0.;
    B[3][3 * k + 2] = 0.;

    Bnl[0][3 * k] = d_S[k][0] * oodx[0];
    Bnl[1][3 * k] = 0.;
    Bnl[2][3 * k] = 0.;
    Bnl[0][3 * k + 1] = 0.;
    Bnl[1][3 * k + 1] = d_S[k][1] * oodx[1];
    Bnl[2][3 * k + 1] = 0.;
    Bnl[0][3 * k + 2] = 0.;
    Bnl[1][3 * k + 2] = 0.;
    Bnl[2][3 * k + 2] = d_S[k][2] * oodx[2];

    // Need to loop over the neighboring patches l2g to get the right
    // dof number.
    l2g_node_num = l2g[ni[k]];
    dof[3 * k] = l2g_node_num;
    dof[3 * k + 1] = l2g_node_num + 1;
    dof[3 * k + 2] = l2g_node_num + 2;
  }
}

