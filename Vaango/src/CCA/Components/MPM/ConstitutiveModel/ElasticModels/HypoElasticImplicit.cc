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

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/HypoElasticImplicit.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/Constants.h>
#include <CCA/Components/MPM/MPMUtils.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Labels/MPMLabel.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/FastMatrix.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <fstream>
#include <iostream>

using std::cerr;
using namespace Uintah;

HypoElasticImplicit::HypoElasticImplicit(ProblemSpecP& ps, MPMFlags* Mflag)
  : HypoElastic(ps, Mflag)
  , ImplicitCM()
{
}

HypoElasticImplicit::HypoElasticImplicit(const HypoElasticImplicit* cm)
  : HypoElastic(cm)
  , ImplicitCM(cm)
{
}

std::unique_ptr<ConstitutiveModel>
HypoElasticImplicit::clone()
{
  return std::make_unique<HypoElasticImplicit>(*this);
}

void 
HypoElasticImplicit::addParticleState(std::vector<const VarLabel*>& from,
                                      std::vector<const VarLabel*>& to)
{
}

void
HypoElasticImplicit::initializeCMData(const Patch* patch,
                                      const MPMMaterial* matl,
                                      DataWarehouse* new_dw)
{
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);
  ParticleVariable<Matrix3> pStress;
  new_dw->allocateAndPut(pStress, lb->pStressLabel, pset);

  for (auto idx : *pset) {
    pStress[idx] = Vaango::Util::Zero;
  }
}

void 
HypoElasticImplicit::computeStableTimestep(const Patch* patch,
                                           const MPMMaterial* matl,
                                           DataWarehouse* new_dw)
{
}

void
HypoElasticImplicit::addComputesAndRequires(Task* task, const MPMMaterial* matl,
                                            const PatchSet*,
                                            const bool /*recurse*/,
                                            const bool SchedParent) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  bool reset = flag->d_doGridReset;
  addSharedCRForImplicitHypo(task, matlset, reset, true, SchedParent);
}

void
HypoElasticImplicit::computeStressTensorImplicit(const PatchSubset* patches,
                                                 const MPMMaterial* matl,
                                                 DataWarehouse* old_dw,
                                                 DataWarehouse* new_dw,
                                                 Solver* solver, const bool)

{
  double G = d_modelParam.G;
  double K = d_modelParam.K;
  double rho_orig = matl->getInitialDensity();

  int matID = matl->getDWIndex();
  DataWarehouse* parent_old_dw =
    new_dw->getOtherDataWarehouse(Task::ParentOldDW);

  for (int pp = 0; pp < patches->size(); pp++) {
    const Patch* patch = patches->get(pp);
    ParticleSubset* pset = parent_old_dw->getParticleSubset(matID, patch);

    Vector dx = patch->dCell();
    double oodx[3] = { 1. / dx.x(), 1. / dx.y(), 1. / dx.z() };

    IntVector lowIndex  = patch->getNodeLowIndex();
    IntVector highIndex = patch->getNodeHighIndex() + IntVector(1, 1, 1);
    Array3<int> l2g(lowIndex, highIndex);
    solver->copyL2G(l2g, patch);

    auto interpolator = flag->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());

    constParticleVariable<Point> pX;
    constParticleVariable<double> pMass, pVolume_new;
    constParticleVariable<Matrix3> pSize, pDefGrad_new, 
                                   pDispGrad, pStress_old;
    parent_old_dw->get(pX, lb->pXLabel, pset);
    parent_old_dw->get(pMass, lb->pMassLabel, pset);
    parent_old_dw->get(pSize, lb->pSizeLabel, pset);
    parent_old_dw->get(pStress_old, lb->pStressLabel, pset);

    new_dw->get(pVolume_new, lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pDispGrad, lb->pDispGradLabel_preReloc, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    ParticleVariable<double> pdTdt;
    ParticleVariable<Matrix3> pStress_new;
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);

    double B[6][24];
    double Bnl[3][24];
    double K_vec[576];
    int dof[24];
    double D[6][6];
    double K_mat[24][24];
    double K_geo[24][24];
    double sig[3][3];

    if (matl->getIsRigid()) {
      for (int idx : *pset) {
        pdTdt[idx] = 0.;
        pStress_new[idx] = Vaango::Util::Zero;
      }
    } else {
      for (int idx : *pset) {
        pdTdt[idx] = 0.;

        interpolator->findCellAndShapeDerivatives(pX[idx], ni, d_S,
                                                  pSize[idx], pDefGrad_new[idx]);
        loadBMats(l2g, dof, B, Bnl, d_S, ni, oodx);

        // Calculate the strain and deviatoric rate 
        Matrix3 e = (pDispGrad[idx] + pDispGrad[idx].Transpose()) * .5;
        Matrix3 ePrime = e - Vaango::Util::Identity * Vaango::Util::one_third * e.Trace();

        pStress_new[idx] = pStress_old[idx] + 
          (ePrime * 2. * G + Vaango::Util::Identity * K * e.Trace());


        Matrix66 D_mat = computeTangentModulus();
        for (int ii = 0; ii < 6; ii++) {
          for (int jj = 0; jj < 6; jj++) {
            D[ii][jj] = D_mat(ii,jj);
          }
        }

        // K_mat = B.transpose()*D*B*volume_old
        BtDB(B, D, K_mat);
        // K_geo = Bnl.transpose*sig*Bnl*volume_new;
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            sig[i][j] = pStress_old[idx](i, j);
          }
        }
        BnltDBnl(Bnl, sig, K_geo);

        double J = pDefGrad_new[idx].Determinant();
        double volume_init = (pMass[idx] / rho_orig);
        double volume_new = volume_init * J;

        for (int ii = 0; ii < 24; ii++) {
          for (int jj = 0; jj < 24; jj++) {
            K_vec[24 * ii + jj] = K_mat[ii][jj] * volume_init + K_geo[ii][jj] * volume_new;
          }
        }

        solver->fillMatrix(24, dof, 24, dof, K_vec);
      }
    }
  }
}

Matrix66
HypoElasticImplicit::computeTangentModulus() const
{
  double G = d_modelParam.G;
  double K = d_modelParam.K;
  double E = 9. * K * G / (3. * K + G);
  double PR = (3. * K - E) / (6. * K);
  double C11 = E * (1. - PR) / ((1. + PR) * (1. - 2. * PR));
  double C12 = E * PR / ((1. + PR) * (1. - 2. * PR));
  double C44 = G;

  Matrix66 D = Matrix66::Zero();
  D(0,0) = C11; D(0,1) = C12; D(0,2) = C12; D(0,3) = 0. ; D(0,4) = 0. ; D(0,5) = 0.;
                D(1,1) = C11; D(1,2) = C12; D(1,3) = 0. ; D(1,4) = 0. ; D(1,5) = 0.; 
                              D(2,2) = C11; D(2,3) = 0. ; D(2,4) = 0. ; D(2,5) = 0.;
                                            D(3,3) = C44; D(3,4) = 0. ; D(3,5) = 0.;
                                                          D(4,4) = C44; D(4,5) = 0.;
                                                                        D(5,5) = C44;
  return D;
}

void
HypoElasticImplicit::addComputesAndRequires(Task* task, const MPMMaterial* matl,
                                            const PatchSet*) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  bool reset = flag->d_doGridReset;

  addSharedCRForImplicitHypo(task, matlset, reset);
  task->computes(lb->StrainEnergyLabel);
}

void
HypoElasticImplicit::computeStressTensorImplicit(const PatchSubset* patches,
                                                 const MPMMaterial* matl,
                                                 DataWarehouse* old_dw,
                                                 DataWarehouse* new_dw)

{
  int matID = matl->getDWIndex();

  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  for (int pp = 0; pp < patches->size(); pp++) {
    const Patch* patch = patches->get(pp);
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    constParticleVariable<Matrix3> pDispGrad, pStress_old;
    constParticleVariable<double> pVolume_new;
    old_dw->get(pStress_old, lb->pStressLabel, pset);
    new_dw->get(pVolume_new, lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pDispGrad, lb->pDispGradLabel_preReloc, pset);

    ParticleVariable<double> pdTdt;
    ParticleVariable<Matrix3> pStress_new;
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);

    double G = d_modelParam.G;
    double bulk = d_modelParam.K;

    double strainEnergy = 0.0;
    if (matl->getIsRigid()) {
      for (int idx : *pset) {
        pStress_new[idx] = Matrix3(0.0);
        pdTdt[idx] = 0.;
      }
    } else {
      for (int idx : *pset) {
        pdTdt[idx] = 0.;

        // Calculate the strain (here called D), and deviatoric rate DPrime
        Matrix3 D = (pDispGrad[idx] + pDispGrad[idx].Transpose()) * .5;
        Matrix3 DPrime = D - Vaango::Util::Identity * Vaango::Util::one_third * D.Trace();

        // This is the (updated) Cauchy stress
        pStress_new[idx] = pStress_old[idx] + 
          (DPrime * 2. * G + Vaango::Util::Identity * bulk * D.Trace());

        // Compute the strain energy for all the particles
        Matrix3 avgStress = (pStress_new[idx] + pStress_old[idx]) * .5;
        double rateOfWork = computeRateOfWork(avgStress, D);
        strainEnergy += (rateOfWork * pVolume_new[idx] * delT);
      }

      if (flag->d_reductionVars->accStrainEnergy ||
          flag->d_reductionVars->strainEnergy) {
        new_dw->put(sum_vartype(strainEnergy), lb->StrainEnergyLabel);
      }
    }
  }
}

