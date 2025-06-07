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

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/TransIsoHyperImplicit.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/Constants.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Ports/DataWarehouse.h>
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
#include <Core/Math/Short27.h>              //for Fracture
#include <Core/Math/TangentModulusTensor.h> //added this for stiffness
#include <Core/ProblemSpec/ProblemSpec.h>
#include <fstream>
#include <iostream>

using namespace Uintah;

TransIsoHyperImplicit::TransIsoHyperImplicit(ProblemSpecP& ps, MPMFlags* Mflag)
  : TransIsoHyper(ps, Mflag)
  , ImplicitCM()
{
}

TransIsoHyperImplicit::TransIsoHyperImplicit(const TransIsoHyperImplicit* cm)
  : TransIsoHyper(cm)
  , ImplicitCM(cm)
{
}

TransIsoHyperImplicit::~TransIsoHyperImplicit()
{
}

std::unique_ptr<ConstitutiveModel>
TransIsoHyperImplicit::clone()
{
  return std::make_unique<TransIsoHyperImplicit>(*this);
}

void
TransIsoHyperImplicit::initializeCMData(const Patch* patch,
                                        const MPMMaterial* matl,
                                        DataWarehouse* new_dw)
{
  Matrix3 Identity, Zero(0.);
  Identity.Identity();

  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);
  ParticleVariable<Matrix3> pStress;
  ParticleVariable<double> pStretch, pFail;

  new_dw->allocateAndPut(pStress, lb->pStressLabel, pset);
  new_dw->allocateAndPut(pStretch, pStretchLabel, pset);
  new_dw->allocateAndPut(pFail, pFailureLabel, pset);

  for (auto idx : *pset) {
    pFail[idx]    = 0.0;
    pStress[idx]  = Zero;
    pStretch[idx] = 1.0;
  }
}

void
TransIsoHyperImplicit::addComputesAndRequires(Task* task,
                                              const MPMMaterial* matl,
                                              const PatchSet*,
                                              const bool /*recurse*/,
                                              const bool SchedParent) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  bool reset                    = flag->d_doGridReset;

  addSharedCRForImplicit(task, matlset, reset, true, SchedParent);

  if (SchedParent) {
    task->needs(Task::ParentOldDW, lb->pFiberDirLabel, matlset, Ghost::None);
    task->needs(Task::ParentOldDW, pFailureLabel, matlset, Ghost::None);
  } else {
    task->needs(Task::OldDW, lb->pFiberDirLabel, matlset, Ghost::None);
    task->needs(Task::OldDW, pFailureLabel, matlset, Ghost::None);
  }

  task->computes(lb->pFiberDirLabel_preReloc, matlset);
  task->computes(pStretchLabel_preReloc, matlset);
}

void
TransIsoHyperImplicit::computeStressTensorImplicit(const PatchSubset* patches,
                                                   const MPMMaterial* matl,
                                                   DataWarehouse* old_dw,
                                                   DataWarehouse* new_dw,
                                                   Solver* solver,
                                                   const bool)
{
  double rho_orig    = matl->getInitialDensity();
  double failure     = d_param.failure;

  DataWarehouse* parent_old_dw =
    new_dw->getOtherDataWarehouse(Task::ParentOldDW);
  int matID = matl->getDWIndex();

  for (int pp = 0; pp < patches->size(); pp++) {
    const Patch* patch   = patches->get(pp);
    ParticleSubset* pset = parent_old_dw->getParticleSubset(matID, patch);

    Vector dx      = patch->dCell();
    double oodx[3] = { 1. / dx.x(), 1. / dx.y(), 1. / dx.z() };

    IntVector lowIndex  = patch->getNodeLowIndex();
    IntVector highIndex = patch->getNodeHighIndex() + IntVector(1, 1, 1);
    Array3<int> l2g(lowIndex, highIndex);
    solver->copyL2G(l2g, patch);

    auto interpolator = flag->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());

    constParticleVariable<Point> pX;
    constParticleVariable<double> pVolume_new, pMass, pFail_old;
    constParticleVariable<Vector> pVelocity, pFiberDir;
    constParticleVariable<Matrix3> pSize, pDefGrad_new;

    parent_old_dw->get(pX, lb->pXLabel, pset);
    parent_old_dw->get(pMass, lb->pMassLabel, pset);
    parent_old_dw->get(pFail_old, pFailureLabel, pset);
    parent_old_dw->get(pFiberDir, lb->pFiberDirLabel, pset);
    parent_old_dw->get(pSize, lb->pSizeLabel, pset);

    new_dw->get(pVolume_new, lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    ParticleVariable<double> pStretch, pFail;
    ParticleVariable<Matrix3> pStress;

    new_dw->allocateAndPut(pStretch, pStretchLabel_preReloc, pset);
    new_dw->allocateAndPut(pFail, pFailureLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress, lb->pStressLabel_preReloc, pset);

    double B[6][24];
    double Bnl[3][24];
    double K_mat[24][24];
    double K_geo[24][24];
    double K_vec[576];
    double D[6][6]; // stiffness matrix

    if (matl->getIsRigid()) {
      for (int idx : *pset) {
        pStress[idx] = Vaango::Util::Zero;
      }
    } else {
      for (int idx : *pset) {

        // Get the node indices that surround the cell
        interpolator->findCellAndShapeDerivatives(
          pX[idx], ni, d_S, pSize[idx], pDefGrad_new[idx]);
        int dof[24];
        loadBMats(l2g, dof, B, Bnl, d_S, ni, oodx);

        // Compute deformation state
        TransIsoHyperState ss =
          computeDeformationState(pDefGrad_new[idx], pFiberDir[idx]);

        pStretch[idx] = std::sqrt(ss.I4);

        // Failure+Stress+Stiffness
        pFail[idx] = 0.;
        Matrix3 pressure, deviatoric_stress, fiber_stress;
        if (failure == 1) {
          computeStressWithFailure(ss,
                                   pStretch[idx],
                                   pFail_old[idx],
                                   pFail[idx],
                                   pressure,
                                   deviatoric_stress,
                                   fiber_stress);
        } else {
          pressure          = computeHydrostaticStress(ss);
          deviatoric_stress = computeDeviatoricStress(ss);
          fiber_stress      = computeFiberStress(ss);
        }

        // Cauchy stress
        pStress[idx] = pressure + deviatoric_stress + fiber_stress;
        double p     = pressure.Trace() / 3.0;

        // Stiffness matrix
        Matrix66 D_mat = computeTangentStiffness(ss, p, pFail[idx]);
        for (int ii = 0; ii < 6; ++ii) {
          for (int jj = 0; jj < 6; ++jj) {
            D[ii][jj] = D_mat(ii, jj);
          }
        }

        // K_mat = B.transpose()*D*B*volume_old
        BtDB(B, D, K_mat);
        // K_geo = Bnl.transpose*sig*Bnl*volume_new;
        double sig[3][3];
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            sig[i][j] = pStress[idx](i, j);
          }
        }
        BnltDBnl(Bnl, sig, K_geo);
        double volume_old = (pMass[idx] / rho_orig);
        double volume_new = pVolume_new[idx];
        for (int ii = 0; ii < 24; ii++) {
          for (int jj = 0; jj < 24; jj++) {
            K_mat[ii][jj] *= volume_old;
            K_geo[ii][jj] *= volume_new;
          }
        }
        for (int I = 0; I < 24; I++) {
          for (int J = 0; J < 24; J++) {
            K_vec[24 * I + J] = K_mat[I][J] + K_geo[I][J];
          }
        }
        solver->fillMatrix(24, dof, 24, dof, K_vec);
      } // end of loop over particles
    }
  }
}

Matrix66
TransIsoHyperImplicit::computeTangentStiffness(const TransIsoHyperState& ss,
                                               double p,
                                               double pFailed) const
{
  Matrix66 C_vol = computeVolTangentStiffness(ss, p);
  Matrix66 C_dev = Matrix66::Zero();
  if (pFailed != 1.0 && pFailed != 3.0) {
    C_dev = computeDevTangentStiffness(ss);
  }
  Matrix66 C_fiber = Matrix66::Zero();
  if (pFailed != 2.0 && pFailed != 3.0) {
    C_fiber = computeFiberTangentStiffness(ss);
  }
  Matrix66 D_mat = C_vol + C_dev + C_fiber;
  
  return D_mat;
}

Matrix66
TransIsoHyperImplicit::computeVolTangentStiffness(const TransIsoHyperState& ss,
                                                  double p) const
{
  double J       = ss.J;
  double K       = d_param.bulkModulus;
  Matrix66 C_vol = Matrix66::Zero();
  C_vol(0, 0) = K * (1. / J) - 2 * p;
  C_vol(0, 1) = K * (1. / J);
  C_vol(0, 2) = K * (1. / J);
  C_vol(1, 1) = K * (1. / J) - 2 * p;
  C_vol(1, 2) = K * (1. / J);
  C_vol(2, 2) = K * (1. / J) - 2 * p;
  C_vol(3, 3) = -p;
  C_vol(4, 4) = -p;
  C_vol(5, 5) = -p;

  return C_vol;
}

Matrix66
TransIsoHyperImplicit::computeDevTangentStiffness(
  const TransIsoHyperState& ss) const
{
  double c1 = d_param.c1;
  double c2 = d_param.c2;

  double J      = ss.J;
  double I1_bar = ss.I1_bar;
  double I2_bar = ss.I2_bar;
  Matrix3 Bbar  = ss.B_bar;
  Matrix3 Bbar2 = Bbar * Bbar;

  double cc2   = (4. / 3.) * (1. / J) * (c1 * I1_bar + 2. * c2 * I2_bar);
  Matrix3 devs = (Bbar * (c1 + c2 * I1_bar) - Bbar2 * c2 +
                  Vaango::Util::Identity * (c1 * I1_bar + 2 * c2 * I2_bar)) *
                 (2. / J) * (-2. / 3.);
  Matrix3 termMR = Bbar * (1. / J) * c2 * I1_bar - Bbar2 * (1. / J) * c2;

  Matrix66 C_dev = Matrix66::Zero();
  C_dev(0, 0) =
    (4. / J) * c2 * Bbar(0, 0) * Bbar(0, 0) -
    (4. / J) * c2 * (Bbar(0, 0) * Bbar(0, 0) + Bbar(0, 0) * Bbar(0, 0)) +
    (2. / 3.) * cc2 + (4. / 9.) * (1. / J) * 2 * c2 * I2_bar + devs(0, 0) +
    devs(0, 0) + (-4. / 3.) * (termMR(0, 0) + termMR(0, 0));
  C_dev(0, 1) =
    (4. / J) * c2 * Bbar(0, 0) * Bbar(1, 1) -
    (4. / J) * c2 * (Bbar(0, 1) * Bbar(0, 1) + Bbar(0, 1) * Bbar(0, 1)) +
    (-1. / 3.) * cc2 + devs(0, 0) + devs(1, 1) +
    (-4. / 3.) * (termMR(0, 0) + termMR(1, 1));
  C_dev(0, 2) =
    (4. / J) * c2 * Bbar(0, 0) * Bbar(2, 2) -
    (4. / J) * c2 * (Bbar(0, 2) * Bbar(0, 2) + Bbar(0, 2) * Bbar(0, 2)) +
    (-1. / 3.) * cc2 + devs(0, 0) + devs(2, 2) +
    (-4. / 3.) * (termMR(0, 0) + termMR(2, 2));
  C_dev(0, 3) =
    (4. / J) * c2 * Bbar(0, 0) * Bbar(0, 1) -
    (4. / J) * c2 * (Bbar(0, 0) * Bbar(0, 1) + Bbar(0, 1) * Bbar(0, 0)) +
    devs(0, 1) + (-4. / 3.) * termMR(0, 1);
  C_dev(0, 4) =
    (4. / J) * c2 * Bbar(0, 0) * Bbar(1, 2) -
    (4. / J) * c2 * (Bbar(0, 1) * Bbar(0, 2) + Bbar(0, 2) * Bbar(0, 1)) +
    devs(1, 2) + (-4. / 3.) * termMR(1, 2);
  C_dev(0, 5) =
    (4. / J) * c2 * Bbar(0, 0) * Bbar(2, 0) -
    (4. / J) * c2 * (Bbar(0, 2) * Bbar(0, 0) + Bbar(0, 0) * Bbar(0, 2)) +
    devs(2, 0) + (-4. / 3.) * termMR(2, 0);

  C_dev(1, 1) =
    (4. / J) * c2 * Bbar(1, 1) * Bbar(1, 1) -
    (4. / J) * c2 * (Bbar(1, 1) * Bbar(1, 1) + Bbar(1, 1) * Bbar(1, 1)) +
    (2. / 3.) * cc2 + (4. / 9.) * (1. / J) * 2 * c2 * I2_bar + devs(1, 1) +
    devs(1, 1) + (-4. / 3.) * (termMR(1, 1) + termMR(1, 1));
  C_dev(1, 2) =
    (4. / J) * c2 * Bbar(1, 1) * Bbar(2, 2) -
    (4. / J) * c2 * (Bbar(1, 2) * Bbar(1, 2) + Bbar(1, 2) * Bbar(1, 2)) +
    (-1. / 3.) * cc2 + devs(1, 1) + devs(2, 2) +
    (-4. / 3.) * (termMR(1, 1) + termMR(2, 2));
  C_dev(1, 3) =
    (4. / J) * c2 * Bbar(1, 1) * Bbar(0, 1) -
    (4. / J) * c2 * (Bbar(1, 0) * Bbar(1, 1) + Bbar(1, 1) * Bbar(1, 0)) +
    devs(0, 1) + (-4. / 3.) * termMR(0, 1);
  C_dev(1, 4) =
    (4. / J) * c2 * Bbar(1, 1) * Bbar(1, 2) -
    (4. / J) * c2 * (Bbar(1, 1) * Bbar(1, 2) + Bbar(1, 2) * Bbar(1, 1)) +
    devs(1, 2) + (-4. / 3.) * termMR(1, 2);
  C_dev(1, 5) =
    (4. / J) * c2 * Bbar(1, 1) * Bbar(2, 0) -
    (4. / J) * c2 * (Bbar(1, 2) * Bbar(1, 0) + Bbar(1, 0) * Bbar(1, 2)) +
    devs(2, 0) + (-4. / 3.) * termMR(2, 0);

  C_dev(2, 2) =
    (4. / J) * c2 * Bbar(2, 2) * Bbar(2, 2) -
    (4. / J) * c2 * (Bbar(2, 2) * Bbar(2, 2) + Bbar(2, 2) * Bbar(2, 2)) +
    (2. / 3.) * cc2 + (4. / 9.) * (1. / J) * 2 * c2 * I2_bar + devs(2, 2) +
    devs(2, 2) + (-4. / 3.) * (termMR(2, 2) + termMR(2, 2));
  C_dev(2, 3) =
    (4. / J) * c2 * Bbar(2, 2) * Bbar(0, 1) -
    (4. / J) * c2 * (Bbar(2, 0) * Bbar(2, 1) + Bbar(2, 1) * Bbar(2, 0)) +
    devs(0, 1) + (-4. / 3.) * termMR(0, 1);
  C_dev(2, 4) =
    (4. / J) * c2 * Bbar(2, 2) * Bbar(1, 2) -
    (4. / J) * c2 * (Bbar(2, 1) * Bbar(2, 2) + Bbar(2, 2) * Bbar(2, 1)) +
    devs(1, 2) + (-4. / 3.) * termMR(1, 2);
  C_dev(2, 5) =
    (4. / J) * c2 * Bbar(2, 2) * Bbar(2, 0) -
    (4. / J) * c2 * (Bbar(2, 2) * Bbar(2, 0) + Bbar(2, 0) * Bbar(2, 2)) +
    devs(2, 0) + (-4. / 3.) * termMR(2, 0);

  C_dev(3, 3) =
    (4. / J) * c2 * Bbar(0, 1) * Bbar(0, 1) -
    (4. / J) * c2 * (Bbar(0, 0) * Bbar(1, 1) + Bbar(0, 1) * Bbar(1, 0)) +
    (1. / 2.) * cc2 + (4. / 9.) * (1. / J) * 2 * c2 * I2_bar;
  C_dev(3, 4) =
    (4. / J) * c2 * Bbar(0, 1) * Bbar(1, 2) -
    (4. / J) * c2 * (Bbar(0, 1) * Bbar(1, 2) + Bbar(0, 2) * Bbar(1, 1));
  C_dev(3, 5) =
    (4. / J) * c2 * Bbar(0, 1) * Bbar(2, 0) -
    (4. / J) * c2 * (Bbar(0, 2) * Bbar(1, 0) + Bbar(0, 0) * Bbar(1, 2));

  C_dev(4, 4) =
    (4. / J) * c2 * Bbar(1, 2) * Bbar(1, 2) -
    (4. / J) * c2 * (Bbar(1, 1) * Bbar(2, 2) + Bbar(1, 2) * Bbar(2, 1)) +
    (1. / 2.) * cc2 + (4. / 9.) * (1. / J) * 2 * c2 * I2_bar;
  C_dev(4, 5) =
    (4. / J) * c2 * Bbar(1, 2) * Bbar(2, 0) -
    (4. / J) * c2 * (Bbar(1, 2) * Bbar(2, 0) + Bbar(1, 0) * Bbar(2, 2));

  C_dev(5, 5) =
    (4. / J) * c2 * Bbar(2, 0) * Bbar(2, 0) -
    (4. / J) * c2 * (Bbar(2, 2) * Bbar(0, 0) + Bbar(2, 0) * Bbar(0, 2)) +
    (1. / 2.) * cc2 + (4. / 9.) * (1. / J) * 2 * c2 * I2_bar;

  return C_dev;
}

Matrix66
TransIsoHyperImplicit::computeFiberTangentStiffness(
  const TransIsoHyperState& ss) const
{
  double J         = ss.J;
  double dWdI4_bar = ss.dWdI4_bar;
  double I4_bar    = ss.I4_bar;

  Matrix3 ff_dyad(ss.fiberDir_new, ss.fiberDir_new);

  double cc1   = ss.d2WdI4_bar2 * I4_bar * I4_bar;
  double cc2   = (4. / 3.) * (1. / J) * dWdI4_bar * I4_bar;
  Matrix3 devs = (ff_dyad - Vaango::Util::Identity * (1. / 3.)) * dWdI4_bar *
                 I4_bar * (2. / J) * (-2. / 3.);

  Matrix66 C_fiber = Matrix66::Zero();
  C_fiber(0, 0) =
    (2. / 3.) * cc2 + (4. / 9.) * (1. / J) * cc1 + devs(0, 0) + devs(0, 0) +
    (-4. / 3.) * (1. / J) * cc1 * (ff_dyad(0, 0) + ff_dyad(0, 0)) +
    (4. / J) * cc1 * ff_dyad(0, 0) * ff_dyad(0, 0);
  C_fiber(0, 1) =
    (-1. / 3.) * cc2 + devs(0, 0) + devs(1, 1) +
    (-4. / 3.) * (1. / J) * cc1 * (ff_dyad(0, 0) + ff_dyad(1, 1)) +
    (4. / J) * cc1 * ff_dyad(0, 0) * ff_dyad(1, 1);
  C_fiber(0, 2) =
    (-1. / 3.) * cc2 + devs(0, 0) + devs(2, 2) +
    (-4. / 3.) * (1. / J) * cc1 * (ff_dyad(0, 0) + ff_dyad(2, 2)) +
    (4. / J) * cc1 * ff_dyad(0, 0) * ff_dyad(2, 2);
  C_fiber(1, 1) =
    (2. / 3.) * cc2 + (4. / 9.) * (1. / J) * cc1 + devs(1, 1) + devs(1, 1) +
    (-4. / 3.) * (1. / J) * cc1 * (ff_dyad(1, 1) + ff_dyad(1, 1)) +
    (4. / J) * cc1 * ff_dyad(1, 1) * ff_dyad(1, 1);
  C_fiber(1, 2) =
    (-1. / 3.) * cc2 + devs(1, 1) + devs(2, 2) +
    (-4. / 3.) * (1. / J) * cc1 * (ff_dyad(1, 1) + ff_dyad(2, 2)) +
    (4. / J) * cc1 * ff_dyad(1, 1) * ff_dyad(2, 2);
  C_fiber(2, 2) =
    (2. / 3.) * cc2 + (4. / 9.) * (1. / J) * cc1 + devs(2, 2) + devs(2, 2) +
    (-4. / 3.) * (1. / J) * cc1 * (ff_dyad(2, 2) + ff_dyad(2, 2)) +
    (4. / J) * cc1 * ff_dyad(2, 2) * ff_dyad(2, 2);
  C_fiber(0, 3) = devs(0, 1) + (-4. / 3.) * (1. / J) * cc1 * ff_dyad(0, 1) +
                  (4. / J) * cc1 * ff_dyad(0, 0) * ff_dyad(0, 1);
  C_fiber(0, 4) = devs(1, 2) + (-4. / 3.) * (1. / J) * cc1 * ff_dyad(1, 2) +
                  (4. / J) * cc1 * ff_dyad(0, 0) * ff_dyad(1, 2);
  C_fiber(0, 5) = devs(2, 0) + (-4. / 3.) * (1. / J) * cc1 * ff_dyad(2, 0) +
                  (4. / J) * cc1 * ff_dyad(0, 0) * ff_dyad(2, 0);
  C_fiber(1, 3) = devs(0, 1) + (-4. / 3.) * (1. / J) * cc1 * ff_dyad(0, 1) +
                  (4. / J) * cc1 * ff_dyad(1, 1) * ff_dyad(0, 1);
  C_fiber(1, 4) = devs(1, 2) + (-4. / 3.) * (1. / J) * cc1 * ff_dyad(1, 2) +
                  (4. / J) * cc1 * ff_dyad(1, 1) * ff_dyad(1, 2);
  C_fiber(1, 5) = devs(2, 0) + (-4. / 3.) * (1. / J) * cc1 * ff_dyad(2, 0) +
                  (4. / J) * cc1 * ff_dyad(1, 1) * ff_dyad(2, 0);
  C_fiber(2, 3) = devs(0, 1) + (-4. / 3.) * (1. / J) * cc1 * ff_dyad(0, 1) +
                  (4. / J) * cc1 * ff_dyad(2, 2) * ff_dyad(0, 1);
  C_fiber(2, 4) = devs(1, 2) + (-4. / 3.) * (1. / J) * cc1 * ff_dyad(1, 2) +
                  (4. / J) * cc1 * ff_dyad(2, 2) * ff_dyad(1, 2);
  C_fiber(2, 5) = devs(2, 0) + (-4. / 3.) * (1. / J) * cc1 * ff_dyad(2, 0) +
                  (4. / J) * cc1 * ff_dyad(2, 2) * ff_dyad(2, 0);
  C_fiber(3, 3) = (1. / 2.) * cc2 + (4. / 9.) * (1. / J) * cc1 +
                  (4. / J) * cc1 * ff_dyad(0, 1) * ff_dyad(0, 1);
  C_fiber(3, 4) = (4. / J) * cc1 * ff_dyad(0, 1) * ff_dyad(1, 2);
  C_fiber(3, 5) = (4. / J) * cc1 * ff_dyad(0, 1) * ff_dyad(2, 0);
  C_fiber(4, 4) = (1. / 2.) * cc2 + (4. / 9.) * (1. / J) * cc1 +
                  (4. / J) * cc1 * ff_dyad(1, 2) * ff_dyad(1, 2);
  C_fiber(4, 5) = (4. / J) * cc1 * ff_dyad(1, 2) * ff_dyad(2, 0);
  C_fiber(5, 5) = (1. / 2.) * cc2 + (4. / 9.) * (1. / J) * cc1 +
                  (4. / J) * cc1 * ff_dyad(2, 0) * ff_dyad(2, 0);

  return C_fiber;
}

void
TransIsoHyperImplicit::addComputesAndRequires(Task* task,
                                              const MPMMaterial* matl,
                                              const PatchSet*) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  bool reset                    = flag->d_doGridReset;

  addSharedCRForImplicit(task, matlset, reset);

  task->needs(Task::OldDW, lb->pFiberDirLabel, matlset, Ghost::None);
  task->needs(Task::OldDW, pFailureLabel, matlset, Ghost::None);

  task->computes(lb->pFiberDirLabel_preReloc, matlset);
  task->computes(pStretchLabel_preReloc, matlset);
  task->computes(pFailureLabel_preReloc, matlset);
}

void
TransIsoHyperImplicit::computeStressTensorImplicit(const PatchSubset* patches,
                                                   const MPMMaterial* matl,
                                                   DataWarehouse* old_dw,
                                                   DataWarehouse* new_dw)
{
  double failure     = d_param.failure;
  int matID = matl->getDWIndex();

  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  for (int pp = 0; pp < patches->size(); pp++) {
    const Patch* patch = patches->get(pp);
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    constParticleVariable<double> pFail_old;
    constParticleVariable<Vector> pFiberDir;
    constParticleVariable<Matrix3> pDefGrad_new;

    old_dw->get(pFail_old, pFailureLabel, pset);
    old_dw->get(pFiberDir, lb->pFiberDirLabel, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    ParticleVariable<double> pStretch, pFail;
    ParticleVariable<Vector> pFiberDir_copy;
    ParticleVariable<Matrix3> pStress;

    new_dw->allocateAndPut(pStretch, pStretchLabel_preReloc, pset);
    new_dw->allocateAndPut(pFail, pFailureLabel_preReloc, pset);
    new_dw->allocateAndPut(pFiberDir_copy, lb->pFiberDirLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress, lb->pStressLabel_preReloc, pset);

    if (matl->getIsRigid()) {
      for (int idx : *pset) {
        pStress[idx] = Vaango::Util::Zero;
      }
    } else {
      for (int idx : *pset) {

        // carry forward fiber direction
        pFiberDir_copy[idx] = pFiberDir[idx];

        // Compute deformation state
        TransIsoHyperState ss =
          computeDeformationState(pDefGrad_new[idx], pFiberDir[idx]);

        pStretch[idx] = std::sqrt(ss.I4);

        // Failure+Stress+Stiffness
        pFail[idx] = 0.;
        Matrix3 pressure, deviatoric_stress, fiber_stress;
        if (failure == 1) {
          computeStressWithFailure(ss,
                                   pStretch[idx],
                                   pFail_old[idx],
                                   pFail[idx],
                                   pressure,
                                   deviatoric_stress,
                                   fiber_stress);
        } else {
          pressure          = computeHydrostaticStress(ss);
          deviatoric_stress = computeDeviatoricStress(ss);
          fiber_stress      = computeFiberStress(ss);
        }

        // Cauchy stress
        pStress[idx] = pressure + deviatoric_stress + fiber_stress;

      } // end loop over particles
    }   // isn't rigid
  }
}
