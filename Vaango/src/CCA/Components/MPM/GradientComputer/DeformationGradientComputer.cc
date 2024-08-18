/*
 * The MIT License
 *
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

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/Constants.h>
#include <CCA/Components/MPM/GradientComputer/DeformationGradientComputer.h>
#include <CCA/Components/MPM/GradientComputer/DisplacementGradientComputer.h>
#include <CCA/Components/MPM/GradientComputer/VelocityGradientComputer.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Util/DebugStream.h>

#include <iostream>

// #define IGNORE_NEGATIVE_JACOBIANS
// #define DEBUG_WITH_PARTICLE_ID

using namespace Uintah;

static DebugStream dbg_doing("DefGrad_doing", false);
static DebugStream dbg("DefGrad", false);

const Matrix3 DeformationGradientComputer::Identity(1.0,
                                                    0.0,
                                                    0.0,
                                                    0.0,
                                                    1.0,
                                                    0.0,
                                                    0.0,
                                                    0.0,
                                                    1.0);
const Matrix3 DeformationGradientComputer::Zero(0.0);

DeformationGradientComputer::DeformationGradientComputer(
  MaterialManagerP& ss,
  const MPMLabel* mpm_labels,
  const MPMFlags* mpm_flags)
{
  lb   = mpm_labels;
  flag = mpm_flags;
  // std::cout << "d_or_27 = " << flag->d_8or27 << "\n";
  if (flag->d_8or27 == 8) {
    NGN = 1;
  } else {
    NGN = 2;
  }
  d_mat_manager = ss;
}

DeformationGradientComputer::DeformationGradientComputer(
  const DeformationGradientComputer* dg)
{
  lb            = dg->lb;
  flag          = dg->flag;
  NGN           = dg->NGN;
  NGP           = dg->NGP;
  d_mat_manager = dg->d_mat_manager;
}

DeformationGradientComputer*
DeformationGradientComputer::clone()
{
  return scinew DeformationGradientComputer(*this);
}

DeformationGradientComputer::~DeformationGradientComputer()
{
}

// Initial velocity/displacement/deformation gradients
void
DeformationGradientComputer::addInitialComputesAndRequires(
  Task* task,
  const MPMMaterial* mpm_matl,
  const PatchSet*)
{
  dbg_doing << "Doing DefGrad::addInitialComputesAndRequires\n";

  const MaterialSubset* matlset = mpm_matl->thisMaterial();

  task->computes(lb->pVelGradLabel, matlset);
  task->computes(lb->pDispGradLabel, matlset);
  task->computes(lb->pDefGradLabel, matlset);

  ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
  if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
    task->computes(lb->pRemoveLabel, matlset);
    task->computes(lb->pPolarDecompRLabel, matlset);
  }
}

void
DeformationGradientComputer::initializeGradient(const Patch* patch,
                                                const MPMMaterial* mpm_matl,
                                                DataWarehouse* new_dw)
{
  dbg_doing << "Doing DefGrad::initializeGradient\n";

  if (flag->d_integrator == MPMFlags::Implicit) {
    initializeGradientImplicit(patch, mpm_matl, new_dw);
  } else {
    initializeGradientExplicit(patch, mpm_matl, new_dw);
  }

  ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
  if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
    ParticleSubset* pset =
      new_dw->getParticleSubset(mpm_matl->getDWIndex(), patch);
    ParticleVariable<int> pRemove;
    ParticleVariable<Matrix3> pPolarDecompR;
    new_dw->allocateAndPut(pRemove, lb->pRemoveLabel, pset);
    new_dw->allocateAndPut(pPolarDecompR, lb->pPolarDecompRLabel, pset);
    for (auto particle : *pset) {
      pRemove[particle]       = 0;
      pPolarDecompR[particle] = Identity;
    }
  }
}

void
DeformationGradientComputer::initializeGradientExplicit(
  const Patch* patch,
  const MPMMaterial* mpm_matl,
  DataWarehouse* new_dw)
{
  dbg_doing << "Doing DefGrad::initializeGradientExplicit\n";

  ParticleSubset* pset =
    new_dw->getParticleSubset(mpm_matl->getDWIndex(), patch);
  ParticleVariable<Matrix3> pVelGrad, pDispGrad, pDefGrad;
  new_dw->allocateAndPut(pVelGrad, lb->pVelGradLabel, pset);
  new_dw->allocateAndPut(pDispGrad, lb->pDispGradLabel, pset);
  new_dw->allocateAndPut(pDefGrad, lb->pDefGradLabel, pset);

  for (auto particle : *pset) {
    pVelGrad[particle]  = Zero;
    pDispGrad[particle] = Zero;
    pDefGrad[particle]  = Identity;
  }
}

void
DeformationGradientComputer::initializeGradientImplicit(
  const Patch* patch,
  const MPMMaterial* mpm_matl,
  DataWarehouse* new_dw)
{
  dbg_doing << "Doing DefGrad::initializeGradientImplicit\n";

  ParticleSubset* pset =
    new_dw->getParticleSubset(mpm_matl->getDWIndex(), patch);
  ParticleVariable<Matrix3> pVelGrad, pDispGrad, pDefGrad;
  new_dw->allocateAndPut(pVelGrad, lb->pVelGradLabel, pset);
  new_dw->allocateAndPut(pDispGrad, lb->pDispGradLabel, pset);
  new_dw->allocateAndPut(pDefGrad, lb->pDefGradLabel, pset);

  for (auto particle : *pset) {
    pVelGrad[particle]  = Zero;
    pDispGrad[particle] = Zero;
    pDefGrad[particle]  = Identity;
  }
}

// Task scheduling for general case
void
DeformationGradientComputer::addComputesAndRequires(Task* task,
                                                    const MPMMaterial* mpm_matl,
                                                    const PatchSet*)
{
  dbg_doing << "Doing DefGrad::addComputesAndRequires\n";

  if (flag->d_integrator == MPMFlags::Implicit) {
    addComputesAndRequiresImplicit(task, mpm_matl);
  } else {
    addComputesAndRequiresExplicit(task, mpm_matl);
  }

  ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
  if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
    const MaterialSubset* matlset = mpm_matl->thisMaterial();
    task->requires(Task::OldDW, lb->pPolarDecompRLabel, matlset, Ghost::None);
    task->computes(lb->pRemoveLabel_preReloc, matlset);
    task->computes(lb->pPolarDecompRLabel_preReloc, matlset);
    task->computes(lb->pPolarDecompRMidLabel, matlset);
  }
}

void
DeformationGradientComputer::addComputesOnly(Task* task,
                                             const MPMMaterial* mpm_matl,
                                             const PatchSet*)
{
  dbg_doing << "Doing DefGrad::addComputesOnly\n";

  std::ostringstream out;
  out << "**ERROR**: addComputesOnly Not implemented "
      << "\n";
  throw InvalidValue(out.str(), __FILE__, __LINE__);
}

void
DeformationGradientComputer::addComputesAndRequiresExplicit(
  Task* task,
  const MPMMaterial* mpm_matl)
{
  dbg_doing << "Doing DefGrad::addComputesAndRequiresExplicit\n";

  Ghost::GhostType gnone = Ghost::None;
  Ghost::GhostType gac   = Ghost::AroundCells;

  const MaterialSubset* matlset = mpm_matl->thisMaterial();

  // Requires (for explicit)
  task->requires(Task::OldDW, lb->delTLabel);
  task->requires(Task::OldDW, lb->pXLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pMassLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pSizeLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pVolumeLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pVelocityLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pVelGradLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pDefGradLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pParticleIDLabel, matlset, gnone);
  if (flag->d_doGridReset) {
    task->requires(Task::NewDW, lb->gVelocityStarLabel, matlset, gac, NGN);
    if (flag->d_fracture) {
      task->requires(Task::NewDW, lb->pgCodeLabel, matlset, gnone);
      task->requires(Task::NewDW, lb->GVelocityStarLabel, matlset, gac, NGN);
    }
  } else {
    task->requires(Task::NewDW, lb->gDisplacementLabel, matlset, gac, NGN);
  }

  // Computes (for explicit)
  task->computes(lb->pVelGradLabel_preReloc, matlset);
  task->computes(lb->pDispGradLabel_preReloc, matlset);
  task->computes(lb->pDefGradLabel_preReloc, matlset);
  task->computes(lb->pVolumeLabel_preReloc, matlset);
}

void
DeformationGradientComputer::addComputesAndRequiresImplicit(
  Task* task,
  const MPMMaterial* mpm_matl)
{
  dbg_doing << "Doing DefGrad::addComputesAndRequiresImplicit\n";

  Ghost::GhostType gnone = Ghost::None;
  Ghost::GhostType gac   = Ghost::AroundCells;

  const MaterialSubset* matlset = mpm_matl->thisMaterial();

  // Requires (for implicit)
  task->requires(Task::OldDW, lb->delTLabel);
  task->requires(Task::OldDW, lb->pXLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pMassLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pSizeLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pVolumeLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pDefGradLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pParticleIDLabel, matlset, gnone);
  if (flag->d_doGridReset) {
    task->requires(Task::NewDW, lb->dispNewLabel, matlset, gac, NGN);
  } else {
    task->requires(Task::NewDW, lb->gDisplacementLabel, matlset, gac, NGN);
  }

  // Computes (for implicit)
  task->computes(lb->pVelGradLabel_preReloc, matlset);
  task->computes(lb->pDispGradLabel_preReloc, matlset);
  task->computes(lb->pDefGradLabel_preReloc, matlset);
  task->computes(lb->pVolumeLabel_preReloc, matlset);
}

//-----------------------------------------------------------------------
//  Actually compute deformation gradient
//-----------------------------------------------------------------------
void
DeformationGradientComputer::computeDeformationGradient(
  const PatchSubset* patches,
  DataWarehouse* old_dw,
  DataWarehouse* new_dw)
{
  dbg_doing << "Doing DefGrad::computeDeformationGradient\n";

  // The explicit code uses the velocity gradient to compute the
  // deformation gradient.  The implicit code uses displacements.
  if (flag->d_integrator == MPMFlags::Implicit) {

    // std::cout << "Compute def grad .. Implicit ..\n";
    // std::cout << "old_dw = " << old_dw << "\n";

    // Get delT
    delt_vartype delT;
    old_dw->get(delT, lb->delTLabel, getLevel(patches));

    // Loop thru patches
    for (int pp = 0; pp < patches->size(); pp++) {
      const Patch* patch = patches->get(pp);

      int numMPMMatls = d_mat_manager->getNumMaterials("MPM");
      for (int m = 0; m < numMPMMatls; m++) {

        // Get particle info and patch info
        MPMMaterial* mpm_matl =
          static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));

        // Compute deformation gradient
        computeDeformationGradientImplicit(patch,
                                           mpm_matl,
                                           delT,
                                           old_dw,
                                           new_dw);
      }
    }

  } else {

    // std::cout << "Compute def grad .. Explicit ..\n";

    // Get delT
    delt_vartype delT;
    old_dw->get(delT, lb->delTLabel, getLevel(patches));

    // Loop thru patches
    for (int pp = 0; pp < patches->size(); pp++) {
      const Patch* patch = patches->get(pp);

      // std::cout << "Compute def grad .. 1 .. pp = " << pp << " patch = " <<
      // patch << "\n";
      int numMPMMatls = d_mat_manager->getNumMaterials("MPM");
      for (int m = 0; m < numMPMMatls; m++) {

        // Get particle info and patch info
        MPMMaterial* mpm_matl =
          static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
        // std::cout << "Compute def grad .. 2 .. pp = " << pp << " m = " << m
        // << " mpm_matl = " << mpm_matl << "\n";

        // Compute deformation gradient
        computeDeformationGradientExplicit(patch,
                                           mpm_matl,
                                           delT,
                                           old_dw,
                                           new_dw);
        // std::cout << "Compute def grad .. 3 .. complete" << "\n";
      }
    }
  }
}

//-----------------------------------------------------------------------
//  Actually compute deformation gradient (for implicit only)
//-----------------------------------------------------------------------
void
DeformationGradientComputer::computeDeformationGradient(
  const PatchSubset* patches,
  DataWarehouse* old_dw,
  DataWarehouse* new_dw,
  bool recurse)
{
  dbg_doing << "Doing DefGrad::computeDeformationGradient implicit recursion\n";

  DataWarehouse* parent_old_dw =
    new_dw->getOtherDataWarehouse(Task::ParentOldDW);
  // std::cout << "parent_old_dw = " << parent_old_dw << " old_dw = " << old_dw
  //           << " new_dw = " << new_dw << "\n";
  //  Get delT
  delt_vartype delT;
  parent_old_dw->get(delT, lb->delTLabel, getLevel(patches));

  for (int pp = 0; pp < patches->size(); pp++) {
    const Patch* patch = patches->get(pp);
    int numMPMMatls    = d_mat_manager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      computeDeformationGradientImplicit(patch,
                                         mpm_matl,
                                         delT,
                                         old_dw,
                                         parent_old_dw,
                                         new_dw);
    }
  }
}

// Compute deformation gradient for explicit computations from velocity gradient
void
DeformationGradientComputer::computeDeformationGradientExplicit(
  const Patch* patch,
  const MPMMaterial* mpm_matl,
  const double& delT,
  DataWarehouse* old_dw,
  DataWarehouse* new_dw)
{
  ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
  int dwi               = mpm_matl->getDWIndex();

  dbg_doing << "Doing DefGrad::computeDeformationGradient explicit: mat = "
            << dwi << " modelType = " << cm->modelType() << "\n";

  // Constants
  // Ghost::GhostType  gnone = Ghost::None;
  Ghost::GhostType gac = Ghost::AroundCells;

  // Get particle info and patch info
  ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);
  Vector dx            = patch->dCell();
  double oodx[3]       = { 1. / dx.x(), 1. / dx.y(), 1. / dx.z() };

  // std::cout << "One . patch = " << patch << " mpm_matl = " << mpm_matl <<
  // "\n";

  // Get initial density
  double rho_orig = mpm_matl->getInitialDensity();

  // Get Interpolator
  // std::cout << "interpolator = " << &(flag->d_interpolator) << "\n";
  auto interpolator = flag->d_interpolator->clone(patch);

  // Set up variables to store old particle and grid data
  // for vel grad and def grad calculation
  constParticleVariable<Short27> pgCode; // For fracture MPM
  constParticleVariable<double> pMass;
  constParticleVariable<Point> px;
  constParticleVariable<Matrix3> pDefGrad_old, pVelGrad_old, pDispGrad_old;
  constParticleVariable<Matrix3> pSize;

  // Set up variables to store new particle and grid data
  // for vel grad and def grad calculation
  ParticleVariable<double> pVolume_new;
  ParticleVariable<Matrix3> pDefGrad_new, pVelGrad_new, pDispGrad_new;
  constNCVariable<Vector> gDisp;
  constNCVariable<Vector> gVelocityStar;
  constNCVariable<Vector> GVelocityStar;

  // Get the old data
  if (flag->d_doGridReset) {
    new_dw->get(gVelocityStar, lb->gVelocityStarLabel, dwi, patch, gac, NGN);
    if (flag->d_fracture) {
      new_dw->get(pgCode, lb->pgCodeLabel, pset);
      new_dw->get(GVelocityStar, lb->GVelocityStarLabel, dwi, patch, gac, NGN);
    }
  } else {
    new_dw->get(gDisp, lb->gDisplacementLabel, dwi, patch, gac, NGN);
  }
  old_dw->get(px, lb->pXLabel, pset);
  old_dw->get(pMass, lb->pMassLabel, pset);
  old_dw->get(pSize, lb->pSizeLabel, pset);
  old_dw->get(pDefGrad_old, lb->pDefGradLabel, pset);
  old_dw->get(pVelGrad_old, lb->pVelGradLabel, pset);

  constParticleVariable<long64> pParticleID;
  old_dw->get(pParticleID, lb->pParticleIDLabel, pset);

  // Allocate new data
  // std::cout << "Two . Before allocate and put" << "\n";
  new_dw->allocateAndPut(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);
  new_dw->allocateAndPut(pVelGrad_new, lb->pVelGradLabel_preReloc, pset);
  new_dw->allocateAndPut(pDispGrad_new, lb->pDispGradLabel_preReloc, pset);
  new_dw->allocateAndPut(pVolume_new, lb->pVolumeLabel_preReloc, pset);
  // std::cout << "Three . After allocate and put" << "\n";

  ParticleVariable<int> pRemove_new;
  ParticleVariable<Matrix3> pPolarDecompR_new, pPolarDecompR_mid;
  if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
    new_dw->allocateAndPut(pRemove_new, lb->pRemoveLabel_preReloc, pset);
    new_dw->allocateAndPut(pPolarDecompR_new,
                           lb->pPolarDecompRLabel_preReloc,
                           pset);
    new_dw->allocateAndPut(pPolarDecompR_mid, lb->pPolarDecompRMidLabel, pset);
  }

  // Loop through particles
  double J = 1.0;
  for (auto particle : *pset) {

    // std::cout << "particle = " << particle
    //           << " px = " << px[particle]
    //           << " pSize = " << pSize[particle] << "\n";
    // Initialize variables
    Matrix3 defGrad_new(0.0); // **WARNING** should be one and not zero
    Matrix3 defGrad_inc(0.0); // **WARNING** should be one and not zero

    if (flag->d_doGridReset) {
      // Compute velocity gradient
      VelocityGradientComputer gradComp(flag);
      Matrix3 velGrad_new(0.0);
      short pgFld[27];
      if (flag->d_fracture) {
        for (int k = 0; k < 27; k++) {
          pgFld[k] = pgCode[particle][k];
        }
      }
      gradComp.computeVelGrad(interpolator.get(),
                              oodx,
                              pgFld,
                              px[particle],
                              pSize[particle],
                              pDefGrad_old[particle],
                              gVelocityStar,
                              GVelocityStar,
                              velGrad_new);
      // std::cout << "Six . After compute vel grad." << "\n";

      // Compute the deformation gradient from velocity
      computeDeformationGradientFromVelocity(pVelGrad_old[particle],
                                             velGrad_new,
                                             pDefGrad_old[particle],
                                             delT,
                                             defGrad_new,
                                             defGrad_inc);

      // std::cout << "Seven . After compute def grad." << "\n";
      //  Update velocity gradient
      pVelGrad_new[particle] = velGrad_new;

      // Update displacement gradient
      pDispGrad_new[particle] = velGrad_new * delT;

      // std::cout << "Eight . After compute disp grad." << "\n";
    } else {
      // Compute displacement gradient
      DisplacementGradientComputer gradComp(flag);
      Matrix3 dispGrad_new(0.0);
      gradComp.computeDispGrad(interpolator.get(),
                               oodx,
                               px[particle],
                               pSize[particle],
                               pDefGrad_old[particle],
                               gDisp,
                               dispGrad_new);

      // Compute the deformation gradient from displacement
      computeDeformationGradientFromTotalDisplacement(dispGrad_new,
                                                      pDefGrad_old[particle],
                                                      defGrad_new,
                                                      defGrad_inc);

      // Update displacement gradient
      pDispGrad_new[particle] = dispGrad_new;

      // Update velocity gradient
      pVelGrad_new[particle] = dispGrad_new / delT;
    }

    // Update deformation gradient
    pDefGrad_new[particle] = defGrad_new;

    // if (pParticleID[particle] == 111670263811) {
    //  std::cout << "Vel grad = " << pVelGrad_new[particle]
    //            << " Def grad = " << pDefGrad_new[particle] << "\n";
    // }

    // std::cout << "Nine . Before jacobian check" << "\n";
    //  Check 1: Look at Jacobian
    double J = defGrad_new.Determinant();
    if (!(J > 0.0)) {
      std::cerr << "matl = " << mpm_matl << " dwi = " << dwi
                << " particle = " << particle
                << " particleID = " << pParticleID[particle] << "\n";
      std::cerr << "velGrad = " << pVelGrad_new[particle] << "\n";
      std::cerr << "F_old = " << pDefGrad_old[particle] << "\n";
      std::cerr << "F_inc = " << defGrad_inc << "\n";
      std::cerr << "F_new = " << pDefGrad_new[particle] << "\n";
      std::cerr << "J = " << J << "\n";
      std::cerr
        << "**ERROR** Negative Jacobian of deformation gradient in material # ="
        << dwi << " and particle " << pParticleID[particle]
        << " which has mass " << pMass[particle] << "\n";
#ifdef IGNORE_NEGATIVE_JACOBIANS
      std::cerr
        << "\t Ignoring new deformation gradient and resetting it back to "
        << " value at the end of the previous timestep.\n";
      J                      = pDefGrad_old[particle].Determinant();
      pDefGrad_new[particle] = pDefGrad_old[particle];
/*
std::cerr << "\t Ignoring new deformation gradient and resetting it back to
I.\n"
          << "\t This action assumes that the particle has fully fractured and
\n"
          << "\t all the accumulated strain energy has been released.\n";
//J = 1;
//pDefGrad_new[particle] = Identity;
*/
#else
      throw InvalidValue("**ERROR**:", __FILE__, __LINE__);
#endif
    }

    //  Compute updated volume
    pVolume_new[particle] = (pMass[particle] / rho_orig) * J;
    // if (pParticleID[particle] == 111670263811) {
    //   std::cout << "mass = " << pMass[particle] << " vol = " <<
    //   pVolume_new[particle]
    //             << " J = " << J << "\n";
    // }

    if (!flag->d_doPressureStabilization) {
      if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
        //---------------------------------------------------------
        // Use polar decomposition to compute the rotation and stretch tensors.
        // These checks prevent
        // failure of the polar decomposition algorithm if [F_new] has some
        // extreme values.
        Matrix3 RR, UU;
        RR.Identity();
        Matrix3 FF_new  = pDefGrad_new[particle];
        double Fmax_new = FF_new.MaxAbsElem();
        double JJ_new   = FF_new.Determinant();
        if ((Fmax_new > 1.0e16) || (JJ_new < 1.0e-16) || (JJ_new > 1.0e16)) {
          pRemove_new[particle] = -999;
          proc0cout << "Deformation gradient component unphysical: [F] = "
                    << FF_new << "\n";
          proc0cout << "Resetting [F]=[I] for this step and deleting particle"
                    << " idx = " << particle
                    << " particleID = " << pParticleID[particle] << "\n";
          Vaango::Util::Identity.polarDecompositionRMB(UU, RR);
        } else {
          pRemove_new[particle] = 0;
          FF_new.polarDecompositionRMB(UU, RR);
        }
        pPolarDecompR_new[particle] = RR;
        pPolarDecompR_mid[particle] = RR;
      }
    }

    // std::cout << "Eight . Particle " << particle<< " : complete. " << "\n";
  } // End of loop over particles

  // The following is used only for pressure stabilization
  if (flag->d_doPressureStabilization) {
    CCVariable<double> J_CC, vol_0_CC, vol_CC;
    new_dw->allocateTemporary(J_CC, patch);
    new_dw->allocateTemporary(vol_0_CC, patch);
    new_dw->allocateTemporary(vol_CC, patch);
    J_CC.initialize(0.);
    vol_0_CC.initialize(0.);
    vol_CC.initialize(0.);

    // Step 1: loop thru particles and cells to compute cell centered J
    for (auto particle : *pset) {

      IntVector cell_index;
      patch->findCell(px[particle], cell_index);

      vol_CC[cell_index] += pVolume_new[particle];
      vol_0_CC[cell_index] += pMass[particle] / rho_orig;
    }

    // Compute cell centered J
    for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      J_CC[c]     = vol_CC[c] / vol_0_CC[c];
    }

    // Step 2: loop thru particles again to compute corrected def grad
    Matrix3 defGrad_inc(0.0);
    for (auto particle : *pset) {
      IntVector cell_index;
      patch->findCell(px[particle], cell_index);

      // get the original volumetric part of the deformation
      J = pDefGrad_new[particle].Determinant();

      // Change F such that the determinant is equal to the average for
      // the cell
      pDefGrad_new[particle] *= cbrt(J_CC[cell_index] / J);
      defGrad_inc = pDefGrad_new[particle] * pDefGrad_old[particle].Inverse();

      // Update the deformed volume
      J                     = pDefGrad_new[particle].Determinant();
      pVolume_new[particle] = (pMass[particle] / rho_orig) * J;

      // if (pParticleID[particle] == 111670263811) {
      //  std::cout << "Vel grad = " << pVelGrad_new[particle]
      //            << " Def grad = " << pDefGrad_new[particle] << "\n";
      // }
      // if (pParticleID[particle] == 111670263811) {
      //   std::cout << "mass = " << pMass[particle] << " vol = " <<
      //   pVolume_new[particle]
      //             << " J = " << J << "\n";
      // }

      // Check 1: Look at Jacobian
      if (!(J > 0.0)) {
        std::cerr << "after pressure stab "
                  << "\n";
        std::cerr << "matl = " << mpm_matl << "\n";
        std::cerr << "F_old = " << pDefGrad_old[particle] << "\n";
        std::cerr << "F_inc = " << defGrad_inc << "\n";
        std::cerr << "F_new = " << pDefGrad_new[particle] << "\n";
        std::cerr << "J = " << J << "\n";
        constParticleVariable<long64> pParticleID;
        old_dw->get(pParticleID, lb->pParticleIDLabel, pset);
        std::cerr << "ParticleID = " << pParticleID[particle] << "\n";
        std::cerr << "**ERROR** Negative Jacobian of deformation gradient"
                  << " in particle " << pParticleID[particle]
                  << " which has mass " << pMass[particle] << "\n";
        // pDefGrad_new[particle] = Identity;
        throw InvalidValue("**ERROR**:Negative Jacobian in UCNH",
                           __FILE__,
                           __LINE__);
      }

      if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
        //---------------------------------------------------------
        // Use polar decomposition to compute the rotation and stretch tensors.
        // These checks prevent
        // failure of the polar decomposition algorithm if [F_new] has some
        // extreme values.
        Matrix3 RR, UU;
        RR.Identity();
        Matrix3 FF_new  = pDefGrad_new[particle];
        double Fmax_new = FF_new.MaxAbsElem();
        double JJ_new   = FF_new.Determinant();
        if ((Fmax_new > 1.0e16) || (JJ_new < 1.0e-16) || (JJ_new > 1.0e16)) {
          pRemove_new[particle] = -999;
          proc0cout << "Deformation gradient component unphysical: [F] = "
                    << FF_new << "\n";
          proc0cout << "Resetting [F]=[I] for this step and deleting particle"
                    << " idx = " << particle
                    << " particleID = " << pParticleID[particle] << "\n";
          Vaango::Util::Identity.polarDecompositionRMB(UU, RR);
        } else {
          pRemove_new[particle] = 0;
          FF_new.polarDecompositionRMB(UU, RR);
        }
        pPolarDecompR_new[particle] = RR;
        pPolarDecompR_mid[particle] = RR;
      }
    }

  } // end of pressureStabilization loop  at the patch level

  return;
}

// Compute deformation gradient for implicit computations from velocity gradient
void
DeformationGradientComputer::computeDeformationGradientImplicit(
  const Patch* patch,
  const MPMMaterial* mpm_matl,
  const double& delT,
  DataWarehouse* old_dw,
  DataWarehouse* new_dw)
{
  dbg_doing << "Doing DefGrad::computeDeformationGradient implicit\n";

  // Constants
  // Ghost::GhostType  gnone = Ghost::None;
  Ghost::GhostType gac = Ghost::AroundCells;

  // Get particle info and patch info
  int dwi               = mpm_matl->getDWIndex();
  ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
  ParticleSubset* pset  = old_dw->getParticleSubset(dwi, patch);

  Vector dx      = patch->dCell();
  double oodx[3] = { 1. / dx.x(), 1. / dx.y(), 1. / dx.z() };

  // Get initial density
  double rho_orig = mpm_matl->getInitialDensity();

  // Get Interpolator
  auto interpolator = flag->d_interpolator->clone(patch);

  // Set up variables to store old particle and grid data
  // for disp grad and def grad calculation
  constParticleVariable<double> pMass;
  constParticleVariable<double> pVolume_old;
  constParticleVariable<Point> px;
  constParticleVariable<Matrix3> pSize;
  constParticleVariable<Matrix3> pDefGrad_old;

  // Set up variables to store new particle and grid data
  // for disp grad and def grad calculation
  ParticleVariable<double> pVolume_new;
  ParticleVariable<Matrix3> pDefGrad_new;
  ParticleVariable<Matrix3> pVelGrad_new;
  ParticleVariable<Matrix3> pDispGrad_new;

  // Get data from old data warehouse
  old_dw->get(pMass, lb->pMassLabel, pset);
  old_dw->get(pVolume_old, lb->pVolumeLabel, pset);
  old_dw->get(px, lb->pXLabel, pset);
  old_dw->get(pSize, lb->pSizeLabel, pset);
  old_dw->get(pDefGrad_old, lb->pDefGradLabel, pset);

  // Allocate data to new data warehouse
  new_dw->allocateAndPut(pVolume_new, lb->pVolumeLabel_preReloc, pset);
  new_dw->allocateAndPut(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);
  new_dw->allocateAndPut(pVelGrad_new, lb->pVelGradLabel_preReloc, pset);
  new_dw->allocateAndPut(pDispGrad_new, lb->pDispGradLabel_preReloc, pset);

  constParticleVariable<Matrix3> pPolarDecompR_old;
  ParticleVariable<Matrix3> pPolarDecompR_new;
  if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
    old_dw->get(pPolarDecompR_old, lb->pPolarDecompRLabel, pset);
    new_dw->allocateAndPut(pPolarDecompR_new,
                           lb->pPolarDecompRLabel_preReloc,
                           pset);
  }

  // Rigid material
  if (mpm_matl->getIsRigid()) {
    for (auto particle : *pset) {
      pVolume_new[particle]   = pVolume_old[particle];
      pDefGrad_new[particle]  = pDefGrad_old[particle];
      pVelGrad_new[particle]  = Zero;
      pDispGrad_new[particle] = Zero;
      if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
        pPolarDecompR_new[particle] = pPolarDecompR_old[particle];
      }
    }
    return;
  }

  // Deformable material
  if (flag->d_doGridReset) {
    constNCVariable<Vector> dispNew;
    new_dw->get(dispNew, lb->dispNewLabel, dwi, patch, gac, 1);

    for (auto particle : *pset) {

      // Compute incremental displacement gradient
      DisplacementGradientComputer gradComp(flag);
      Matrix3 dispGrad_new(0.0);
      gradComp.computeDispGrad(interpolator.get(),
                               oodx,
                               px[particle],
                               pSize[particle],
                               pDefGrad_old[particle],
                               dispNew,
                               dispGrad_new);
      pDispGrad_new[particle] = dispGrad_new;

      // Compute the deformation gradient from displacement
      Matrix3 defGrad_new(0.0);
      Matrix3 defGrad_inc(0.0);
      computeDeformationGradientFromIncrementalDisplacement(
        dispGrad_new,
        pDefGrad_old[particle],
        defGrad_new,
        defGrad_inc);

      // Update deformation gradient
      /*
      std::cout << "pDefGrad_old = " << pDefGrad_old[particle] << "\n"
                << "pDefGrad_new = " << defGrad_new << "\n"
                << "pDispGrad_new = " << dispGrad_new << "\n";
      */
      pDefGrad_new[particle] = defGrad_new;

      // Update velocity gradient
      Matrix3 Fdot           = (defGrad_new - pDefGrad_old[particle]) / delT;
      Matrix3 Finv           = defGrad_new.Inverse();
      pVelGrad_new[particle] = Fdot * Finv;

      //  Compute updated volume
      double J              = defGrad_new.Determinant();
      pVolume_new[particle] = (pMass[particle] / rho_orig) * J;

      if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
        Matrix3 RR, UU;
        pDefGrad_new[particle].polarDecompositionRMB(UU, RR);
        pPolarDecompR_new[particle] = RR;
      }
    }
  } else {
    constNCVariable<Vector> gDisp;
    new_dw->get(gDisp, lb->gDisplacementLabel, dwi, patch, gac, 1);

    for (auto particle : *pset) {

      // Compute total displacement gradient
      DisplacementGradientComputer gradComp(flag);
      Matrix3 dispGrad_new(0.0);
      gradComp.computeDispGrad(interpolator.get(),
                               oodx,
                               px[particle],
                               pSize[particle],
                               pDefGrad_old[particle],
                               gDisp,
                               dispGrad_new);

      pDispGrad_new[particle] = dispGrad_new;

      // Compute the deformation gradient from displacement
      Matrix3 defGrad_new(0.0);
      Matrix3 defGrad_inc(0.0);
      computeDeformationGradientFromTotalDisplacement(dispGrad_new,
                                                      pDefGrad_old[particle],
                                                      defGrad_new,
                                                      defGrad_inc);

      // Update deformation gradient
      pDefGrad_new[particle] = defGrad_new;

      // Update velocity gradient
      Matrix3 Fdot           = (defGrad_new - pDefGrad_old[particle]) / delT;
      Matrix3 Finv           = defGrad_new.Inverse();
      pVelGrad_new[particle] = Fdot * Finv;

      //  Compute updated volume
      double J              = defGrad_new.Determinant();
      pVolume_new[particle] = (pMass[particle] / rho_orig) * J;

      if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
        Matrix3 RR, UU;
        pDefGrad_new[particle].polarDecompositionRMB(UU, RR);
        pPolarDecompR_new[particle] = RR;
      }
    }
  }

  return;
}

// Task scheduling for implicit only
void
DeformationGradientComputer::addComputesAndRequires(
  Task* task,
  const MPMMaterial* matl,
  const PatchSet* patches,
  const bool /*recurse*/,
  const bool SchedParent) const
{
  dbg_doing << "Doing DefGrad::addComputesAndRequires implicit\n";

  ConstitutiveModel* cm         = matl->getConstitutiveModel();
  const MaterialSubset* matlset = matl->thisMaterial();

  Ghost::GhostType gnone = Ghost::None;
  Ghost::GhostType gac   = Ghost::AroundCells;

  if (SchedParent) {
    // For subscheduler

    printSchedule(patches,
                  dbg_doing,
                  "SchedParent::DefGrad::scheduleComputeDeformationGradient");

    task->requires(Task::ParentOldDW, lb->delTLabel);
    task->requires(Task::ParentOldDW, lb->pXLabel, matlset, gnone);
    task->requires(Task::ParentOldDW, lb->pSizeLabel, matlset, gnone);
    task->requires(Task::ParentOldDW, lb->pMassLabel, matlset, gnone);
    task->requires(Task::ParentOldDW, lb->pVolumeLabel, matlset, gnone);
    task->requires(Task::ParentOldDW, lb->pDefGradLabel, matlset, gnone);
    task->requires(Task::ParentOldDW, lb->pParticleIDLabel, matlset, gnone);

    task->computes(lb->pVelGradLabel_preReloc, matlset);
    task->computes(lb->pDispGradLabel_preReloc, matlset);
    task->computes(lb->pDefGradLabel_preReloc, matlset);
    task->computes(lb->pVolumeLabel_preReloc, matlset);
    if (flag->d_doGridReset) {
      task->requires(Task::OldDW, lb->dispNewLabel, matlset, gac, 1);
    } else {
      task->requires(Task::OldDW, lb->gDisplacementLabel, matlset, gac, 1);
    }

    if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
      task->requires(Task::ParentOldDW, lb->pPolarDecompRLabel, matlset, gnone);
      task->computes(lb->pRemoveLabel_preReloc, matlset);
      task->computes(lb->pPolarDecompRLabel_preReloc, matlset);
      task->computes(lb->pPolarDecompRMidLabel, matlset);
    }

  } else {
    printSchedule(patches,
                  dbg_doing,
                  "Recurse::DefGrad::scheduleComputeDeformationGradient");

    task->requires(Task::OldDW, lb->delTLabel);
    task->requires(Task::OldDW, lb->pXLabel, matlset, gnone);
    task->requires(Task::OldDW, lb->pSizeLabel, matlset, gnone);
    task->requires(Task::OldDW, lb->pMassLabel, matlset, gnone);
    task->requires(Task::OldDW, lb->pVolumeLabel, matlset, gnone);
    task->requires(Task::OldDW, lb->pDefGradLabel, matlset, gnone);
    task->requires(Task::OldDW, lb->pParticleIDLabel, matlset, gnone);

    if (flag->d_doGridReset) {
      task->requires(Task::NewDW, lb->dispNewLabel, matlset, gac, 1);
    } else {
      task->requires(Task::OldDW, lb->gDisplacementLabel, matlset, gac, 1);
    }

    if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
      task->requires(Task::OldDW, lb->pPolarDecompRLabel, matlset, gnone);
    }
  }
}

// Compute deformation gradient for implicit computations from velocity gradient
void
DeformationGradientComputer::computeDeformationGradientImplicit(
  const Patch* patch,
  const MPMMaterial* mpm_matl,
  const double& delT,
  DataWarehouse* old_dw,
  DataWarehouse* parent_old_dw,
  DataWarehouse* new_dw)
{
  dbg_doing << "Doing DefGrad::computeDeformationGradient parent_dw implicit\n";

  // Constants
  // Ghost::GhostType  gnone = Ghost::None;
  Ghost::GhostType gac = Ghost::AroundCells;
  // DataWarehouse* parent_new_dw =
  // new_dw->getOtherDataWarehouse(Task::ParentNewDW);

  // Get particle info and patch info
  int dwi               = mpm_matl->getDWIndex();
  ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
  ParticleSubset* pset  = parent_old_dw->getParticleSubset(dwi, patch);
  Vector dx             = patch->dCell();
  double oodx[3]        = { 1. / dx.x(), 1. / dx.y(), 1. / dx.z() };

  // Get initial density
  double rho_orig = mpm_matl->getInitialDensity();

  // Get Interpolator
  auto interpolator = flag->d_interpolator->clone(patch);

  // Set up variables to store old particle and grid data
  // for disp grad and def grad calculation
  constParticleVariable<double> pMass;
  constParticleVariable<double> pVolume_old;
  constParticleVariable<Point> px;
  constParticleVariable<Matrix3> pSize;
  constParticleVariable<Matrix3> pDefGrad_old;

  // Set up variables to store new particle and grid data
  // for disp grad and def grad calculation
  ParticleVariable<double> pVolume_new;
  ParticleVariable<Matrix3> pDefGrad_new;
  ParticleVariable<Matrix3> pVelGrad_new;
  ParticleVariable<Matrix3> pDispGrad_new;

  // Get data from parent data warehouse
  parent_old_dw->get(pMass, lb->pMassLabel, pset);
  parent_old_dw->get(pVolume_old, lb->pVolumeLabel, pset);
  parent_old_dw->get(px, lb->pXLabel, pset);
  parent_old_dw->get(pSize, lb->pSizeLabel, pset);
  parent_old_dw->get(pDefGrad_old, lb->pDefGradLabel, pset);

  // Allocate data to new data warehouse
  new_dw->allocateAndPut(pVolume_new, lb->pVolumeLabel_preReloc, pset);
  new_dw->allocateAndPut(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);
  new_dw->allocateAndPut(pVelGrad_new, lb->pVelGradLabel_preReloc, pset);
  new_dw->allocateAndPut(pDispGrad_new, lb->pDispGradLabel_preReloc, pset);

  constParticleVariable<Matrix3> pPolarDecompR_old;
  ParticleVariable<Matrix3> pPolarDecompR_new;
  if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
    parent_old_dw->get(pPolarDecompR_old, lb->pPolarDecompRLabel, pset);
    new_dw->allocateAndPut(pPolarDecompR_new,
                           lb->pPolarDecompRLabel_preReloc,
                           pset);
  }

  // Rigid material
  if (mpm_matl->getIsRigid()) {
    for (auto particle : *pset) {
      pVolume_new[particle]  = pVolume_old[particle];
      pDefGrad_new[particle] = pDefGrad_old[particle];
      // Set to zero for temporrary convenience. Fix if ever needed or remove.
      pVelGrad_new[particle]  = Zero;
      pDispGrad_new[particle] = Zero;
      if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
        pPolarDecompR_new[particle] = pPolarDecompR_old[particle];
      }
    }
    return;
  }

  // Deformable material
  if (flag->d_doGridReset) {
    constNCVariable<Vector> gDisp;
    old_dw->get(gDisp, lb->dispNewLabel, dwi, patch, gac, 1);
    for (auto particle : *pset) {

      // Compute incremental displacement gradient
      DisplacementGradientComputer gradComp(flag);
      Matrix3 dispGrad_new(0.0);
      gradComp.computeDispGrad(interpolator.get(),
                               oodx,
                               px[particle],
                               pSize[particle],
                               pDefGrad_old[particle],
                               gDisp,
                               dispGrad_new);

      pDispGrad_new[particle] = dispGrad_new;

      // Compute the deformation gradient from displacement
      Matrix3 defGrad_new(0.0);
      Matrix3 defGrad_inc(0.0);
      computeDeformationGradientFromIncrementalDisplacement(
        dispGrad_new,
        pDefGrad_old[particle],
        defGrad_new,
        defGrad_inc);

      // Update deformation gradient
      pDefGrad_new[particle] = defGrad_new;

      // Update velocity gradient
      Matrix3 Fdot           = (defGrad_new - pDefGrad_old[particle]) / delT;
      Matrix3 Finv           = defGrad_new.Inverse();
      pVelGrad_new[particle] = Fdot * Finv;

      //  Compute updated volume
      double J              = defGrad_new.Determinant();
      pVolume_new[particle] = (pMass[particle] / rho_orig) * J;

      if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
        Matrix3 RR, UU;
        pDefGrad_new[particle].polarDecompositionRMB(UU, RR);
        pPolarDecompR_new[particle] = RR;
      }
    }
  } else {
    constNCVariable<Vector> gDisp;
    old_dw->get(gDisp, lb->gDisplacementLabel, dwi, patch, gac, 1);
    for (auto particle : *pset) {

      // Compute total displacement gradient
      DisplacementGradientComputer gradComp(flag);
      Matrix3 dispGrad_new(0.0);
      gradComp.computeDispGrad(interpolator.get(),
                               oodx,
                               px[particle],
                               pSize[particle],
                               pDefGrad_old[particle],
                               gDisp,
                               dispGrad_new);

      pDispGrad_new[particle] = dispGrad_new;

      // Compute the deformation gradient from displacement
      Matrix3 defGrad_new(0.0);
      Matrix3 defGrad_inc(0.0);
      computeDeformationGradientFromTotalDisplacement(dispGrad_new,
                                                      pDefGrad_old[particle],
                                                      defGrad_new,
                                                      defGrad_inc);

      // Update deformation gradient
      pDefGrad_new[particle] = defGrad_new;

      // Update velocity gradient
      Matrix3 Fdot           = (defGrad_new - pDefGrad_old[particle]) / delT;
      Matrix3 Finv           = defGrad_new.Inverse();
      pVelGrad_new[particle] = Fdot * Finv;

      //  Compute updated volume
      double J              = defGrad_new.Determinant();
      pVolume_new[particle] = (pMass[particle] / rho_orig) * J;

      if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
        Matrix3 RR, UU;
        pDefGrad_new[particle].polarDecompositionRMB(UU, RR);
        pPolarDecompR_new[particle] = RR;
      }
    }
  }

  return;
}

void
DeformationGradientComputer::addRequiresForConvert(Task* task,
                                                   const MPMMaterial* mpm_matl)
{
  Ghost::GhostType gnone        = Ghost::None;
  const MaterialSubset* matlset = mpm_matl->thisMaterial();

  if (flag->d_integrator == MPMFlags::Implicit) {
    // Requires (for implicit)
    task->requires(Task::NewDW, lb->pVolumeLabel, matlset, gnone);
    task->requires(Task::NewDW, lb->pVelGradLabel, matlset, gnone);
    task->requires(Task::NewDW, lb->pDispGradLabel, matlset, gnone);
    task->requires(Task::NewDW, lb->pDefGradLabel, matlset, gnone);
  } else {
    // Requires (for explicit)
    task->requires(Task::NewDW, lb->pVolumeLabel, matlset, gnone);
    task->requires(Task::NewDW, lb->pVelGradLabel, matlset, gnone);
    task->requires(Task::NewDW, lb->pDispGradLabel, matlset, gnone);
    task->requires(Task::NewDW, lb->pDefGradLabel, matlset, gnone);
  }

  // **WARNING and TODO** Will not work for INCREMENTAL models unless we
  // figure out a way of determining the constitutive models for the addset and
  // delset.  Needs the sharedState to be passed.
  // ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
  // if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
  //  task->requires(Task::NewDW, lb->pRemoveLabel,                matlset,
  //  gnone); task->requires(Task::NewDW, lb->pPolarDecompRLabel_preReloc,
  //  matlset, gnone);
  //}
}

void
DeformationGradientComputer::copyAndDeleteForConvert(
  DataWarehouse* new_dw,
  ParticleSubset* addset,
  std::map<const VarLabel*, ParticleVariableBase*>* newState,
  ParticleSubset* delset,
  DataWarehouse* old_dw)
{
  // **WARNING and TODO** Will not work for INCREMENTAL models unless we
  // figure out a way of determining the constitutive models for the addset and
  // delset.  Needs the sharedState to be passed.
  // For example:
  //    auto matID_add = addset->getMatlIndex();
  //    auto matID_del = delset->getMatlIndex();
  //    MPMMaterial* mpm_matl_add = d_mat_manager->getMaterial("MPM",
  //    matID_add); ConstitutiveModel* cm_add =
  //    mpm_matl_add->getConstitutiveModel();
  //    .....

  // Copy the data common to all constitutive models from the particle to be
  // deleted to the particle to be added.

  constParticleVariable<double> o_pVolume;
  constParticleVariable<Matrix3> o_pVelGrad, o_pDispGrad, o_pDefGrad;
  new_dw->get(o_pVolume, lb->pVolumeLabel_preReloc, delset);
  new_dw->get(o_pVelGrad, lb->pVelGradLabel_preReloc, delset);
  new_dw->get(o_pDispGrad, lb->pDispGradLabel_preReloc, delset);
  new_dw->get(o_pDefGrad, lb->pDefGradLabel_preReloc, delset);

  ParticleVariable<double> pVolume;
  ParticleVariable<Matrix3> pVelGrad, pDispGrad, pDefGrad;
  new_dw->allocateTemporary(pVolume, addset);
  new_dw->allocateTemporary(pVelGrad, addset);
  new_dw->allocateTemporary(pDispGrad, addset);
  new_dw->allocateTemporary(pDefGrad, addset);

  ParticleSubset::iterator o, n = addset->begin();
  for (o = delset->begin(); o != delset->end(); o++, n++) {
    pVolume[*n]   = o_pVolume[*o];
    pVelGrad[*n]  = o_pVelGrad[*o];
    pDispGrad[*n] = o_pDispGrad[*o];
    pDefGrad[*n]  = o_pDefGrad[*o];
  }

  (*newState)[lb->pVolumeLabel]   = pVolume.clone();
  (*newState)[lb->pVelGradLabel]  = pVelGrad.clone();
  (*newState)[lb->pDispGradLabel] = pDispGrad.clone();
  (*newState)[lb->pDefGradLabel]  = pDefGrad.clone();
}

void
DeformationGradientComputer::computeDeformationGradientFromVelocity(
  const Matrix3& velGrad_old,
  const Matrix3& velGrad_new,
  const Matrix3& defGrad_old,
  const double& delT,
  Matrix3& defGrad_new,
  Matrix3& defGrad_inc)
{
  // Compute deformation gradient
  // **NOTE** Use function pointers in next iteration (TO DO)
  if (flag->d_defGradAlgorithm == "first_order") {
    seriesUpdateConstantVelGrad(velGrad_new,
                                defGrad_old,
                                delT,
                                defGrad_new,
                                defGrad_inc);
  } else if (flag->d_defGradAlgorithm == "subcycling") {
    subcycleUpdateConstantVelGrad(velGrad_new,
                                  defGrad_old,
                                  delT,
                                  defGrad_new,
                                  defGrad_inc);
  } else if (flag->d_defGradAlgorithm == "taylor_series") {
    seriesUpdateConstantVelGrad(velGrad_new,
                                defGrad_old,
                                delT,
                                defGrad_new,
                                defGrad_inc);
  } else if (flag->d_defGradAlgorithm == "cayley_hamilton") {
    cayleyUpdateConstantVelGrad(velGrad_new,
                                defGrad_old,
                                delT,
                                defGrad_new,
                                defGrad_inc);
  } else {
    std::ostringstream out;
    out << "**ERROR** Deformation gradient algorithm"
        << flag->d_defGradAlgorithm << " not implemented."
        << "\n";
    throw ParameterNotFound(out.str(), __FILE__, __LINE__);
  }
  return;
}

void
DeformationGradientComputer::computeDeformationGradientFromTotalDisplacement(
  const Matrix3& dispGrad_new,
  const Matrix3& defGrad_old,
  Matrix3& defGrad_new,
  Matrix3& defGrad_inc)
{
  // Update the deformation gradient tensor to its time n+1 value.
  // Compute the deformation gradient from the displacement gradient
  defGrad_new = Identity + dispGrad_new;
  defGrad_inc = defGrad_new * defGrad_old.Inverse();

  return;
}

// Use Taylor series expansion of exact solution
// Assume constant velocity gradient over timestep
void
DeformationGradientComputer::seriesUpdateConstantVelGrad(
  const Matrix3& velGrad_new,
  const Matrix3& defGrad_old,
  const double& delT,
  Matrix3& defGrad_new,
  Matrix3& defGrad_inc)
{
  Matrix3 Amat = velGrad_new * delT;
  defGrad_inc  = Amat.Exponential(flag->d_numTermsSeriesDefGrad);
  defGrad_new  = defGrad_inc * defGrad_old;
  // std::cout << " Amat = " << Amat << " defGrad_inc = " << defGrad_inc <<
  // "\n";
  return;
}

// Use Cayley-Hamilton theorem for series series expansion of exact solution
// Assume constant velocity gradient over timestep
void
DeformationGradientComputer::cayleyUpdateConstantVelGrad(
  const Matrix3& velGrad_new,
  const Matrix3& defGrad_old,
  const double& delT,
  Matrix3& defGrad_new,
  Matrix3& defGrad_inc)
{
  Matrix3 Id;
  Id.Identity();
  Matrix3 gradu    = velGrad_new * delT;
  Matrix3 graduSq  = gradu * gradu;
  double c0        = gradu.Determinant();
  double c2        = gradu.Trace();
  double c3        = graduSq.Trace();
  double c1        = 0.5 * (c3 - c2 * c2);
  double alpha0    = 1.0;
  double alpha1    = 1.0;
  double alpha2    = 0.5;
  double factorial = 1.0;
  // double beta0 = 0.0;
  double beta1 = 0.0;
  double beta2 = 1.0;
  for (int kk = 2; kk < flag->d_numTermsSeriesDefGrad; kk++) {
    factorial /= (double)(kk + 1);
    double beta0_new = c0 * beta2;
    double beta1_new = c1 * beta2;
    double beta2_new = c2 * beta2 + beta1;
    alpha0 += beta0_new * factorial;
    alpha1 += beta1_new * factorial;
    alpha2 += beta2_new * factorial;
    // beta0 = beta0_new;
    beta1 = beta1_new;
    beta2 = beta2_new;
  }
  defGrad_inc = Id * alpha0 + gradu * alpha1 + graduSq * alpha2;
  defGrad_new = defGrad_inc * defGrad_old;
  return;
}

// Use first term of series expansion of exact solution
// and subcycling.
// Assume constant velocity gradient over timestep
void
DeformationGradientComputer::subcycleUpdateConstantVelGrad(
  const Matrix3& velGrad_new,
  const Matrix3& defGrad_old,
  const double& delT,
  Matrix3& defGrad_new,
  Matrix3& defGrad_inc)
{
  Matrix3 Identity;
  Identity.Identity();
  defGrad_new       = defGrad_old;
  double Lnorm_dt   = velGrad_new.Norm() * delT;
  int MAX_SUBCYCLES = 1000;
  int num_scs       = std::min(std::max(1, 2 * ((int)Lnorm_dt)), MAX_SUBCYCLES);
  double dtsc       = delT / (double(num_scs));
  Matrix3 OP_tensorL_DT = Identity + velGrad_new * dtsc;
  for (int n = 0; n < num_scs; n++) {
    defGrad_new = OP_tensorL_DT * defGrad_new;
    // if(num_scs >1000){
    //   std::cerr <<  "n = " << n << "\n";
    //   std::cerr <<  "F = " << defGrad_new << "\n";
    //   std::cerr <<  "J = " << defGrad_new.Determinant() << "\n" << "\n";
    // }
  }
  defGrad_inc = defGrad_new * defGrad_old.Inverse();
  return;
}

void
DeformationGradientComputer::
  computeDeformationGradientFromIncrementalDisplacement(
    const Matrix3& dispGrad_new,
    const Matrix3& defGrad_old,
    Matrix3& defGrad_new,
    Matrix3& defGrad_inc)
{
  // Update the deformation gradient tensor to its time n+1 value.
  // Compute the deformation gradient increment
  defGrad_inc = Identity + dispGrad_new;
  defGrad_new = defGrad_inc * defGrad_old;
  return;
}
