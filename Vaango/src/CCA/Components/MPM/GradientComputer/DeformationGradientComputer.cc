/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <CCA/Components/MPM/GradientComputer/DeformationGradientComputer.h>
#include <CCA/Components/MPM/GradientComputer/VelocityGradientComputer.h>
#include <CCA/Components/MPM/GradientComputer/DisplacementGradientComputer.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Patch.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Malloc/Allocator.h>
#include <iostream>

using namespace Uintah;

const Matrix3 DeformationGradientComputer::Identity(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
const Matrix3 DeformationGradientComputer::Zero(0.0);

DeformationGradientComputer::DeformationGradientComputer(MPMFlags* Mflag, SimulationStateP& ss)
{
  lb = scinew MPMLabel();
  flag = Mflag;
  if(flag->d_8or27==8){
    NGN=1;
  } else{ 
    NGN=2;
  }
  d_sharedState = ss;
}

DeformationGradientComputer::DeformationGradientComputer(const DeformationGradientComputer* dg)
{
  lb = scinew MPMLabel();
  flag = dg->flag;
  NGN = dg->NGN;
  NGP = dg->NGP;
  d_sharedState = dg->d_sharedState;
}

DeformationGradientComputer* DeformationGradientComputer::clone()
{
  return scinew DeformationGradientComputer(*this);
}

DeformationGradientComputer::~DeformationGradientComputer()
{
  delete lb;
}

void
DeformationGradientComputer::addInitialComputesAndRequires(Task* task,
                                                           const MPMMaterial* mpm_matl,
                                                           const PatchSet*)
{
  const MaterialSubset* matlset = mpm_matl->thisMaterial();

  // Computes (for explicit)
  task->computes(lb->pVelGradLabel,  matlset);
  task->computes(lb->pDispGradLabel, matlset);
  task->computes(lb->pDefGradLabel,  matlset);
}

void
DeformationGradientComputer::addComputesAndRequires(Task* task,
                                                    const MPMMaterial* mpm_matl,
                                                    const PatchSet*)
{
  if (flag->d_integrator == MPMFlags::Implicit) {
    addComputesAndRequiresImplicit(task, mpm_matl);
  } else {
    addComputesAndRequiresExplicit(task, mpm_matl);
  }
}

void
DeformationGradientComputer::addComputesOnly(Task* task,
                                             const MPMMaterial* mpm_matl,
                                             const PatchSet*)
{
  if (flag->d_integrator == MPMFlags::Implicit) {
    std::ostringstream out;
    out << "**ERROR**: Not implemented for implicit integration" << std::endl;
    throw InvalidValue(out.str(), __FILE__, __LINE__);
  } else {
    // Computes (for explicit)
    const MaterialSubset* matlset = mpm_matl->thisMaterial();
    task->computes(lb->pVelGradLabel,  matlset);
    task->computes(lb->pDispGradLabel, matlset);
    task->computes(lb->pDefGradLabel,  matlset);
    //task->computes(lb->pVolumeLabel,   matlset);
  }
}

void 
DeformationGradientComputer::addComputesAndRequires(Task* task,
                                                    const MPMMaterial* matl,
                                                    const PatchSet* patches,
                                                    const bool /*recurse*/,
                                                    const bool SchedParent) const
{
  const MaterialSubset* matlset = matl->thisMaterial();

  Ghost::GhostType  gnone = Ghost::None;
  Ghost::GhostType  gac   = Ghost::AroundCells;

  if(SchedParent){
    // For subscheduler
    task->requires(Task::ParentOldDW, lb->pXLabel,           matlset, gnone);
    task->requires(Task::ParentOldDW, lb->pSizeLabel,        matlset, gnone);
    task->requires(Task::ParentOldDW, lb->pMassLabel,        matlset, gnone);
    task->requires(Task::ParentOldDW, lb->pVolumeLabel,      matlset, gnone);
    task->requires(Task::ParentOldDW, lb->pDefGradLabel,     matlset, gnone);
    task->requires(Task::ParentOldDW, lb->pParticleIDLabel,  matlset, gnone);

    task->computes(lb->pDefGradLabel_preReloc,     matlset);
    task->computes(lb->pVolumeLabel_preReloc,       matlset);
    if (flag->d_doGridReset) {
      task->requires(Task::OldDW, lb->dispNewLabel,       matlset, gac,1);
    } else {
      task->requires(Task::OldDW, lb->gDisplacementLabel, matlset, gac,1);
    }
  } else {
    task->requires(Task::OldDW, lb->pXLabel,           matlset, gnone);
    task->requires(Task::OldDW, lb->pSizeLabel,        matlset, gnone);
    task->requires(Task::OldDW, lb->pMassLabel,        matlset, gnone);
    task->requires(Task::OldDW, lb->pVolumeLabel,      matlset, gnone);
    task->requires(Task::OldDW, lb->pDefGradLabel,     matlset, gnone);
    task->requires(Task::OldDW, lb->pParticleIDLabel,  matlset, gnone);
  }
}

void
DeformationGradientComputer::addComputesAndRequiresExplicit(Task* task,
                                                            const MPMMaterial* mpm_matl)
{
  Ghost::GhostType  gnone = Ghost::None;
  Ghost::GhostType  gac   = Ghost::AroundCells;

  const MaterialSubset* matlset = mpm_matl->thisMaterial();

  // Requires (for explicit)
  task->requires(Task::OldDW, lb->delTLabel);
  task->requires(Task::OldDW, lb->pXLabel,                  matlset, gnone);
  task->requires(Task::OldDW, lb->pMassLabel,               matlset, gnone);
  task->requires(Task::OldDW, lb->pSizeLabel,               matlset, gnone);
  task->requires(Task::OldDW, lb->pVolumeLabel,             matlset, gnone);
  task->requires(Task::OldDW, lb->pVelocityLabel,           matlset, gnone);
  task->requires(Task::OldDW, lb->pVelGradLabel,            matlset, gnone);
  task->requires(Task::OldDW, lb->pDefGradLabel,            matlset, gnone);
  task->requires(Task::OldDW, lb->pParticleIDLabel,         matlset, gnone);
  if(flag->d_doGridReset){
    task->requires(Task::NewDW, lb->gVelocityStarLabel,     matlset, gac, NGN);
    if (flag->d_fracture) {
      task->requires(Task::NewDW, lb->pgCodeLabel,          matlset, gnone);
      task->requires(Task::NewDW, lb->GVelocityStarLabel,   matlset, gac, NGN);
    }
  } else {
    task->requires(Task::NewDW, lb->gDisplacementLabel,     matlset, gac, NGN);
  }

  // Computes (for explicit)
  task->computes(lb->pVelGradLabel_preReloc,            matlset);
  task->computes(lb->pDispGradLabel_preReloc,           matlset);
  task->computes(lb->pDefGradLabel_preReloc,            matlset);
  task->computes(lb->pVolumeLabel_preReloc,             matlset);
}

void
DeformationGradientComputer::addComputesAndRequiresImplicit(Task* task,
                                                            const MPMMaterial* mpm_matl)
{
  Ghost::GhostType  gnone = Ghost::None;
  Ghost::GhostType  gac   = Ghost::AroundCells;

  const MaterialSubset* matlset = mpm_matl->thisMaterial();

  // Requires (for implicit)
  task->requires(Task::OldDW, lb->delTLabel);
  task->requires(Task::OldDW, lb->pXLabel,                  matlset, gnone);
  task->requires(Task::OldDW, lb->pMassLabel,               matlset, gnone);
  task->requires(Task::OldDW, lb->pSizeLabel,               matlset, gnone);
  task->requires(Task::OldDW, lb->pVolumeLabel,             matlset, gnone);
  task->requires(Task::OldDW, lb->pDefGradLabel,            matlset, gnone);
  task->requires(Task::OldDW, lb->pParticleIDLabel,         matlset, gnone);
  if(flag->d_doGridReset){
    task->requires(Task::NewDW, lb->dispNewLabel,           matlset, gac, NGN);
  } else {
    task->requires(Task::NewDW, lb->gDisplacementLabel,     matlset, gac, NGN);
  }

  // Computes (for implicit)
  task->computes(lb->pDispGradLabel_preReloc,           matlset);
  task->computes(lb->pDefGradLabel_preReloc,            matlset);
  task->computes(lb->pVolumeLabel_preReloc,             matlset);
}

void
DeformationGradientComputer::addRequiresForConvert(Task* task,
                                                   const MPMMaterial* mpm_matl)
{
  Ghost::GhostType  gnone = Ghost::None;
  const MaterialSubset* matlset = mpm_matl->thisMaterial();

  if (flag->d_integrator == MPMFlags::Implicit) {
    // Requires (for implicit)
    task->requires(Task::NewDW, lb->pVolumeLabel,   matlset, gnone);
    task->requires(Task::NewDW, lb->pDefGradLabel,  matlset, gnone);
    task->requires(Task::NewDW, lb->pDispGradLabel, matlset, gnone);
  } else {
    // Requires (for explicit)
    task->requires(Task::NewDW, lb->pVolumeLabel,   matlset, gnone);
    task->requires(Task::NewDW, lb->pVelGradLabel,  matlset, gnone);
    task->requires(Task::NewDW, lb->pDefGradLabel,  matlset, gnone);
    task->requires(Task::NewDW, lb->pDispGradLabel, matlset, gnone);
  }
}

void DeformationGradientComputer::copyAndDeleteForConvert(DataWarehouse* new_dw,
                                                          ParticleSubset* addset,
                                                          map<const VarLabel*,
                                                          ParticleVariableBase*>* newState,
                                                          ParticleSubset* delset,
                                                          DataWarehouse* old_dw )
{
  // Copy the data common to all constitutive models from the particle to be 
  // deleted to the particle to be added. 
  if(flag->d_integrator != MPMFlags::Implicit){
    constParticleVariable<Matrix3>     o_pVelGrad;
    ParticleVariable<Matrix3>     pVelGrad;
    new_dw->get(o_pVelGrad, lb->pVelGradLabel_preReloc, delset);
    new_dw->allocateTemporary(pVelGrad,   addset);
    ParticleSubset::iterator o,n = addset->begin();
    for (o=delset->begin(); o != delset->end(); o++, n++) {
      pVelGrad[*n]  = o_pVelGrad[*o];
    }
    (*newState)[lb->pVelGradLabel]  = pVelGrad.clone();
  }

  constParticleVariable<double>      o_pVolume;
  constParticleVariable<Matrix3>     o_pDispGrad, o_pDefGrad;
  new_dw->get(o_pVolume,   lb->pVolumeLabel_preReloc,          delset);
  new_dw->get(o_pDefGrad,  lb->pDefGradLabel_preReloc,        delset);
  new_dw->get(o_pDispGrad, lb->pDispGradLabel_preReloc,       delset);

  ParticleVariable<double>      pVolume;
  ParticleVariable<Matrix3>     pDispGrad, pDefGrad;
  new_dw->allocateTemporary(pVolume,   addset);
  new_dw->allocateTemporary(pDispGrad, addset);
  new_dw->allocateTemporary(pDefGrad,   addset);

  ParticleSubset::iterator o,n = addset->begin();
  for (o=delset->begin(); o != delset->end(); o++, n++) {
    pVolume[*n]   = o_pVolume[*o];
    pDispGrad[*n] = o_pDispGrad[*o];
    pDefGrad[*n]  = o_pDefGrad[*o];
  }

  (*newState)[lb->pVolumeLabel]   = pVolume.clone();
  (*newState)[lb->pDefGradLabel]  = pDefGrad.clone();
  (*newState)[lb->pDispGradLabel] = pDispGrad.clone();
}

// Initial velocity/displacement/deformation gradients
void 
DeformationGradientComputer::initializeGradient(const Patch* patch,
                                                const MPMMaterial* mpm_matl,
                                                DataWarehouse* new_dw)
{
  if (flag->d_integrator == MPMFlags::Implicit) {
    initializeGradientImplicit(patch, mpm_matl, new_dw);
  } else {
    initializeGradientExplicit(patch, mpm_matl, new_dw);
  }
}

void 
DeformationGradientComputer::initializeGradientExplicit(const Patch* patch,
                                                        const MPMMaterial* mpm_matl,
                                                        DataWarehouse* new_dw)
{
  ParticleSubset* pset = new_dw->getParticleSubset(mpm_matl->getDWIndex(), patch);
  ParticleVariable<Matrix3> pVelGrad, pDispGrad, pDefGrad;
  new_dw->allocateAndPut(pVelGrad,  lb->pVelGradLabel,  pset);
  new_dw->allocateAndPut(pDispGrad, lb->pDispGradLabel, pset);
  new_dw->allocateAndPut(pDefGrad,  lb->pDefGradLabel,  pset);

  for(ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){
    particleIndex idx = *iter;
    pVelGrad[idx] = Zero;
    pDispGrad[idx] = Zero;
    pDefGrad[idx] = Identity;
  }
}

void 
DeformationGradientComputer::initializeGradientImplicit(const Patch* patch,
                                                        const MPMMaterial* mpm_matl,
                                                        DataWarehouse* new_dw)
{
  ParticleSubset* pset = new_dw->getParticleSubset(mpm_matl->getDWIndex(), patch);
  ParticleVariable<Matrix3> pDispGrad, pDefGrad;
  new_dw->allocateAndPut(pDispGrad, lb->pDispGradLabel, pset);
  new_dw->allocateAndPut(pDefGrad,  lb->pDefGradLabel,  pset);

  for(ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){
    particleIndex idx = *iter;
    pDispGrad[idx] = Zero;
    pDefGrad[idx] = Identity;
  }
}

//-----------------------------------------------------------------------
//  Actually compute deformation gradient
//-----------------------------------------------------------------------
void 
DeformationGradientComputer::computeDeformationGradient(const PatchSubset* patches,
                                                        DataWarehouse* old_dw,
                                                        DataWarehouse* new_dw)
{
  //std::cout << "Compute def grad .. 1 .." << std::endl; 
  // Get delT
  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  // The explicit code uses the velocity gradient to compute the
  // deformation gradient.  The implicit code uses displacements.
  if (flag->d_integrator == MPMFlags::Implicit) {

    // Loop thru patches
    for (int pp = 0; pp < patches->size(); pp++) {
      const Patch* patch = patches->get(pp);

      int numMPMMatls=d_sharedState->getNumMPMMatls();
      for(int m = 0; m < numMPMMatls; m++){

        // Get particle info and patch info
        MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );

        // Compute deformation gradient
        computeDeformationGradientImplicit(patch, mpm_matl, old_dw, new_dw);
      }
    }

  } else {

    // Loop thru patches
    for (int pp = 0; pp < patches->size(); pp++) {
      const Patch* patch = patches->get(pp);

      //std::cout << "Compute def grad .. 1 .. pp = " << pp << " patch = " << patch << std::endl; 
      int numMPMMatls=d_sharedState->getNumMPMMatls();
      for(int m = 0; m < numMPMMatls; m++){

        // Get particle info and patch info
        MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
        //std::cout << "Compute def grad .. 2 .. pp = " << pp << " m = " << m << " mpm_matl = " << mpm_matl << std::endl; 

        // Compute deformation gradient
        computeDeformationGradientExplicit(patch, mpm_matl, delT, old_dw, new_dw);
        //std::cout << "Compute def grad .. 3 .. complete" << std::endl; 
      }
    }

  }
}

//-----------------------------------------------------------------------
//  Actually compute deformation gradient (for implicit only)
//-----------------------------------------------------------------------
void 
DeformationGradientComputer::computeDeformationGradient(const PatchSubset* patches,
                                                        const MPMMaterial* mpm_matl,
                                                        DataWarehouse* old_dw,
                                                        DataWarehouse* new_dw,
                                                        bool recurse)
{
  std::cerr << "Fully implicit computeDeformationGradient not implemented yet." << endl;
  exit(1);
}

// Compute deformation gradient for explicit computations from velocity gradient
void
DeformationGradientComputer::computeDeformationGradientExplicit(const Patch* patch,
                                                                const MPMMaterial* mpm_matl,
                                                                const double& delT,
                                                                DataWarehouse* old_dw,
                                                                DataWarehouse* new_dw)
{
  // Constants
  //Ghost::GhostType  gnone = Ghost::None;
  Ghost::GhostType  gac   = Ghost::AroundCells;

  // Get particle info and patch info
  int dwi = mpm_matl->getDWIndex();
  ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);
  Vector dx = patch->dCell();
  double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};

  //std::cout << "One . patch = " << patch << " mpm_matl = " << mpm_matl << std::endl;

  // Get initial density
  double rho_orig = mpm_matl->getInitialDensity();

  // Get Interpolator
  ParticleInterpolator* interpolator = flag->d_interpolator->clone(patch);

  // Set up variables to store old particle and grid data 
  // for vel grad and def grad calculation
  constParticleVariable<Short27> pgCode; // For fracture MPM
  constParticleVariable<double>  pMass;
  constParticleVariable<long64>  pParticleID;
  constParticleVariable<Point>   px;
  constParticleVariable<Matrix3> pDefGrad_old, pVelGrad_old, pDispGrad_old;
  constParticleVariable<Matrix3> pSize;

  // Set up variables to store new particle and grid data 
  // for vel grad and def grad calculation
  ParticleVariable<double>     pVolume_new;
  ParticleVariable<Matrix3>    pDefGrad_new, pVelGrad_new, pDispGrad_new;
  constNCVariable<Vector>      gDisp;
  constNCVariable<Vector>      gVelocityStar;
  constNCVariable<Vector>      GVelocityStar;

  // Get the old data
  if (flag->d_doGridReset) {
    new_dw->get(gVelocityStar, lb->gVelocityStarLabel, dwi, patch, gac, NGN);
    if (flag->d_fracture) {
      new_dw->get(pgCode,        lb->pgCodeLabel, pset);
      new_dw->get(GVelocityStar, lb->GVelocityStarLabel, dwi, patch, gac, NGN);
    }
  } else {
    new_dw->get(gDisp, lb->gDisplacementLabel, dwi, patch, gac, NGN);
  }
  old_dw->get(px,                lb->pXLabel,                  pset);
  old_dw->get(pMass,             lb->pMassLabel,               pset);
  old_dw->get(pSize,             lb->pSizeLabel,               pset);
  old_dw->get(pDefGrad_old,      lb->pDefGradLabel,            pset);
  old_dw->get(pVelGrad_old,      lb->pVelGradLabel,            pset);

  // Allocate new data
  //std::cout << "Two . Before allocate and put" << std::endl;
  new_dw->allocateAndPut(pDefGrad_new,  lb->pDefGradLabel_preReloc,  pset);
  new_dw->allocateAndPut(pVelGrad_new,  lb->pVelGradLabel_preReloc,  pset);
  new_dw->allocateAndPut(pDispGrad_new, lb->pDispGradLabel_preReloc,  pset);
  new_dw->allocateAndPut(pVolume_new,   lb->pVolumeLabel_preReloc, pset);
  //std::cout << "Three . After allocate and put" << std::endl;
      
  // Loop through particles
  double J = 1.0;
  for(auto iter = pset->begin(); iter != pset->end(); iter++){
    particleIndex idx = *iter;

    // Initialize variables
    Matrix3 defGrad_new(0.0); // **WARNING** should be one and not zero
    Matrix3 defGrad_inc(0.0); // **WARNING** should be one and not zero

    if (flag->d_doGridReset) {
      // Compute velocity gradient
      VelocityGradientComputer gradComp(flag);
      Matrix3 velGrad_new(0.0);
      short pgFld[27];
      if (flag->d_fracture) {
        for(int k=0; k<27; k++) pgFld[k]=pgCode[idx][k];
      }
      gradComp.computeVelGrad(interpolator, oodx, pgFld, px[idx], pSize[idx], pDefGrad_old[idx], 
                              gVelocityStar, GVelocityStar, velGrad_new);
      //std::cout << "Six . After compute vel grad." << std::endl;

      // Compute the deformation gradient from velocity
      computeDeformationGradientFromVelocity(pVelGrad_old[idx], velGrad_new, pDefGrad_old[idx], delT, defGrad_new,
                                             defGrad_inc);

      //std::cout << "Seven . After compute def grad." << std::endl;
      // Update velocity gradient
      pVelGrad_new[idx] = velGrad_new;

      // Update displacement gradient
      pDispGrad_new[idx] = velGrad_new*delT;

      //std::cout << "Eight . After compute disp grad." << std::endl;
    } else {
      // Compute displacement gradient
      DisplacementGradientComputer gradComp(flag);
      Matrix3 dispGrad_new(0.0);
      gradComp.computeDispGrad(interpolator, oodx, px[idx], pSize[idx], 
                               pDefGrad_old[idx], gDisp, dispGrad_new);

      // Compute the deformation gradient from displacement
      computeDeformationGradientFromTotalDisplacement(dispGrad_new, pDefGrad_old[idx], defGrad_new, defGrad_inc);

      // Update displacement gradient
      pDispGrad_new[idx] = dispGrad_new;

      // Update velocity gradient
      pVelGrad_new[idx] = dispGrad_new/delT;
    }

    // Update deformation gradient
    pDefGrad_new[idx] = defGrad_new;

    // std::cout << "Vel grad = " << pVelGrad_new[idx]
    //           << " Def grad = " << pDefGrad_new[idx] << std::endl;

    //std::cout << "Nine . Before jacobian check" << std::endl;
    // Check 1: Look at Jacobian
    double J = defGrad_new.Determinant();
    if (!(J > 0.0)) {
      constParticleVariable<long64> pParticleID;
      old_dw->get(pParticleID, lb->pParticleIDLabel, pset);
      std::cerr << "matl = "  << mpm_matl << " dwi = " << dwi << " particle = " << idx
           << " particleID = " << pParticleID[idx] << endl;
      std::cerr << "velGrad = " << pVelGrad_new[idx] << endl;
      std::cerr << "F_old = " << pDefGrad_old[idx]     << endl;
      std::cerr << "F_inc = " << defGrad_inc       << endl;
      std::cerr << "F_new = " << pDefGrad_new[idx] << endl;
      std::cerr << "J = "     << J                 << endl;
      std::cerr << "**ERROR** Negative Jacobian of deformation gradient in material # ="
           << mpm_matl << " and particle " << pParticleID[idx]  << " which has mass "
           << pMass[idx] << endl;
      throw InvalidValue("**ERROR**:", __FILE__, __LINE__);
    }
  
    //  Compute updated volume
    pVolume_new[idx]=(pMass[idx]/rho_orig)*J;

    //std::cout << "Eight . Particle " << idx << " : complete. " << std::endl;
  } // End of loop over particles

  // The following is used only for pressure stabilization
  if(flag->d_doPressureStabilization) {
    CCVariable<double> J_CC, vol_0_CC, vol_CC;
    new_dw->allocateTemporary(J_CC,     patch);
    new_dw->allocateTemporary(vol_0_CC, patch);
    new_dw->allocateTemporary(vol_CC,   patch);
    J_CC.initialize(0.);
    vol_0_CC.initialize(0.);
    vol_CC.initialize(0.);
  
    // Step 1: loop thru particles and cells to compute cell centered J
    for(ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){
      particleIndex idx = *iter;

      IntVector cell_index;
      patch->findCell(px[idx],cell_index);

      vol_CC[cell_index]  +=pVolume_new[idx];
      vol_0_CC[cell_index]+=pMass[idx]/rho_orig;
    }

    // Compute cell centered J
    for(CellIterator iter=patch->getCellIterator(); !iter.done();iter++){
      IntVector c = *iter;
      J_CC[c]=vol_CC[c]/vol_0_CC[c];
    }

    // Step 2: loop thru particles again to compute corrected def grad
    Matrix3 defGrad_inc(0.0);
    for(ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){
      particleIndex idx = *iter;
      IntVector cell_index;
      patch->findCell(px[idx],cell_index);

      // get the original volumetric part of the deformation
      J = pDefGrad_new[idx].Determinant();

      // Change F such that the determinant is equal to the average for
      // the cell
      pDefGrad_new[idx]*=cbrt(J_CC[cell_index]/J);
      defGrad_inc = pDefGrad_new[idx]*pDefGrad_old[idx].Inverse();

      // Update the deformed volume
      J = pDefGrad_new[idx].Determinant();
      pVolume_new[idx]= (pMass[idx]/rho_orig)*J;

      // Check 1: Look at Jacobian
      if (!(J > 0.0)) {
        std::cerr << "after pressure stab "          << endl;
        std::cerr << "matl = "  << mpm_matl          << endl;
        std::cerr << "F_old = " << pDefGrad_old[idx]     << endl;
        std::cerr << "F_inc = " << defGrad_inc       << endl;
        std::cerr << "F_new = " << pDefGrad_new[idx] << endl;
        std::cerr << "J = "     << J                 << endl;
        constParticleVariable<long64> pParticleID;
        old_dw->get(pParticleID, lb->pParticleIDLabel, pset);
        std::cerr << "ParticleID = " << pParticleID[idx] << endl;
        std::cerr << "**ERROR** Negative Jacobian of deformation gradient"
             << " in particle " << pParticleID[idx]  << " which has mass "
             << pMass[idx] << endl;
        //pDefGrad_new[idx] = Identity;
        throw InvalidValue("**ERROR**:Negative Jacobian in UCNH", __FILE__, __LINE__);
      }
    }

  } //end of pressureStabilization loop  at the patch level

  return;
}

// Compute deformation gradient for implicit computations from velocity gradient
void
DeformationGradientComputer::computeDeformationGradientImplicit(const Patch* patch,
                                                                const MPMMaterial* mpm_matl,
                                                                DataWarehouse* old_dw,
                                                                DataWarehouse* new_dw)
{
  // Constants
  //Ghost::GhostType  gnone = Ghost::None;
  Ghost::GhostType  gac   = Ghost::AroundCells;

  // Get particle info and patch info
  int dwi = mpm_matl->getDWIndex();
  ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);
  Vector dx = patch->dCell();
  double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};

  // Get initial density
  double rho_orig = mpm_matl->getInitialDensity();

  // Get Interpolator
  ParticleInterpolator* interpolator = flag->d_interpolator->clone(patch);

  // Set up variables to store old particle and grid data 
  // for disp grad and def grad calculation
  constParticleVariable<double>  pMass;
  constParticleVariable<double>  pVolume_old;
  constParticleVariable<Point>   px;
  constParticleVariable<Matrix3> pSize;
  constParticleVariable<Matrix3> pDefGrad_old;

  // Set up variables to store new particle and grid data 
  // for disp grad and def grad calculation
  ParticleVariable<double>     pVolume_new;
  ParticleVariable<Matrix3>    pDefGrad_new;

  // Get data from old data warehouse
  old_dw->get(pMass,        lb->pMassLabel,    pset);
  old_dw->get(pVolume_old,  lb->pVolumeLabel,  pset);
  old_dw->get(px,           lb->pXLabel,       pset);
  old_dw->get(pSize,        lb->pSizeLabel,    pset);
  old_dw->get(pDefGrad_old, lb->pDefGradLabel, pset);

  // Allocate data to new data warehouse
  new_dw->allocateAndPut(pVolume_new, lb->pVolumeLabel_preReloc,  pset);
  new_dw->allocateTemporary(pDefGrad_new, pset);

  // Rigid material
  if (mpm_matl->getIsRigid()) {
    for(ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){
      particleIndex idx = *iter;
      pVolume_new[idx] = pVolume_old[idx];
      pDefGrad_new[idx] = pDefGrad_old[idx];
    }
    return;
  }

  // Deformable material
  if (flag->d_doGridReset) {
    constNCVariable<Vector> gDisp;
    old_dw->get(gDisp, lb->dispNewLabel, dwi, patch, gac, 1);
    for(ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){

      particleIndex idx = *iter;

      // Compute incremental displacement gradient
      DisplacementGradientComputer gradComp(flag);
      Matrix3 dispGrad_new(0.0);
      gradComp.computeDispGrad(interpolator, oodx, px[idx], pSize[idx], 
                               pDefGrad_old[idx], gDisp, dispGrad_new);

      // Compute the deformation gradient from displacement
      Matrix3 defGrad_new(0.0);
      Matrix3 defGrad_inc(0.0);
      computeDeformationGradientFromIncrementalDisplacement(dispGrad_new, pDefGrad_old[idx], defGrad_new, defGrad_inc);

      // Update deformation gradient
      pDefGrad_new[idx] = defGrad_new;

      //  Compute updated volume
      double J = defGrad_new.Determinant();
      pVolume_new[idx]=(pMass[idx]/rho_orig)*J;
    }
  } else {
    constNCVariable<Vector> gDisp;
    old_dw->get(gDisp, lb->gDisplacementLabel, dwi, patch, gac, 1);
    for(ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){

      particleIndex idx = *iter;

      // Compute total displacement gradient
      DisplacementGradientComputer gradComp(flag);
      Matrix3 dispGrad_new(0.0);
      gradComp.computeDispGrad(interpolator, oodx, px[idx], pSize[idx], 
                               pDefGrad_old[idx], gDisp, dispGrad_new);

      // Compute the deformation gradient from displacement
      Matrix3 defGrad_new(0.0);
      Matrix3 defGrad_inc(0.0);
      computeDeformationGradientFromTotalDisplacement(dispGrad_new, pDefGrad_old[idx], defGrad_new, defGrad_inc);

      // Update deformation gradient
      pDefGrad_new[idx] = defGrad_new;

      //  Compute updated volume
      double J = defGrad_new.Determinant();
      pVolume_new[idx]=(pMass[idx]/rho_orig)*J;
    }
  }

  return;
}

// Compute deformation gradient for implicit computations from velocity gradient
void
DeformationGradientComputer::computeDeformationGradientImplicit(const Patch* patch,
                                                                const MPMMaterial* mpm_matl,
                                                                DataWarehouse* old_dw,
                                                                DataWarehouse* parent_old_dw,
                                                                DataWarehouse* new_dw)
{
  // Constants
  //Ghost::GhostType  gnone = Ghost::None;
  Ghost::GhostType  gac   = Ghost::AroundCells;

  // Get particle info and patch info
  int dwi = mpm_matl->getDWIndex();
  ParticleSubset* pset = parent_old_dw->getParticleSubset(dwi, patch);
  Vector dx = patch->dCell();
  double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};

  // Get initial density
  double rho_orig = mpm_matl->getInitialDensity();

  // Get Interpolator
  ParticleInterpolator* interpolator = flag->d_interpolator->clone(patch);

  // Set up variables to store old particle and grid data 
  // for disp grad and def grad calculation
  constParticleVariable<double>  pMass;
  constParticleVariable<double>  pVolume_old;
  constParticleVariable<Point>   px;
  constParticleVariable<Matrix3> pSize;
  constParticleVariable<Matrix3> pDefGrad_old;

  // Set up variables to store new particle and grid data 
  // for disp grad and def grad calculation
  ParticleVariable<double>     pVolume_new;
  ParticleVariable<Matrix3>    pDefGrad_new;

  // Get data from parent data warehouse
  parent_old_dw->get(pMass,        lb->pMassLabel,    pset);
  parent_old_dw->get(pVolume_old,  lb->pVolumeLabel,  pset);
  parent_old_dw->get(px,           lb->pXLabel,       pset);
  parent_old_dw->get(pSize,        lb->pSizeLabel,    pset);
  parent_old_dw->get(pDefGrad_old, lb->pDefGradLabel, pset);

  // Allocate data to new data warehouse
  new_dw->allocateAndPut(pVolume_new, lb->pVolumeLabel_preReloc,  pset);
  new_dw->allocateTemporary(pDefGrad_new, pset);

  // Rigid material
  if (mpm_matl->getIsRigid()) {
    for(ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){
      particleIndex idx = *iter;
      pVolume_new[idx] = pVolume_old[idx];
      pDefGrad_new[idx] = pDefGrad_old[idx];
    }
    return;
  }

  // Deformable material
  if (flag->d_doGridReset) {
    constNCVariable<Vector> gDisp;
    old_dw->get(gDisp, lb->dispNewLabel, dwi, patch, gac, 1);
    for(ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){

      particleIndex idx = *iter;

      // Compute incremental displacement gradient
      DisplacementGradientComputer gradComp(flag);
      Matrix3 dispGrad_new(0.0);
      gradComp.computeDispGrad(interpolator, oodx, px[idx], pSize[idx], 
                               pDefGrad_old[idx], gDisp, dispGrad_new);

      // Compute the deformation gradient from displacement
      Matrix3 defGrad_new(0.0);
      Matrix3 defGrad_inc(0.0);
      computeDeformationGradientFromIncrementalDisplacement(dispGrad_new, pDefGrad_old[idx], defGrad_new, defGrad_inc);

      // Update deformation gradient
      pDefGrad_new[idx] = defGrad_new;

      //  Compute updated volume
      double J = defGrad_new.Determinant();
      pVolume_new[idx]=(pMass[idx]/rho_orig)*J;
    }
  } else {
    constNCVariable<Vector> gDisp;
    old_dw->get(gDisp, lb->gDisplacementLabel, dwi, patch, gac, 1);
    for(ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){

      particleIndex idx = *iter;

      // Compute total displacement gradient
      DisplacementGradientComputer gradComp(flag);
      Matrix3 dispGrad_new(0.0);
      gradComp.computeDispGrad(interpolator, oodx, px[idx], pSize[idx], 
                               pDefGrad_old[idx], gDisp, dispGrad_new);

      // Compute the deformation gradient from displacement
      Matrix3 defGrad_new(0.0);
      Matrix3 defGrad_inc(0.0);
      computeDeformationGradientFromTotalDisplacement(dispGrad_new, pDefGrad_old[idx], defGrad_new, defGrad_inc);

      // Update deformation gradient
      pDefGrad_new[idx] = defGrad_new;

      //  Compute updated volume
      double J = defGrad_new.Determinant();
      pVolume_new[idx]=(pMass[idx]/rho_orig)*J;
    }
  }

  return;
}

void 
DeformationGradientComputer::computeDeformationGradientFromVelocity(const Matrix3& velGrad_old,
                                                                    const Matrix3& velGrad_new,
                                                                    const Matrix3& defGrad_old,
                                                                    const double& delT,
                                                                    Matrix3& defGrad_new,
                                                                    Matrix3& defGrad_inc)
{
  // Compute deformation gradient
  // **NOTE** Use function pointers in next iteration (TO DO)
  if (flag->d_defgrad_algorithm == "first_order") {
    flag->d_numTermsSeriesDefGrad = 1;
    seriesUpdateConstantVelGrad(velGrad_new, defGrad_old, delT, defGrad_new, defGrad_inc);
  } else if (flag->d_defgrad_algorithm == "subcycle") {
    flag->d_numTermsSeriesDefGrad = 1;
    subcycleUpdateConstantVelGrad(velGrad_new, defGrad_old, delT, defGrad_new, defGrad_inc);
  } else if (flag->d_defgrad_algorithm == "taylor_series") {
    seriesUpdateConstantVelGrad(velGrad_new, defGrad_old, delT, defGrad_new, defGrad_inc);
  } else if (flag->d_defgrad_algorithm == "cayley_hamilton") {
    cayleyUpdateConstantVelGrad(velGrad_new, defGrad_old, delT, defGrad_new, defGrad_inc);
  } else {
    std::ostringstream out;
    out << "**ERROR** Deformation gradient algorithm"
        << flag->d_defgrad_algorithm
        << " not implemented." << std::endl;
    throw ParameterNotFound(out.str(), __FILE__, __LINE__);
  }
  return;
}

void
DeformationGradientComputer::computeDeformationGradientFromTotalDisplacement(const Matrix3& dispGrad_new,
                                                                             const Matrix3& defGrad_old,
                                                                             Matrix3& defGrad_new,
                                                                             Matrix3& defGrad_inc)
{
  // Update the deformation gradient tensor to its time n+1 value.
  // Compute the deformation gradient from the displacement gradient
  defGrad_new = Identity + dispGrad_new;
  defGrad_inc = defGrad_new*defGrad_old.Inverse();

  return;
}

// Use Taylor series expansion of exact solution
// Assume constant velocity gradient over timestep
void
DeformationGradientComputer::seriesUpdateConstantVelGrad(const Matrix3& velGrad_new,
                                                         const Matrix3& defGrad_old,
                                                         const double& delT,
                                                         Matrix3& defGrad_new,
                                                         Matrix3& defGrad_inc)
{
  Matrix3 Amat = velGrad_new*delT;
  defGrad_inc = Amat.Exponential(flag->d_numTermsSeriesDefGrad);
  defGrad_new = defGrad_inc*defGrad_old;
  // std::cout << " Amat = " << Amat << " defGrad_inc = " << defGrad_inc << std::endl;
  return;
}

// Use Cayley-Hamilton theorem for series series expansion of exact solution
// Assume constant velocity gradient over timestep
void
DeformationGradientComputer::cayleyUpdateConstantVelGrad(const Matrix3& velGrad_new,
                                                         const Matrix3& defGrad_old,
                                                         const double& delT,
                                                         Matrix3& defGrad_new,
                                                         Matrix3& defGrad_inc)
{
  Matrix3 Id; Id.Identity();
  Matrix3 gradu = velGrad_new*delT;
  Matrix3 graduSq = gradu*gradu;
  double c0 = gradu.Determinant();
  double c2 = gradu.Trace();
  double c3 = graduSq.Trace();
  double c1 = 0.5*(c3 - c2*c2);
  double alpha0 = 1.0;
  double alpha1 = 1.0;
  double alpha2 = 0.5;
  double factorial = 1.0;
  double beta0 = 0.0;
  double beta1 = 0.0;
  double beta2 = 1.0;
  for (int kk = 2; kk < flag->d_numTermsSeriesDefGrad; kk++) {
    factorial /= (double) (kk+1);
    double beta0_new = c0*beta2;
    double beta1_new = c1*beta2;
    double beta2_new = c2*beta2 + beta1;
    alpha0 += beta0_new*factorial;
    alpha1 += beta1_new*factorial;
    alpha2 += beta2_new*factorial;
    beta0 = beta0_new;
    beta1 = beta1_new;
    beta2 = beta2_new;
  }
  defGrad_inc = Id*alpha0 + gradu*alpha1 + graduSq*alpha2;
  defGrad_new = defGrad_inc*defGrad_old;
  return;
}

// Use first term of series expansion of exact solution
// and subcycling.
// Assume constant velocity gradient over timestep
void
DeformationGradientComputer::subcycleUpdateConstantVelGrad(const Matrix3& velGrad_new,
                                                           const Matrix3& defGrad_old,
                                                           const double& delT,
                                                           Matrix3& defGrad_new,
                                                           Matrix3& defGrad_inc)
{
  Matrix3 Identity; Identity.Identity();
  defGrad_new = defGrad_old;
  double Lnorm_dt = velGrad_new.Norm()*delT;
  int MAX_SUBCYCLES = 1000;
  int num_scs = std::min(std::max(1, 2*((int) Lnorm_dt)), MAX_SUBCYCLES);
  double dtsc = delT/(double (num_scs));
  Matrix3 OP_tensorL_DT = Identity + velGrad_new*dtsc;
  for(int n=0;n<num_scs;n++){
    defGrad_new = OP_tensorL_DT*defGrad_new;
    // if(num_scs >1000){
    //   cerr << "n = " << n << endl;
    //   cerr << "F = " << defGrad_new << endl;
    //   cerr << "J = " << defGrad_new.Determinant() << endl << endl;
    // }
  }
  defGrad_inc = defGrad_new*defGrad_old.Inverse();
  return;
}

void
DeformationGradientComputer::computeDeformationGradientFromIncrementalDisplacement(const Matrix3& dispGrad_new,
                                                                                   const Matrix3& defGrad_old,
                                                                                   Matrix3& defGrad_new,
                                                                                   Matrix3& defGrad_inc)
{
  // Update the deformation gradient tensor to its time n+1 value.
  // Compute the deformation gradient increment
  defGrad_inc = Identity + dispGrad_new;
  defGrad_new = defGrad_inc*defGrad_old;
  return;
}
