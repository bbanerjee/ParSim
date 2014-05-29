#include <CCA/Components/Peridynamics/InternalForceComputer/BondInternalForceComputer.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>

#include <CCA/Ports/DataWarehouse.h>

#include <Core/Grid/Patch.h>

#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/NeighborList.h>
#include <Core/Grid/Variables/NeighborConnectivity.h>

#include <Core/Util/DebugStream.h>

#include <fstream>
#include <iostream>

using namespace Vaango;

using Uintah::Task;
using Uintah::Patch;
using Uintah::PatchSet;
using Uintah::PatchSubset;
using Uintah::DataWarehouse;
using Uintah::ParticleSubset;
using Uintah::MaterialSubset;
using Uintah::Ghost;
using Uintah::ParticleVariable;
using Uintah::constParticleVariable;
using Uintah::long64;
using Uintah::Matrix3;
using Uintah::NeighborList;
using Uintah::NeighborConnectivity;
using Uintah::NeighborBondInternalForce;
using SCIRun::Point;
using SCIRun::IntVector;
using SCIRun::Vector;

//__________________________________
//  To turn on debug flags
//  csh/tcsh : setenv SCI_DEBUG "PDIntForceDoing:+,PDIntForceDebug:+".....
//  bash     : export SCI_DEBUG="PDIntForceDoing:+,PDIntForceDebug:+" )
//  default is OFF
using Uintah::DebugStream;
static DebugStream cout_doing("PDIntForceDoing", false);
static DebugStream cout_dbg("PDIntForceDebug", false);


BondInternalForceComputer::BondInternalForceComputer(PeridynamicsFlags* flags,
                                                     PeridynamicsLabel* labels)
{
  d_flags = flags;
  d_label = labels;
}

BondInternalForceComputer::~BondInternalForceComputer()
{
}

/*! Computes and requires for the internal force computer */
void 
BondInternalForceComputer::addComputesAndRequires(Task* task,
                                                  const PeridynamicsMaterial* matl,
                                                  const PatchSet* patches) const
{
  cout_doing << "\t Scheduling task variables in bond internal force computer: Peridynamics: " 
             << __FILE__ << ":" << __LINE__ << std::endl;

  // Constants
  Ghost::GhostType gac = Ghost::AroundCells;
  int numGhostCells = d_flags->d_numCellsInHorizon;
  
  // Get the current material
  const MaterialSubset* matlset = matl->thisMaterial();
  
  // List the variables needed for this task to execute
  task->requires(Task::OldDW, d_label->pParticleIDLabel,         matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_label->pPositionLabel,           matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_label->pDisplacementLabel,       matlset, gac, numGhostCells);
  task->requires(Task::NewDW, d_label->pPK1StressLabel_preReloc,      matlset, gac, numGhostCells);
  task->requires(Task::NewDW, d_label->pShapeTensorInvLabel_preReloc, matlset, gac, numGhostCells);

  task->requires(Task::OldDW, d_label->pNeighborListLabel,       matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_label->pNeighborConnLabel,       matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_label->pNeighborCountLabel,      matlset, gac, numGhostCells);

  // List the variables computed by this task
  task->computes(d_label->pNeighborBondForceLabel_preReloc, matlset);

}

void
BondInternalForceComputer::computeInternalForce(const PatchSubset* patches,
                                                const PeridynamicsMaterial* matl,
                                                DataWarehouse* old_dw,
                                                DataWarehouse* new_dw)
{
  cout_doing << "\t Computing bond internal force: Peridynamics: " 
             << __FILE__ << ":" << __LINE__ << std::endl;

  // Assume influence function is always 1 (**TODO** Compute it using a factory.)
  double pInfluence = 1.0;

  // Loop through patches
  for (int pp = 0; pp < patches->size(); pp++) {

    // Get the current patch
    const Patch* patch = patches->get(pp);

    // Get the index of current material in the data warehouse
    int matlIndex = matl->getDWIndex();

    // Get the particle subset for this material in this patch
    ParticleSubset* pset = old_dw->getParticleSubset(matlIndex, patch);

    // Get the particle subset for the neighbor particles of the same material in this patch + ghost regions
    ParticleSubset* familySet = old_dw->getParticleSubset(matlIndex, patch, Ghost::AroundCells,
                                                          d_flags->d_numCellsInHorizon,
                                                          d_label->pPositionLabel);

    // Create a map that takes particle IDs to the array index in the larger particle subset
    Uintah::ParticleIDMap familyIdMap;
    old_dw->createParticleIDMap(familySet, d_label->pParticleIDLabel, familyIdMap);

    // Get the particle data needed for damage computation
    constParticleVariable<Point> pPosition, pPosition_family;
    old_dw->get(pPosition,        d_label->pPositionLabel, pset);
    old_dw->get(pPosition_family, d_label->pPositionLabel, familySet);

    constParticleVariable<Vector> pDisp, pDisp_family;
    old_dw->get(pDisp,        d_label->pDisplacementLabel, pset);
    old_dw->get(pDisp_family, d_label->pDisplacementLabel, familySet);

    constParticleVariable<Matrix3> pPK1Stress_new;
    new_dw->get(pPK1Stress_new, d_label->pPK1StressLabel_preReloc, pset);

    constParticleVariable<Matrix3> pShapeInv_new;
    new_dw->get(pShapeInv_new, d_label->pShapeTensorInvLabel_preReloc, pset);

    constParticleVariable<int> pNeighborCount;
    old_dw->get(pNeighborCount, d_label->pNeighborCountLabel, pset);

    constParticleVariable<NeighborList> pNeighborList;
    old_dw->get(pNeighborList, d_label->pNeighborListLabel, pset);

    constParticleVariable<NeighborConnectivity> pNeighborConn;
    old_dw->get(pNeighborConn, d_label->pNeighborConnLabel, pset);

    // Initialize variables that will be updated in this task
    ParticleVariable<NeighborBondInternalForce> pNeighborBondForce_new;
    new_dw->allocateAndPut(pNeighborBondForce_new, d_label->pNeighborBondForceLabel_preReloc, pset);

    // Loop through particles
    for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++) {

      // Get particle index
      particleIndex idx = *iter;

      // Find the reference position of the current particle
      Point cur_pos = pPosition[idx];
      Vector cur_disp = pDisp[idx];
      Point cur_ref_pos = cur_pos - cur_disp;

      // Find the first PK stress and shape tensor inverse
      Matrix3 cur_PK1_stress = pPK1Stress_new[idx];
      Matrix3 cur_shape_inv = pShapeInv_new[idx];

      // Compute the internal force state
      Matrix3 cur_force_state = (cur_PK1_stress*cur_shape_inv)*pInfluence;

      // Get the neighbor data
      NeighborList family = pNeighborList[idx];
      NeighborConnectivity connected = pNeighborConn[idx];

      // Loop through the neighbor list
      int neighborCount = pNeighborCount[idx];
      for (int ii = 0; ii < neighborCount; ii++) {
       
        // If the bond exists
        if (connected[ii]) {

          // Find the idx of the neighbor
	  Uintah::ParticleID familyID = family[ii];
          particleIndex family_idx;
          old_dw->getParticleIndex(familyIdMap, familyID, family_idx);

          // Find the reference position of the neighbor
          Point family_pos = pPosition_family[family_idx];
          Vector family_disp = pDisp_family[family_idx];
          Point family_ref_pos = family_pos - family_disp;

          // Compute relative reference position xi = x_family - x
          Vector xi = family_ref_pos - cur_ref_pos;

          // Update the bond force
          Vector bond_force_new = cur_force_state*xi;
          (pNeighborBondForce_new[idx])[ii] = bond_force_new;

        } // Endif bond connected/broken
      }  // End neighbor particle loop

    }  // End particle loop
  }  // End patch loop
}

