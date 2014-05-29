#include <CCA/Components/Peridynamics/InternalForceComputer/ParticleInternalForceComputer.h>
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
//  csh/tcsh : setenv SCI_DEBUG "PDPIntForceDoing:+,PDPIntForceDebug:+".....
//  bash     : export SCI_DEBUG="PDPIntForceDoing:+,PDPIntForceDebug:+" )
//  default is OFF
using Uintah::DebugStream;
static DebugStream cout_doing("PDPIntForceDoing", false);
static DebugStream cout_dbg("PDPIntForceDebug", false);


ParticleInternalForceComputer::ParticleInternalForceComputer(PeridynamicsFlags* flags,
                                                             PeridynamicsLabel* labels)
{
  d_flags = flags;
  d_label = labels;
}

ParticleInternalForceComputer::~ParticleInternalForceComputer()
{
}

/*! Computes and requires for the internal force computer */
void 
ParticleInternalForceComputer::addComputesAndRequires(Task* task,
                                                      const PeridynamicsMaterial* matl,
                                                      const PatchSet* patches) const
{
  cout_doing << "\t Scheduling task variables in particle internal force computer: Peridynamics: " 
             << __FILE__ << ":" << __LINE__ << std::endl;

  // Constants
  Ghost::GhostType gac = Ghost::AroundCells;
  int numGhostCells = d_flags->d_numCellsInHorizon;
  
  // Get the current material
  const MaterialSubset* matlset = matl->thisMaterial();
  
  // List the variables needed for this task to execute
  task->requires(Task::OldDW, d_label->pParticleIDLabel,         matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_label->pVolumeLabel,             matlset, gac, numGhostCells);

  task->requires(Task::OldDW, d_label->pNeighborListLabel,               matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_label->pNeighborConnLabel,               matlset, gac, numGhostCells);
  task->requires(Task::OldDW, d_label->pNeighborCountLabel,              matlset, gac, numGhostCells);
  task->requires(Task::NewDW, d_label->pNeighborBondForceLabel_preReloc, matlset, gac, numGhostCells);

  // List the variables computed by this task
  task->computes(d_label->pInternalForceLabel_preReloc, matlset);

}

void
ParticleInternalForceComputer::computeInternalForce(const PatchSubset* patches,
                                                    const PeridynamicsMaterial* matl,
                                                    DataWarehouse* old_dw,
                                                    DataWarehouse* new_dw)
{
  cout_doing << "\t Computing particle internal force: Peridynamics: " 
             << __FILE__ << ":" << __LINE__ << std::endl;

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
    constParticleVariable<long64> pParticleID;
    old_dw->get(pParticleID, d_label->pParticleIDLabel, pset);

    constParticleVariable<double> pVolume_family;
    old_dw->get(pVolume_family, d_label->pVolumeLabel, familySet);

    constParticleVariable<int> pNeighborCount, pNeighborCount_family;
    old_dw->get(pNeighborCount, d_label->pNeighborCountLabel, pset);
    old_dw->get(pNeighborCount_family, d_label->pNeighborCountLabel, familySet);

    constParticleVariable<NeighborList> pNeighborList, pNeighborList_family;
    old_dw->get(pNeighborList, d_label->pNeighborListLabel, pset);
    old_dw->get(pNeighborList_family, d_label->pNeighborListLabel, familySet);

    constParticleVariable<NeighborConnectivity> pNeighborConn;
    old_dw->get(pNeighborConn, d_label->pNeighborConnLabel, pset);

    constParticleVariable<NeighborBondInternalForce> pNeighborBondForce_new, pNeighborBondForce_new_family;
    new_dw->get(pNeighborBondForce_new, d_label->pNeighborBondForceLabel_preReloc, pset);
    new_dw->get(pNeighborBondForce_new_family, d_label->pNeighborBondForceLabel_preReloc, familySet);

    // Initialize variables that will be updated in this task
    ParticleVariable<Vector> pInternalForce_new;
    new_dw->allocateAndPut(pInternalForce_new, d_label->pInternalForceLabel_preReloc, pset);

    // Loop through particles
    for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++) {

      // Get particle index
      particleIndex idx = *iter;

      // Initialize the internal force
      Vector internal_force(0.0);

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

          // Get the volume of the neighbor
          double volume = pVolume_family[family_idx];

          // Get the bond internal force
          Vector bond_int_force = (pNeighborBondForce_new[idx])[ii];

          // Add force*vol to internal force 
          internal_force += (bond_int_force*volume);

          // Now loop through the neighbor's family and find the bond to
          // this particle
          int neighborFamilyCount = pNeighborCount_family[family_idx];
          for (int jj = 0; jj < neighborFamilyCount; jj++) {

            // If the current particle is in the list
            if ((pNeighborList_family[family_idx])[jj] == pParticleID[idx]) {

              // Get the bond internal force
              Vector reverse_bond_int_force = (pNeighborBondForce_new_family[family_idx])[jj];

              // Subtract force*vol from internal force 
              internal_force -= (reverse_bond_int_force*volume);
             
              // Break out of the jj loop
              break;
            }

          } // end loop through neighbor's family

        } // Endif bond connected/broken
      }  // End neighbor particle loop

      // Update particle internal force
      pInternalForce_new[idx] = internal_force;

    }  // End particle loop
  }  // End patch loop
}

