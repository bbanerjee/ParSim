#include <CCA/Components/Peridynamics/FamilyComputer/FamilyComputer.h>
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

//__________________________________
//  To turn on debug flags
//  csh/tcsh : setenv SCI_DEBUG "PDFamDoing:+,PDFamDebug:+".....
//  bash     : export SCI_DEBUG="PDFamDoing:+,PDFamDebug:+" )
//  default is OFF
using Uintah::DebugStream;
static DebugStream cout_doing("PDFamDoing", false);
static DebugStream cout_dbg("PDFamDebug", false);


FamilyComputer::FamilyComputer(PeridynamicsFlags* flags,
                               PeridynamicsLabel* labels)
{
  d_flags = flags;
  d_label = labels;
}

FamilyComputer::~FamilyComputer()
{
}

/*! Initial computes and requires for the family computer */
void 
FamilyComputer::addInitialComputesAndRequires(Uintah::Task* task,
                                              const PeridynamicsMaterial* matl,
                                              const Uintah::PatchSet* patches) const
{
  cout_doing << "\t Scheduling task variables in family computer: Peridynamics: " 
             << __FILE__ << ":" << __LINE__ << std::endl;

  // Identify this material
  const Uintah::MaterialSubset* matlset = matl->thisMaterial();

  // The quantities that are required by this task
  Uintah::Ghost::GhostType  gn = Uintah::Ghost::None;
  task->requires(Uintah::Task::OldDW, d_label->pPositionLabel,    gn);
  task->requires(Uintah::Task::OldDW, d_label->pParticleIDLabel,  gn);
  task->requires(Uintah::Task::OldDW, d_label->pHorizonLabel,     gn);

  // The quantities that are computed in this task
  task->computes(d_label->pNeighborListLabel, matlset);
  task->computes(d_label->pNeighborConnLabel, matlset);
  task->computes(d_label->pNeighborCountLabel, matlset);
  task->computes(d_label->pNeighborBondEnergyLabel, matlset);
}

void
FamilyComputer::createNeighborList(PeridynamicsMaterial* matl,
                                   const Uintah::Patch* patch,
                                   Uintah::DataWarehouse* new_dw)
{
  cout_doing << "\t Creating Neighbor list: Peridynamics: " << __FILE__ << ":" << __LINE__ << std::endl;

  // Get the material index in the Datawarehouse
  int matlIndex = matl->getDWIndex();

  // Get the subset of particles in this patch
  // This subset should contain all particles within the horizon of a particle.  
  // The number of extra cells to be carried around will therefore depend on the horizon
  Uintah::ParticleSubset* pset = new_dw->getParticleSubset(matlIndex, patch);

  // Get the positions and horizons of the particles from the new data warehouse
  Uintah::constParticleVariable<SCIRun::Point> pPosition;
  new_dw->get(pPosition, d_label->pPositionLabel, pset);

  Uintah::constParticleVariable<double> pHorizon;
  new_dw->get(pHorizon, d_label->pHorizonLabel, pset);

  // Get the particle IDs from the data warehouse
  Uintah::constParticleVariable<Uintah::long64> pParticleID;
  new_dw->get(pParticleID, d_label->pParticleIDLabel, pset);

  // Create allocation for the family of each particle
  Uintah::ParticleVariable<Uintah::NeighborList> pNeighborList;
  Uintah::ParticleVariable<Uintah::NeighborConnectivity> pNeighborConn;
  Uintah::ParticleVariable<int> pNeighborCount;
  Uintah::ParticleVariable<Uintah::NeighborBondEnergy> pNeighborBondEnergy;
  new_dw->allocateAndPut(pNeighborList, d_label->pNeighborListLabel, pset);
  new_dw->allocateAndPut(pNeighborConn, d_label->pNeighborConnLabel, pset);
  new_dw->allocateAndPut(pNeighborCount, d_label->pNeighborCountLabel, pset);
  new_dw->allocateAndPut(pNeighborBondEnergy, d_label->pNeighborBondEnergyLabel, pset);

  // Create a map that takes particle IDs to the array index in the particle subset
  Uintah::ParticleIDMap idMap;
  new_dw->createParticleIDMap(pset, d_label->pParticleIDLabel, idMap);

  if (cout_dbg.active()) {
    cout_dbg << "\t" << "ParticleID <-> Particle index Map for patch" << patch << std::endl;
    for (auto iter = idMap.begin(); iter != idMap.end(); iter++) {
      cout_dbg << "\t\t" << " ID = " << iter->first << " index = " << iter->second << std::endl;
    }
  }

  // Loop through the particle list
  Uintah::ParticleSubset::iterator iter = pset->begin();
  for (; iter != pset->end(); iter++) {
    particleIndex idx = *iter;

    // Find the cell in which this particle sits and its neighbors within the
    // horizon
    SCIRun::IntVector cellLow, cellHigh;
    findCellsInHorizon(patch, pPosition[idx], pHorizon[idx], cellLow, cellHigh);

    // Now that the cells have been found, create a list of particles that reside within these cells
    Uintah::ParticleSubset* pFamilySet = new_dw->getParticleSubset(matlIndex, patch, cellLow, cellHigh);

    // Get the positions and particle IDs of the particles from the new data warehouse
    Uintah::constParticleVariable<SCIRun::Point> pFamilyPos;
    Uintah::constParticleVariable<Uintah::long64> pFamilyPID;
    new_dw->get(pFamilyPos, d_label->pPositionLabel, pFamilySet);
    new_dw->get(pFamilyPID, d_label->pParticleIDLabel, pFamilySet);

    // Create a family particle vector
    std::vector<Uintah::ParticleID> family;

    // Loop through the (potential) family particle set
    Uintah::ParticleSubset::iterator pFamilyIter = pFamilySet->begin();
    for (; pFamilyIter != pFamilySet->end(); pFamilyIter++) {
      particleIndex pFamilyIdx = *pFamilyIter;

      // Calculate distance from particle to family particle
      if (pParticleID[idx] != pFamilyPID[pFamilyIdx]) {

	SCIRun::Vector bondVector = pFamilyPos[pFamilyIdx] - pPosition[idx];
        double horizon = pHorizon[idx];
        if (!(bondVector.length2() > horizon*horizon)) {
          family.push_back(pFamilyPID[pFamilyIdx]);
        }
      }
    }

    // Fill the neighbor information
    Uintah::NeighborList neighbors(family);
    Uintah::NeighborConnectivity intact;
    Uintah::NeighborBondEnergy bondEnergy;
    pNeighborList[idx] = neighbors;
    pNeighborConn[idx] = intact;
    pNeighborBondEnergy[idx] = bondEnergy;
    pNeighborCount[idx] = (int) family.size();

    // Check whether get particle index works
    // particleIndex idx_test;
    // new_dw->getParticleIndex(idMap, pParticleID[idx], idx_test);

  } // end loop over particles
}

// Compute the min max cells which intersect the horizon ball
// Algorithm:  Let p be the point and let d = (dx, dy, dz) be the diagonal of a grid cell
//             Then the unit vectors along the three coordinate directions are
//             \hat{v}_x = (1,0,0) ,  \hat{v}_y = (0,1,0) ,  \hat{v}_z = (0,0,1)
//             The vectors going from the point to the max and min cells are then
//             v_x = h  \hat{v}_x, v_y = h  \hat{v}_y, v_z = h  \hat{v}_z, 
//             The points in the min max cells are 
//             p_x_max = p + v_x , p_x_min = p - v_x 
//             p_y_max = p + v_y , p_y_min = p - v_y 
//             p_z_max = p + v_z , p_z_min = p - v_z 
// Returns: The lowest and highest cell indices in the range
void
FamilyComputer::findCellsInHorizon(const Uintah::Patch* patch,
                                   const SCIRun::Point& pos,
                                   const double& horizon,
                                   SCIRun::IntVector& cellLow,
                                   SCIRun::IntVector& cellHigh)
{
  SCIRun::Vector vx(horizon, 0.0, 0.0), vy(0.0, horizon, 0.0), vz(0.0, 0.0, horizon);

  SCIRun::Point pxMin = pos - vx;
  SCIRun::Point pyMin = pos - vy;
  SCIRun::Point pzMin = pos - vz;
  SCIRun::IntVector cellPxMin = patch->getLevel()->getCellIndex(pxMin);
  SCIRun::IntVector cellPyMin = patch->getLevel()->getCellIndex(pyMin);
  SCIRun::IntVector cellPzMin = patch->getLevel()->getCellIndex(pzMin);
  cellLow.x(cellPxMin.x());
  cellLow.y(cellPyMin.y());
  cellLow.z(cellPzMin.z());

  SCIRun::Point pxMax = pos + vx;
  SCIRun::Point pyMax = pos + vy;
  SCIRun::Point pzMax = pos + vz;
  SCIRun::IntVector cellPxMax = patch->getLevel()->getCellIndex(pxMax);
  SCIRun::IntVector cellPyMax = patch->getLevel()->getCellIndex(pyMax);
  SCIRun::IntVector cellPzMax = patch->getLevel()->getCellIndex(pzMax);
  cellHigh.x(cellPxMax.x());
  cellHigh.y(cellPyMax.y());
  cellHigh.z(cellPzMax.z());
}
