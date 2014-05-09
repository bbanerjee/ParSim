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

#include <fstream>
#include <iostream>

using namespace Vaango;

FamilyComputer::FamilyComputer(PeridynamicsFlags* flags,
                               PeridynamicsLabel* labels)
{
  d_flags = flags;
  d_labels = labels;
}

FamilyComputer::~FamilyComputer()
{
}

void
FamilyComputer::createNeighborList(PeridynamicsMaterial* matl,
                                   const Uintah::Patch* patch,
                                   Uintah::DataWarehouse* new_dw)
{
  // Get the material index in the Datawarehouse
  int dwi = matl->getDWIndex();

  // Get the subset of particles in this patch
  // This subset should contain all particles within the horizon of a particle.  
  // The number of extra cells to be carried around will therefore depend on the horizon
  Uintah::ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);

  // Get the positions and horizons of the particles from the new data warehouse
  Uintah::constParticleVariable<SCIRun::Point> pPosition;
  new_dw->get(pPosition, d_labels->pPositionLabel, pset);

  Uintah::constParticleVariable<double> pHorizon;
  new_dw->get(pHorizon, d_labels->pHorizonLabel, pset);

  // Create allocation for the family of each particle
  Uintah::ParticleVariable<Uintah::NeighborList> pNeighorList;
  Uintah::ParticleVariable<Uintah::NeighborConnectivity> pNeighorConn;
  Uintah::ParticleVariable<int> pNeighorCount;
  new_dw->allocateAndPut(pNeighorList, d_labels->pNeighborListLabel, pset);
  new_dw->allocateAndPut(pNeighorConn, d_labels->pNeighborConnLabel, pset);
  new_dw->allocateAndPut(pNeighorCount, d_labels->pNeighborCountLabel, pset);

  // Loop through the particle list
  Uintah::ParticleSubset::iterator iter = pset->begin();
  for (; iter != pset->end(); iter++) {
    particleIndex idx = *iter;

    // Find the cell in which this particle sits and its neighbors within the
    // horizon
    std::vector<SCIRun::IntVector> cells;
    findCellsInHorizon(patch, pPosition[idx], pHorizon[idx], cells);
  }
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
// Returns: a vector of cells in the following order:
//            cell(p), cell(p_x_max), cell(p_y_max), cell(p_z_max),
//            cell(p_x_min), cell(p_y_min), cell(p_z_min)
void
FamilyComputer::findCellsInHorizon(const Uintah::Patch* patch,
                                   const SCIRun::Point& pos,
                                   const double& horizon,
                                   std::vector<SCIRun::IntVector>& cells)
{
  // Find the cell in which the particle resides
  SCIRun::IntVector cellIndex = patch->getLevel()->getCellIndex(pos);
  cells.push_back(cellIndex);

  SCIRun::Vector vx(horizon, 0.0, 0.0), vy(0.0, horizon, 0.0), vz(0.0, 0.0, horizon);
  SCIRun::Point pxMax = pos + vx;
  SCIRun::Point pyMax = pos + vy;
  SCIRun::Point pzMax = pos + vz;
  SCIRun::IntVector cellPxMax = patch->getLevel()->getCellIndex(pxMax);
  SCIRun::IntVector cellPyMax = patch->getLevel()->getCellIndex(pyMax);
  SCIRun::IntVector cellPzMax = patch->getLevel()->getCellIndex(pzMax);
  cells.push_back(cellPxMax);
  cells.push_back(cellPyMax);
  cells.push_back(cellPzMax);
  SCIRun::Point pxMin = pos - vx;
  SCIRun::Point pyMin = pos - vy;
  SCIRun::Point pzMin = pos - vz;
  SCIRun::IntVector cellPxMin = patch->getLevel()->getCellIndex(pxMin);
  SCIRun::IntVector cellPyMin = patch->getLevel()->getCellIndex(pyMin);
  SCIRun::IntVector cellPzMin = patch->getLevel()->getCellIndex(pzMin);
  cells.push_back(cellPxMin);
  cells.push_back(cellPyMin);
  cells.push_back(cellPzMin);
}
