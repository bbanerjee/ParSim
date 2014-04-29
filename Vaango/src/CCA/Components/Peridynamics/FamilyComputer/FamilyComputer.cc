#include <CCA/Components/Peridynamics/FamilyComputer/FamilyComputer.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>
#include <CCA/Components/Peridynamics/FamilyComputer/NeighborList.h>

#include <CCA/Ports/DataWarehouse.h>

#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/VarLabel.h>

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
  Uintah::ParticleVariable<NeighborList> pFamily;
  new_dw->allocateAndPut(pFamily, d_labels->pFamilyLabel, pset);

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

void
FamilyComputer::findCellsInHorizon(const Uintah::Patch* patch,
                                   const SCIRun::Point& pos,
                                   const double& horizon,
                                   std::vector<SCIRun::IntVector>& cells)
{
  SCIRun::Point cellpos = d_patch->getLevel()->positionToIndex(pos);
}
