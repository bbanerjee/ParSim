/*
 * The MIT License
 *
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

#include <CCA/Components/Peridynamics/FamilyComputer/FamilyComputer.h>

#include <CCA/Components/Peridynamics/Core/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsMaterial.h>

#include <Core/Grid/Variables/VarLabel.h>

#include <Core/Grid/Variables/NeighborConnectivity.h>
#include <Core/Grid/Variables/NeighborList.h>

#include <Core/Util/DebugStream.h>

#include <fstream>
#include <iostream>

using namespace Vaango;

using Uintah::constParticleVariable;
using Uintah::DataWarehouse;
using Uintah::Ghost;
using Uintah::IntVector;
using Uintah::long64;
using Uintah::MaterialSubset;
using Uintah::NeighborBondEnergy;
using Uintah::NeighborBondInternalForce;
using Uintah::NeighborConnectivity;
using Uintah::NeighborList;
using Uintah::ParticleSubset;
using Uintah::ParticleVariable;
using Uintah::Patch;
using Uintah::PatchSet;
using Uintah::Point;
using Uintah::Task;
using Uintah::Vector;

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

FamilyComputer::~FamilyComputer() {}

/*! Initial computes and requires for the family computer */
void
FamilyComputer::addInitialComputesAndRequires(Task* task,
                                              const PeridynamicsMaterial* matl,
                                              [[maybe_unused]] const PatchSet* patches) const
{
  cout_doing
    << "\t Scheduling task variables in family computer: Peridynamics: "
    << __FILE__ << ":" << __LINE__ << std::endl;

  // Identify this material
  const MaterialSubset* matlset = matl->thisMaterial();

  // The quantities that are required by this task
  Ghost::GhostType gac = Ghost::AroundCells;
  task->needs(
    Task::NewDW, d_label->pPositionLabel, gac, d_flags->d_numCellsInHorizon);
  task->needs(
    Task::NewDW, d_label->pParticleIDLabel, gac, d_flags->d_numCellsInHorizon);
  task->needs(
    Task::NewDW, d_label->pHorizonLabel, gac, d_flags->d_numCellsInHorizon);

  // The quantities that are computed in this task
  task->computes(d_label->pNeighborListLabel, matlset);
  task->computes(d_label->pNeighborConnLabel, matlset);
  task->computes(d_label->pNeighborCountLabel, matlset);
  task->computes(d_label->pNeighborBondEnergyLabel, matlset);
  task->computes(d_label->pNeighborBondForceLabel, matlset);
}

void
FamilyComputer::createNeighborList(PeridynamicsMaterial* matl,
                                   const Patch* patch,
                                   DataWarehouse* new_dw)
{
  cout_doing << "\t Creating Neighbor list: Peridynamics: " << __FILE__ << ":"
             << __LINE__ << std::endl;

  // Get the material index in the Datawarehouse
  int matlIndex = matl->getDWIndex();

  // Get the subset of particles in this patch
  // This subset should contain all particles within the horizon of a particle.
  // The number of extra cells to be carried around will therefore depend on the
  // horizon
  ParticleSubset* pset = new_dw->getParticleSubset(matlIndex, patch);

  // Get the positions and horizons of the particles from the new data warehouse
  constParticleVariable<Point> pPosition;
  new_dw->get(pPosition, d_label->pPositionLabel, pset);

  constParticleVariable<double> pHorizon;
  new_dw->get(pHorizon, d_label->pHorizonLabel, pset);

  // Get the particle IDs from the data warehouse
  constParticleVariable<long64> pParticleID;
  new_dw->get(pParticleID, d_label->pParticleIDLabel, pset);

  // Create allocation for the family of each particle
  ParticleVariable<NeighborList> pNeighborList;
  ParticleVariable<NeighborConnectivity> pNeighborConn;
  ParticleVariable<int> pNeighborCount;
  ParticleVariable<NeighborBondEnergy> pNeighborBondEnergy;
  ParticleVariable<NeighborBondInternalForce> pNeighborBondForce;
  new_dw->allocateAndPut(pNeighborList, d_label->pNeighborListLabel, pset);
  new_dw->allocateAndPut(pNeighborConn, d_label->pNeighborConnLabel, pset);
  new_dw->allocateAndPut(pNeighborCount, d_label->pNeighborCountLabel, pset);
  new_dw->allocateAndPut(
    pNeighborBondEnergy, d_label->pNeighborBondEnergyLabel, pset);
  new_dw->allocateAndPut(
    pNeighborBondForce, d_label->pNeighborBondForceLabel, pset);

  // Create a map that takes particle IDs to the array index in the particle
  // subset
  Uintah::ParticleIDMap idMap;
  new_dw->createParticleIDMap(pset, d_label->pParticleIDLabel, idMap);

  if (cout_dbg.active()) {
    cout_dbg << "\t"
             << "ParticleID <-> Particle index Map for patch " << patch
             << std::endl;
    for (auto iter = idMap.begin(); iter != idMap.end(); iter++) {
      cout_dbg << "\t\t"
               << " ID = " << iter->first << " index = " << iter->second
               << std::endl;
    }
  }

  // Loop through the particle list
  for (auto idx : *pset) {

    // Find the cell in which this particle sits and its neighbors within the
    // horizon
    IntVector cellLow, cellHigh;
    findCellsInHorizon(patch, pPosition[idx], pHorizon[idx], cellLow, cellHigh);

    if (cout_dbg.active()) {
      cout_dbg << "\t Particle index = " << idx
               << " position = " << pPosition[idx]
               << " horizon = " << pHorizon[idx] << std::endl;
      cout_dbg << "\t Cells in horizon: low = " << cellLow
               << " high = " << cellHigh << std::endl;
    }

    // Now that the cells have been found, create a list of particles that
    // reside within these cells
    ParticleSubset* pFamilySet =
      new_dw->getParticleSubset(matlIndex,
                                patch,
                                cellLow,
                                cellHigh,
                                Ghost::AroundCells,
                                d_flags->d_numCellsInHorizon,
                                d_label->pPositionLabel);

    cout_doing << "\t Got family particle set " << pFamilySet << " " << __FILE__
               << ":" << __LINE__ << std::endl;

    // Get the positions and particle IDs of the particles from the new data
    // warehouse
    constParticleVariable<Point> pFamilyPos;
    constParticleVariable<long64> pFamilyPID;
    new_dw->get(pFamilyPos, d_label->pPositionLabel, pFamilySet);
    new_dw->get(pFamilyPID, d_label->pParticleIDLabel, pFamilySet);

    // Create a family particle vector
    std::vector<Uintah::ParticleID> family;
    std::vector<bool> connected;
    std::vector<double> energy;
    std::vector<Uintah::Vector> force;

    // Loop through the (potential) family particle set
    for (auto pFamilyIdx : *pFamilySet) {

      // Calculate distance from particle to family particle
      if (pParticleID[idx] != pFamilyPID[pFamilyIdx]) {

        Vector bondVector = pFamilyPos[pFamilyIdx] - pPosition[idx];
        double horizon    = pHorizon[idx];
        if (!(bondVector.length2() > horizon * horizon)) {
          family.push_back(pFamilyPID[pFamilyIdx]);
          connected.push_back(true);
          energy.push_back(0.0);
          force.push_back(Vector(0.0, 0.0, 0.0));
        }
      }
      cout_doing << "\t\t Family particle  " << pFamilyIdx << " " << __FILE__
                 << ":" << __LINE__ << std::endl;
    }

    // Fill the neighbor information
    NeighborList neighbors;
    NeighborConnectivity intact;
    NeighborBondEnergy bondEnergy;
    NeighborBondInternalForce bondForce;

    int ii = 0;
    for (auto iter = family.begin(); iter != family.end(); iter++) {
      neighbors[ii]  = family[ii];
      intact[ii]     = connected[ii];
      bondEnergy[ii] = energy[ii];
      bondForce[ii]  = force[ii];
      ii++;
    }

    pNeighborList[idx]       = neighbors;
    pNeighborConn[idx]       = intact;
    pNeighborBondEnergy[idx] = bondEnergy;
    pNeighborBondForce[idx]  = bondForce;
    pNeighborCount[idx]      = (int)family.size();

    cout_doing << "\t\t Completed neighbor fill for particle " << idx << " "
               << __FILE__ << ":" << __LINE__ << std::endl;

    if (cout_dbg.active()) {
      cout_dbg << "\t\t\t  Neighbor count for particle  " << idx << " = "
               << pNeighborCount[idx] << std::endl;
      cout_dbg << "\t\t\t  Neighbor list for particle  " << idx << " = "
               << pNeighborList[idx] << std::endl;
      cout_dbg << "\t\t\t  Neighbor conn for particle  " << idx << " = "
               << pNeighborConn[idx] << std::endl;
      cout_dbg << "\t\t\t  Neighbor bond energy for particle  " << idx << " = "
               << pNeighborBondEnergy[idx] << std::endl;
      cout_dbg << "\t\t\t  Neighbor bond force for particle  " << idx << " = "
               << pNeighborBondForce[idx] << std::endl;
    }

    // Check whether get particle index works
    // particleIndex idx_test;
    // new_dw->getParticleIndex(idMap, pParticleID[idx], idx_test);

  } // end loop over particles
}

// Compute the min max cells which intersect the horizon ball
// Algorithm:  Let p be the point and let d = (dx, dy, dz) be the diagonal of a
// grid cell
//             Then the unit vectors along the three coordinate directions are
//             \hat{v}_x = (1,0,0) ,  \hat{v}_y = (0,1,0) ,  \hat{v}_z = (0,0,1)
//             The vectors going from the point to the max and min cells are
//             then v_x = h  \hat{v}_x, v_y = h  \hat{v}_y, v_z = h  \hat{v}_z,
//             The points in the min max cells are
//             p_x_max = p + v_x , p_x_min = p - v_x
//             p_y_max = p + v_y , p_y_min = p - v_y
//             p_z_max = p + v_z , p_z_min = p - v_z
// Returns: The lowest and highest cell indices in the range
void
FamilyComputer::findCellsInHorizon(const Patch* patch,
                                   const Point& pos,
                                   const double& horizon,
                                   IntVector& cellLow,
                                   IntVector& cellHigh)
{
  // Get the index range for this patch
  IntVector lowIndex =
    patch->getExtraCellLowIndex(d_flags->d_numCellsInHorizon);
  IntVector highIndex =
    patch->getExtraCellHighIndex(d_flags->d_numCellsInHorizon);
  if (cout_dbg.active()) {
    cout_dbg << "\t Patch: Low index = " << lowIndex
             << " High index = " << highIndex << std::endl;
  }

  // Set the horizon range vectors
  Vector vx(horizon, 0.0, 0.0), vy(0.0, horizon, 0.0), vz(0.0, 0.0, horizon);

  Point pxMin         = pos - vx;
  Point pyMin         = pos - vy;
  Point pzMin         = pos - vz;
  IntVector cellPxMin = patch->getLevel()->getCellIndex(pxMin);
  IntVector cellPyMin = patch->getLevel()->getCellIndex(pyMin);
  IntVector cellPzMin = patch->getLevel()->getCellIndex(pzMin);

  // If the computed min are outside the patch limits, set them to the
  // limits
  cellLow.x(cellPxMin.x() < lowIndex.x() ? lowIndex.x() : cellPxMin.x());
  cellLow.y(cellPyMin.y() < lowIndex.y() ? lowIndex.y() : cellPxMin.y());
  cellLow.z(cellPzMin.z() < lowIndex.z() ? lowIndex.z() : cellPxMin.z());

  Point pxMax         = pos + vx;
  Point pyMax         = pos + vy;
  Point pzMax         = pos + vz;
  IntVector cellPxMax = patch->getLevel()->getCellIndex(pxMax);
  IntVector cellPyMax = patch->getLevel()->getCellIndex(pyMax);
  IntVector cellPzMax = patch->getLevel()->getCellIndex(pzMax);

  // If the computed max are outside the patch limits, set them to the
  // limits
  cellHigh.x(cellPxMax.x() > highIndex.x() ? highIndex.x() : cellPxMax.x());
  cellHigh.y(cellPyMax.y() > highIndex.y() ? highIndex.y() : cellPyMax.y());
  cellHigh.z(cellPzMax.z() > highIndex.z() ? highIndex.z() : cellPzMax.z());
}
