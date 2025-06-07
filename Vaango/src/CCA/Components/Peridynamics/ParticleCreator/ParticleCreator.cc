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

#include <CCA/Components/Peridynamics/ParticleCreator/ParticleCreator.h>

#include <CCA/Components/Peridynamics/Core/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsMaterial.h>
#include <CCA/Components/Peridynamics/DamageModels/PeridynamicsDamageModel.h>
#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModel.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticleForceBC.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticleLoadBCFactory.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticleNormalForceBC.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticlePressureBC.h>

#include <CCA/Ports/DataWarehouse.h>

#include <Core/GeometryPiece/FileGeometryPiece.h>
#include <Core/GeometryPiece/GeometryObject.h>
#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/GeometryPiece/SpecialGeomPiece.h>

#include <Core/Grid/Box.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/VarLabel.h>

#include <Core/Util/DebugStream.h>

#include <fstream>
#include <iostream>

using namespace Vaango;

//__________________________________
//  To turn on debug flags
//  csh/tcsh : setenv SCI_DEBUG "PDPartDoing:+,PDPartDebug:+".....
//  bash     : export SCI_DEBUG="PDPartDoing:+,PDPartDebug:+" )
//  default is OFF
using Uintah::DebugStream;
static DebugStream cout_doing("PDPartDoing", false);
static DebugStream cout_dbg("PDPartDebug", false);

ParticleCreator::ParticleCreator(PeridynamicsMaterial* matl,
                                 const PeridynamicsFlags* flags)
{
  d_varLabel = std::make_unique<PeridynamicsLabel>();
  d_flags    = flags;
  registerPermanentParticleState(matl);
}

ParticleCreator::~ParticleCreator() {}

particleIndex
ParticleCreator::createParticles(PeridynamicsMaterial* matl,
                                 Uintah::CCVariable<short int>& cellNAPID,
                                 const Uintah::Patch* patch,
                                 Uintah::DataWarehouse* new_dw,
                                 const VecGeometryObjectSP& geom_objs)
{
  ObjectVars obj_vars;
  particleIndex numParticles = 0;
  for (const auto& geom : geom_objs) {
    numParticles += countAndCreateParticles(patch, geom.get(), obj_vars);
  }

  // Allocate the space for particle variables associated with this material
  ParticleVars pvars;
  int matlIndex = matl->getDWIndex();
  [[maybe_unused]] Uintah::ParticleSubset* subset =
    allocateVariables(numParticles, matlIndex, patch, new_dw, pvars);

  // Loop through geometry objects
  particleIndex start = 0;

  for (auto& obj : geom_objs) {
    particleIndex count = 0;

    // Get the geometry piece and its bounding box
    Uintah::GeometryPieceP piece = obj->getPiece();
    Uintah::Box b1               = piece->getBoundingBox();
    Uintah::Box b2               = patch->getExtraBox();
    Uintah::Box b                = b1.intersect(b2);
    if (b.degenerate()) {
      count = 0;
      continue;
    }

    // Get the distance from the actual boundary of the object to the
    // particle nearest the boundary.  This is needed later for locating
    // boundary particles
    Uintah::Vector dxpp = patch->dCell() / obj->getInitialData_IntVector("res");

    // Special case exception for SpecialGeomPieces and FileGeometryPieces
    // FileGeometryPieces are derived from SpecialGeomPiece and contain the
    // particle data in a file while smooth geometry pieces generate these same
    // data.
    // **WARNING** Not sure what the effect is on Abaqus type input files.
    Uintah::SpecialGeomPiece* sgp =
      dynamic_cast<Uintah::SpecialGeomPiece*>(piece.get());

    // Set up pointers to SpecialGeometryPiece particle data
    const std::vector<double>* pVolArray           = nullptr;
    const std::vector<Uintah::Vector>* pForceArray = nullptr;
    const std::vector<Uintah::Vector>* pVelArray   = nullptr;

    // Create pairs for the particle variables created by
    // the SpecialGeometryPieces
    auto vol_pair   = std::make_pair("p.volume", obj.get());
    auto force_pair = std::make_pair("p.externalforce", obj.get());
    auto vel_pair   = std::make_pair("p.velocity", obj.get());

    // Set up iterators for SpecialGeometryObject data
    std::vector<double>::const_iterator sgp_vol_iter;
    std::vector<Uintah::Vector>::const_iterator sgp_force_iter;
    std::vector<Uintah::Vector>::const_iterator sgp_velocity_iter;

    if (sgp) {
      if ((pVolArray = sgp->getScalar("p.volume"))) {
        sgp_vol_iter = obj_vars.scalars.at(vol_pair).begin();
      }
      if ((pForceArray = sgp->getVector("p.externalforce"))) {
        sgp_force_iter = obj_vars.vectors.at(force_pair).begin();
      }
      if ((pVelArray = sgp->getVector("p.velocity"))) {
        sgp_velocity_iter = obj_vars.vectors.at(vel_pair).begin();
      }
    }

    // Loop through object points
    for (auto pPosition : obj_vars.points.at(obj.get())) {

      Uintah::IntVector cell_idx;
      if (!patch->findCell(pPosition, cell_idx)) {
        continue;
      }
      if (!patch->containsPoint(pPosition)) {
        continue;
      }

      particleIndex pidx = start + count;
      cout_dbg << "\t\t\t CreateParticles: Point[" << pidx << "]=" << pPosition
               << " Cell = " << cell_idx << std::endl;

      initializeParticle(
        patch, obj.get(), matl, pPosition, cell_idx, pidx, cellNAPID, pvars);

      if (sgp) {
        if (sgp_vol_iter != obj_vars.scalars.at(vol_pair).end()) {
          pvars.pVolume[pidx] = *sgp_vol_iter;
          pvars.pMass[pidx]   = matl->getInitialDensity() * pvars.pVolume[pidx];
          ++sgp_vol_iter;
        }

        if (sgp_velocity_iter != obj_vars.vectors.at(vel_pair).end()) {
          pvars.pVelocity[pidx] = *sgp_velocity_iter;
          ++sgp_velocity_iter;
        }

        if (sgp_force_iter != obj_vars.vectors.at(force_pair).end()) {
          pvars.pExternalForce[pidx] = *sgp_force_iter;
          ++sgp_force_iter;
        }
      }

      // If the particle is on the surface and if there is
      // a physical BC attached to it then mark with the
      // physical BC pointer
      if (checkForSurface(piece, pPosition, dxpp)) {
        pvars.pLoadCurveID[pidx] = getLoadCurveID(pPosition, dxpp);
      } else {
        pvars.pLoadCurveID[pidx] = 0;
      }

      count++;
    }
    start += count;
  }
  return numParticles;
}

particleIndex
ParticleCreator::countAndCreateParticles(const Uintah::Patch* patch,
                                         Uintah::GeometryObject* obj,
                                         ObjectVars& obj_vars)
{
  Uintah::GeometryPieceP piece = obj->getPiece();
  Uintah::Box b1               = piece->getBoundingBox();
  Uintah::Box b2               = patch->getExtraBox();
  Uintah::Box b                = b1.intersect(b2);
  if (b.degenerate()) {
    return 0;
  }

  // Create pairs for the particle variables created by
  // the SpecialGeometryPieces
  auto vol_pair   = std::make_pair("p.volume", obj);
  auto force_pair = std::make_pair("p.externalforce", obj);
  auto vel_pair   = std::make_pair("p.velocity", obj);

  // If the object is a SpecialGeometryPiece (e.g. FileGeometryPiece or
  // SmoothCylGeomPiece) then use the particle creators in that
  // class to do the counting
  Uintah::SpecialGeomPiece* sgp =
    dynamic_cast<Uintah::SpecialGeomPiece*>(piece.get());
  if (sgp) {
    int numPts = 0;

    // Check whether it is actually a FileGeometryPiece
    Uintah::FileGeometryPiece* fgp =
      dynamic_cast<Uintah::FileGeometryPiece*>(sgp);
    if (fgp) {
      fgp->readPoints(patch->getID());
      numPts = fgp->returnPointCount();
    } else {
      Uintah::Vector dxpp =
        patch->dCell() / obj->getInitialData_IntVector("res");
      double dx = std::min(std::min(dxpp.x(), dxpp.y()), dxpp.z());
      sgp->setParticleSpacing(dx);
      sgp->setCellSize(patch->dCell());
      numPts = sgp->createPoints();
      std::cout << "Smooth Geom Piece: Number of points created = " << numPts
                << std::endl;
    }

    auto points      = sgp->getPoints();
    auto vols        = sgp->getScalar("p.volume");
    auto pvelocities = sgp->getVector("p.velocity");
    auto pforces     = sgp->getVector("p.externalforce");

    Uintah::Point pp;
    Uintah::IntVector cell_idx;

    for (int ii = 0; ii < numPts; ++ii) {
      pp = points->at(ii);
      if (patch->findCell(pp, cell_idx)) {
        if (patch->containsPoint(pp)) {
          obj_vars.points[obj].push_back(pp);

          if (vols) {
            double vol = vols->at(ii);
            obj_vars.scalars[vol_pair].push_back(vol);
          }
          if (pvelocities) {
            Uintah::Vector pvel = pvelocities->at(ii);
            obj_vars.vectors[vel_pair].push_back(pvel);
          }
          if (pforces) {
            Uintah::Vector pforce = pforces->at(ii);
            obj_vars.vectors[force_pair].push_back(pforce);
          }
        } // patch contains cell
      }
    }
  } else { // Not a SpecialGeomPiece or FileGeometryPiece
    createPoints(patch, obj, obj_vars);
  }

  return static_cast<particleIndex>(obj_vars.points[obj].size());
}

void
ParticleCreator::createPoints(const Uintah::Patch* patch,
                              Uintah::GeometryObject* obj,
                              ObjectVars& obj_vars)
{
  Uintah::GeometryPieceP geom_piece = obj->getPiece();
  Uintah::Box box_with_extra_cells  = patch->getExtraBox();

  Uintah::IntVector particles_per_cell = obj->getInitialData_IntVector("res");
  Uintah::Vector dxpp                  = patch->dCell() / particles_per_cell;
  Uintah::Vector dcorner               = dxpp * 0.5;

  // Iterate through cells in patch
  for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {

    Uintah::IntVector cell = *iter;
    Uintah::Point lower    = patch->nodePosition(cell) + dcorner;

    for (int ix = 0; ix < particles_per_cell.x(); ix++) {
      for (int iy = 0; iy < particles_per_cell.y(); iy++) {
        for (int iz = 0; iz < particles_per_cell.z(); iz++) {

          Uintah::IntVector idx(ix, iy, iz);
          Uintah::Point point = lower + dxpp * idx;

          if (!box_with_extra_cells.contains(point)) {
            throw Uintah::InternalError(
              "Particle created outside of patch ?", __FILE__, __LINE__);
          }

          if (geom_piece->inside(point)) {
            obj_vars.points[obj].push_back(point);
          }
        }
      }
    }
  } // end cell iterator

  cout_dbg << "\t\t Number of points created in patch " << patch << " = "
           << obj_vars.points[obj].size() << std::endl;
}

Uintah::ParticleSubset*
ParticleCreator::allocateVariables(particleIndex numParticles,
                                   int matlIndex,
                                   const Uintah::Patch* patch,
                                   Uintah::DataWarehouse* new_dw,
                                   ParticleVars& pvars)
{
  Uintah::ParticleSubset* subset =
    new_dw->createParticleSubset(numParticles, matlIndex, patch);

  new_dw->allocateAndPut(
    pvars.pParticleID, d_varLabel->pParticleIDLabel, subset);
  new_dw->allocateAndPut(pvars.position, d_varLabel->pPositionLabel, subset);
  new_dw->allocateAndPut(pvars.pMass, d_varLabel->pMassLabel, subset);
  new_dw->allocateAndPut(pvars.pSize, d_varLabel->pSizeLabel, subset);
  new_dw->allocateAndPut(pvars.pVolume, d_varLabel->pVolumeLabel, subset);
  new_dw->allocateAndPut(
    pvars.pDisplacement, d_varLabel->pDisplacementLabel, subset);
  new_dw->allocateAndPut(pvars.pVelocity, d_varLabel->pVelocityLabel, subset);
  new_dw->allocateAndPut(
    pvars.pExternalForce, d_varLabel->pExternalForceLabel, subset);
  new_dw->allocateAndPut(pvars.pHorizon, d_varLabel->pHorizonLabel, subset);
  new_dw->allocateAndPut(
    pvars.pLoadCurveID, d_varLabel->pLoadCurveIDLabel, subset);

  return subset;
}

// Get the LoadCurveID applicable for this material point
// WARNING : Should be called only once per particle during a simulation
// because it updates the number of particles to which a BC is applied.
int
ParticleCreator::getLoadCurveID(const Uintah::Point& pp,
                                const Uintah::Vector& dxpp)
{
  int ret = -1; // Default load curve ID for particles without an associate
                // particle load BC
  if (ParticleLoadBCFactory::particleLoadBCs.size() == 0) {
    return ret;
  }

  for (auto iter = ParticleLoadBCFactory::particleLoadBCs.begin();
       iter != ParticleLoadBCFactory::particleLoadBCs.end();
       iter++) {

    std::string bcs_type = (*iter)->getType();

    // cerr << " BC Type = " << bcs_type << std::endl;
    if (bcs_type == "Pressure") {
      ParticlePressureBC* bc = dynamic_cast<ParticlePressureBC*>(*iter);
      if (bc->flagSurfaceParticle(pp, dxpp)) {
        ret = bc->loadCurveID();
      }
    } else if (bcs_type == "Force") {
      ParticleForceBC* bc = dynamic_cast<ParticleForceBC*>(*iter);
      if (bc->flagSurfaceParticle(pp, dxpp)) {
        ret = bc->loadCurveID();
      }
    } else if (bcs_type == "NormalForce") {
      ParticleNormalForceBC* bc = dynamic_cast<ParticleNormalForceBC*>(*iter);
      if (bc->flagSurfaceParticle(pp, dxpp)) {
        ret = bc->loadCurveID();
      }
    }
  }
  return ret;
}

void
ParticleCreator::initializeParticle(const Uintah::Patch* patch,
                                    Uintah::GeometryObject* obj,
                                    PeridynamicsMaterial* matl,
                                    Uintah::Point pPosition,
                                    Uintah::IntVector cell_idx,
                                    particleIndex pidx,
                                    Uintah::CCVariable<short int>& cellNAPID,
                                    ParticleVars& pvars)
{
  Uintah::IntVector ppc = obj->getInitialData_IntVector("res");
  //Uintah::Vector dxpp   = patch->dCell() / obj->getInitialData_IntVector("res");
  Uintah::Vector dxcc   = patch->dCell();

  Uintah::Matrix3 size(1. / ((double)ppc.x()),
                       0.,
                       0.,
                       0.,
                       1. / ((double)ppc.y()),
                       0.,
                       0.,
                       0.,
                       1. / ((double)ppc.z()));
  cout_dbg << "\t\t\t\t Particle size = " << size << std::endl;
  cout_dbg << "\t\t\t\t Particle location = " << pPosition << std::endl;
  cout_dbg << "\t\t\t\t Particle index = " << pidx << std::endl;

  pvars.position[pidx]  = pPosition;
  pvars.pVolume[pidx]   = size.Determinant() * dxcc.x() * dxcc.y() * dxcc.z();
  pvars.pSize[pidx]     = size;
  pvars.pVelocity[pidx] = obj->getInitialData_Vector("velocity");
  pvars.pMass[pidx]     = matl->getInitialDensity() * pvars.pVolume[pidx];
  pvars.pDisplacement[pidx] = Uintah::Vector(0., 0., 0.);

  // Compute the max length of the side of a cell and set the horizon
  // accordingly
  double maxCellEdge   = std::max(std::max(dxcc.x(), dxcc.y()), dxcc.z());
  pvars.pHorizon[pidx] = maxCellEdge * d_flags->d_numCellsInHorizon;

  // Initialize external force (updated later using ParticleLoadBCs)
  Uintah::Vector pExtForce(0.0, 0.0, 0.0);
  pvars.pExternalForce[pidx] = pExtForce;

  ASSERT(cell_idx.x() <= 0xffff && cell_idx.y() <= 0xffff &&
         cell_idx.z() <= 0xffff);

  Uintah::long64 cellID = ((Uintah::long64)cell_idx.x() << 16) |
                          ((Uintah::long64)cell_idx.y() << 32) |
                          ((Uintah::long64)cell_idx.z() << 48);

  short int& myCellNAPID = cellNAPID[cell_idx];
  pvars.pParticleID[pidx]    = (cellID | (Uintah::long64)myCellNAPID);
  ASSERT(myCellNAPID < 0x7fff);
  myCellNAPID++;
}

int
ParticleCreator::checkForSurface(const Uintah::GeometryPieceP piece,
                                 const Uintah::Point p,
                                 const Uintah::Vector dxpp)
{

  //  Check the candidate points which surround the point just passed
  //   in.  If any of those points are not also inside the object
  //  the current point is on the surface

  int ss = 0;
  // Check to the left (-x)
  if (!piece->inside(p - Uintah::Vector(dxpp.x(), 0., 0.))) {
    ss++;
  }
  // Check to the right (+x)
  if (!piece->inside(p + Uintah::Vector(dxpp.x(), 0., 0.))) {
    ss++;
  }
  // Check behind (-y)
  if (!piece->inside(p - Uintah::Vector(0., dxpp.y(), 0.))) {
    ss++;
  }
  // Check in front (+y)
  if (!piece->inside(p + Uintah::Vector(0., dxpp.y(), 0.))) {
    ss++;
  }
  // Check below (-z)
  if (!piece->inside(p - Uintah::Vector(0., 0., dxpp.z()))) {
    ss++;
  }
  // Check above (+z)
  if (!piece->inside(p + Uintah::Vector(0., 0., dxpp.z()))) {
    ss++;
  }

  if (ss > 0) {
    return 1;
  } else {
    return 0;
  }
}

std::vector<const Uintah::VarLabel*>
ParticleCreator::returnParticleState()
{
  return particle_state;
}

std::vector<const Uintah::VarLabel*>
ParticleCreator::returnParticleStatePreReloc()
{
  return particle_state_preReloc;
}

void
ParticleCreator::registerPermanentParticleState(PeridynamicsMaterial* matl)
{
  particle_state.push_back(d_varLabel->pDisplacementLabel);
  particle_state_preReloc.push_back(d_varLabel->pDisplacementLabel_preReloc);

  particle_state.push_back(d_varLabel->pVelocityLabel);
  particle_state_preReloc.push_back(d_varLabel->pVelocityLabel_preReloc);

  particle_state.push_back(d_varLabel->pExternalForceLabel);
  particle_state_preReloc.push_back(d_varLabel->pExternalForceLabel_preReloc);

  particle_state.push_back(d_varLabel->pMassLabel);
  particle_state_preReloc.push_back(d_varLabel->pMassLabel_preReloc);

  particle_state.push_back(d_varLabel->pVolumeLabel);
  particle_state_preReloc.push_back(d_varLabel->pVolumeLabel_preReloc);

  particle_state.push_back(d_varLabel->pParticleIDLabel);
  particle_state_preReloc.push_back(d_varLabel->pParticleIDLabel_preReloc);

  particle_state.push_back(d_varLabel->pSizeLabel);
  particle_state_preReloc.push_back(d_varLabel->pSizeLabel_preReloc);

  particle_state.push_back(d_varLabel->pDefGradLabel);
  particle_state_preReloc.push_back(d_varLabel->pDefGradLabel_preReloc);

  particle_state.push_back(d_varLabel->pStressLabel);
  particle_state_preReloc.push_back(d_varLabel->pStressLabel_preReloc);

  particle_state.push_back(d_varLabel->pPK1StressLabel);
  particle_state_preReloc.push_back(d_varLabel->pPK1StressLabel_preReloc);

  particle_state.push_back(d_varLabel->pHorizonLabel);
  particle_state_preReloc.push_back(d_varLabel->pHorizonLabel_preReloc);

  particle_state.push_back(d_varLabel->pNeighborListLabel);
  particle_state_preReloc.push_back(d_varLabel->pNeighborListLabel_preReloc);

  particle_state.push_back(d_varLabel->pNeighborConnLabel);
  particle_state_preReloc.push_back(d_varLabel->pNeighborConnLabel_preReloc);

  particle_state.push_back(d_varLabel->pNeighborCountLabel);
  particle_state_preReloc.push_back(d_varLabel->pNeighborCountLabel_preReloc);

  particle_state.push_back(d_varLabel->pNeighborBondForceLabel);
  particle_state_preReloc.push_back(
    d_varLabel->pNeighborBondForceLabel_preReloc);

  particle_state.push_back(d_varLabel->pNeighborBondEnergyLabel);
  particle_state_preReloc.push_back(
    d_varLabel->pNeighborBondEnergyLabel_preReloc);

  particle_state.push_back(d_varLabel->pInternalForceLabel);
  particle_state_preReloc.push_back(d_varLabel->pInternalForceLabel_preReloc);

  particle_state.push_back(d_varLabel->pDamageLabel);
  particle_state_preReloc.push_back(d_varLabel->pDamageLabel_preReloc);

  particle_state.push_back(d_varLabel->pLoadCurveIDLabel);
  particle_state_preReloc.push_back(d_varLabel->pLoadCurveIDLabel_preReloc);

  matl->getMaterialModel()->addParticleState(particle_state,
                                             particle_state_preReloc);

  matl->getDamageModel()->addParticleState(particle_state,
                                           particle_state_preReloc);
}
