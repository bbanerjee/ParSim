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

#include <CCA/Components/MPM/ParticleCreator/ShellParticleCreator.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/AMRMPMLabel.h>
#include <CCA/Components/MPM/Core/HydroMPMLabel.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/GeometryPiece/GeometryObject.h>
#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/GeometryPiece/ShellGeometryPiece.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <iostream>

using namespace Uintah;
using std::cerr;
using std::vector;

/////////////////////////////////////////////////////////////////////////
//
// Constructor
//
ShellParticleCreator::ShellParticleCreator(MPMMaterial* matl, MPMFlags* flags)
  : ParticleCreator(matl, flags)
{
}

/////////////////////////////////////////////////////////////////////////
//
// Destructor
//
ShellParticleCreator::~ShellParticleCreator() {}

/////////////////////////////////////////////////////////////////////////
//
// Actually create particles using geometry
//
particleIndex
ShellParticleCreator::createParticles(MPMMaterial* matl,
                                      CCVariable<short int>& cellNAPID,
                                      const Patch* patch,
                                      DataWarehouse* new_dw,
                                      std::vector<GeometryObject*>& d_geom_objs)
{
  // Print the physical boundary conditions
  printPhysicalBCs();

  ObjectVars vars;
  particleIndex numParticles = 0;
  std::vector<GeometryObject*>::const_iterator geom;
  for (geom = d_geom_objs.begin(); geom != d_geom_objs.end(); ++geom) {
    numParticles += countAndCreateParticles(patch, *geom, vars);
  }

  // Get datawarehouse index
  int dwi = matl->getDWIndex();

  // Create a particle subset for the patch
  ParticleVars pvars;
  ParticleSubset* subset =
    ParticleCreator::allocateVariables(numParticles, dwi, patch, new_dw, pvars);
  // Create the variables that go with each shell particle
  ParticleVariable<double> pThickTop0, pThickBot0, pThickTop, pThickBot;
  ParticleVariable<Vector> pNormal0, pNormal;
  new_dw->allocateAndPut(pThickTop, d_mpm_labels->pThickTopLabel, subset);
  new_dw->allocateAndPut(pThickTop0,
                         d_mpm_labels->pInitialThickTopLabel,
                         subset);
  new_dw->allocateAndPut(pThickBot, d_mpm_labels->pThickBotLabel, subset);
  new_dw->allocateAndPut(pThickBot0,
                         d_mpm_labels->pInitialThickBotLabel,
                         subset);
  new_dw->allocateAndPut(pNormal, d_mpm_labels->pNormalLabel, subset);
  new_dw->allocateAndPut(pNormal0, d_mpm_labels->pInitialNormalLabel, subset);

  // Initialize the global particle index
  particleIndex start = 0;

  // Loop thru the geometry objects
  std::vector<GeometryObject*>::const_iterator obj;
  for (obj = d_geom_objs.begin(); obj != d_geom_objs.end(); ++obj) {

    // Initialize the per geometryObject particle count
    particleIndex count = 0;

    // If the geometry piece is outside the patch, look
    // for the next geometry piece
    GeometryPieceP piece = (*obj)->getPiece();
    Box b = (piece->getBoundingBox()).intersect(patch->getExtraBox());
    if (b.degenerate()) {
      count = 0;
      continue;
    }

    // Find volume of influence of each particle as a
    // fraction of the cell size
    IntVector ppc  = (*obj)->getInitialData_IntVector("res");
    Vector dxpp    = patch->dCell() / (*obj)->getInitialData_IntVector("res");
    Vector dcorner = dxpp * 0.5;
    Matrix3 size(1. / ((double)ppc.x()),
                 0.,
                 0.,
                 0.,
                 1. / ((double)ppc.y()),
                 0.,
                 0.,
                 0.,
                 1. / ((double)ppc.z()));

    // If the geometry object is a shell perform special
    // operations else just treat the geom object in the standard
    // way
    ShellGeometryPiece* shell = dynamic_cast<ShellGeometryPiece*>(piece.get());

    // Create the appropriate particles
    if (shell) {

      // The position, volume and size variables are from the
      // ParticleCreator class
      int numP = shell->createParticles(patch,
                                        pvars.position,
                                        pvars.pVolume,
                                        pThickTop,
                                        pThickBot,
                                        pNormal,
                                        pvars.pSize,
                                        start);

      // Update the other variables that are attached to each particle
      // (declared in the ParticleCreator class)
      for (int idx = 0; idx < numP; idx++) {
        particleIndex pidx       = start + idx;
        pvars.pVelocity[pidx]    = (*obj)->getInitialData_Vector("velocity");
        pvars.pTemperature[pidx] = (*obj)->getInitialData_double("temperature");
        pvars.pDisp[pidx]        = Vector(0., 0., 0.);
        pvars.pFiberDir[pidx]    = Vector(0.0, 0.0, 0.0);

        // Calculate particle mass
        double partMass   = matl->getInitialDensity() * pvars.pVolume[pidx];
        pvars.pMass[pidx] = partMass;

        // The particle can be tagged as a surface particle.
        // If there is a physical BC attached to it then mark with the
        // physical BC pointer
        if (d_useLoadCurves) {
          pvars.pLoadCurveID[pidx] = getLoadCurveID(pvars.position[pidx], dxpp);
        }

        // Apply the force BC if applicable
        Vector pExtForce(0, 0, 0);
        ParticleCreator::applyForceBC(dxpp,
                                      pvars.position[pidx],
                                      partMass,
                                      pExtForce);
        pvars.pExternalForce[pidx] = pExtForce;

        // Assign a particle id
        IntVector cell_idx;
        if (patch->findCell(pvars.position[pidx], cell_idx)) {
          long64 cellID = ((long64)cell_idx.x() << 16) |
                          ((long64)cell_idx.y() << 32) |
                          ((long64)cell_idx.z() << 48);
          short int& myCellNAPID = cellNAPID[cell_idx];
          ASSERT(myCellNAPID < 0x7fff);
          myCellNAPID++;
          pvars.pParticleID[pidx] = cellID | (long64)myCellNAPID;
        } else {
          double x = pvars.position[pidx].x();
          double y = pvars.position[pidx].y();
          double z = pvars.position[pidx].z();
          if (std::abs(x) < 1.0e-15) {
            x = 0.0;
          }
          if (std::abs(y) < 1.0e-15) {
            y = 0.0;
          }
          if (std::abs(z) < 1.0e-15) {
            z = 0.0;
          }
          double px = patch->getExtraBox().upper().x();
          double py = patch->getExtraBox().upper().y();
          double pz = patch->getExtraBox().upper().z();
          if (fabs(px) < 1.0e-15) {
            px = 0.0;
          }
          if (fabs(py) < 1.0e-15) {
            py = 0.0;
          }
          if (fabs(pz) < 1.0e-15) {
            pz = 0.0;
          }
          if (x == px) {
            x -= 1.0e-10;
          }
          if (y == py) {
            y -= 1.0e-10;
          }
          if (z == pz) {
            z -= 1.0e-10;
          }
          pvars.position[pidx] = Point(x, y, z);
          if (!patch->findCell(pvars.position[pidx], cell_idx)) {
            std::cerr <<  "Pidx = " << pidx << " Pos = " << pvars.position[pidx]
                 << " patch BBox = " << patch->getExtraBox()
                 << " cell_idx = " << cell_idx
                 << " low = " << patch->getExtraCellLowIndex()
                 << " high = " << patch->getExtraCellHighIndex()
                 << " : Particle not in any cell." << std::endl;
            pvars.pParticleID[pidx] = 0;
          } else {
            long64 cellID = ((long64)cell_idx.x() << 16) |
                            ((long64)cell_idx.y() << 32) |
                            ((long64)cell_idx.z() << 48);
            short int& myCellNAPID = cellNAPID[cell_idx];
            ASSERT(myCellNAPID < 0x7fff);
            myCellNAPID++;
            pvars.pParticleID[pidx] = cellID | (long64)myCellNAPID;
          }
        }

        // The shell specific variables
        pThickTop0[pidx] = pThickTop[pidx];
        pThickBot0[pidx] = pThickBot[pidx];
        pNormal0[pidx]   = pNormal[pidx];

      } // End of loop thry particles per geom-object

    } else {

      // Loop thru cells and assign particles
      for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
        Point lower = patch->nodePosition(*iter) + dcorner;
        for (int ix = 0; ix < ppc.x(); ix++) {
          for (int iy = 0; iy < ppc.y(); iy++) {
            for (int iz = 0; iz < ppc.z(); iz++) {
              IntVector idx(ix, iy, iz);
              Point p            = lower + dxpp * idx;
              IntVector cell_idx = *iter;
              // If the assertion fails then we may just need to change
              // the format of particle ids such that the cell indices
              // have more bits.
              ASSERT(cell_idx.x() <= 0xffff && cell_idx.y() <= 0xffff &&
                     cell_idx.z() <= 0xffff);
              long64 cellID = ((long64)cell_idx.x() << 16) |
                              ((long64)cell_idx.y() << 32) |
                              ((long64)cell_idx.z() << 48);
              if (piece->inside(p)) {
                particleIndex pidx   = start + count;
                pvars.position[pidx] = p;
                pvars.pDisp[pidx]    = Vector(0., 0., 0.);
                pvars.pVolume[pidx]  = dxpp.x() * dxpp.y() * dxpp.z();
                pvars.pVelocity[pidx] =
                  (*obj)->getInitialData_Vector("velocity");
                pvars.pTemperature[pidx] =
                  (*obj)->getInitialData_double("temperature");
                pvars.pSpecificVolume[pidx] = 1.0 / matl->getInitialDensity();
                pvars.pFiberDir[pidx]       = Vector(0.0, 0.0, 0.0);

                // Calculate particle mass
                double partMass =
                  matl->getInitialDensity() * pvars.pVolume[pidx];
                pvars.pMass[pidx] = partMass;

                // If the particle is on the surface and if there is
                // a physical BC attached to it then mark with the
                // physical BC pointer
                if (d_useLoadCurves) {
                  if (checkForSurface(piece, p, dxpp)) {
                    pvars.pLoadCurveID[pidx] = getLoadCurveID(p, dxpp);
                  } else {
                    pvars.pLoadCurveID[pidx] = 0;
                  }
                }

                // Apply the force BC if applicable
                Vector pExtForce(0, 0, 0);
                ParticleCreator::applyForceBC(dxpp, p, partMass, pExtForce);
                pvars.pExternalForce[pidx] = pExtForce;

                // Assign particle id
                short int& myCellNAPID  = cellNAPID[cell_idx];
                pvars.pParticleID[pidx] = cellID | (long64)myCellNAPID;
                pvars.pSize[pidx]       = size;
                ASSERT(myCellNAPID < 0x7fff);
                myCellNAPID++;
                count++;

                // Assign dummy values to shell-specific variables
                pThickTop[pidx]  = 1.0;
                pThickTop0[pidx] = 1.0;
                pThickBot[pidx]  = 1.0;
                pThickBot0[pidx] = 1.0;
                pNormal[pidx]    = Vector(0, 1, 0);
                pNormal0[pidx]   = Vector(0, 1, 0);
              } // if inside
            }   // loop in z
          }     // loop in y
        }       // loop in x
      }         // for
    }           // end of else
    start += count;
  }

  return numParticles;
}

/////////////////////////////////////////////////////////////////////////
//
// Return number of particles
//
particleIndex
ShellParticleCreator::countAndCreateParticles(const Patch* patch,
                                              GeometryObject* obj,
                                              ObjectVars& vars)
{

  GeometryPieceP piece      = obj->getPiece();
  ShellGeometryPiece* shell = dynamic_cast<ShellGeometryPiece*>(piece.get());
  if (shell) {
    return shell->returnParticleCount(patch);
  }
  return ParticleCreator::countAndCreateParticles(patch, obj, vars);
}

/////////////////////////////////////////////////////////////////////////
//
// Register variables for crossing patches
//
void
ShellParticleCreator::registerPermanentParticleState(MPMMaterial*)
{
  particle_state.push_back(d_mpm_labels->pThickTopLabel);
  particle_state.push_back(d_mpm_labels->pInitialThickTopLabel);
  particle_state.push_back(d_mpm_labels->pThickBotLabel);
  particle_state.push_back(d_mpm_labels->pInitialThickBotLabel);
  particle_state.push_back(d_mpm_labels->pNormalLabel);
  particle_state.push_back(d_mpm_labels->pInitialNormalLabel);

  particle_state_preReloc.push_back(d_mpm_labels->pThickTopLabel_preReloc);
  particle_state_preReloc.push_back(
    d_mpm_labels->pInitialThickTopLabel_preReloc);
  particle_state_preReloc.push_back(d_mpm_labels->pThickBotLabel_preReloc);
  particle_state_preReloc.push_back(
    d_mpm_labels->pInitialThickBotLabel_preReloc);
  particle_state_preReloc.push_back(d_mpm_labels->pNormalLabel_preReloc);
  particle_state_preReloc.push_back(d_mpm_labels->pInitialNormalLabel_preReloc);
}
