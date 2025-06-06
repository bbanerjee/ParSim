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

#include <CCA/Components/MPM/ParticleCreator/MembraneParticleCreator.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/AMRMPMLabel.h>
#include <CCA/Components/MPM/Core/HydroMPMLabel.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/GeometryPiece/GeometryObject.h>
#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/GeometryPiece/SphereMembraneGeometryPiece.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <iostream>

using namespace Uintah;

MembraneParticleCreator::MembraneParticleCreator(MPMMaterial* matl,
                                                 MPMFlags* flags)

  : ParticleCreator(matl, flags)
{
  registerPermanentParticleState(matl);
  pTang1Label =
    VarLabel::create("p.tang1", ParticleVariable<Vector>::getTypeDescription());

  pTang2Label =
    VarLabel::create("p.tang2", ParticleVariable<Vector>::getTypeDescription());

  pNormLabel =
    VarLabel::create("p.norm", ParticleVariable<Vector>::getTypeDescription());

  pTang1Label_preReloc =
    VarLabel::create("p.tang1+",
                     ParticleVariable<Vector>::getTypeDescription());

  pTang2Label_preReloc =
    VarLabel::create("p.tang2+",
                     ParticleVariable<Vector>::getTypeDescription());

  pNormLabel_preReloc =
    VarLabel::create("p.norm+", ParticleVariable<Vector>::getTypeDescription());
}

MembraneParticleCreator::~MembraneParticleCreator()
{
  VarLabel::destroy(pTang1Label);
  VarLabel::destroy(pTang1Label_preReloc);
  VarLabel::destroy(pTang2Label);
  VarLabel::destroy(pTang2Label_preReloc);
  VarLabel::destroy(pNormLabel);
  VarLabel::destroy(pNormLabel_preReloc);
}

particleIndex
MembraneParticleCreator::createParticles(
  MPMMaterial* matl,
  CCVariable<short int>& cellNAPID,
  const Patch* patch,
  DataWarehouse* new_dw,
  std::vector<GeometryObject*>& d_geom_objs)
{
  ObjectVars vars;
  particleIndex numParticles = 0;
  std::vector<GeometryObject*>::const_iterator geom;
  for (geom = d_geom_objs.begin(); geom != d_geom_objs.end(); ++geom) {
    numParticles += countAndCreateParticles(patch, *geom, vars);
  }
  int dwi = matl->getDWIndex();

  ParticleVars pvars;
  allocateVariables(numParticles, dwi, patch, new_dw, pvars);

  particleIndex start = 0;

  std::vector<GeometryObject*>::const_iterator obj;
  for (obj = d_geom_objs.begin(); obj != d_geom_objs.end(); ++obj) {
    particleIndex count  = 0;
    GeometryPieceP piece = (*obj)->getPiece();
    Box b1               = piece->getBoundingBox();
    Box b2               = patch->getExtraBox();
    Box b                = b1.intersect(b2);
    if (b.degenerate()) {
      count = 0;
    }

    IntVector ppc  = (*obj)->getInitialData_IntVector("res");
    Vector dxpp    = patch->dCell() / (*obj)->getInitialData_IntVector("res");
    Vector dcorner = dxpp * 0.5;
    // Size as a fraction of the cell size
    Matrix3 size(1. / ((double)ppc.x()),
                 0.,
                 0.,
                 0.,
                 1. / ((double)ppc.y()),
                 0.,
                 0.,
                 0.,
                 1. / ((double)ppc.z()));

    SphereMembraneGeometryPiece* SMGP =
      dynamic_cast<SphereMembraneGeometryPiece*>(piece.get());
    if (SMGP) {
      int numP = SMGP->createParticles(patch,
                                       pvars.position,
                                       pvars.pVolume,
                                       pvars.pTang1,
                                       pvars.pTang2,
                                       pvars.pNorm,
                                       pvars.pSize,
                                       start); // CPTI
      for (int idx = 0; idx < (start + numP); idx++) {
        pvars.pVelocity[start + idx] =
          (*obj)->getInitialData_Vector("velocity");
        pvars.pTemperature[start + idx] =
          (*obj)->getInitialData_double("temperature");
        pvars.pSpecificVolume[start + idx] = 1.0 / matl->getInitialDensity();
        pvars.pMass[start + idx] =
          matl->getInitialDensity() * pvars.pVolume[start + idx];
        // Determine if particle is on the surface
        pvars.pExternalForce[start + idx] = Vector(0, 0, 0); // for now
        IntVector cell_idx;
        if (patch->findCell(pvars.position[start + idx], cell_idx)) {
          long64 cellID = ((long64)cell_idx.x() << 16) |
                          ((long64)cell_idx.y() << 32) |
                          ((long64)cell_idx.z() << 48);
          short int& myCellNAPID = cellNAPID[cell_idx];
          ASSERT(myCellNAPID < 0x7fff);
          myCellNAPID++;
          pvars.pParticleID[start + idx] = cellID | (long64)myCellNAPID;
        } else {
          std::cerr << "cellID is not right" << std::endl;
        }
      }
    } else {
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
                pvars.position[start + count] = p;
                pvars.pVolume[start + count]  = dxpp.x() * dxpp.y() * dxpp.z();
                pvars.pVelocity[start + count] =
                  (*obj)->getInitialData_Vector("velocity");
                pvars.pTemperature[start + count] =
                  (*obj)->getInitialData_double("temperature");
                pvars.pSpecificVolume[start + count] =
                  1.0 / matl->getInitialDensity();
                // Calculate particle mass
                double partMass =
                  matl->getInitialDensity() * pvars.pVolume[start + count];
                pvars.pMass[start + count] = partMass;

                // Apply the force BC if applicable
                Vector pExtForce(0, 0, 0);
                ParticleCreator::applyForceBC(dxpp, p, partMass, pExtForce);
                pvars.pExternalForce[start + count] = pExtForce;

                // Determine if particle is on the surface
                pvars.pSize[start + count]       = size;
                pvars.pTang1[start + count]      = Vector(1, 0, 0);
                pvars.pTang2[start + count]      = Vector(0, 0, 1);
                pvars.pNorm[start + count]       = Vector(0, 1, 0);
                short int& myCellNAPID           = cellNAPID[cell_idx];
                pvars.pParticleID[start + count] = cellID | (long64)myCellNAPID;
                ASSERT(myCellNAPID < 0x7fff);
                myCellNAPID++;

                count++;

              } // if inside
            }   // loop in z
          }     // loop in y
        }       // loop in x
      }         // for
    }           // else
    start += count;
  }

  return numParticles;
}

particleIndex
MembraneParticleCreator::countAndCreateParticles(const Patch* patch,
                                                 GeometryObject* obj,
                                                 ObjectVars& vars)
{

  GeometryPieceP piece = obj->getPiece();

  SphereMembraneGeometryPiece* SMGP =
    dynamic_cast<SphereMembraneGeometryPiece*>(piece.get());

  if (SMGP) {
    return SMGP->returnParticleCount(patch);
  } else {
    return ParticleCreator::countAndCreateParticles(patch, obj, vars);
  }
}

ParticleSubset*
MembraneParticleCreator::allocateVariables(particleIndex numParticles,
                                           int dwi,
                                           const Patch* patch,
                                           DataWarehouse* new_dw,
                                           ParticleVars& pvars)
{

  ParticleSubset* subset =
    ParticleCreator::allocateVariables(numParticles, dwi, patch, new_dw, pvars);

  new_dw->allocateAndPut(pvars.pTang1, pTang1Label, subset);
  new_dw->allocateAndPut(pvars.pTang2, pTang2Label, subset);
  new_dw->allocateAndPut(pvars.pNorm, pNormLabel, subset);

  return subset;
}

void
MembraneParticleCreator::registerPermanentParticleState(MPMMaterial*)
{
  particle_state.push_back(pTang1Label);
  particle_state_preReloc.push_back(pTang1Label_preReloc);

  particle_state.push_back(pTang2Label);
  particle_state_preReloc.push_back(pTang2Label_preReloc);

  particle_state.push_back(pNormLabel);
  particle_state_preReloc.push_back(pNormLabel_preReloc);
}
