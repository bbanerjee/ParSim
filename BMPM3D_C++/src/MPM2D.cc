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

/*
 * MPM2D.cc
 *
 *  Created on: 11/10/2013
 *      Author: banerjee
 */

#include <MPM2D.h>
#include <MPMDatawarehouse.h>
#include <MPMPatch.h>
#include <MPMMaterial.h>
#include <Contact/MPMContact.h>
#include <Contact/MPMContactP.h>

using namespace BrMPM;

MPM2D::MPM2D()
{
	// TODO Auto-generated constructor stub

}

MPM2D::~MPM2D()
{
	// TODO Auto-generated destructor stub
}

void MPM2D::timeAdvance(MPMDatawarehouseP& dw, MPMPatchP& patch,
                        MPMMaterialsList& mats, MPMContactPList& contacts)
{
    // Advance timestep
    updateMats( dw, patch, mats );
    applyExternalLoads( dw, patch, mats );
    interpolateParticlesToGrid( dw, patch, mats );
    exchMomentumInterpolated( dw, contacts );
    computeStressTensor( dw, patch, mats );
    computeInternalForce( dw, patch, mats );
    exchForceInterpolated( dw, contacts );
    computeAndIntegrateAcceleration( dw, patch, mats );
    exchMomentumIntegrated( dw, contacts );
    setGridBoundaryConditions( dw, patch );
    interpolateToParticlesAndUpdate( dw, patch, mats );
}

void MPM2D::updateMats(MPMDatawarehouseP& dw, MPMPatchP& patch,
                       MPMMaterialsList& mats)
{
  for (auto iter = mats.begin(); iter != mats.end(); ++iter) {
    MPMMaterial mat = *iter;
    mat.updateContributions(dw, patch);
  }
}

void MPM2D::applyExternalLoads(MPMDatawarehouseP& dw, MPMPatchP& patch,
                               MPMMaterialsList& mats)
{
  for (auto iter = mats.begin(); iter != mats.end(); ++iter) {
    MPMMaterial mat = *iter;
    mat.applyExternalLoads(dw, patch);
  }
}

void MPM2D::interpolateParticlesToGrid(MPMDatawarehouseP& dw, MPMPatchP& patch,
                                       MPMMaterialsList& mats)
{
  for (auto iter = mats.begin(); iter != mats.end(); ++iter) {
    MPMMaterial mat = *iter;
    mat.interpolateParticlesToGrid(dw, patch);
  }
}

void MPM2D::exchMomentumInterpolated(MPMDatawarehouseP dw,
                                     MPMContactPList& contacts)
{
  for (auto iter = contacts.begin(); iter != contacts.end(); ++iter) {
    MPMContactP contact = *iter;
    contact->exchMomentumInterpolated(dw);
  }
}

void MPM2D::computeStressTensor(MPMDatawarehouseP& dw, MPMPatchP& patch,
                                MPMMaterialsList& mats)
{
  for (auto iter = mats.begin(); iter != mats.end(); ++iter) {
    MPMMaterial mat = *iter;
    mat.computeStressTensor(dw, patch);
  }
}

void MPM2D::computeInternalForce(MPMDatawarehouseP& dw, MPMPatchP& patch,
                                 MPMMaterialsList& mats)
{
  for (auto iter = mats.begin(); iter != mats.end(); ++iter) {
    MPMMaterial mat = *iter;
    mat.computeInternalForce(dw, patch);
  }
}

void MPM2D::exchForceInterpolated(MPMDatawarehouseP& dw,
                                  MPMContactPList& contacts)
{
  for (auto iter = contacts.begin(); iter != contacts.end(); ++iter) {
    MPMContactP contact = *iter;
    contact->exchForceInterpolated(dw);
  }
}

void MPM2D::computeAndIntegrateAcceleration(MPMDatawarehouseP& dw,
                                            MPMPatchP& patch, MPMMaterialsList& mats)
{
  double tol = patch->getTolerance();
  for (auto iter = mats.begin(); iter != mats.end(); ++iter) {
    MPMMaterial mat = *iter;
    mat.computeAndIntegrateAcceleration(dw, patch, tol);
  }
}

void MPM2D::exchMomentumIntegrated(MPMDatawarehouseP& dw,
                                   MPMContactPList& contacts)
{
  for (auto iter = contacts.begin(); iter != contacts.end(); ++iter) {
    MPMContactP contact = *iter;
    contact->exchMomentumInterpolated(dw);
  }
}

void MPM2D::setGridBoundaryConditions(MPMDatawarehouseP& dw, MPMPatchP& patch)
{
  /*  TODO: Implement MPMBCs and MPMBC */
  /*
  MPMBCs patchBCs = patch->getBCs();
  double tol = patch.getTolerance();
  for (auto iter = patchBCs.begin(); iter != patchBCs.end(); ++iter) {
    MPMBC bc = *iter;
    bc.setBoundCond(dw, patch, tol);
  }
  */
}

void MPM2D::interpolateToParticlesAndUpdate(MPMDatawarehouseP& dw,
                                            MPMPatchP& patch, MPMMaterialsList& mats)
{
  for (auto iter = mats.begin(); iter != mats.end(); ++iter) {
    MPMMaterial mat = *iter;
    mat.interpolateToParticlesAndUpdate(dw, patch);
  }
  patch->stepTime();
}
