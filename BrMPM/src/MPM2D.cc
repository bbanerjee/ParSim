/*
 * MPM2D.cc
 *
 *  Created on: 11/10/2013
 *      Author: banerjee
 */

#include <MPM2D.h>

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
                        MPMMaterialsList& mats, MPMContactList& contacts)
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
                                     MPMContactList& contacts)
{
  for (auto iter = contacts.begin(); iter != contacts.end(); ++iter) {
    MPMContact contact = *iter;
    contact.exchangeMomentumInterpolated(dw);
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
                                  MPMContactList& contacts)
{
  for (auto iter = contacts.begin(); iter != contacts.end(); ++iter) {
    MPMContact contact = *iter;
    contact.exchForceInterpolated(dw);
  }
}

void MPM2D::computeAndIntegrateAcceleration(MPMDatawarehouseP& dw,
                                            MPMPatchP& patch, MPMMaterialsList& mats)
{
  double tol = patch.getTolerance();
  for (auto iter = mats.begin(); iter != mats.end(); ++iter) {
    MPMMaterial mat = *iter;
    mat.computeAndIntegrateAcceleration(dw, patch, tol);
  }
}

void MPM2D::exchMomentumIntegrated(MPMDatawarehouseP& dw,
                                   MPMContactList& contacts)
{
  for (auto iter = contacts.begin(); iter != contacts.end(); ++iter) {
    MPMContact contact = *iter;
    contact.exchMomentumInterpolated(dw);
  }
}

void MPM2D::setGridBoundaryConditions(MPMDatawarehouseP& dw, MPMPatchP& patch)
{
  MPMBCs patchBCs = patch.getBCs();
  double tol = patch.getTolerance();
  for (auto iter = patchBCs.begin(); iter != patchBCs.end(); ++iter) {
    MPMBC bc = *iter;
    bc.setBoundCond(dw, patch, tol);
  }
}

void MPM2D::interpolateToParticlesAndUpdate(MPMDatawarehouseP& dw,
                                            MPMPatchP& patch, MPMMaterialsList& mats)
{
  for (auto iter = mats.begin(); iter != mats.end(); ++iter) {
    MPMMaterial mat = *iter;
    mat.interpolateToParticlesAndUpdate(dw, patch, tol);
  }
  patch.stepTime();
}
