/*
 * MPMVelocityContact.cc
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#include <MPMVelocityContact.h>
#include <MPMDatawarehouse.h>
#include <MPMDataTypes.h>
#include <MPMUtils.h>

using namespace BrMPM;

MPMVelocityContact::MPMVelocityContact(std::vector<int>& dwis,
                                       MPMPatchP& patch)
   : MPMFrictionlessContact(dwis, patch)
{
}

MPMVelocityContact::~MPMVelocityContact()
{
}

void
MPMVelocityContact::findIntersection(MPMDatawarehouseP& dw)
{
  for (auto iter = d_dwis.begin(); iter != d_dwis.end(); ++iter) {
    int dwi = *iter;

    // Get the particle interpolation information
    VectorIntParticleData cIdx;
    VectorDoubleParticleData cGradx;
    VectorDoubleParticleData cGrady;
    VectorDoubleParticleData cGradz;
    dw->get("cIdx", dwi, cIdx);
    dw->get("cGradx", dwi, cGradx);
    dw->get("cGrady", dwi, cGrady);
    dw->get("cGradz", dwi, cGradz);

    // Get the particle state information
    DoubleParticleData pm;
    DoubleParticleData pVol;
    dw->get("pm", dwi, pm);
    dw->get("pVol", dwi, pVol);

    // Get the node data
    Vector3DNodeData gGm;
    dw->get("gGm", dwi, gGm);

    // Compute gradient
    MPMUtils::gradscalar(cIdx, cGradx, cGrady, cGradz, pm, gGm);
  }
}

void
MPMVelocityContact::exchMomentumInterpolated(MPMDatawarehouseP& dw) {
}

void
MPMVelocityContact::exchForceInterpolated(MPMDatawarehouseP& dw) {
}


void
MPMVelocityContact::exchMomentumIntegrated(MPMDatawarehouseP& dw) {
}
/* namespace BrMPM */
