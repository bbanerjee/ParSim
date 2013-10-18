/*
 * MPMVelocityContact.cc
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#include <MPMVelocityContact.h>
#include <MPMParticleData.h>

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
    dw.get("cIdx", cIdx, dwi);
    dw.get("cGrad", cGrad, dwi);

    // Get the particle state information
    DoubleParticleData pm;
    DoubleParticleData pVol;
    DoubleParticleData gGm;
    dw.get("pm", pm, dwi);
    dw.get("pVol", pVol, dwi);
    dw.get("gGm", gGm, dwi);

    // Compute gradient
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
