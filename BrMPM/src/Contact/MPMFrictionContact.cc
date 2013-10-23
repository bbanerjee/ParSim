/*
 * MPMFrictionContact.cc
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#include <Contact/MPMFrictionContact.h>

using namespace BrMPM;

MPMFrictionContact::MPMFrictionContact(std::vector<int>& dwis, MPMPatchP& patch)
   : MPMFrictionlessContact(dwis, patch)
{
}

MPMFrictionContact::~MPMFrictionContact() {
  // TODO Auto-generated destructor stub
}

void MPMFrictionContact::exchMomentumInterpolated(MPMDatawarehouseP& dw) {
}

void MPMFrictionContact::exchForceInterpolated(MPMDatawarehouseP& dw) {
}

void MPMFrictionContact::exchMomentumIntegrated(MPMDatawarehouseP& dw) {
}
