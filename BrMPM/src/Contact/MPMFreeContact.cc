/*
 * MPMFreeContact.cc
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#include <Contact/MPMFreeContact.h>

using namespace BrMPM;

MPMFreeContact::MPMFreeContact(std::vector<int>& dwis, MPMPatchP& patch)
   : MPMContact(dwis, patch)
{
}

MPMFreeContact::~MPMFreeContact() {
  // TODO Auto-generated destructor stub
}

void
MPMFreeContact::exchMomentumInterpolated(MPMDatawarehouseP& dw) {
  // TODO
}

void
MPMFreeContact::exchForceInterpolated(MPMDatawarehouseP& dw) {
  // TODO
}

void
MPMFreeContact::exchMomentumIntegrated(MPMDatawarehouseP& dw) {
  // TODO
}
