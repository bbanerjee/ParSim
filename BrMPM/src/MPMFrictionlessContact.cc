/*
 * MPMFrictionlessContact.cc
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#include <MPMFrictionlessContact.h>

using namespace BrMPM;

MPMFrictionlessContact::MPMFrictionlessContact(std::vector<int>& dwis,
                                               MPMPatchP& patch)
  : MPMContact(dwis, patch)
{
  // TODO Auto-generated constructor stub

}

MPMFrictionlessContact::~MPMFrictionlessContact() {
  // TODO Auto-generated destructor stub
}

void
MPMFrictionlessContact::exchMomentumInterpolated(MPMDatawarehouseP& dw) {
}

void
MPMFrictionlessContact::exchForceInterpolated(MPMDatawarehouseP& dw) {
}

void
MPMFrictionlessContact::exchMomentumIntegrated(MPMDatawarehouseP& dw) {
}
