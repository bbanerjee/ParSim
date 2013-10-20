/*
 * MPMVelocityContact.h
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#ifndef MPMVELOCITYCONTACT_H_
#define MPMVELOCITYCONTACT_H_

#include <MPMFrictionlessContact.h>

namespace BrMPM {

  class MPMVelocityContact: public MPMFrictionlessContact
  {

  public:

    MPMVelocityContact(std::vector<int>& dwis, MPMPatchP& patch);

    virtual ~MPMVelocityContact();

    void exchForceInterpolated(MPMDatawarehouseP& dw);
    void exchMomentumIntegrated(MPMDatawarehouseP& dw);

  };

} /* namespace BrMPM */
#endif /* MPMVELOCITYCONTACT_H_ */
