/*
 * MPMFrictionlessContact.h
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#ifndef MPMFRICTIONLESSCONTACT_H_
#define MPMFRICTIONLESSCONTACT_H_

#include <MPMContact.h>

namespace BrMPM {

  class MPMFrictionlessContact: public MPMContact {

    public:

      MPMFrictionlessContact(std::vector<int>& dwis,
                             MPMPatchP& patch);

      virtual ~MPMFrictionlessContact();

      virtual void exchMomentumInterpolated(MPMDatawarehouseP& dw);
      virtual void exchForceInterpolated(MPMDatawarehouseP& dw);
      virtual void exchMomentumIntegrated(MPMDatawarehouseP& dw);
  };

} /* namespace BrMPM */
#endif /* MPMFRICTIONLESSCONTACT_H_ */
