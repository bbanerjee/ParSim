/*
 * MPMFrictionContact.h
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#ifndef MPMFRICTIONCONTACT_H_
#define MPMFRICTIONCONTACT_H_

#include <MPMFrictionlessContact.h>

namespace BrMPM {

  class MPMFrictionContact: public MPMFrictionlessContact {

    public:

      MPMFrictionContact();
      virtual ~MPMFrictionContact();

      virtual void exchMomentumInterpolated(MPMDatawarehouseP& dw);
      virtual void exchForceInterpolated(MPMDatawarehouseP& dw);
      virtual void exchMomentumIntegrated(MPMDatawarehouseP& dw);
  };

} /* namespace BrMPM */

#endif /* MPMFRICTIONCONTACT_H_ */
