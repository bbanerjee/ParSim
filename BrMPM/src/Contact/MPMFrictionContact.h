/*
 * MPMFrictionContact.h
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#ifndef MPMFRICTIONCONTACT_H_
#define MPMFRICTIONCONTACT_H_

#include <Contact/MPMFrictionlessContact.h>

namespace BrMPM {

  class MPMFrictionContact: public MPMFrictionlessContact {

    public:

      MPMFrictionContact(std::vector<int>& dwis, MPMPatchP& patch);
      virtual ~MPMFrictionContact();

      virtual void exchMomentumInterpolated(MPMDatawarehouseP& dw);
      virtual void exchForceInterpolated(MPMDatawarehouseP& dw);
      virtual void exchMomentumIntegrated(MPMDatawarehouseP& dw);

    private:

      MPMFrictionContact();
  };

} /* namespace BrMPM */

#endif /* MPMFRICTIONCONTACT_H_ */
