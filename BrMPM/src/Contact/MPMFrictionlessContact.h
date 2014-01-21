/*
 * MPMFrictionlessContact.h
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#ifndef MPMFRICTIONLESSCONTACT_H_
#define MPMFRICTIONLESSCONTACT_H_

#include <Contact/MPMContact.h>

namespace BrMPM {

  class MPMFrictionlessContact: public MPMContact {

    public:

      MPMFrictionlessContact(std::vector<int>& dwis,
                             MPMPatchP& patch);

      virtual ~MPMFrictionlessContact();

      void findIntersection(MPMDatawarehouseP& dw);
      void exchMomentumInterpolated(MPMDatawarehouseP& dw);
      virtual void exchForceInterpolated(MPMDatawarehouseP& dw);
      virtual void exchMomentumIntegrated(MPMDatawarehouseP& dw) {}

    private:

      MPMFrictionlessContact();
  };

} /* namespace BrMPM */
#endif /* MPMFRICTIONLESSCONTACT_H_ */
