/*
 * MPMFreeContact.h
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#ifndef MPMFREECONTACT_H_
#define MPMFREECONTACT_H_

#include <Contact/MPMContact.h>

namespace BrMPM {

  class MPMFreeContact: public MPMContact {

  public:
    MPMFreeContact(std::vector<int>& dwis, MPMPatchP& patch);
    virtual ~MPMFreeContact();

    void exchMomentumInterpolated(MPMDatawarehouseP& dw);
    void exchForceInterpolated(MPMDatawarehouseP& dw);

  private:

    // Don't allow construction without arguments
    MPMFreeContact();

  };

} /* namespace BrMPM */

#endif /* MPMFREECONTACT_H_ */
