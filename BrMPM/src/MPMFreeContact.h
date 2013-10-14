/*
 * MPMFreeContact.h
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#ifndef MPMFREECONTACT_H_
#define MPMFREECONTACT_H_

#include <MPMContact.h>

namespace BrMPM {

  class MPMFreeContact: public MPMContact {

  public:
    MPMFreeContact();
    virtual ~MPMFreeContact();

    void exchMomentumInterpolated(MPMDatawarehouseP& dw);
    void exchForceInterpolated(MPMDatawarehouseP& dw);
    void exchMomentumIntegrated(MPMDatawarehouseP& dw);

  };

} /* namespace BrMPM */

#endif /* MPMFREECONTACT_H_ */
