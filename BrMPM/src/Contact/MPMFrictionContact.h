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

      MPMFrictionContact(std::vector<int>& dwis, MPMPatchP& patch, double mu);
      virtual ~MPMFrictionContact();

      void exchForceInterpolated(MPMDatawarehouseP& dw);

    private:

      double d_mu; // Coefficient of friction
      double d_dt; // Time increment

      MPMFrictionContact();
  };

} /* namespace BrMPM */

#endif /* MPMFRICTIONCONTACT_H_ */
