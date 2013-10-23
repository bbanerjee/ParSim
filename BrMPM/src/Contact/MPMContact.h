/*
 * MPMContact.h
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#ifndef MPMCONTACT_H_
#define MPMCONTACT_H_

#include <MPMDatawarehouseP.h>
#include <MPMPatchP.h>
#include <vector>

namespace BrMPM {

  class MPMContact {

    public:

      MPMContact();
      MPMContact(std::vector<int>& dwis, MPMPatchP& patch);
      virtual ~MPMContact();

      virtual void findIntersection(MPMDatawarehouseP& dw);

      virtual void exchMomentumInterpolated(MPMDatawarehouseP& dw);
      virtual void exchForceInterpolated(MPMDatawarehouseP& dw);
      virtual void exchMomentumIntegrated(MPMDatawarehouseP& dw);

    protected:

      void findIntersectionSimple(MPMDatawarehouseP& dw);

      std::vector<int> d_dwis;
      MPMPatchP d_patch;
      std::vector<int> d_nodes;
      double d_mtol;

  }; // end class

} /* namespace BrMPM */

#endif /* MPMCONTACT_H_ */
