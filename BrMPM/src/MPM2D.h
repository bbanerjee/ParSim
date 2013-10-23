/*
 * MPM2D.h
 *
 *  Created on: 11/10/2013
 *      Author: banerjee
 */

#ifndef MPM2D_H_
#define MPM2D_H_

#include <MPMDatawarehouseP.h>
#include <MPMPatchP.h>
#include <MPMMaterialsList.h>
#include <Contact/MPMContactPList.h>

namespace BrMPM {

class MPM2D {
public:
  MPM2D();
  virtual ~MPM2D();
  void timeAdvance(MPMDatawarehouseP& dw,
                   MPMPatchP& patch,
                   MPMMaterialsList& mats,
                   MPMContactPList& contacts);

protected:

  void updateMats(MPMDatawarehouseP& dw,
                  MPMPatchP& patch,
                  MPMMaterialsList& mats);

  void applyExternalLoads(MPMDatawarehouseP& dw,
                          MPMPatchP& patch,
                          MPMMaterialsList& mats);

  void interpolateParticlesToGrid(MPMDatawarehouseP& dw,
                                  MPMPatchP& patch,
                                  MPMMaterialsList& mats );

  void exchMomentumInterpolated(MPMDatawarehouseP dw,
                                MPMContactPList& contacts );

  void computeStressTensor(MPMDatawarehouseP& dw,
                           MPMPatchP& patch,
                           MPMMaterialsList& mats );

  void computeInternalForce(MPMDatawarehouseP& dw,
                            MPMPatchP& patch,
                            MPMMaterialsList& mats );

  void exchForceInterpolated(MPMDatawarehouseP& dw,
                             MPMContactPList& contacts );

  void computeAndIntegrateAcceleration(MPMDatawarehouseP& dw,
                                       MPMPatchP& patch,
                                       MPMMaterialsList& mats );

  void exchMomentumIntegrated(MPMDatawarehouseP& dw,
                              MPMContactPList& contacts );

  void setGridBoundaryConditions(MPMDatawarehouseP& dw,
                                 MPMPatchP& patch );

  void interpolateToParticlesAndUpdate(MPMDatawarehouseP& dw,
                                       MPMPatchP& patch,
                                       MPMMaterialsList& mats );

};

} /* namespace BrMPM */
#endif /* MPM2D_H_ */
