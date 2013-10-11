/*
 * MPMMaterial.h
 *
 *  Created on: 11/10/2013
 *      Author: banerjee
 */

#ifndef MPMMATERIAL_H_
#define MPMMATERIAL_H_

#include <MPMDatawarehouseP.h>
#include <MPMPatchP.h>

namespace BrMPM {

class MPMMaterial {
public:
  MPMMaterial();
  virtual ~MPMMaterial();

  inline int getdwi() {return d_dwi;}

  void updateContributions(MPMDatawarehouseP& dw,
                           MPMPatchP& patch );

  void setVelocity(MPMDatawarehouseP& dw,
                   MPMParticleVar<Vector>& v );

  void setExternalLoad(MPMDatawarehouseP& dw,
                       MPMParticleVar<Vector>fe );

  void setExternalAcceleration(MPMDatawarehouseP& dw,
                               MPMParticleVar<Vector>acc );

  void applyExternalLoads(MPMDatawarehouseP& dw,
                          MPMPatchP& patch );

  void interpolateParticlesToGrid(MPMDatawarehouseP& dw,
                                  MPMPatchP& patch );

  void computeStressTensor(MPMDatawarehouseP& dw,
                           MPMPatchP& patch );

  void computeInternalForce(MPMDatawarehouseP& dw,
                            MPMPatchP& patch );

  void computeAndIntegrateAcceleration(MPMDatawarehouseP& dw,
                                       MPMPatchP& patch,
                                       double& tol );

  void interpolateToParticlesAndUpdate(MPMDatawarehouseP& dw,
                                       MPMPatchP& patch );
private:

  int d_dwi;

};

} /* namespace BrMPM */
#endif /* MPMMATERIAL_H_ */
