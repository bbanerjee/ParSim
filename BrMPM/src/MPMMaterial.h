/*
 * MPMMaterial.h
 *
 *  Created on: 11/10/2013
 *      Author: banerjee
 */

#ifndef MPMMATERIAL_H_
#define MPMMATERIAL_H_

#include <MPMDatawarehouseP.h>
#include <MPMConstitutiveModelP.h>
#include <MPMShapeFunctionP.h>
#include <MPMPatchP.h>
#include <MPMDataTypes.h>

namespace BrMPM {

class MPMMaterial {
public:

  MPMMaterial(const int dwi,
		      const MPMConstitutiveModelP& model,
		      const MPMShapeFunctionP& shape);

  virtual ~MPMMaterial();

  inline int getdwi() {return d_dwi;}

  void updateContributions(MPMDatawarehouseP& dw,
                           MPMPatchP& patch);

  void setVelocity(MPMDatawarehouseP& dw,
                   Vector3DParticleData& velocity);

  void setExternalLoad(MPMDatawarehouseP& dw,
                       Vector3DParticleData& externalForce);

  void setExternalAcceleration(MPMDatawarehouseP& dw,
                               Vector3DParticleData& acceleration);

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
  MPMConstitutiveModelP& d_model;
  MPMShapeFunctionP& d_shape;

  // Don't allow default constructor
  MPMMaterial();

};

} /* namespace BrMPM */
#endif /* MPMMATERIAL_H_ */
