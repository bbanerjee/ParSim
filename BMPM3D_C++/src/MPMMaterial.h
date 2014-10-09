/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

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
#include <ShapeFunctions/MPMShapeFunctionP.h>
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
  MPMConstitutiveModelP d_model;
  MPMShapeFunctionP d_shape;

  // Don't allow default constructor
  MPMMaterial();

};

} /* namespace BrMPM */
#endif /* MPMMATERIAL_H_ */
