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
