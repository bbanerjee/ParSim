/*
 * MPMShapeFunction.h
 *
 *  Created on: 22/10/2013
 *      Author: banerjee
 */

#ifndef MPMSHAPEFUNCTION_H_
#define MPMSHAPEFUNCTION_H_

#include <MPMDatawarehouseP.h>
#include <MPMPatchP.h>

namespace BrMPM
{

class MPMShapeFunction
{
public:
  MPMShapeFunction() {}
  virtual ~MPMShapeFunction() {}

  virtual void updateContribList(MPMDatawarehouseP& dw,
                                 const MPMPatchP& patch,
                                 int dwi) = 0;
};

} /* namespace BrMPM */
#endif /* MPMSHAPEFUNCTION_H_ */
