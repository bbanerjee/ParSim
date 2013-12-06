/*
 * GIMPShapeFunction.h
 *
 *  Created on: 3/12/2013
 *      Author: banerjee
 */

#ifndef GIMPSHAPEFUNCTION_H_
#define GIMPSHAPEFUNCTION_H_

#include <ShapeFunctions/MPMShapeFunction.h>
#include <GeometryMath/Point3D.h>

namespace BrMPM {

class GIMPShapeFunction: public MPMShapeFunction
{
public:
  GIMPShapeFunction();
  virtual ~GIMPShapeFunction();

  void updateContribList(MPMDatawarehouseP& dw,
                         const MPMPatchP& patch,
                         int dwi);

private:

  int getCell(const MPMPatchP& patch, const Point3D& pos);

  void getWeightAndDerivative(const double& xx, const double& hh, const double& pSize,
                              double& weight, double& derivative);
};

} // end namespace

#endif /* GIMPSHAPEFUNCTION_H_ */
