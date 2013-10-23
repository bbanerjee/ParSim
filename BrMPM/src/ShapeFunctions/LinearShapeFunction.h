/*
 * LinearShapeFunction.h
 *
 *  Created on: 22/10/2013
 *      Author: banerjee
 */

#ifndef LINEARSHAPEFUNCTION_H_
#define LINEARSHAPEFUNCTION_H_

#include <ShapeFunctions/MPMShapeFunction.h>
#include <GeometryMath/Point3D.h>

namespace BrMPM
{

class LinearShapeFunction : public MPMShapeFunction
{
public:
  LinearShapeFunction();
  virtual ~LinearShapeFunction();

  void updateContribList(MPMDatawarehouseP& dw,
                         const MPMPatchP& patch,
                         int dwi);

private:

  int getCell(const MPMPatchP& patch, const Point3D& pos);

  void getWeightAndDerivative(const double& xx, const double& hh,
                              double& weight, double& derivative);


};

} /* namespace BrMPM */
#endif /* LINEARSHAPEFUNCTION_H_ */
