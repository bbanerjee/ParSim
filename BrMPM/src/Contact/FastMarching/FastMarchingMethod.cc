/*
 * FastMarchingMethod.cc
 *
 *  Created on: 5/12/2013
 *      Author: banerjee
 */

#include <Contact/FastMarching/FastMarchingMethod.h>
#include <Exception.h>

using namespace BrMPM;

FastMarchingMethod::FastMarchingMethod()
{
}

FastMarchingMethod::~FastMarchingMethod()
{
}

void
FastMarchingMethod::distance(Double3DArray& phi, const Vector3D& dx,
                             bool self_test, int order,
                             Double3DArray& distance)
{
  // Check order is correct
  if (! (order == 1 || order == 2)) {
    throw Exception("Distance fast marching algorithm order must be 1 or 2", __FILE__, __LINE__);
  }

  // Check that dx is greater than 0
  if (!(dx > 0.0)) {
    throw Exception("Distance fast marching algorithm: dx must be greater than zero.", __FILE__, __LINE__);
  }

  // A flag array to match skfmm
  Double3DArray::extent_gen extents;
  Int3DArray flag(extents[phi.shape()[0]][phi.shape()[1]][phi.shape()[2]]);
  std::fill(flag.origin(), flag.origin()+flag.size(), 0);

  // Create a level set object to do the calculation
  BaseMarcher* marcher = new DistanceMarcher(phi.data(), dx, flag.data(), distance.data(), 3,
                                             phi.shape(), self_test, order);
  marcher->march();
  int error = marcher->getError();
  delete marcher;

  switch(error) {
  case 0:
    break;
  case 1:
    throw Exception("Unknown error is fast marching algorithm.", __FILE__, __LINE__);
  case 2:
    throw Exception("Fast marching algorithm: Array phi contains no zero contour.", __FILE__, __LINE__);
  }

}
