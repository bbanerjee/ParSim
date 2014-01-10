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
FastMarchingMethod::distance(Double3DArray& phi, Double3DArray& speed, Vector3D& dx,
                             Int3DArray& flag, bool self_test, int order,
                             Double3DArray& distance)
{
  // Check order is correct
  if (! (order == 1 || order ==1)) {
    throw Exception("Distance fast marching algorithm order must be 1 or 2", __FILE__, __LINE__);
  }

  // Check that phi and speed have the same dimensions
  if (phi.shape() != speed.shape()) {
    throw Exception("Distance fast marching algorithm phi and speed must have the same dimensions", __FILE__, __LINE__);
  }

  // Check that dx is greater than 0
  if (!(dx > 0.0)) {
    throw Exception("Distance fast marching algorithm: dx must be greater than zero.", __FILE__, __LINE__);
  }

  // Create the distance array and initialize to 0.0
  //Double3DArray distance(phi);
  std::fill(distance.origin(), distance.origin() + distance.size(), 0.0);

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
