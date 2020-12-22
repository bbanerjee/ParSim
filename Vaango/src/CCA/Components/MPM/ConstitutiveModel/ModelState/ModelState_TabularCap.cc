/*
 * The MIT License
 *
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_TabularCap.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/Constants.h>
#include <Core/Exceptions/InternalError.h>
#include <iostream>

using namespace Vaango;

ModelState_TabularCap::ModelState_TabularCap() 
 : ModelState_Tabular()
 , capX(0)
 , I1_min(0)
 , I1_max(0)
 , sqrtJ2_max(0)
{
}

ModelState_TabularCap::ModelState_TabularCap(const ModelState_TabularCap* state)
{
  *this = *state;
}

ModelState_TabularCap*
ModelState_TabularCap::operator=(const ModelState_TabularCap* state)
{
  if (this == state)
    return this;
  
  *this = *state;

  return this;
}

void 
ModelState_TabularCap::updateYieldSurface(const Polyline& yield_poly)
{

  // Copy the polyline
  yield_f_pts = yield_poly;

  // Convert tabular data to z-rprime coordinates
  convertToZRprime();

  // Save as a point cloud
  z_r_cloud = std::make_shared<Util::PolylinePointCloud>(z_r_table);

  // Build index
  z_r_index = std::make_shared<Util::PolylineKDTree>(2 /*dim*/, *z_r_cloud,
                nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  z_r_index->buildIndex();
}

/* Convert yield function data to z_rprime coordinates */
void
ModelState_TabularCap::convertToZRprime()
{
  // Get the bulk and shear moduli and compute sqrt(3/2 K/G)
  double sqrtKG = std::sqrt(1.5 * bulkModulus / shearModulus);

  // Compute z and r' for the yield surface points
  for (const auto& pt : yield_f_pts) {
    double p_bar   = pt.x();
    double sqrt_J2 = pt.y();
    double z = -Util::sqrt_three * p_bar;
    double r_prime = Util::sqrt_two * sqrt_J2 * sqrtKG;
    z_r_table.emplace_back(Uintah::Point(z, r_prime, 0));
  }
}

