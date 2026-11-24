/*
 * The MIT License
 *
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#include <Core/GeometryPiece/LocalSDF.h>
#include <Core/Math/Matrix3.h>
#include <cmath>
#include <algorithm>

namespace Uintah {

LocalSDF::LocalSDF(const Box& box, const IntVector& res)
  : d_box(box), d_res(res)
{
  d_distances.resize(res.x() * res.y() * res.z(), 1.0e10);
  Vector range = box.upper() - box.lower();
  d_spacing = Vector(range.x() / (res.x() - 1), 
                     range.y() / (res.y() - 1), 
                     range.z() / (res.z() - 1));
}

std::unique_ptr<LocalSDF> LocalSDF::clone() const {
  auto newSDF = std::make_unique<LocalSDF>(d_box, d_res);
  newSDF->d_distances = d_distances;
  newSDF->d_spacing = d_spacing;
  return newSDF;
}

void LocalSDF::setDistance(const IntVector& idx, double dist) {
  d_distances[getIndex(idx)] = dist;
}

double LocalSDF::getDistance(const Point& p) const {
  Vector localP = p - d_box.lower();
  
  // Grid indices
  double fx = localP.x() / d_spacing.x();
  double fy = localP.y() / d_spacing.y();
  double fz = localP.z() / d_spacing.z();

  int i = std::clamp((int)std::floor(fx), 0, d_res.x() - 2);
  int j = std::clamp((int)std::floor(fy), 0, d_res.y() - 2);
  int k = std::clamp((int)std::floor(fz), 0, d_res.z() - 2);

  double tx = fx - i;
  double ty = fy - j;
  double tz = fz - k;

  // Tri-linear interpolation
  auto d = [&](int ii, int jj, int kk) {
    return d_distances[getIndex(IntVector(ii, jj, kk))];
  };

  double c00 = d(i, j, k) * (1 - tx) + d(i + 1, j, k) * tx;
  double c01 = d(i, j, k + 1) * (1 - tx) + d(i + 1, j, k + 1) * tx;
  double c10 = d(i, j + 1, k) * (1 - tx) + d(i + 1, j + 1, k) * tx;
  double c11 = d(i, j + 1, k + 1) * (1 - tx) + d(i + 1, j + 1, k + 1) * tx;

  double c0 = c00 * (1 - ty) + c10 * ty;
  double c1 = c01 * (1 - ty) + c11 * ty;

  return c0 * (1 - tz) + c1 * tz;
}

Vector LocalSDF::getGradient(const Point& p) const {
  double eps = std::min({d_spacing.x(), d_spacing.y(), d_spacing.z()}) * 0.1;
  double dx = (getDistance(p + Vector(eps, 0, 0)) - getDistance(p - Vector(eps, 0, 0))) / (2 * eps);
  double dy = (getDistance(p + Vector(0, eps, 0)) - getDistance(p - Vector(0, eps, 0))) / (2 * eps);
  double dz = (getDistance(p + Vector(0, 0, eps)) - getDistance(p - Vector(0, 0, eps))) / (2 * eps);
  Vector grad(dx, dy, dz);
  if (grad.length2() > 1e-12) {
    grad.normalize();
  }
  return grad;
}

Point LocalSDF::worldToLocal(const Point& worldP, const Point& center, const Matrix3& orientation) const {
  Vector rel = worldP - center;
  // orientation is R, we need R^T * rel
  // Assuming orientation stores columns? Matrix3 in Uintah is usually row-major access but let's check.
  // Actually, we can use the Transpose.
  Matrix3 rotT = orientation.Transpose();
  return Point(0,0,0) + rotT * rel;
}

Vector LocalSDF::localToWorld(const Vector& localV, const Matrix3& orientation) const {
  return orientation * localV;
}

} // namespace Uintah
