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

#ifndef __CORE_GEOMETRYPIECE_LOCALSDF_H__
#define __CORE_GEOMETRYPIECE_LOCALSDF_H__

#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Math/Matrix3.h>
#include <Core/Grid/Box.h>
#include <vector>
#include <memory>

namespace Uintah {

class LocalSDF {
public:
  LocalSDF(const Box& box, const IntVector& res);
  ~LocalSDF() = default;

  std::unique_ptr<LocalSDF> clone() const;

  // Set the distance value at a grid index
  void setDistance(const IntVector& idx, double dist);

  // Get interpolated distance at a local point
  double getDistance(const Point& p) const;

  // Get interpolated gradient at a local point
  Vector getGradient(const Point& p) const;

  // Coordinate transformations
  Point worldToLocal(const Point& worldP, const Point& center, const Matrix3& orientation) const;
  Vector localToWorld(const Vector& localV, const Matrix3& orientation) const;

  // Getters
  const Box& getBoundingBox() const { return d_box; }
  const IntVector& getResolution() const { return d_res; }

private:
  Box d_box;
  IntVector d_res;
  std::vector<double> d_distances;
  Vector d_spacing;

  inline size_t getIndex(const IntVector& idx) const {
    return (size_t)idx.x() + (size_t)d_res.x() * ((size_t)idx.y() + (size_t)d_res.y() * (size_t)idx.z());
  }
};

} // namespace Uintah

#endif
