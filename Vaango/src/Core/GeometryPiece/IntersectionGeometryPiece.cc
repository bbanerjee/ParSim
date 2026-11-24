/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <Core/Geometry/Point.h>
#include <Core/GeometryPiece/GeometryPieceFactory.h>
#include <Core/GeometryPiece/IntersectionGeometryPiece.h>
#include <Core/Grid/Box.h>
#include <Core/Malloc/Allocator.h>

#include <memory>

namespace Uintah {

const std::string IntersectionGeometryPiece::TYPE_NAME = "intersection";

IntersectionGeometryPiece::IntersectionGeometryPiece(ProblemSpecP& ps) {
  d_name = "Unnamed Intersection";
  GeometryPieceFactory::create(ps, d_children);
}

IntersectionGeometryPiece::IntersectionGeometryPiece(
    const IntersectionGeometryPiece& rhs) {
  for (auto& child : rhs.d_children) {
    d_children.emplace_back(child->clone());
  }
}

IntersectionGeometryPiece&
IntersectionGeometryPiece::operator=(const IntersectionGeometryPiece& rhs) {
  if (this == &rhs) return *this;
  d_children.clear();
  for (auto& child : rhs.d_children) {
    d_children.emplace_back(child->clone());
  }
  return *this;
}

void
IntersectionGeometryPiece::outputHelper(ProblemSpecP& ps) const {
  for (auto& child : d_children) {
    child->outputProblemSpec(ps);
  }
}

GeometryPieceP
IntersectionGeometryPiece::clone() const {
  return std::make_shared<IntersectionGeometryPiece>(*this);
}

bool
IntersectionGeometryPiece::inside(const Point& p) const {
  for (const auto& child : d_children) {
    if (!child->inside(p)) return false;
  }
  return true;
}

double
IntersectionGeometryPiece::getSDF(const Point& p) const {
  double max_sdf = -1.0e30;
  for (const auto& child : d_children) {
    max_sdf = std::max(max_sdf, child->getSDF(p));
  }
  return max_sdf;
}

Vector
IntersectionGeometryPiece::getSDFGradient(const Point& p) const {
  double max_sdf = -1.0e30;
  size_t max_idx = 0;
  for (size_t i = 0; i < d_children.size(); ++i) {
    double s = d_children[i]->getSDF(p);
    if (s > max_sdf) {
      max_sdf = s;
      max_idx = i;
    }
  }
  return d_children[max_idx]->getSDFGradient(p);
}

Box
IntersectionGeometryPiece::getBoundingBox() const {
  Point lo, hi;

  // Initialize the lo and hi points to the first element
  lo = d_children[0]->getBoundingBox().lower();
  hi = d_children[0]->getBoundingBox().upper();

  for (const auto& child : d_children) {
    Box box = child->getBoundingBox();
    lo      = Min(lo, box.lower());
    hi      = Max(hi, box.upper());
  }

  return Box(lo, hi);
}

double
IntersectionGeometryPiece::volume() const {
  double min_vol = 1.0e30;
  for (const auto& child : d_children) {
    min_vol = std::min(min_vol, child->volume());
  }
  return min_vol;
}

Point
IntersectionGeometryPiece::getCenter() const {
  double min_vol = 1.0e30;
  size_t min_idx = 0;
  for (size_t i = 0; i < d_children.size(); ++i) {
    double vol = d_children[i]->volume();
    if (vol < min_vol) {
      min_vol = vol;
      min_idx = i;
    }
  }
  return d_children[min_idx]->getCenter();
}

Matrix3
IntersectionGeometryPiece::getInertiaTensor() const {
  double min_vol = 1.0e30;
  size_t min_idx = 0;
  for (size_t i = 0; i < d_children.size(); ++i) {
    double vol = d_children[i]->volume();
    if (vol < min_vol) {
      min_vol = vol;
      min_idx = i;
    }
  }
  return d_children[min_idx]->getInertiaTensor();
}

} // end namespace Uintah