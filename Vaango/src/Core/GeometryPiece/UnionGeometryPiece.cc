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
#include <Core/GeometryPiece/UnionGeometryPiece.h>
#include <Core/Grid/Box.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <memory>

namespace Uintah {

const std::string UnionGeometryPiece::TYPE_NAME = "union";

UnionGeometryPiece::UnionGeometryPiece(ProblemSpecP& ps) {
  d_name = "Unnamed " + TYPE_NAME + " from PS";
  // Need to loop through all the geometry pieces
  GeometryPieceFactory::create(ps, d_children);
}

UnionGeometryPiece::UnionGeometryPiece(const std::vector<GeometryPieceP>& child)
    : d_children(child) {
  d_name = "Unnamed " + TYPE_NAME + " from vector";
}

UnionGeometryPiece&
UnionGeometryPiece::operator=(const UnionGeometryPiece& rhs) {
  if (this == &rhs) return *this;

  d_children.clear();

  // Copy in the new values
  for (const auto& child : rhs.d_children) {
    d_children.push_back(child->clone());
  }

  return *this;
}

void
UnionGeometryPiece::outputHelper(ProblemSpecP& ps) const {

  // If this is a named object, then only output the children the first time.
  for (const auto& child : d_children) {
    child->outputProblemSpec(ps);
  }
}

GeometryPieceP
UnionGeometryPiece::clone() const {
  return std::make_shared<UnionGeometryPiece>(*this);
}

bool
UnionGeometryPiece::inside(const Point& p) const {
  for (const auto& child : d_children) {
    if (child->inside(p)) {
      return true;
    }
  }
  return false;
}

double
UnionGeometryPiece::getSDF(const Point& p) const {
  double min_sdf = 1.0e30;
  for (const auto& child : d_children) {
    min_sdf = std::min(min_sdf, child->getSDF(p));
  }
  return min_sdf;
}

Vector
UnionGeometryPiece::getSDFGradient(const Point& p) const {
  double min_sdf = 1.0e30;
  size_t min_idx = 0;
  for (size_t i = 0; i < d_children.size(); ++i) {
    double s = d_children[i]->getSDF(p);
    if (s < min_sdf) {
      min_sdf = s;
      min_idx = i;
    }
  }
  return d_children[min_idx]->getSDFGradient(p);
}

Box
UnionGeometryPiece::getBoundingBox() const {
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
UnionGeometryPiece::volume() const {
  double vol = 0;
  for (const auto& child : d_children) {
    vol += child->volume();
  }
  return vol;
}

Point
UnionGeometryPiece::getCenter() const {
  Point center(0, 0, 0);
  double total_vol = 0;
  for (const auto& child : d_children) {
    double vol = child->volume();
    center    = center + child->getCenter().asVector() * vol;
    total_vol += vol;
  }
  if (total_vol > 1e-12) {
    center = center / total_vol;
  }
  return center;
}

Matrix3
UnionGeometryPiece::getInertiaTensor() const {
  Matrix3 I(0.0);
  Point center = getCenter();
  for (const auto& child : d_children) {
    double vol = child->volume();
    Point c_i = child->getCenter();
    Vector d = c_i - center;
    Matrix3 I_i = child->getInertiaTensor();
    
    // Parallel axis theorem: I = I_cm + M * (d.d * I - d cross d)
    // Here I is the identity matrix, and d cross d is dyadic product d*d^T
    double d_dot_d = d.length2();
    Matrix3 d_dyad_d(d, d);
    Matrix3 identity; identity.Identity();
    Matrix3 I_shifted = I_i + (identity * d_dot_d - d_dyad_d) * vol;
    I += I_shifted;
  }
  return I;
}

}  // end namespace Uintah