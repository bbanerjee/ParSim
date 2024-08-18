/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

} // end namespace Uintah