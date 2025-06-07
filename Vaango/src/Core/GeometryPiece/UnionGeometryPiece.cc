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

}  // end namespace Uintah