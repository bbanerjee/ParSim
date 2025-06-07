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

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Geometry/Point.h>
#include <Core/GeometryPiece/DifferenceGeometryPiece.h>
#include <Core/GeometryPiece/GeometryPieceFactory.h>
#include <Core/Grid/Box.h>
#include <Core/Malloc/Allocator.h>

#include <sstream>
#include <vector>
#include <memory>

namespace Uintah {

const string DifferenceGeometryPiece::TYPE_NAME = "difference";

DifferenceGeometryPiece::DifferenceGeometryPiece(ProblemSpecP& ps) {
  d_name = "Unnamed " + TYPE_NAME + " from PS";
  std::vector<GeometryPieceP> objs;
  GeometryPieceFactory::create(ps, objs);

  if (objs.size() != 2) {
    std::ostringstream warn;
    warn << "\nERROR:Input File:Geom_object:Difference:  You need two "
            "geom_objects in order to take the difference of them";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  d_left  = objs[0];
  d_right = objs[1];
}

DifferenceGeometryPiece::DifferenceGeometryPiece(GeometryPieceP p1,
                                                 GeometryPieceP p2)
    : d_left(p1), d_right(p2) {
  d_name = "Unnamed " + TYPE_NAME + " from pieces";
}

DifferenceGeometryPiece::DifferenceGeometryPiece(
    const DifferenceGeometryPiece& rhs) {
  d_name = "Unnamed " + TYPE_NAME + " from CpyCnstr";

  d_left  = rhs.d_left->clone();
  d_right = rhs.d_right->clone();
}

DifferenceGeometryPiece&
DifferenceGeometryPiece::operator=(const DifferenceGeometryPiece& rhs) {
  if (this == &rhs) return *this;

  d_left  = rhs.d_left->clone();
  d_right = rhs.d_right->clone();

  return *this;
}

void
DifferenceGeometryPiece::outputHelper(ProblemSpecP& ps) const {
  d_left->outputProblemSpec(ps);
  d_right->outputProblemSpec(ps);
}

GeometryPieceP
DifferenceGeometryPiece::clone() const {
  return std::make_shared<DifferenceGeometryPiece>(*this);
}

bool
DifferenceGeometryPiece::inside(const Point& p) const {
  return (d_left->inside(p) && !d_right->inside(p));
}

Box
DifferenceGeometryPiece::getBoundingBox() const {
  // Initialize the lo and hi points to the left element
  Point d_leftlo  = d_left->getBoundingBox().lower();
  Point d_lefthi  = d_left->getBoundingBox().upper();
  Point d_rightlo = d_right->getBoundingBox().lower();
  Point d_righthi = d_right->getBoundingBox().upper();

  Point lo = Min(d_leftlo, d_rightlo);
  Point hi = Max(d_lefthi, d_righthi);

  return Box(lo, hi);
}

}  // end namespace Uintah