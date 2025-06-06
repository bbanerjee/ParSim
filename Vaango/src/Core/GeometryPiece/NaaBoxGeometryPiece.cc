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
#include <Core/GeometryPiece/NaaBoxGeometryPiece.h>
#include <Core/Grid/Box.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Util/DebugStream.h>

#include <iostream>
#include <memory>

namespace Uintah {

static DebugStream dbg("GeometryPiece", false);

const std::string NaaBoxGeometryPiece::TYPE_NAME = "parallelepiped";

NaaBoxGeometryPiece::NaaBoxGeometryPiece(ProblemSpecP& ps) {
  string gp_label = "Unamed";

  if (!ps->getAttribute("label", gp_label)) {
    // "label" and "name" are both used... so check for "label" first, and if it
    // isn't found, then check for "name".
    ps->getAttribute("name", gp_label);
  }

  d_name = gp_label + " " + TYPE_NAME + " from PS";

  Point p1, p2, p3, p4;
  ps->require("p1", p1);
  ps->require("p2", p2);
  ps->require("p3", p3);
  ps->require("p4", p4);

  init(p1, p2, p3, p4);
}

NaaBoxGeometryPiece::NaaBoxGeometryPiece(const Point& p1,
                                         const Point& p2,
                                         const Point& p3,
                                         const Point& p4) {
  d_name = "Unnamed " + TYPE_NAME + " from points";
  init(p1, p2, p3, p4);
}

void
NaaBoxGeometryPiece::init(const Point& p1,
                          const Point& p2,
                          const Point& p3,
                          const Point& p4) {
  p1_ = p1;
  p2_ = p2;
  p3_ = p3;
  p4_ = p4;

  Vector p2minusP1, p3minusP1, p4minusP1;
  p2minusP1 = p2 - p1;
  p3minusP1 = p3 - p1;
  p4minusP1 = p4 - p1;

  // p5 is the opposite corner to p1 and is used for the bounding box.
  // Point p5 = p1 + (p2minusP1 + p3minusP1 + p4minusP1);

  // Find the bounding box with the following gross code
  double lowX = std::min(
      std::min(std::min(p1.x(), p2.x()), std::min(p2.x(), p3.x())), p4.x());
  double lowY = std::min(
      std::min(std::min(p1.y(), p2.y()), std::min(p2.y(), p3.y())), p4.y());
  double lowZ = std::min(
      std::min(std::min(p1.z(), p2.z()), std::min(p2.z(), p3.z())), p4.z());
  double highX = std::max(
      std::max(std::max(p1.x(), p2.x()), std::max(p2.x(), p3.x())), p4.x());
  double highY = std::max(
      std::max(std::max(p1.y(), p2.y()), std::max(p2.y(), p3.y())), p4.y());
  double highZ = std::max(
      std::max(std::max(p1.z(), p2.z()), std::max(p2.z(), p3.z())), p4.z());

  Point blow  = Point(lowX, lowY, lowZ);
  Point bhigh = Point(highX, highY, highZ);

  d_boundingBox = Box(blow, bhigh);

  if (d_boundingBox.degenerate()) {
    // 1st point must be '<' second point, so flip them.
    d_boundingBox.fixBoundingBox();
    if (d_boundingBox.degenerate()) {
      // If there are still problems, throw an exception...

      std::ostringstream error;
      error << "NaaBoxGeometryPiece.cc: d_boundingBox for '" + d_name +
                   "' is degenerate..."
            << d_boundingBox << "\n";
      error << "See src/Core/GeometryPiece/NaaBoxGeometryPiece.h or the Users "
               "Guide for details\n";

      throw ProblemSetupException(error.str(), __FILE__, __LINE__);
    }
  }

  dbg << "Creating NaaBoxx with BBox of: " << d_boundingBox << "\n";

  // Map the arbitrary box to a unit cube...
  Matrix3 mat(p2minusP1.x(),
              p3minusP1.x(),
              p4minusP1.x(),
              p2minusP1.y(),
              p3minusP1.y(),
              p4minusP1.y(),
              p2minusP1.z(),
              p3minusP1.z(),
              p4minusP1.z());

  d_toUnitCube = mat.Inverse();
}

void
NaaBoxGeometryPiece::outputHelper(ProblemSpecP& ps) const {
  ps->appendElement("p1", p1_);
  ps->appendElement("p2", p2_);
  ps->appendElement("p3", p3_);
  ps->appendElement("p4", p4_);
}

GeometryPieceP
NaaBoxGeometryPiece::clone() const {
  return std::make_shared<NaaBoxGeometryPiece>(*this);
}

//********************************************
//                                          //
//             *-------------*              //
//            / .           / \             //
//           /   .         /   \            //
//          P4-------------*    \           //
//           \    .         \    \          //
//            \   P2.........\....*         //
//             \ .            \  /          //
//             P1--------------P3           //
//
//  Returns true if the point is inside (or on) the parallelepiped.
//  (The order of p2, p3, and p4 don't really matter.)
//
//  The arbitrary box has been transformed into a unit cube... we take
//  the Point to check and transform it the same way, then just check
//  to see if the Pt is in the unit cube.
//
bool
NaaBoxGeometryPiece::inside(const Point& pt) const {
  Vector result = d_toUnitCube * (pt - p1_);

  if ((result.minComponent() > 0) && (result.maxComponent() <= 1.0))
    return true;
  else
    return false;
}

Box
NaaBoxGeometryPiece::getBoundingBox() const {
  return d_boundingBox;
}

} // end namespace Uintah