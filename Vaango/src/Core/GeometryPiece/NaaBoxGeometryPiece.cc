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

#include <algorithm>
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

double
NaaBoxGeometryPiece::getSDF(const Point& pt) const {
  Vector result = d_toUnitCube * (pt - p1_);
  
  // result is in [0,1]^3 if inside
  double d_x = std::max(-result.x(), result.x() - 1.0);
  double d_y = std::max(-result.y(), result.y() - 1.0);
  double d_z = std::max(-result.z(), result.z() - 1.0);
  
  // This is in unit cube space. To get real distance, we'd need to scale by 
  // the Jacobian or plane normals. 
  // Simplified: use max of these as a rough indicator, 
  // then multiply by a characteristic length.
  double max_d = std::max({d_x, d_y, d_z});
  Vector diag = d_boundingBox.upper() - d_boundingBox.lower();
  double volume = diag.x() * diag.y() * diag.z();
  double scale = std::pow(volume, 1.0/3.0);
  return max_d * scale;
}

Vector
NaaBoxGeometryPiece::getSDFGradient(const Point& pt) const {
  Vector result = d_toUnitCube * (pt - p1_);
  
  double d_x_low = -result.x();
  double d_x_high = result.x() - 1.0;
  double d_y_low = -result.y();
  double d_y_high = result.y() - 1.0;
  double d_z_low = -result.z();
  double d_z_high = result.z() - 1.0;

  double d_x = std::max(d_x_low, d_x_high);
  double d_y = std::max(d_y_low, d_y_high);
  double d_z = std::max(d_z_low, d_z_high);

  Vector grad_unit(0, 0, 0);
  if (d_x > d_y && d_x > d_z) {
    grad_unit = Vector(d_x_high > d_x_low ? 1 : -1, 0, 0);
  } else if (d_y > d_z) {
    grad_unit = Vector(0, d_y_high > d_y_low ? 1 : -1, 0);
  } else {
    grad_unit = Vector(0, 0, d_z_high > d_z_low ? 1 : -1);
  }

  // Transform unit normal back to world space
  // This is a simplification; for a general parallelepiped, 
  // the normal is the plane normal.
  Matrix3 inv = d_toUnitCube.Inverse();
  Vector grad = inv * grad_unit;
  if (grad.length2() > 1e-12) grad.normalize();
  return grad;
}

Box
NaaBoxGeometryPiece::getBoundingBox() const {
  return d_boundingBox;
}

double
NaaBoxGeometryPiece::volume() const {
  Vector v1 = p2_ - p1_;
  Vector v2 = p3_ - p1_;
  Vector v3 = p4_ - p1_;
  return std::abs(Dot(v1, Cross(v2, v3)));
}

Point
NaaBoxGeometryPiece::getCenter() const {
  Vector v1 = p2_ - p1_;
  Vector v2 = p3_ - p1_;
  Vector v3 = p4_ - p1_;
  return p1_ + (v1 + v2 + v3) * 0.5;
}

Matrix3
NaaBoxGeometryPiece::getInertiaTensor() const {
  double mass = volume();
  Vector v1 = p2_ - p1_;
  Vector v2 = p3_ - p1_;
  Vector v3 = p4_ - p1_;
  
  Matrix3 identity; identity.Identity();
  Matrix3 I = (identity * v1.length2() - Matrix3(v1, v1)) +
              (identity * v2.length2() - Matrix3(v2, v2)) +
              (identity * v3.length2() - Matrix3(v3, v3));
  return I * (mass / 12.0);
}

} // end namespace Uintah