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

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Geometry/Vector.h>
#include <Core/GeometryPiece/EllipsoidGeometryPiece.h>
#include <Core/Grid/Box.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/Matrix3.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <cmath>
#include <memory>

namespace Uintah {

const std::string EllipsoidGeometryPiece::TYPE_NAME = "ellipsoid";

EllipsoidGeometryPiece::EllipsoidGeometryPiece(ProblemSpecP& ps) {
  d_name = "Unnamed " + TYPE_NAME + " from PS";

  ps->require("origin", d_origin);

  // Get Vector axes
  Vector x_axis = Vector(1.0, 0.0, 0.0);
  Vector y_axis = Vector(0.0, 1.0, 0.0);
  ps->getWithDefault("v1", d_v1, x_axis);
  ps->getWithDefault("v2", d_v2, y_axis);

  // Get orthogonal axes
  ps->getWithDefault("r1", d_r1, 1.0);
  ps->getWithDefault("r2", d_r2, 1.0);
  ps->getWithDefault("r3", d_r3, 1.0);

  // Compute v3
  d_v3 = Cross(d_v1, d_v2);

  // Check the input
  checkEllipsoidData();

  // Normalize
  d_v1 = d_v1 / d_v1.length();
  d_v2 = d_v2 / d_v2.length();
  d_v3 = d_v3 / d_v3.length();

  // Rescale
  d_v1 = d_v1 * d_r1;
  d_v2 = d_v2 * d_r2;
  d_v3 = d_v3 * d_r3;
  std::cout << "v1 = " << d_v1 << " v2 = " << d_v2 << " v3 = " << d_v3 << "\n";
  std::cout << "r1 = " << d_r1 << " r2 = " << d_r2 << " r3 = " << d_r3 << "\n";

  // Run helper function to determine if inputs are correct
  initializeEllipsoidData();
}

EllipsoidGeometryPiece::EllipsoidGeometryPiece(const Point& origin,
                                               double radx,
                                               double rady,
                                               double radz) {
  d_origin = origin;
  d_r1     = radx;
  d_r2     = rady;
  d_r3     = radz;

  d_v1 = Vector(d_r1, 0.0, 0.0);
  d_v2 = Vector(0.0, d_r2, 0.0);
  d_v3 = Vector(0.0, 0.0, d_r3);

  checkEllipsoidData();
  initializeEllipsoidData();
}

EllipsoidGeometryPiece::EllipsoidGeometryPiece(const Point& origin,
                                               Vector one,
                                               Vector two,
                                               Vector three) {
  d_origin = origin;
  d_v1     = one;
  d_v2     = two;
  d_v3     = three;

  d_r1 = d_v1.length();
  d_r2 = d_v2.length();
  d_r3 = d_v3.length();

  checkEllipsoidData();
  initializeEllipsoidData();
}

void
EllipsoidGeometryPiece::checkEllipsoidData() {
  // Check lengths
  if (std::abs(d_v1.length()) < 1.0e-12 || std::abs(d_v2.length()) < 1.0e-12 ||
      std::abs(d_v3.length()) < 1.0e-12) {
    std::ostringstream err;
    err << "**ERROR** Input ellipsoid axis vectors have zero length.  Please "
           "check input file";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }

  // Check for orthogonality
  if ((std::abs(Dot(d_v1, d_v2)) >= 1e-12) ||
      (std::abs(Dot(d_v2, d_v3)) >= 1e-12) ||
      (std::abs(Dot(d_v3, d_v1)) >= 1e-12)) {
    std::ostringstream err;
    err << "**ERROR** Input File Error: (Ellipsoid initialization) input "
           "vectors (v1,v2,v3) \n"
           "are not orthogonal to within 1e-12 or each other.\n";
    err << "v1 . v2 = " << Dot(d_v1, d_v2) << " v2 . v3 = " << Dot(d_v2, d_v3)
        << " v3 . v1 = " << Dot(d_v3, d_v1) << "\n";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }
}

void
EllipsoidGeometryPiece::initializeEllipsoidData() {
  Vector x_axis(1.0, 0.0, 0.0);
  Vector y_axis(0.0, 1.0, 0.0);

  // First find rotation matrix that takes v1 to the x-axis
  Vector x_rot_axis  = Cross(d_v1 / d_r1, x_axis);
  double x_rot_angle = std::acos(Dot(d_v1 / d_r1, x_axis));
  Matrix3 R_x(x_rot_angle, x_rot_axis);

  // Next find rotation matrix that takes v2_rot to the y-axis
  Vector v2_rot      = R_x.Transpose() * d_v2;
  v2_rot             = v2_rot / v2_rot.length();
  Vector y_rot_axis  = Cross(v2_rot, y_axis);
  double y_rot_angle = std::acos(Dot(v2_rot, y_axis));
  Matrix3 R_y(y_rot_angle, y_rot_axis);

  // Total rotation to make the ellipsoid axes align with
  // x, y, z axes
  d_rotation = R_y.Transpose() * R_x.Transpose();
}

void
EllipsoidGeometryPiece::outputHelper(ProblemSpecP& ps) const {
  ps->appendElement("origin", d_origin);
  ps->appendElement("v1", d_v1);
  ps->appendElement("v2", d_v2);
  ps->appendElement("v3", d_v3);
  ps->appendElement("r1", d_r1);
  ps->appendElement("r2", d_r1);
  ps->appendElement("r3", d_r1);
}

GeometryPieceP
EllipsoidGeometryPiece::clone() const {
  return std::make_shared<EllipsoidGeometryPiece>(*this);
}

bool
EllipsoidGeometryPiece::inside(const Point& p) const {
  // Translate point
  Vector pTransformed = p - d_origin;

  // Rotate point
  pTransformed = d_rotation * pTransformed;

  // Check if inside unit distance from sphere center after scaling
  if (std::sqrt(pTransformed.x() * pTransformed.x() / (d_r1 * d_r1) +
                pTransformed.y() * pTransformed.y() / (d_r2 * d_r2) +
                pTransformed.z() * pTransformed.z() / (d_r3 * d_r3)) <= 1.0) {
    // std::cout << "p = " << p << " p_rot = " << pTransformed << "\n";
    return true;
  }
  return false;
}

Box
EllipsoidGeometryPiece::getBoundingBox() const {
  double highX = 0.0;
  double highY = 0.0;
  double highZ = 0.0;

  // Use vectors to find highest xyz
  // X
  highX = d_v1.x();
  if (std::abs(d_v2.x()) > highX) {
    highX = d_v2.x();
  }
  if (std::abs(d_v3.x()) > highX) {
    highX = d_v3.x();
  }
  // Y
  highY = d_v1.y();
  if (std::abs(d_v2.y()) > highY) {
    highY = d_v2.y();
  }
  if (std::abs(d_v3.y()) > highY) {
    highY = d_v3.y();
  }
  // X
  highZ = d_v1.z();
  if (std::abs(d_v2.z()) > highZ) {
    highZ = d_v2.z();
  }
  if (std::abs(d_v3.z()) > highZ) {
    highZ = d_v3.z();
  }

  Point low(d_origin.x() - std::abs(highX),
            d_origin.y() - std::abs(highY),
            d_origin.z() - std::abs(highZ));

  Point high(d_origin.x() + std::abs(highX),
             d_origin.y() + std::abs(highY),
             d_origin.z() + std::abs(highZ));

  return Box(low, high);
}

} // end namespace Uintah