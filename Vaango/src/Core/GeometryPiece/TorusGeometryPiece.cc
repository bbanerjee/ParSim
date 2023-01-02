/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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
#include <Core/GeometryPiece/TorusGeometryPiece.h>
#include <Core/Grid/Box.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <memory>

namespace Uintah {

const string TorusGeometryPiece::TYPE_NAME = "torus";

TorusGeometryPiece::TorusGeometryPiece() {
  d_name = "Unnamed " + TYPE_NAME + " from BasicCtor";

  d_center       = Point(0., 0., 0.);
  d_major_radius = 1.0;
  d_minor_radius = 0.0;
  d_axis_vec     = Vector(0, 0, 1.0);
  checkInput();
  computeRotation();
}

TorusGeometryPiece::TorusGeometryPiece(ProblemSpecP& ps) {
  d_name = "Unnamed " + TYPE_NAME + " from PS";

  ps->require("center", d_center);
  ps->require("major_radius", d_major_radius);
  ps->require("minor_radius", d_minor_radius);
  ps->require("axis_vector", d_axis_vec);

  checkInput();
  d_axis_vec /= d_axis_vec.length();
  computeRotation();
}

TorusGeometryPiece::TorusGeometryPiece(const Point& center,
                                       const Vector& axis_vec,
                                       const double major,
                                       const double minor) {
  d_name = "Unnamed " + TYPE_NAME + " from center/major/minor";

  d_center       = center;
  d_major_radius = major;
  d_minor_radius = minor;
  d_axis_vec     = axis_vec;

  checkInput();
  d_axis_vec /= d_axis_vec.length();
  computeRotation();
}

void
TorusGeometryPiece::checkInput() const {
  if (d_minor_radius <= 0.0) {
    std::ostringstream err;
    err << "**ERROR** Input File: Torus minor_radius must be > 0.0.";
    SCI_THROW(ProblemSetupException(err.str(), __FILE__, __LINE__));
  }
  if (d_major_radius <= 0.0) {
    std::ostringstream err;
    err << "**ERROR** Input File: Torus major_radius must be > 0.0.";
    SCI_THROW(ProblemSetupException(err.str(), __FILE__, __LINE__));
  }
  if (d_major_radius <= d_minor_radius) {
    std::ostringstream err;
    err << "**ERROR** Input File: Torus major_radius must be greater than "
           "torus minor radius.";
    SCI_THROW(ProblemSetupException(err.str(), __FILE__, __LINE__));
  }
  if (std::abs(d_axis_vec.length()) < 1.0e-12) {
    std::ostringstream err;
    err << "**ERROR** Input File: Torus axis vector has zero length.  Please "
           "check input file.";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }
}

void
TorusGeometryPiece::computeRotation() {
  // Find rotation matrix that takes the axis vector to the z-axis
  Vector z_axis(0.0, 0.0, 1.0);
  Vector z_rot_axis  = Cross(d_axis_vec, z_axis);
  double z_rot_angle = std::acos(Dot(d_axis_vec, z_axis));
  d_rotation         = Matrix3(z_rot_angle, z_rot_axis);
  d_rotation         = d_rotation.Transpose();
}

void
TorusGeometryPiece::outputHelper(ProblemSpecP& ps) const {
  ps->appendElement("center", d_center);
  ps->appendElement("major_radius", d_major_radius);
  ps->appendElement("minor_radius", d_minor_radius);
  ps->appendElement("axis_vector", d_axis_vec);
}

GeometryPieceP
TorusGeometryPiece::clone() const {
  return std::make_shared<TorusGeometryPiece>(*this);
}

bool
TorusGeometryPiece::inside(const Point& p) const {
  // Translate point
  Vector pTransformed = p - d_center;

  // Rotate point so that torus axis is along "z"
  pTransformed = d_rotation * pTransformed;

  double x = pTransformed.x();
  double y = pTransformed.y();
  double z = pTransformed.z();
  if ((d_major_radius - sqrt(x * x + y * y)) *
              (d_major_radius - sqrt(x * x + y * y)) +
          z * z <
      d_minor_radius * d_minor_radius) {
    return true;
  }
  return false;
}

Box
TorusGeometryPiece::getBoundingBox() const {
  // This is an overly generous bounding box.
  double R    = d_major_radius + d_minor_radius;
  Point minBB = Point(d_center.x() - R, d_center.y() - R, d_center.z() - R);
  Point maxBB = Point(d_center.x() + R, d_center.y() + R, d_center.z() + R);

  return Box(minBB, maxBB);
}

//////////
// Calculate the unit normal vector to axis from point
Vector
TorusGeometryPiece::radialDirection(const Point& pt) const {
  // Translate point
  Vector pNormal = pt - d_center;

  // Rotate point so that torus axis is along "z"
  pNormal = d_rotation * pNormal;

  // Normal vector has z = 0
  pNormal[2] = 0.0;

  // Rotate back
  pNormal = d_rotation.Transpose() * pNormal;
  pNormal /= pNormal.length();

  return pNormal;
}

}  // end namespace Uintah