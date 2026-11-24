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
#include <Core/Geometry/Vector.h>
#include <Core/GeometryPiece/CylinderGeometryPiece.h>
#include <Core/Grid/Box.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <memory>

namespace Uintah {

const string CylinderGeometryPiece::TYPE_NAME = "cylinder";

CylinderGeometryPiece::CylinderGeometryPiece() {
  d_name = "Unnamed " + TYPE_NAME + " from BasicCtor";

  Point top, bottom;
  d_bottom            = bottom;
  d_top               = top;
  d_radius            = 0.0;
  d_cylinder_end      = false;
  d_axisymmetric_end  = false;
  d_axisymmetric_side = false;
}

CylinderGeometryPiece::CylinderGeometryPiece(ProblemSpecP& ps) {
  d_name = "Unnamed " + TYPE_NAME + " from PS";
  Point top, bottom;
  double rad;

  ps->require("bottom", bottom);
  ps->require("top", top);
  ps->require("radius", rad);
  ps->getWithDefault("cylinder_end", d_cylinder_end, false);
  ps->getWithDefault("axisymmetric_end", d_axisymmetric_end, false);
  ps->getWithDefault("axisymmetric_side", d_axisymmetric_side, false);

  double near_zero = 1e-100;
  Vector axis      = top - bottom;

  if (axis.length() < near_zero) {
    SCI_THROW(ProblemSetupException(
        "Input File Error: Cylinder axes has zero length", __FILE__, __LINE__));
  }
  if (rad <= 0.0) {
    SCI_THROW(ProblemSetupException(
        "Input File Error: Cylinder radius must be > 0.0", __FILE__, __LINE__));
  }
  d_bottom = bottom;
  d_top    = top;
  d_radius = rad;
}

CylinderGeometryPiece::CylinderGeometryPiece(const Point& top,
                                             const Point& bottom,
                                             double radius) {
  d_name = "Unnamed " + TYPE_NAME + " from top/bottom/radius";

  double near_zero = 1e-100;
  Vector axis      = top - bottom;

  if (axis.length() < near_zero) {
    SCI_THROW(ProblemSetupException(
        "Input File Error: Cylinder axes has zero length", __FILE__, __LINE__));
  }
  if (radius <= 0.0) {
    SCI_THROW(ProblemSetupException(
        "Input File Error: Cylinder radius must be > 0.0", __FILE__, __LINE__));
  }
  d_bottom            = bottom;
  d_top               = top;
  d_radius            = radius;
  d_cylinder_end      = false;
  d_axisymmetric_end  = false;
  d_axisymmetric_side = false;
}

void
CylinderGeometryPiece::outputHelper(ProblemSpecP& ps) const {
  ps->appendElement("bottom", d_bottom);
  ps->appendElement("top", d_top);
  ps->appendElement("radius", d_radius);
  if (d_cylinder_end || d_axisymmetric_end || d_axisymmetric_side) {
    ps->appendElement("cylinder_end", d_cylinder_end);
    ps->appendElement("axisymmetric_end", d_axisymmetric_end);
    ps->appendElement("axisymmetric_side", d_axisymmetric_side);
  }
}

GeometryPieceP
CylinderGeometryPiece::clone() const {
  return std::make_shared<CylinderGeometryPiece>(*this);
}

bool
CylinderGeometryPiece::inside(const Point& p) const {
  Vector axis    = d_top - d_bottom;
  double height2 = axis.length2();

  Vector tobot = p - d_bottom;

  // pt is the "test" point
  double h = Dot(tobot, axis);
  if (h < 0.0 || h > height2) return false;  // Above or below the cylinder

  double area = Cross(axis, tobot).length2();
  double d    = area / height2;
  if (d > d_radius * d_radius) return false;
  return true;
}

double
CylinderGeometryPiece::getSDF(const Point& p) const {
  Vector ba = d_top - d_bottom;
  Vector pa = p - d_bottom;
  double l_ba = ba.length();
  double h = Dot(pa, ba) / (l_ba * l_ba);
  Vector rel_proj = pa - ba * h;
  double dist_x = rel_proj.length() - d_radius;
  double dist_z = std::abs(h - 0.5) * l_ba - 0.5 * l_ba;
  
  double outside_dist = Vector(std::max(dist_x, 0.0), std::max(dist_z, 0.0), 0.0).length();
  double inside_dist = std::min(std::max(dist_x, dist_z), 0.0);
  
  return outside_dist + inside_dist;
}

Vector
CylinderGeometryPiece::getSDFGradient(const Point& p) const {
  Vector ba = d_top - d_bottom;
  Vector pa = p - d_bottom;
  double l_ba = ba.length();
  double h = Dot(pa, ba) / (l_ba * l_ba);
  Vector rel_proj = pa - ba * h;
  double r = rel_proj.length();
  
  double dist_x = r - d_radius;
  double dist_z = std::abs(h - 0.5) * l_ba - 0.5 * l_ba;
  
  if (std::max(dist_x, dist_z) > 0) {
    // Outside
    Vector grad_radial = (r > 1e-12) ? (rel_proj / r) : Vector(1, 0, 0);
    Vector grad_axial = (h > 0.5) ? (ba / l_ba) : (-ba / l_ba);
    Vector total_grad(0, 0, 0);
    if (dist_x > 0) total_grad = total_grad + grad_radial * dist_x;
    if (dist_z > 0) total_grad = total_grad + grad_axial * dist_z;
    if (total_grad.length2() > 1e-12) total_grad.normalize();
    return total_grad;
  } else {
    // Inside
    if (dist_x > dist_z) return (r > 1e-12) ? (rel_proj / r) : Vector(1, 0, 0);
    return (h > 0.5) ? (ba / l_ba) : (-ba / l_ba);
  }
}

Box
CylinderGeometryPiece::getBoundingBox() const {
  Point minBB = Min(d_bottom, d_top);
  Point maxBB = Max(d_bottom, d_top);

  double x_sqrd   = pow((d_bottom.x() - d_top.x()), 2);
  double y_sqrd   = pow((d_bottom.y() - d_top.y()), 2);
  double z_sqrd   = pow((d_bottom.z() - d_top.z()), 2);
  double all_sqrd = x_sqrd + y_sqrd + z_sqrd;

  double kx = sqrt((y_sqrd + z_sqrd) / all_sqrd);
  double ky = sqrt((x_sqrd + z_sqrd) / all_sqrd);
  double kz = sqrt((x_sqrd + y_sqrd) / all_sqrd);

  Vector tmp(kx * d_radius, ky * d_radius, kz * d_radius);
  minBB -= tmp;
  maxBB += tmp;

  return Box(minBB, maxBB);
}

Point
CylinderGeometryPiece::getCenter() const {
  return d_bottom + (d_top - d_bottom) * 0.5;
}

Matrix3
CylinderGeometryPiece::getInertiaTensor() const {
  double r = d_radius;
  double h = height();
  double mass = volume(); // Unit density
  
  double izz_local = 0.5 * mass * r * r;
  double ixx_local = (1.0 / 12.0) * mass * (3.0 * r * r + h * h);
  double iyy_local = ixx_local;
  
  Matrix3 I_local(ixx_local, 0, 0, 0, iyy_local, 0, 0, 0, izz_local);
  
  Vector e3 = (d_top - d_bottom);
  if (h > 1e-12) {
    e3 /= h;
  } else {
    return I_local; // Degenerate case
  }
  
  Vector e1, e2;
  e3.findOrthogonal(e1, e2);
  
  Matrix3 R(e1[0], e2[0], e3[0],
            e1[1], e2[1], e3[1],
            e1[2], e2[2], e3[2]);
            
  return rotateInertiaTensor(I_local, R);
}

//////////
// Calculate the unit normal vector to axis from point
Vector
CylinderGeometryPiece::radialDirection(const Point& pt) const {
  Vector axis       = d_top - d_bottom;
  double height2    = axis.length();
  Vector pbot       = pt - d_bottom;
  double tt         = Dot(pbot, axis) / height2;
  Vector projOnAxis = tt * axis / height2;
  Vector normal     = pbot - projOnAxis;
  return (normal / normal.length());
}

} // end namespace Uintah