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
#include <Core/GeometryPiece/BoxGeometryPiece.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <memory>

#ifndef DMIN
#define DMIN(a, b) (a < b ? a : b)
#endif

namespace Uintah {

const std::string BoxGeometryPiece::TYPE_NAME = "box";

BoxGeometryPiece::BoxGeometryPiece(ProblemSpecP& ps) {
  d_name = "Unnamed " + TYPE_NAME + " from PS";

  Point min, max;
  ps->require("min", min);
  ps->require("max", max);

  double near_zero = 1e-100;
  double xdiff     = max.x() - min.x();
  double ydiff     = max.y() - min.y();
  double zdiff     = max.z() - min.z();

  if (xdiff < near_zero || ydiff < near_zero || zdiff < near_zero) {
    std::ostringstream warn;
    warn << "Input File Error: box max " << max << " <= min " << min
         << " coordinates";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  d_box = Box(min, max);
}

BoxGeometryPiece::BoxGeometryPiece(const Point& p1, const Point& p2)
    : d_box(Min(p1, p2), Max(p1, p2)) {
  d_name = "Unnamed " + TYPE_NAME + " from p1,p2";

  if (d_box.degenerate())
    throw ProblemSetupException("degenerate box", __FILE__, __LINE__);
}

void
BoxGeometryPiece::outputHelper(ProblemSpecP& ps) const {
  ps->appendElement("min", d_box.lower());
  ps->appendElement("max", d_box.upper());
}

GeometryPieceP
BoxGeometryPiece::clone() const {
  return std::make_shared<BoxGeometryPiece>(*this);
}

bool
BoxGeometryPiece::inside(const Point& p) const {
  // Check p with the lower coordinates

  if (p == Max(p, d_box.lower()) && p == Min(p, d_box.upper()))
    return true;
  else
    return false;
}

double
BoxGeometryPiece::getSDF(const Point& p) const {
  Point center = d_box.lower() + (d_box.upper() - d_box.lower()) * 0.5;
  Vector b = (d_box.upper() - d_box.lower()) * 0.5;
  Vector p_v = p - center;
  Vector q = Vector(std::abs(p_v.x()), std::abs(p_v.y()), std::abs(p_v.z())) - b;
  
  Vector q_max = Vector(std::max(q.x(), 0.0), std::max(q.y(), 0.0), std::max(q.z(), 0.0));
  double outside_dist = q_max.length();
  double inside_dist = std::min(std::max({q.x(), q.y(), q.z()}), 0.0);
  
  return outside_dist + inside_dist;
}

Vector
BoxGeometryPiece::getSDFGradient(const Point& p) const {
  Point center = d_box.lower() + (d_box.upper() - d_box.lower()) * 0.5;
  Vector b = (d_box.upper() - d_box.lower()) * 0.5;
  Vector p_v = p - center;
  Vector q = Vector(std::abs(p_v.x()), std::abs(p_v.y()), std::abs(p_v.z())) - b;
  
  if (std::max({q.x(), q.y(), q.z()}) > 0) {
    // Outside
    Vector grad(
      (p_v.x() > 0 ? 1 : -1) * std::max(q.x(), 0.0),
      (p_v.y() > 0 ? 1 : -1) * std::max(q.y(), 0.0),
      (p_v.z() > 0 ? 1 : -1) * std::max(q.z(), 0.0)
    );
    if (grad.length2() > 1e-12) grad.normalize();
    return grad;
  } else {
    // Inside: normal to closest face
    if (q.x() > q.y() && q.x() > q.z()) return Vector(p_v.x() > 0 ? 1 : -1, 0, 0);
    if (q.y() > q.z()) return Vector(0, p_v.y() > 0 ? 1 : -1, 0);
    return Vector(0, 0, p_v.z() > 0 ? 1 : -1);
  }
}

Box
BoxGeometryPiece::getBoundingBox() const {
  return d_box;
}

Point
BoxGeometryPiece::getCenter() const {
  return d_box.lower() + (d_box.upper() - d_box.lower()) * 0.5;
}

Matrix3
BoxGeometryPiece::getInertiaTensor() const {
  double w = (d_box.upper()).x() - (d_box.lower()).x();
  double h = (d_box.upper()).y() - (d_box.lower()).y();
  double d = (d_box.upper()).z() - (d_box.lower()).z();
  double mass = w * h * d; // Unit density
  double ixx = (1.0 / 12.0) * mass * (h * h + d * d);
  double iyy = (1.0 / 12.0) * mass * (w * w + d * d);
  double izz = (1.0 / 12.0) * mass * (w * w + h * h);
  return Matrix3(ixx, 0, 0, 0, iyy, 0, 0, 0, izz);
}

double
BoxGeometryPiece::volume() const {
  double dx = (d_box.upper()).x() - (d_box.lower()).x();
  double dy = (d_box.upper()).y() - (d_box.lower()).y();
  double dz = (d_box.upper()).z() - (d_box.lower()).z();
  return (dx * dy * dz);
}

double
BoxGeometryPiece::smallestSide() const {
  double dx = (d_box.upper()).x() - (d_box.lower()).x();
  double dy = (d_box.upper()).y() - (d_box.lower()).y();
  double dz = (d_box.upper()).z() - (d_box.lower()).z();
  return DMIN(DMIN(dx, dy), dz);
}

unsigned int
BoxGeometryPiece::thicknessDirection() const {
  double dx  = (d_box.upper()).x() - (d_box.lower()).x();
  double dy  = (d_box.upper()).y() - (d_box.lower()).y();
  double dz  = (d_box.upper()).z() - (d_box.lower()).z();
  double min = DMIN(DMIN(dx, dy), dz);
  if (dx == min)
    return 0;
  else if (dy == min)
    return 1;
  return 2;
}

} // end namespace Uintah