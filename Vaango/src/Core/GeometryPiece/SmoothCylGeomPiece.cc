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
#include <Core/GeometryPiece/SmoothCylGeomPiece.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Patch.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/Matrix3.h>
#include <Core/Parallel/Parallel.h>  // for proc0cout
#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>

namespace Uintah {

const std::string SmoothCylGeomPiece::TYPE_NAME = "smoothcyl";

SmoothCylGeomPiece::SmoothCylGeomPiece(ProblemSpecP& ps) {
  ps->require("bottom", d_bottom);
  ps->require("top", d_top);
  ps->require("radius", d_radius);
  d_numRadial = 5;
  d_numAxial  = 5;
  ps->getWithDefault("thickness", d_thickness, d_radius);
  ps->getWithDefault("endcap_thickness", d_capThick, 0.0);

  d_arcStart      = 0.0;
  double arcStart = 0.0;
  ps->getWithDefault("arc_start_angle_degree", arcStart, 0.0);
  if (arcStart > 0.0) d_arcStart = (M_PI / 180.0) * arcStart;

  d_angle      = 2.0 * M_PI;
  double angle = 360.0;
  ps->getWithDefault("arc_angle_degree", angle, 360.0);
  if (angle > 0.0) d_angle = (M_PI / 180.0) * angle;

  d_fileName = "none";
  ps->get("output_file", d_fileName);

  /* Save the domain limits.  This is needed so that points are not
     created outside the domain, leading to errors in the application of
     boundary conditions */
  // BBox domain;
  // grid->getSpatialRange(domain);
  // d_domainMin = domain.min();
  // d_domainMax = domain.max();

  checkInput();
  computeRotation();
}

void
SmoothCylGeomPiece::checkInput() const {
  std::ostringstream msg;
  if ((d_top - d_bottom).length2() <= 1.0e-12) {
    msg << "**ERROR** SmoothCylGeom: Top and bottom of cylinder are the same";
    throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
  }

  if (d_radius <= 0.0) {
    msg << "**ERROR** SmoothCylGeom: Radius < 0";
    throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
  }

  if (d_thickness > d_radius) {
    msg << "**ERROR** SmoothCylGeom: Thickness > Radius";
    throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
  }

  if (d_capThick < 0.0) {
    msg << "**ERROR** SmoothCylGeom: Cap Thickness < 0.0";
    throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
  }

  if (d_arcStart < 0.0 || d_arcStart > 2.0 * M_PI) {
    msg << "**ERROR** SmoothCylGeom: Arc Start Angle < 0.0 || > 360 degrees";
    throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
  }

  if (d_angle < 0.0 || d_angle > 2.0 * M_PI) {
    msg << "**ERROR** SmoothCylGeom: Total arc angle < 0.0 || > 360 degrees";
    throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
  }
}

void
SmoothCylGeomPiece::computeRotation() {
  // Compute axis vector
  d_axis   = d_top - d_bottom;
  d_height = d_axis.length();
  d_axis /= d_height;

  // Find rotation matrix that takes the axis vector to the z-axis
  Vector z_axis(0.0, 0.0, 1.0);
  Vector z_rot_axis  = Cross(d_axis, z_axis);
  double z_rot_angle = std::acos(Dot(d_axis, z_axis));
  d_rotation         = Matrix3(z_rot_angle, z_rot_axis);
}

void
SmoothCylGeomPiece::outputHelper(ProblemSpecP& ps) const {
  ps->appendElement("bottom", d_bottom);
  ps->appendElement("top", d_top);
  ps->appendElement("radius", d_radius);
  ps->appendElement("num_radial", d_numRadial);
  ps->appendElement("num_axial", d_numAxial);
  ps->appendElement("thickness", d_thickness);
  ps->appendElement("endcap_thickness", d_capThick);
  ps->appendElement("arc_start_angle", d_arcStart);
  ps->appendElement("arc_angle", d_angle);
  ps->appendElement("output_file", d_fileName);
}

GeometryPieceP
SmoothCylGeomPiece::clone() const {
  return std::make_shared<SmoothCylGeomPiece>(*this);
}

/////////////////////////////////////////////////////////////////////////////
/*! Find if a point is inside the cylinder or end caps */
/////////////////////////////////////////////////////////////////////////////
bool
SmoothCylGeomPiece::inside(const Point& p) const {
  // Translate and rotate point so that cylinder axis is along "z"
  Vector pTransformed = p - d_bottom;
  pTransformed        = d_rotation.Transpose() * pTransformed;
  double x            = pTransformed.x();
  double y            = pTransformed.y();
  double z            = pTransformed.z();

  // a) Check is the point is outside the solid composite cylinder
  if (z < d_capThick || z > d_height + d_capThick) {
    return false;
  }
  double r_sq = x * x + y * y;
  if (r_sq > d_radius * d_radius) {
    return false;
  }

  // b) Find if the point is inside the inner cylinder
  double innerRad = d_radius - d_thickness;
  if (r_sq < innerRad * innerRad) {
    return false;
  }

  return true;
}

/////////////////////////////////////////////////////////////////////////////
/*! Find the bounding box for the cylinder */
/////////////////////////////////////////////////////////////////////////////
Box
SmoothCylGeomPiece::getBoundingBox() const {
  // Find the vector along the axis of the cylinder
  Vector capAxis = d_axis * (d_capThick / d_height);

  Vector bot = d_bottom.asVector() - capAxis;
  Vector top = d_top.asVector() + capAxis;
  Point lo(bot.x() - d_radius, bot.y() - d_radius, bot.z() - d_radius);
  Point hi(top.x() + d_radius, top.y() + d_radius, top.z() + d_radius);

  return Box(lo, hi);
}

//////////////////////////////////////////////////////////////////////////
/* Create particles */
//////////////////////////////////////////////////////////////////////////
unsigned int
SmoothCylGeomPiece::createPoints() {
  if (!(d_dx > 0.0)) {
    proc0cout << "**WARNING** smooth_cyl: Using default values of"
              << " num_radial = 5 and num_axial = 5.  Not creating points based"
              << " on particles per cell input."
              << "\n";
  }

  double axislen = d_height;
  d_numAxial     = std::ceil(axislen / d_dx);
  d_numRadial    = std::ceil(d_radius / d_dx);

  int totCount = 0;
  if (d_capThick > 0.0) {
    int count = createEndCapPoints();
    totCount += count;
  }
  if (d_thickness < d_radius) {
    int count = createHollowCylPoints();
    totCount += count;
  } else {
    int count = createSolidCylPoints();
    totCount += count;
  }

  // Write the output if requested
  if (d_fileName != "none") {
    writePoints(d_fileName, "vol");
  }

  return totCount;
}

//////////////////////////////////////////////////////////////////////////
/*! Create the particles on a circle on the x-y plane and then
  rotate them to the correct position. First particle is located
  at the center. */
//////////////////////////////////////////////////////////////////////////
int
SmoothCylGeomPiece::createEndCapPoints() {
  proc0cout << "Creating particles for the End Caps"
            << "\n";

  // Initialize count of the number of particles
  int count = 0;

  // Calculate the radial and axial material point spacing
  double axisInc  = d_height / (double)d_numAxial;
  int numCapAxial = int(d_capThick / axisInc) - 1;
  double radInc   = d_radius / (double)d_numRadial;

  // Create particles for the bottom end cap
  double currZ = 0.5 * axisInc;
  for (int kk = 0; kk < numCapAxial; ++kk) {
    Vector currCenter = d_bottom.asVector() - d_axis * currZ;

    for (int ii = 0; ii < d_numRadial; ++ii) {
      double prevRadius = ii * radInc;
      double currRadius = prevRadius + 0.5 * radInc;
      double nextRadius = (ii + 1) * radInc;
      int numCircum     = (int)(d_angle * currRadius / radInc);
      double phiInc     = d_angle / (double)numCircum;
      double area =
          0.5 * phiInc * (nextRadius * nextRadius - prevRadius * prevRadius);
      for (int jj = 0; jj < numCircum; ++jj) {
        double phi    = d_arcStart + 0.5 * phiInc + (double)jj * phiInc;
        double cosphi = cos(phi);
        double sinphi = sin(phi);

        // Create points on xy plane
        double x = currRadius * cosphi;
        double y = currRadius * sinphi;
        double z = 0;

        // Rotate points to correct orientation and
        // Translate to correct position
        Vector pp(x, y, z);
        pp = d_rotation * pp + currCenter;

        // Create points for particle "size" (radial, circum, axial order)
        Vector pr((currRadius + 0.5 * radInc) * cosphi,
                  (currRadius + 0.5 * radInc) * sinphi,
                  0);
        Vector ptheta(currRadius * cos(phi + 0.5 * phiInc),
                      currRadius * sin(phi + 0.5 * phiInc),
                      0);
        Vector pz(currRadius * cosphi, currRadius * sinphi, 0.5 * axisInc);
        pr        = d_rotation * pr + currCenter;
        ptheta    = d_rotation * ptheta + currCenter;
        pz        = d_rotation * pz + currCenter;
        Vector r1 = (pr - pp) * 2.0;
        Vector r2 = (ptheta - pp) * 2.0;
        Vector r3 = (pz - pp) * 2.0;

        Matrix3 size;
        size(0, 0) = r1[0];
        size(1, 0) = r1[1];
        size(2, 0) = r1[2];
        size(0, 1) = r2[0];
        size(1, 1) = r2[1];
        size(2, 1) = r2[2];
        size(0, 2) = r3[0];
        size(1, 2) = r3[1];
        size(2, 2) = r3[2];

        Point p(pp);
        if (insideComputationalDomain(p)) {
          d_points.push_back(p);
          d_scalars.at("p.volume").push_back(axisInc * area);
          d_tensors.at("p.size").push_back(size);
          d_vectors.at("p,rvec1").push_back(r1);
          d_vectors.at("p.rvec2").push_back(r2);
          d_vectors.at("p.rvec3").push_back(r3);
          // std::cout << "Point["<<count<<"]="<<p<<"\n";
          count++;
        }
      }
    }
    currZ -= axisInc;
  }

  // Create particles for the top end cap
  currZ = 0.5 * axisInc;
  for (int kk = 0; kk < numCapAxial; ++kk) {
    Vector currCenter = d_top.asVector() + d_axis * currZ;

    for (int ii = 0; ii < d_numRadial; ++ii) {
      double prevRadius = ii * radInc;
      double currRadius = prevRadius + 0.5 * radInc;
      double nextRadius = (ii + 1) * radInc;
      int numCircum     = (int)(d_angle * currRadius / radInc);
      double phiInc     = d_angle / (double)numCircum;
      double area =
          0.5 * phiInc * (nextRadius * nextRadius - prevRadius * prevRadius);
      for (int jj = 0; jj < numCircum; ++jj) {
        double phi    = d_arcStart + 0.5 * phiInc + (double)jj * phiInc;
        double cosphi = cos(phi);
        double sinphi = sin(phi);

        // Create points on xy plane
        double x = currRadius * cosphi;
        double y = currRadius * sinphi;
        double z = 0;

        // Rotate points to correct orientation and
        // Translate to correct position
        Vector pp(x, y, z);
        pp = d_rotation * pp + currCenter;
        Point p(pp);

        // Create points for particle "size" (radial, circum, axial order)
        Vector pr((currRadius + 0.5 * radInc) * cosphi,
                  (currRadius + 0.5 * radInc) * sinphi,
                  0);
        Vector ptheta(currRadius * cos(phi + 0.5 * phiInc),
                      currRadius * sin(phi + 0.5 * phiInc),
                      0);
        Vector pz(currRadius * cosphi, currRadius * sinphi, 0.5 * axisInc);
        pr        = d_rotation * pr + currCenter;
        ptheta    = d_rotation * ptheta + currCenter;
        pz        = d_rotation * pz + currCenter;
        Vector r1 = (pr - pp) * 2.0;
        Vector r2 = (ptheta - pp) * 2.0;
        Vector r3 = (pz - pp) * 2.0;

        Matrix3 size;
        size(0, 0) = r1[0];
        size(1, 0) = r1[1];
        size(2, 0) = r1[2];
        size(0, 1) = r2[0];
        size(1, 1) = r2[1];
        size(2, 1) = r2[2];
        size(0, 2) = r3[0];
        size(1, 2) = r3[1];
        size(2, 2) = r3[2];

        if (insideComputationalDomain(p)) {
          d_points.push_back(p);
          d_scalars.at("p.volume").push_back(axisInc * area);
          d_tensors.at("p.size").push_back(size);
          d_vectors.at("p,rvec1").push_back(r1);
          d_vectors.at("p.rvec2").push_back(r2);
          d_vectors.at("p.rvec3").push_back(r3);
          // std::cout << "Point["<<count<<"]="<<p<<"\n";
          count++;
        }
      }
    }
    currZ += axisInc;
  }

  return count;
}

//////////////////////////////////////////////////////////////////////////
/*! Create the particles on a circle on the x-y plane and then
  rotate them to the correct position.
  First particle is located at the center. */
//////////////////////////////////////////////////////////////////////////
int
SmoothCylGeomPiece::createSolidCylPoints() {
  proc0cout << "Creating particles for the Solid Cylinder"
            << "\n";

  // Initialize count of the number of particles
  int count = 0;

  // Calculate the radial and axial material point spacing
  double axisInc = d_height / (double)d_numAxial;
  double radInc  = d_radius / (double)d_numRadial;

  // Create particles for the solid cylinder
  double currZ = 0.5 * axisInc;
  for (int kk = 0; kk < d_numAxial; ++kk) {
    Vector currCenter = d_bottom.asVector() + d_axis * currZ;

    for (int ii = 0; ii < d_numRadial; ++ii) {
      double prevRadius = ii * radInc;
      double currRadius = prevRadius + 0.5 * radInc;
      double nextRadius = (ii + 1) * radInc;
      int numCircum     = (int)(d_angle * currRadius / radInc);
      double phiInc     = d_angle / (double)numCircum;
      double area =
          0.5 * phiInc * (nextRadius * nextRadius - prevRadius * prevRadius);
      for (int jj = 0; jj < numCircum; ++jj) {
        double phi    = d_arcStart + 0.5 * phiInc + (double)jj * phiInc;
        double cosphi = cos(phi);
        double sinphi = sin(phi);

        // Create points on xy plane
        double x = currRadius * cosphi;
        double y = currRadius * sinphi;
        double z = 0;
        Vector pp(x, y, z);

        // Rotate points to correct orientation and
        // Translate to correct position
        pp = d_rotation * pp + currCenter;
        Point p(pp);

        // Create points for particle "size" (radial, circum, axial order)
        Vector pr((currRadius + 0.5 * radInc) * cosphi,
                  (currRadius + 0.5 * radInc) * sinphi,
                  0);
        Vector ptheta(currRadius * cos(phi + 0.5 * phiInc),
                      currRadius * sin(phi + 0.5 * phiInc),
                      0);
        Vector pz(currRadius * cosphi, currRadius * sinphi, 0.5 * axisInc);
        pr        = d_rotation * pr + currCenter;
        ptheta    = d_rotation * ptheta + currCenter;
        pz        = d_rotation * pz + currCenter;
        Vector r1 = (pr - pp) * 2.0;
        Vector r2 = (ptheta - pp) * 2.0;
        Vector r3 = (pz - pp) * 2.0;

        Matrix3 size;
        size(0, 0) = r1[0];
        size(1, 0) = r1[1];
        size(2, 0) = r1[2];
        size(0, 1) = r2[0];
        size(1, 1) = r2[1];
        size(2, 1) = r2[2];
        size(0, 2) = r3[0];
        size(1, 2) = r3[1];
        size(2, 2) = r3[2];

        if (insideComputationalDomain(p)) {
          d_points.push_back(p);
          d_scalars.at("p.volume").push_back(axisInc * area);
          d_tensors.at("p.size").push_back(size);
          d_vectors.at("p,rvec1").push_back(r1);
          d_vectors.at("p.rvec2").push_back(r2);
          d_vectors.at("p.rvec3").push_back(r3);
          // std::cout << "Point["<<count<<"]="<<p<<"\n";
          count++;
        }
      }
    }
    currZ += axisInc;
  }

  return count;
}

//////////////////////////////////////////////////////////////////////////
/*! Create the particles on a circle on the x-y plane and then
  rotate them to the correct position */
//////////////////////////////////////////////////////////////////////////
int
SmoothCylGeomPiece::createHollowCylPoints() {
  proc0cout << "Creating particles for the Hollow Cylinder"
            << "\n";

  // Initialize count of the number of particles
  int count = 0;

  // Calculate the radial and axial material point spacing
  double axisInc  = d_height / (double)d_numAxial;
  double radInc   = d_radius / (double)d_numRadial;
  int numThick    = (int)(d_thickness / radInc);
  double innerRad = d_radius - d_thickness;

  // Create particles for the hollow cylinder
  double currZ = 0.5 * axisInc;
  for (int kk = 0; kk < d_numAxial; ++kk) {
    Vector currCenter = d_bottom.asVector() + d_axis * currZ;
    for (int ii = 0; ii < numThick; ++ii) {
      double prevRadius = innerRad + ii * radInc;
      double currRadius = prevRadius + radInc * 0.5;
      double nextRadius = innerRad + (ii + 1) * radInc;
      int numCircum     = (int)(d_angle * currRadius / radInc);
      double phiInc     = d_angle / (double)numCircum;
      double area =
          0.5 * phiInc * (nextRadius * nextRadius - prevRadius * prevRadius);
      for (int jj = 0; jj < numCircum; ++jj) {
        double phi    = d_arcStart + 0.5 * phiInc + (double)jj * phiInc;
        double cosphi = cos(phi);
        double sinphi = sin(phi);

        // Create points on xy plane
        double x = currRadius * cosphi;
        double y = currRadius * sinphi;
        double z = 0;

        // Rotate points to correct orientation and
        // Translate to correct position
        Vector pp(x, y, z);
        pp = d_rotation * pp + currCenter;
        Point p(pp);

        // Create points for particle "size" (radial, circum, axial order)
        Vector pr((currRadius + 0.5 * radInc) * cosphi,
                  (currRadius + 0.5 * radInc) * sinphi,
                  0);
        Vector ptheta(currRadius * cos(phi + 0.5 * phiInc),
                      currRadius * sin(phi + 0.5 * phiInc),
                      0);
        Vector pz(currRadius * cosphi, currRadius * sinphi, 0.5 * axisInc);
        pr        = d_rotation * pr + currCenter;
        ptheta    = d_rotation * ptheta + currCenter;
        pz        = d_rotation * pz + currCenter;
        Vector r1 = (pr - pp) * 2.0;
        Vector r2 = (ptheta - pp) * 2.0;
        Vector r3 = (pz - pp) * 2.0;

        Matrix3 size;
        size(0, 0) = r1[0];
        size(1, 0) = r1[1];
        size(2, 0) = r1[2];
        size(0, 1) = r2[0];
        size(1, 1) = r2[1];
        size(2, 1) = r2[2];
        size(0, 2) = r3[0];
        size(1, 2) = r3[1];
        size(2, 2) = r3[2];

        if (insideComputationalDomain(p)) {
          d_points.push_back(p);
          d_scalars.at("p.volume").push_back(axisInc * area);
          d_tensors.at("p.size").push_back(size);
          d_vectors.at("p,rvec1").push_back(r1);
          d_vectors.at("p.rvec2").push_back(r2);
          d_vectors.at("p.rvec3").push_back(r3);
          // std::cout << "Point["<<count<<"]="<<p<<"\n";
          count++;
        }
      }
    }
    currZ += axisInc;
  }

  return count;
}

/*! Test whether the created point is inside the computational domain
    or not */
bool
SmoothCylGeomPiece::insideComputationalDomain(const Point& pt) {
  if (!(pt.x() > d_domainMin.x())) return false;
  if (!(pt.y() > d_domainMin.y())) return false;
  if (!(pt.z() > d_domainMin.z())) return false;
  if (!(pt.x() < d_domainMax.x())) return false;
  if (!(pt.y() < d_domainMax.y())) return false;
  if (!(pt.z() < d_domainMax.z())) return false;

  return true;
}

} // end namespace Uintah