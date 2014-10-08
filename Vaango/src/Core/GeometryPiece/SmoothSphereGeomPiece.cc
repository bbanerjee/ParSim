/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <Core/GeometryPiece/SmoothSphereGeomPiece.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Box.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Grid/Patch.h>
#include <Core/Math/Matrix3.h>
#include <Core/Malloc/Allocator.h>
#include <iostream>

#include <boost/math/special_functions/ellint_2.hpp>
#include <cmath>

using namespace std;
using namespace Uintah;
using namespace SCIRun;

const string SmoothSphereGeomPiece::TYPE_NAME = "smooth_sphere";

//////////
// Constructor : Initialize stuff
SmoothSphereGeomPiece::SmoothSphereGeomPiece(ProblemSpecP& ps)
{
  ps->require("center", d_center);

  ps->require("outer_radius", d_outerRadius);
  if (d_outerRadius <= 0.0)
    SCI_THROW(ProblemSetupException("SmoothSphereGeom: Radius <= 0", __FILE__, __LINE__));

  d_innerRadius = 0.0;
  ps->get("inner_radius", d_innerRadius);
  if (d_innerRadius > d_outerRadius)
    SCI_THROW(ProblemSetupException("SmoothSphereGeom: inner radius > outer radius", __FILE__, __LINE__));

  ps->require("num_radial_pts", d_numRadial);
  if (d_numRadial < 1)
    SCI_THROW(ProblemSetupException("SmoothSphereGeom: Radial Divs < 1", __FILE__, __LINE__));

  d_algorithm = "spiral";
  ps->get("algorithm", d_algorithm);

  d_fileName = "none";
  ps->get("output_file", d_fileName);

}

//////////
// Destructor
SmoothSphereGeomPiece::~SmoothSphereGeomPiece()
{
}

void
SmoothSphereGeomPiece::outputHelper( ProblemSpecP & ps ) const
{
  ps->appendElement("center", d_center);
  ps->appendElement("outer_radius", d_outerRadius);
  ps->appendElement("inner_radius", d_innerRadius);
  ps->appendElement("num_radial_pts", d_numRadial);
  ps->appendElement("algorithm", d_algorithm);
  ps->appendElement("output_file", d_fileName);
}

GeometryPieceP
SmoothSphereGeomPiece::clone() const
{
  return scinew SmoothSphereGeomPiece(*this);
}

/////////////////////////////////////////////////////////////////////////////
/*! Find if a point is inside the _sphere */
/////////////////////////////////////////////////////////////////////////////
bool 
SmoothSphereGeomPiece::inside(const Point& pt) const
{
  // Find the distance of the point to the origin
  Vector vec = pt - d_center;
  double dist = vec.length();

  if (dist > d_innerRadius && dist < d_outerRadius) return true;

  return false;
}

/////////////////////////////////////////////////////////////////////////////
/*! Find the bounding box for the _sphere */
/////////////////////////////////////////////////////////////////////////////
Box 
SmoothSphereGeomPiece::getBoundingBox() const
{
  // Find the vector along the axis of the _sphere
  Point lo(d_center.x() - d_outerRadius, d_center.y() - d_outerRadius,
           d_center.z() - d_outerRadius);
  Point hi(d_center.x() + d_outerRadius, d_center.y() + d_outerRadius,
           d_center.z() + d_outerRadius);

  return Box(lo,hi);
}

//////////////////////////////////////////////////////////////////////////
/* Create particles */
//////////////////////////////////////////////////////////////////////////
unsigned int 
SmoothSphereGeomPiece::createPoints()
{
  int totCount = 0;
  if (d_algorithm == "spiral") {
    totCount = createSpherePointsSpiral();
  } else if (d_algorithm == "equal_area") {
    totCount = createSpherePointsEqualArea();
  } else {
    SCI_THROW(InvalidValue("SmoothSphereGeom: Unknown algorithm for point generation.", 
              __FILE__, __LINE__));
  }

  // Write the output if requested
  if (d_fileName != "none") {
    writePoints(d_fileName, "pts");
    writePoints(d_fileName, "vol");
  }

  return totCount;
}

//////////////////////////////////////////////////////////////////////////
/*! Create the particles using a spiral assignement with each point assigned an
    equal volume */
//////////////////////////////////////////////////////////////////////////
int 
SmoothSphereGeomPiece::createSpherePointsSpiral()
{
  cout << "Creating particles for the Solid Sphere" << endl;

  // Find the characteristic distance between points
  double thickness = d_outerRadius - d_innerRadius; 
  double char_dist = thickness/d_numRadial;
 
  // Create points for each shell
  for (int ii = 0; ii < d_numRadial; ++ii) {
    double shell_inner_radius = d_innerRadius + ii*char_dist; 
    double shell_outer_radius = d_innerRadius + (ii+1)*char_dist; 
    createPointSetSpiral(shell_outer_radius, shell_inner_radius, char_dist);
  }
  
  return d_points.size();
}

//////////////////////////////////////////////////////////////////////////
/*! Create the spiral point set for ecah shell */
//////////////////////////////////////////////////////////////////////////
void
SmoothSphereGeomPiece::createPointSetSpiral(double outer_radius,
                                            double inner_radius,
                                            double char_dist)
{
  // Find the radius of the middle of the annulus
  double mid_radius = 0.5*(outer_radius+inner_radius);

  // This is the point at the center of the sphere that is being discretized
  if (inner_radius < 1.0e-16) {
    d_points.push_back(Point(0.0, 0.0, 0.0));
    d_volume.push_back(volumeOfSphere(mid_radius));
    return;
  }

  // Find the length of a spiral that covers the surface of the sphere with
  // equally spaced points separated from each other by the characteristic distance
  double phi_max = 3.0*M_PI*M_PI*mid_radius/(2.0*char_dist);
  double mm = -(phi_max*phi_max)/(M_PI*M_PI);
  double len_max = elliptic2E(mm);
  len_max *= (2.0*mid_radius);
  int num_points = std::ceil(len_max/char_dist);
  
  // Compute volume of sphere/annulus
  double outer_volume = volumeOfSphere(outer_radius);
  double inner_volume = volumeOfSphere(inner_radius);
  double point_volume = (outer_volume - inner_volume)/(double) num_points;

  // Loop thru set of points
  double theta = 0.0;
  for (int ii = 0; ii < num_points; ++ii) {
    double cosphi = ((double)(ii + 1 - num_points) + (double)(ii))/(double)(num_points-1);
    double sinphi = std::sqrt(1.0 - cosphi*cosphi);
    if (ii == 0 || ii == num_points -1) {
      theta = 0.0;
    } else {
      theta += 3.6/(sinphi*std::sqrt((double) num_points));
      theta = std::fmod(theta, 2.0*M_PI);
    }
    double x = mid_radius*sinphi*cos(theta) + d_center.x();
    double y = mid_radius*sinphi*sin(theta) + d_center.y();
    double z = mid_radius*cosphi + d_center.z();
    d_points.push_back(Point(x, y, z));
    d_volume.push_back(point_volume);
  }
}

//////////////////////////////////////////////////////////////////////////
/*! Create the particles using the equal area algorithm of Paul Leopardi
    Paul Leopardi, "A partition of the unit sphere into regions of equal area and small 
    diameter", Electronic Transactions on Numerical Analysis, Volume 25, 2006, pp. 309-327.
    MR 2280380, Preprint: UNSW Applied Mathematics Report AMR05/18, May 2005, revised June 2006.

    The algorithm is the one used in the EQSP software package, which partitions a finite 
    dimensional unit sphere into regions of equal area and small diameter. */
//////////////////////////////////////////////////////////////////////////
int 
SmoothSphereGeomPiece::createSpherePointsEqualArea()
{
  cout << "Creating particles for the Solid Sphere" << endl;

  // Find the characteristic distance between points
  double thickness = d_outerRadius - d_innerRadius; 
  double char_dist = thickness/d_numRadial;
 
  // Create points for each shell
  for (int ii = 0; ii < d_numRadial; ++ii) {
    double shell_inner_radius = d_innerRadius + ii*char_dist; 
    double shell_outer_radius = d_innerRadius + (ii+1)*char_dist; 
    createPointSetPolar2D(shell_outer_radius, shell_inner_radius, char_dist);
  }
  
  return d_points.size();
}

// Create the point set on a unit 2-sphere with origin at (0.0, 0.0, 0.0) 
// using Leopardi's recursive algorithm
void
SmoothSphereGeomPiece::createPointSetPolar2D(double outer_radius,
                                             double inner_radius,
                                             double char_dist)
{
  // Find the radius of the middle of the annulus
  double mid_radius = 0.5*(outer_radius+inner_radius);

  // This is the point at the center of the sphere that is being discretized
  if (mid_radius < char_dist) {
    d_points.push_back(Point(0.0, 0.0, 0.0));
    d_volume.push_back(volumeOfSphere(mid_radius));
    return;
  }

  // Find the length of a spiral that covers the surface of the sphere with
  // equally spaced points separated from each other by the characteristic distance
  double phi_max = 3.0*M_PI*M_PI*mid_radius/(2.0*char_dist);
  double mm = -(phi_max*phi_max)/(M_PI*M_PI);
  double len_max = elliptic2E(mm);
  len_max *= (2.0*mid_radius);
  int num_points = std::ceil(len_max/char_dist);

  //std::cout << "radius = " << mid_radius << " char_dist = " << char_dist 
  //          << " num_points = " << num_points << std::endl;
  
  // Compute volume of sphere/annulus
  double outer_volume = volumeOfSphere(outer_radius);
  double inner_volume = volumeOfSphere(inner_radius);
  double point_volume = (outer_volume - inner_volume)/(double) num_points;

  // Partition the sphere into caps 
  std::vector<double> cap_colatitudes;
  std::vector<int> int_regions;
  createCaps2D(num_points, cap_colatitudes, int_regions);

  // Create a vector for storing polar angles
  std::vector<double> points_polar_theta;
  std::vector<double> points_polar_phi;

  // Add center point of top polar cap
  points_polar_theta.push_back(0.0);
  points_polar_phi.push_back(0.0);
  d_volume.push_back(point_volume);

  // Loop thru partitions
  int num_collars = int_regions.size() - 2;
  double offset = 0.0;
  for (int collar_id = 0; collar_id < num_collars; ++collar_id) {
    double colatitude_start = cap_colatitudes[collar_id];
    double colatitude_end = cap_colatitudes[collar_id+1];

    int num_regions_in_collar = int_regions[collar_id+1];
    
    //std::cout << "a_top = " << colatitude_start
    //          << " a_bot = " << colatitude_end
    //          << " n_in_collar = " << num_regions_in_collar << std::endl;

    // Create lower dimensional point set
    std::vector<double> points_1D;
    createPointSetPolar1D(num_regions_in_collar, points_1D);

    //std::cout << "points_1 = " ;
    //for (auto iter = points_1D.begin(); iter != points_1D.end(); iter++) {
    //   std::cout << *iter << " ";
    //}
    //std::cout << std::endl;

    // Determine the center points of the lower-dimensional sphere
    double colatitude_mid = 0.5*(colatitude_start + colatitude_end);
    for (auto iter = points_1D.begin(); iter != points_1D.end(); ++iter) {
      double point_angle = *iter + 2.0*M_PI*offset;
      point_angle = std::fmod(point_angle, 2.0*M_PI);
      points_polar_theta.push_back(point_angle);
      points_polar_phi.push_back(colatitude_mid);
      d_volume.push_back(point_volume);

      //std::cout << "  point_angle = " << point_angle << " offset = " << offset << std::endl;
    }

    //std::cout << "n_in_collar = " << num_regions_in_collar 
    //          << " 2+collar_n = " << collar_id+2 << std::endl;
    //std::cout << "n_regions = " ;
    //for (auto iter = int_regions.begin(); iter != int_regions.end(); iter++) {
    //   std::cout << *iter << " ";
    //}
    //std::cout << std::endl;

    // Compute offset
    double circle_offset = circleOffset(num_regions_in_collar, int_regions[collar_id+2]);
    offset += circle_offset;
    offset -= std::floor(offset);

    //std::cout << "points_s(1,:) = " ;
    //for (auto iter = points_polar_theta.begin(); iter != points_polar_theta.end(); iter++) {
    //   std::cout << *iter << " ";
    //}
    //std::cout << std::endl;
    

  }

  // Add center point of bottom polar cap
  points_polar_theta.push_back(0.0);
  points_polar_phi.push_back(M_PI);
  d_volume.push_back(point_volume);

  // Convert polar angles to cartesian coordinates
  for (int ii = 0; ii <  (int) points_polar_theta.size(); ++ii) {
    double theta = points_polar_theta[ii];
    double phi = points_polar_phi[ii];
    double x = mid_radius*sin(phi)*cos(theta) + d_center.x();
    double y = mid_radius*sin(phi)*sin(theta) + d_center.y();
    double z = mid_radius*cos(phi) + d_center.z();
    d_points.push_back(Point(x, y, z));
  }
}

// Create point sets for the collars (1D)
void
SmoothSphereGeomPiece::createPointSetPolar1D(int num_points,
                                             std::vector<double>& points)
{
  if (num_points == 1) {
    points.push_back(0.0);
    return;
  }

  // Partition the 1D-sphere into caps 
  std::vector<double> cap_colatitudes;
  createCaps1D(num_points, cap_colatitudes);

  // We have a circle and cap_colatitudes is an increasing list of angles of sectors,
  // with capColatitude(k) being the cumulative arc length 2*pi/k.  The points are placed half 
  // way along each sector.
  for (auto iter = cap_colatitudes.begin(); iter != cap_colatitudes.end(); iter++) {
    points.push_back(*iter - M_PI/(double) num_points);
  }
}

// Create the set of nested caps in 2D
void
SmoothSphereGeomPiece::createCaps2D(int num_points,
                                    std::vector<double>& cap_colatitudes,
                                    std::vector<int>& int_regions)
{
  // Find the colatitude of the North polar spherical cap.
  double polar_colatitude = M_PI/2.0;
  if (num_points > 2) {
    polar_colatitude = polarColatitude(num_points);
  }
  
  // Determine the ideal angle for spherical collars.
  double ideal_collar_angle = idealCollarAngle(num_points);

  // Find the number of collars between the polar caps.
  int num_collars = numCollars(num_points, polar_colatitude, ideal_collar_angle);

  // Create a list of the ideal real number of regions in each collar,
  // plus the polar caps.
  std::vector<double> real_regions;
  idealRegionList(num_points, polar_colatitude, num_collars, real_regions);
   
  // Create  a list of the natural number of regions in each collar and
  // the polar caps.
  roundToNaturals(num_points, real_regions, int_regions);

  // Compute s_cap, an increasing list of colatitudes of spherical caps 
  // which enclose the same area  as that given by the cumulative sum of regions.
  capColatitudes(num_points, polar_colatitude, int_regions, cap_colatitudes);
}

// Create the set of nested caps in 1D
void
SmoothSphereGeomPiece::createCaps1D(int num_points,
                                    std::vector<double>& cap_colatitudes)
{
  for (int ii = 0; ii < num_points; ++ii) {
    double sector = (double) (ii+1);
    double s_cap = (sector*2.0*M_PI)/(double) num_points;
    //cap_colatitudes.push_back((double)(ii+1)*(2.0*M_PI/(double) num_points));
    cap_colatitudes.push_back(s_cap);
  }
}

// Maximize the minimum distance of center points of collars
double
SmoothSphereGeomPiece::circleOffset(int num_start, int num_end)
{
  double offset = (1.0/(double) num_end - 1.0/(double) num_start)*0.5
                  + (double) gcd(num_start, num_end)/(double) (2.0*num_start*num_end);
  //std::cout << "n_bot = " << num_end << " n_top = " << num_start
  //          << " gcd = " << gcd(num_start, num_end) << " offset = " << offset  << std::endl;
  return offset;
}

// Compute the colatitide of the north polar spherical cap
double
SmoothSphereGeomPiece::polarColatitude(int num_points)
{
  double area = areaOfIdealRegion(num_points);
  return sRadiusOfCap(area);
}

// Compute the numbers of collars between the polar caps
int 
SmoothSphereGeomPiece::numCollars(int num_points,
                                  double polar_colatitude,
                                  double ideal_collar_angle)
{
  return std::max(1, (int) std::round((M_PI - 2.0*polar_colatitude)/ideal_collar_angle));
}

// Compute ideal real number of regions in each collar
void
SmoothSphereGeomPiece::idealRegionList(int num_points,
                                       double polar_colatitude,
                                       int num_collars,
                                       std::vector<double>& real_regions)
{
  real_regions.push_back(1);
  if (num_collars > 0) {
    double angle_fitting = (M_PI - 2.0*polar_colatitude)/(double) num_collars;
    double ideal_region_area = areaOfIdealRegion(num_points);
    for (int collar_id = 0; collar_id < num_collars; ++collar_id) {
      double collar_end_radius = polar_colatitude + (collar_id+1.0)*angle_fitting;
      double collar_start_radius = collar_end_radius - angle_fitting;
      double ideal_collar_area = areaOfCollar(collar_start_radius, collar_end_radius);
      real_regions.push_back(ideal_collar_area/ideal_region_area);
    }
  }
  real_regions.push_back(1);
}

// Round to nearest int
void
SmoothSphereGeomPiece::roundToNaturals(int num_points,
                                       std::vector<double>& real_regions,
                                       std::vector<int>& int_regions)
{
  double discrepancy = 0.0;
  for (int zone_n = 0; zone_n < (int) real_regions.size(); ++zone_n) {
    int_regions.push_back(std::round(real_regions[zone_n] + discrepancy));
    discrepancy += real_regions[zone_n] - int_regions[zone_n];
  }
}

// Compute colatitudes of spherical caps enclosing cumulative sum of regions
void
SmoothSphereGeomPiece::capColatitudes(int num_points, 
                                      double polar_colatitude,
                                      std::vector<int>& int_regions,
                                      std::vector<double>& cap_colatitudes)
{
  cap_colatitudes.push_back(polar_colatitude);
  double ideal_region_area = areaOfIdealRegion(num_points);
  int num_collars = int_regions.size() - 2;
  int subtotal_int_regions = 1;
  for (int collar_id = 0; collar_id < num_collars; ++collar_id) {
    subtotal_int_regions += int_regions[collar_id+1];
    double sph_radius_of_cap = sRadiusOfCap(subtotal_int_regions*ideal_region_area);
    cap_colatitudes.push_back(sph_radius_of_cap);
  }
  cap_colatitudes.push_back(M_PI);
}

// Compute the spherical radius of a cap
double
SmoothSphereGeomPiece::sRadiusOfCap(double area)
{
  return 2.0*std::asin(std::sqrt(area/M_PI)*0.5);
}

// Compute angle for spherecal collars for an equal area partition
double
SmoothSphereGeomPiece::idealCollarAngle(int num_points)
{
  double area = areaOfIdealRegion(num_points);
  return std::sqrt(area); 
}

// Compute area of one region of a equal area partition
double
SmoothSphereGeomPiece::areaOfIdealRegion(int num_points)
{
  double area = surfaceAreaOfSphere(1.0);
  area /= (double) num_points;
  return area;
}

// Compute the area of a spherical collar
double
SmoothSphereGeomPiece::areaOfCollar(double collar_start_radius, double collar_end_radius)
{
  double cap_area_start = areaOfCap(collar_start_radius);
  double cap_area_end = areaOfCap(collar_end_radius);
  return cap_area_end - cap_area_start;
}

// Compute the area of a spherical cap
double 
SmoothSphereGeomPiece::areaOfCap(double collar_radius)
{
  double ss = std::sin(0.5*collar_radius);
  return 4.0*M_PI*ss*ss;
}

// Compute the surface area of a sphere 
double
SmoothSphereGeomPiece::surfaceAreaOfSphere(double radius)
{
  return (4.0*M_PI*radius*radius); 
}

// Compute the volume of a sphere
double
SmoothSphereGeomPiece::volumeOfSphere(double radius)
{
  return (4.0/3.0)*M_PI*(radius*radius*radius);
}

// Compute greatest common divisor of two ints
int 
SmoothSphereGeomPiece::gcd(int u, int v)
{
  int shift;
 
  /* GCD(0,v) == v; GCD(u,0) == u, GCD(0,0) == 0 */
  if (u == 0) return v;
  if (v == 0) return u;
 
  /* Let shift := lg K, where K is the greatest power of 2
        dividing both u and v. */
  for (shift = 0; ((u | v) & 1) == 0; ++shift) {
         u >>= 1;
         v >>= 1;
  }
 
  while ((u & 1) == 0)
    u >>= 1;
 
  /* From here on, u is always odd. */
  do {
       /* remove all factors of 2 in v -- they are not common */
       /*   note: v is not zero, so while will terminate */
       while ((v & 1) == 0)  /* Loop X */
           v >>= 1;
 
       /* Now u and v are both odd. Swap if necessary so u <= v,
          then set v = v - u (which is even). For bignums, the
          swapping is just pointer movement, and the subtraction
          can be done in-place. */
       if (u > v) {
         unsigned int t = v; v = u; u = t;}  // Swap u and v.
       v = v - u;                       // Here v >= u.
     } while (v != 0);
 
  /* restore common factors of 2 */
  return u << shift;
}

// Compute complete elliptic integral of second kind
// In this case mm < 1 so imaginary modulus transformation is needed
double
SmoothSphereGeomPiece::elliptic2E(double mm)
{
  double m_inv = -mm/(1.0-mm);
  double kk = std::sqrt(m_inv);
  double EE = boost::math::ellint_2(kk);
  return std::sqrt(1.0-mm)*EE;
}
