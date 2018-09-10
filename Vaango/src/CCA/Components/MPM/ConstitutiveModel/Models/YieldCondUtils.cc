#include <CCA/Components/MPM/ConstitutiveModel/Models/YieldCondUtils.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/TensorUtils.h>
#include <Core/Geometry/Vector.h>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iterator>

#include <Eigen/Dense>

namespace Vaango {

namespace Util {

/* Compute outward normal at each point of a polyline */
std::vector<Uintah::Vector>
computeNormals(const std::vector<Uintah::Point>& polyline)
{
  const double tminus = 0.99; 
  const double tplus = 0.01; 
  std::vector<Uintah::Vector> tangents;
  for (auto ii = 1u; ii < polyline.size()-1; ii++) {
    Uintah::Vector minus = polyline[ii-1]*(1 - tminus) + polyline[ii]*tminus;
    Uintah::Vector plus = polyline[ii]*(1 - tplus) + polyline[ii+1]*tplus;
    Uintah::Vector tangent = (plus - minus)/2.0;

    tangent.normalize();
    tangents.push_back(tangent);
  }
  //std::cout << "Tangents:";
  //std::copy(tangents.begin(), tangents.end(),
  //          std::ostream_iterator<Uintah::Vector>(std::cout, " "));
  //std::cout << std::endl;

  std::vector<Uintah::Vector> normals;
  for (auto ii = 1u; ii < tangents.size()-1; ii++) {
    Uintah::Vector minus = tangents[ii-1]*(1 - tminus) + tangents[ii]*tminus;
    Uintah::Vector plus = tangents[ii]*(1 - tplus) + tangents[ii+1]*tplus;
    Uintah::Vector normal = -(plus - minus)/2.0;
    normal.normalize();
    normals.push_back(normal);
  }

  return normals;
}

/* Checks whethere three points are in counter-clockwise order */
double
checkOrder(const Uintah::Point& p1, const Uintah::Point& p2,
           const Uintah::Point& p3)
{
  return (p2.x() - p1.x())*(p3.y() - p1.y()) -
         (p2.y() - p1.y())*(p3.x() - p1.x());
}

/* Compute the convex hull of a set of points in 2D */
std::vector<Uintah::Point> 
convexHull2D(const std::vector<Uintah::Point>& points)
{
  std::vector<Uintah::Point> lower;
  if (points.size() < 3) {
    for (const auto& point : points) {
      lower.push_back(point);
    }
    return lower;
  }

  for (const auto& point : points) {
    auto k = lower.size();
    while (k > 1 && 
           Vaango::Util::checkOrder(lower[k-2], lower[k-1], point) <= 0) {
      lower.pop_back();
      k = lower.size();
    }
    lower.push_back(point);
  }
  lower.pop_back();

  std::vector<Uintah::Point> upper;
  for (const auto& point : Vaango::Util::reverse(points)) {
    auto k = upper.size();
    while (k > 1 && 
           Vaango::Util::checkOrder(upper[k-2], upper[k-1], point) <= 0) {
      upper.pop_back();
      k = upper.size();
    }
    upper.push_back(point);
  }
  upper.pop_back();

  lower.insert(lower.end(), 
               std::make_move_iterator(upper.begin()), 
               std::make_move_iterator(upper.end()));
  return lower;
}

/* Get the closest point on the yield surface */
double
findClosestPoint(const Uintah::Point& p, const std::vector<Uintah::Point>& poly,
                 Uintah::Point& min_p)
{
  double TOLERANCE_MIN = 1.0e-12;
  std::vector<Uintah::Point> XP;

  // Loop through the segments of the polyline
  auto iterStart = poly.begin();
  auto iterEnd = poly.end();
  auto iterNext = iterStart;
  ++iterNext;
  for (; iterNext != iterEnd; ++iterStart, ++iterNext) {
    Uintah::Point start = *iterStart;
    Uintah::Point next = *iterNext;

    // Find shortest distance from point to the polyline line
    Uintah::Vector m = next - start;
    Uintah::Vector n = p - start;
    if (m.length2() < TOLERANCE_MIN) {
      XP.push_back(start);
    } else {
      const double t0 = Dot(m, n) / Dot(m, m);
      if (t0 <= 0.0) {
        XP.push_back(start);
      } else if (t0 >= 1.0) {
        XP.push_back(next);
      } else {
        // Shortest distance is inside segment; this is the closest point
        min_p = m * t0 + start;
        XP.push_back(min_p);
        // std::cout << "Closest: " << min_p << std::endl;
        // return;
      }
    }
  }

  double min_d = std::numeric_limits<double>::max();
  for (const auto& xp : XP) {
    // Compute distance sq
    double dSq = (p - xp).length2();
    if (dSq < min_d) {
      min_d = dSq;
      min_p = xp;
    }
  }

  // std::cout << "Closest: " << min_p << std::endl
  //          << "At: " << min_d << std::endl;
  return min_d;
}

/* Find two yield surface segments that are closest to an input point */
std::size_t
getClosestSegments(const Uintah::Point& pt,
                   const std::vector<Uintah::Point>& poly,
                   std::vector<Uintah::Point>& segments)
{
  if (poly.size() < 4) {
    segments = poly;
    return 0;
  }

  segments.clear();
  std::size_t min_index = 0;
  double min_dSq = std::numeric_limits<double>::max();
  for (auto index = 0u; index < poly.size()-1; index++) {
    std::vector<Uintah::Point> segment = { poly[index], poly[index+1] };
    auto segmentLength = (segment[1] - segment[0]).length();
    Uintah::Point closest(0, 0, 0);
    Vaango::Util::findClosestPoint(pt, segment, closest);

    double dSq = (pt - closest).length2();
    if (dSq < min_dSq) {
      min_dSq = dSq;
      if ((closest - segment[0]).length() < 0.5*segmentLength) {
        min_index = index;
      } else {
        min_index = index+1;
      }
    }
  } 

  auto close_index = min_index;
  if (min_index == 0) {
    segments.push_back(poly[min_index]);
    segments.push_back(poly[min_index+1]);
    segments.push_back(poly[min_index+2]);
    close_index = min_index+1;
  } else if (min_index == poly.size()-1) {
    segments.push_back(poly[min_index-2]);
    segments.push_back(poly[min_index-1]);
    segments.push_back(poly[min_index]);
  } else {
    segments.push_back(poly[min_index-1]);
    segments.push_back(poly[min_index]);
    segments.push_back(poly[min_index+1]);
  }

  return close_index;
}

/**
 *  Find three sequential points (two sequential yield surface segments) 
 *  that are closest to an input point (using a binary search)
 *
 *  Modifies:
 *   vector<Point> : 
 *      Point  : first_point
 *      Point  : second_point
 *      Point  : third_point
 *
 *  Returns:
 *   size_t : index of the closest point in the sequence
 */
std::size_t
getClosestSegmentsBinarySearch(const Uintah::Point& pt,
                               const std::vector<Uintah::Point>& polyline,
                               std::vector<Uintah::Point>& segments)
{
  if (polyline.size() < 4) {
    segments = polyline;
    return 0;
  }

  Uintah::Point closest_pt;
  std::size_t closest_index;
  double closest_dist;

  std::tie(closest_pt, closest_index, closest_dist) = 
    Vaango::Util::closestPointBinarySearch(pt, polyline.begin(), polyline.end(),
                                           polyline);

  if (closest_index == 0) {
    segments.push_back(polyline[closest_index]);
    segments.push_back(polyline[closest_index+1]);
    segments.push_back(polyline[closest_index+2]);
  } else if (closest_index == polyline.size()-1) {
    segments.push_back(polyline[closest_index-2]);
    segments.push_back(polyline[closest_index-1]);
    segments.push_back(polyline[closest_index]);
  } else {
    segments.push_back(polyline[closest_index-1]);
    segments.push_back(polyline[closest_index]);
    segments.push_back(polyline[closest_index+1]);
  }

  return closest_index;
}

std::tuple<Uintah::Point, std::size_t, double>
closestPointBinarySearch(const Uintah::Point P0,
                         std::vector<Uintah::Point>::const_iterator poly_begin,
                         std::vector<Uintah::Point>::const_iterator poly_end,
                         const std::vector<Uintah::Point>& poly_orig)
{
  auto num_pts = poly_end - poly_begin;
  auto mid_index = std::round(num_pts/2);

  auto P_left  = *poly_begin;
  auto P_right = *(poly_end-1);
  auto P_mid = *(poly_begin + mid_index);

  auto dist_left = (P_left - P0).length2();
  auto dist_right = (P_right - P0).length2();
  auto dist_mid = (P_mid - P0).length2();

  //std::cout << "P0 = " << P0 << " P_left = " << P_left << " P_right = " << P_right
  //          << " P_mid = " << P_mid << "\n";
  //std::cout << "d_left = " << dist_left << " d_right = " << dist_right
  //          << " d_mid = " << dist_mid << "\n";

  Uintah::Point P_near;
  std::size_t index_near;
  double dist_near = dist_mid;

  if (num_pts == 2) {
    if (dist_left < dist_right) {
      P_near = P_left;
      dist_near = dist_left;
      index_near = poly_begin - poly_orig.begin();
    } else {
      P_near = P_right;
      dist_near = dist_right;
      index_near = poly_end - poly_orig.begin() - 1;
    }
    //std::cout << "final: near = " << P_near << " index = " << index_near 
    //          << " dist = " << dist_near << "\n";
    return std::make_tuple(P_near, index_near, dist_near);
  }

  if (dist_left < dist_mid) {
    std::tie(P_near, index_near, dist_near) =
      closestPointBinarySearch(P0, poly_begin, poly_begin+mid_index+1, poly_orig);
    //std::cout << "left < right: near = " << P_near << " index = " << index_near 
    //          << " dist = " << dist_near << "\n";
  } else if (dist_right < dist_mid) {
    std::tie(P_near, index_near, dist_near) =
      closestPointBinarySearch(P0, poly_begin+mid_index, poly_end, poly_orig);
    //std::cout << "right < left: near = " << P_near << " index = " << index_near 
    //          << " dist = " << dist_near << "\n";
  } else {
    Uintah::Point P_near_left;
    std::size_t index_near_left;
    double dist_near_left;
    std::tie(P_near_left, index_near_left, dist_near_left) =
      closestPointBinarySearch(P0, poly_begin, poly_begin+mid_index+1, poly_orig);

    Uintah::Point P_near_right;
    std::size_t index_near_right;
    double dist_near_right;
    std::tie(P_near_right, index_near_right, dist_near_right) =
      closestPointBinarySearch(P0, poly_begin+mid_index, poly_end, poly_orig);

    if (dist_near_left < dist_near_right) {
      P_near = P_near_left;
      dist_near = dist_near_left;
      index_near = index_near_left;
      //std::cout << "mid right > left: near = " << P_near << " index = " << index_near 
      //        << " dist = " << dist_near << "\n";
    } else {
      P_near = P_near_right;
      dist_near = dist_near_right;
      index_near = index_near_right;
      //std::cout << "mid right <= left: near = " << P_near << " index = " << index_near 
      //        << " dist = " << dist_near << "\n";
    }
  }

  return std::make_tuple(P_near, index_near, dist_near);
}

/* linspace function */
std::vector<double>
linspace(double start, double end, int num)
{
  std::vector<double> linspaced;
  if (num > 0) {
    double delta = (end - start) / (double)num;
    for (int i = 0; i < num + 1; ++i) {
      linspaced.push_back(start + delta * (double)i);
    }
  }
  return linspaced;
}

/* linspace function */
void
linspace(const double& start, const double& end, const int& num,
         std::vector<double>& linspaced)
{
  linspaced = linspace(start, end, num);
  /*
  linspaced.clear();
  if (num > 0) {
    double delta = (end - start) / (double)num;
    for (int i = 0; i < num + 1; ++i) {
      linspaced.push_back(start + delta * (double)i);
    }
  }
  */
  return;
}

/* Create open quadratic uniform B-spline between approximating a polyline */
void
computeOpenUniformQuadraticBSpline(const std::vector<Uintah::Point>& polyline,
                                   size_t ptsPerSegment,
                                   std::vector<Uintah::Point>& spline)
{
  computeOpenUniformQuadraticBSpline(polyline, 0, polyline.size()-1,
                                     ptsPerSegment, spline);
}

/* Create open quadratic uniform B-spline approximating a segment of a polyline */
/* WARNING : Creates duplicate points at shared segement end points */
void
computeOpenUniformQuadraticBSpline(const std::vector<Uintah::Point>& polyline,
                                   size_t segmentStartIndex,
                                   size_t segmentEndIndex,
                                   size_t ptsPerSegment,
                                   std::vector<Uintah::Point>& spline)
{
  auto n = polyline.size() - 1;
  if (n < 2) {
    spline = polyline;
    return;
  }

  auto k = 2u;
  auto tvals = Vaango::Util::linspace(0.0, 1.0, ptsPerSegment);

  Uintah::Point splinePoint(0.0, 0.0, 0.0);
  //for (auto jj = 0u; jj < n - k + 1; jj++) {
  for (auto jj = segmentStartIndex; jj < segmentEndIndex - k + 1; jj++) {
    for (const auto t : tvals) {
      if (jj == 0) {
        splinePoint = 
          computeOpenUniformQuadraticBSpline(t, quadBSplineLo, polyline[jj],
                                             polyline[jj+1], polyline[jj+2]);
      } else if (jj == n-2) {
        splinePoint = 
          computeOpenUniformQuadraticBSpline(t, quadBSplineHi, polyline[jj],
                                             polyline[jj+1], polyline[jj+2]);
      } else {
        splinePoint = 
          computeOpenUniformQuadraticBSpline(t, quadBSpline, polyline[jj],
                                             polyline[jj+1], polyline[jj+2]);
      }
      spline.push_back(splinePoint);
    }
  }
}

/* Get a single point on open quadratic uniform B-spline between three points */
Uintah::Point
computeOpenUniformQuadraticBSpline(const double& t,
                                   const Uintah::Matrix3& splineMatrix,
                                   const Uintah::Point& point_k,
                                   const Uintah::Point& point_k1,
                                   const Uintah::Point& point_k2)
{
  Uintah::Vector A(1, t, t*t);
  Uintah::Vector Px(point_k.x(), point_k1.x(), point_k2.x());
  Uintah::Vector Py(point_k.y(), point_k1.y(), point_k2.y());
  double sx = Uintah::Dot(A, (splineMatrix*Px*0.5));
  double sy = Uintah::Dot(A, (splineMatrix*Py*0.5));
  return Uintah::Point(sx, sy, 0.0);
}

/* Find closest point to open quadratic uniform B-spline approximating a 
   segment of a polyline */
std::tuple<Uintah::Point, Uintah::Vector, Uintah::Vector>
computeClosestPointQuadraticBSpline(const Uintah::Point pt,
                                    const std::vector<Uintah::Point>& polyline,
                                    size_t segmentStartIndex,
                                    size_t segmentEndIndex)
{
  double t;
  Uintah::Point closest(0, 0, 0);
  Uintah::Vector T(0, 0, 0);
  Uintah::Vector N(0, 0, 0);

  auto n = polyline.size() - 1;
  if (n < 2) {
    if (n == 0) {
      closest = polyline[0];
    } else {
      findClosestPoint(pt, polyline, closest);
      T.x(polyline[1].x() - polyline[0].x());
      T.y(polyline[1].y() - polyline[0].y());
      N.x(0);
      N.y(0);
    }
    return std::make_tuple(closest, T, N);
  }

  auto k = 2u;
  for (auto jj = segmentStartIndex; jj < segmentEndIndex - k + 1; jj++) {
    if (jj == 0) {
      std::tie(t, closest, T, N) = 
        findClosestPointToQuadraticBSplineNewton(pt, quadBSplineLo, polyline[jj],
                                                 polyline[jj+1], polyline[jj+2]);
    } else if (jj == n-2) {
      std::tie(t, closest, T, N) = 
        findClosestPointToQuadraticBSplineNewton(pt, quadBSplineHi, polyline[jj],
                                                 polyline[jj+1], polyline[jj+2]);
    } else {
      std::tie(t, closest, T, N) = 
        findClosestPointToQuadraticBSplineNewton(pt, quadBSpline, polyline[jj],
                                                 polyline[jj+1], polyline[jj+2]);
    }
  }

  return std::make_tuple(closest, T, N);
}

/* Find closest point on a quadratic B-spline using Newton iterations */
std::tuple<double, Uintah::Point, Uintah::Vector, Uintah::Vector>
findClosestPointToQuadraticBSplineNewton(const Uintah::Point pt, 
                                         const Uintah::Matrix3& splineMatrix,
                                         const Uintah::Point& point_k,
                                         const Uintah::Point& point_k1,
                                         const Uintah::Point& point_k2)
{
  double t = 0.5;
  double f = 1.0;
  int n = 0;

  Uintah::Point B(0, 0, 0);
  Uintah::Vector T(0, 0, 0);
  Uintah::Vector N(0, 0, 0);

  while (std::abs(f) > 1.0e-10) {
    computeOpenUniformQuadraticBSpline(t, splineMatrix, point_k, point_k1,
                                       point_k2, B, T, N);
    f = Uintah::Dot(pt - B, T);
    double fp = -Uintah::Dot(T, T) + Uintah::Dot(pt - B, N);
    t -= f/fp;
    ++n;
    if (n > 20) break;
  }

  if (!isInBounds<double>(t, 0, 1)) {
    std::vector<Uintah::Point> poly = {{point_k, point_k1, point_k2}};
    findClosestPoint(pt, poly, B);
    T.x(poly[2].x() - poly[0].x());
    T.y(poly[2].y() - poly[0].y());
    T.z(poly[2].z() - poly[0].z());
  }

  return std::make_tuple(t, B, T, N);
}

/* Get a point, tangent, and second derivative of open quadratic uniform B-spline 
   between three points */
void
computeOpenUniformQuadraticBSpline(const double& t,
                                   const Uintah::Matrix3& splineMatrix,
                                   const Uintah::Point& point_k,
                                   const Uintah::Point& point_k1,
                                   const Uintah::Point& point_k2,
                                   Uintah::Point& B,
                                   Uintah::Vector& T,
                                   Uintah::Vector& N)
{
  Uintah::Vector A(1, t, t*t);
  Uintah::Vector dA(0, 1, 2*t);
  Uintah::Vector ddA(0, 0, 2);
  Uintah::Vector Px(point_k.x(), point_k1.x(), point_k2.x());
  Uintah::Vector Py(point_k.y(), point_k1.y(), point_k2.y());
  double bx = Uintah::Dot(A, (splineMatrix*Px*0.5));
  double by = Uintah::Dot(A, (splineMatrix*Py*0.5));
  double tx = Uintah::Dot(dA, (splineMatrix*Px*0.5));
  double ty = Uintah::Dot(dA, (splineMatrix*Py*0.5));
  double nx = Uintah::Dot(ddA, (splineMatrix*Px*0.5));
  double ny = Uintah::Dot(ddA, (splineMatrix*Py*0.5));
  B.x(bx); B.y(by);
  T.x(tx); T.y(ty);
  N.x(nx); N.y(ny);
}

/* 
 * Find point of intersection between two line segments 
 * Returns:
 *  bool  = true  if the point of intersection is inside the two line segments
 *        = false otherwise
 *  double = value of interpolation parameter t1 at the point of intersection
 *  double = value of interpolation parameter t2 at the point of intersection
 *  Point = point of intersection of the two straight lines
 */
std::tuple<bool, double, double, Uintah::Point> 
intersectionPoint(const Uintah::Point& P0,
                  const Uintah::Point& P1,
                  const Uintah::Point& Q0,
                  const Uintah::Point& Q1)
{
  auto P10 = P1 - P0;
  auto Q10 = Q1 - Q0;
  auto QP0 = Q0 - P0;

  Uintah::Matrix3 A(P10.x(), -Q10.x(), 0, P10.y(), -Q10.y(), 0, 0, 0, 1);

  // If the segments are parallel or collinear
  auto detA = A.Determinant();
  if (detA == 0.0) {
    return std::make_tuple(false, 0.0, 0.0, P0);
  }

  auto Ainv = A.Inverse();
  auto t = Ainv * QP0;

  auto intersect = P0 * (1 - t.x()) + P1 * t.x();

  // Point of intersection is outside segments
  if (!(isInBounds<double>(t.x(), 0, 1) && isInBounds<double>(t.y(), 0, 1))) {
    //std::cout << "No line line intersection " << t.x() << " " << t.y() << "\n";
    return std::make_tuple(false, t.x(), t.y(), intersect.asPoint());
  }

  return std::make_tuple(true, t.x(), t.y(), intersect.asPoint());
}

/*
 * Find side of line Q0-Q1 that a point P falls on 
 * 
 * Returns:
 *   0  if P is on the line through Q0-Q1
 *   1  if P is to the left of Q0-Q1
 *  -1  if P is to the right of Q0-Q1
 */
double findSide(const Uintah::Point& P,
                const Uintah::Point& Q0,
                const Uintah::Point& Q1)
{
  auto q0 = Q0 - P;
  auto q1 = Q1 - P;
  auto det = q0.x() * q1.y() - q1.x() * q0.y();
  if (det == 0) {
    return 0.0;
  }
  return std::copysign(1.0, det);
}

/* 
 * Find point of intersection between a line segment P0-P1 and a segment Q0-Q1 
 * of a polyline
 * and the side/s of the polyline segment that the two points of the line segment
 * fall on
 *
 * Returns:
 *  bool  = true  if the point of intersection is inside the two line segments
 *        = false otherwise
 *  double = s0 : side of Q0-Q1 that P0 falls on
 *  double = s1 : side of Q0-Q1 that P1 falls on
 *  double = value of interpolation parameter t_p at the point of intersection
 *           where P = (1 - t_p) P0 + t_p P1
 *  double = value of interpolation parameter t_q at the point of intersection
 *           where Q = (1 - t_q) Q0 + t_p Q1
 *  Point = point of intersection (Q)
 */
std::tuple<bool, double, double, double, double, Uintah::Point> 
intersectionPointAndSide(const Uintah::Point& P0,
                         const Uintah::Point& P1,
                         const Uintah::Point& Q0,
                         const Uintah::Point& Q1)
{
  // Find the position of the ends of the line segment relative to the polyline
  // segment
  auto s0 = findSide(P0, Q0, Q1);
  auto s1 = findSide(P1, Q0, Q1);

  // Find the intersection point
  auto P1P0 = P1 - P0;
  auto Q1Q0 = Q1 - Q0;
  auto Q0P0 = Q0 - P0;

  Uintah::Matrix3 J(P1P0.x(), -Q1Q0.x(), 0, P1P0.y(), -Q1Q0.y(), 0, 0, 0, 1);

  // If the segments are parallel or collinear
  auto detJ = J.Determinant();
  if (detJ == 0.0) {
    auto P0Q0 = P0 - Q0;
    auto P1Q0 = P1 - Q0;
    double cos_angle = 0.0;
    if (P0Q0.length() == 0.0) {
      cos_angle = Uintah::Dot(P1Q0/P1Q0.length(), Q1Q0/Q1Q0.length());
    } else {
      cos_angle = Uintah::Dot(P0Q0/P0Q0.length(), Q1Q0/Q1Q0.length());
    }
    if (std::abs(std::abs(cos_angle)-1) < 1.0e-6) {
      // Collinear
      auto t_q_p0 = P0Q0.x()/Q1Q0.x();
      auto t_q_p1 = P1Q0.x()/Q1Q0.x();
      //std::cout << "cos_angle = " << cos_angle << " t_q_p0 = " << t_q_p0 << " t_q_p1 = " << t_q_p1 << "\n";
      if (!(isInBounds<double>(t_q_p0, 0, 1)) && !(isInBounds<double>(t_q_p1, 0, 1))) {
        // No overlap between segments
        auto t_p = -1.0;
        auto t_q = -1.0;
        auto Q = Q0 * (1 - t_q) + Q1 * t_q;
        return std::make_tuple(false, s0, s1, t_p, t_q, Q.asPoint());
      } else {
        // Segments overlap (take midpoint of overlap)
        t_q_p0 = std::max(0.0, t_q_p0);
        t_q_p0 = std::min(1.0, t_q_p0);
        t_q_p1 = std::max(0.0, t_q_p1);
        t_q_p1 = std::min(1.0, t_q_p1);
        auto t_q = 0.5 * (t_q_p0 + t_q_p1);
        auto Q = Q0 * (1 - t_q) + Q1 * t_q;
        auto t_p = (Q.x() - P0.x())/(P1.x() - P0.x());
        return std::make_tuple(false, s0, s1, t_p, t_q, Q.asPoint());
      }
    } else {
      // Parallel but not collinear
      auto t_p = -1.0;
      auto t_q = -1.0;
      auto Q = Q0 * (1 - t_q) + Q1 * t_q;
      return std::make_tuple(false, s0, s1, t_p, t_q, Q.asPoint());
    }
  }

  // Normal case
  auto Jinv = J.Inverse();
  auto t = Jinv * Q0P0;
  auto t_p = t.x();
  auto t_q = t.y();

  auto Q = Q0 * (1 - t_q) + Q1 * t_q;

  // Point of intersection is outside segments
  if (!(isInBounds<double>(t_p, 0, 1) && isInBounds<double>(t_q, 0, 1))) {
    //std::cout << "No line line intersection " << t_p << " " << t_q << "\n";
    return std::make_tuple(false, s0, s1, t_p, t_q, Q.asPoint());
  }

  return std::make_tuple(true, s0, s1, t_p, t_q, Q.asPoint());
}

/* 
 * Find point of intersection between a polyline and a line segment 
 * using a linear search 
 * Returns:
 *  bool   = true  if the point of intersection is inside the 
 *                 polyline and the line segment
 *         = false otherwise
 *  size_t = index of first point on polyline segment
 *  double = value of interpolation parameter t at the point of intersection
 *  Point  = point of intersection of the polyline and the line segment
 */
std::tuple<bool, std::size_t, double, Uintah::Point> 
intersectionPointLinearSearch(const std::vector<Uintah::Point>& polyline,
                              const Uintah::Point& start_seg,
                              const Uintah::Point& end_seg)
{

  // Find intersection of a box containing the line segment with the polyline
  bool status = false;
  double t1, t2;
  Uintah::Point intersect;

  std::size_t index = 0;
  for (auto iter = polyline.begin(); iter < polyline.end()-1; ++iter) {
    auto start_poly = *iter;
    auto end_poly = *(iter+1);

    //std::cout << start_poly << " " << end_poly << " ";
    std::tie(status, t1, t2, intersect) = 
      intersectionPoint(start_poly, end_poly, start_seg, end_seg);
    //std::cout << status << " " << t1 << " " << t2 << " " << intersect << "\n";
    if (status) {
      return std::make_tuple(true, index, t2, intersect);
    }
    ++index;
  }

  return std::make_tuple(false, index, t2, intersect);
}

/* 
 * Find point of intersection between a polyline and a line segment 
 * using a binary search 
 * Returns:
 *  bool   = true  if the point of intersection is inside the 
 *                 polyline and the line segment
 *         = false otherwise
 *  size_t = index of first point on polyline segment
 *  double = value of interpolation parameter t_p at the point of intersection
 *  double = value of interpolation parameter t_q at the point of intersection
 *  Point  = point of intersection of the polyline and the line segment
 */
std::tuple<bool, std::size_t, double, double, Uintah::Point> 
intersectionPointBinarySearch(const std::vector<Uintah::Point>& polyline,
                              const Uintah::Point& start_seg,
                              const Uintah::Point& end_seg)
{
  bool status_left, status_right;
  double t_p_left, t_p_right;
  double t_q_left, t_q_right;
  Uintah::Point intersection_left, intersection_right;
  std::size_t index_left, index_right;

  std::tie(status_left, index_left, t_p_left, t_q_left, intersection_left,
           status_right, index_right, t_p_right, t_q_right, intersection_right) =
   Vaango::Util::intersectionLinePolyBinary(start_seg, end_seg, polyline.begin(), polyline.end(),
                                            polyline);

  if (status_left) {
    if (isInBounds<double>(t_p_left, 0, 1) && isInBounds<double>(t_q_left, 0, 1)) {
      return std::make_tuple(true, index_left, t_p_left, t_q_left, 
                             intersection_left);
    }
  }

  return std::make_tuple(false, index_left, t_p_left, t_q_left,
                         intersection_left);
}

std::tuple<bool, std::size_t, double, double, Uintah::Point, 
           bool, std::size_t, double, double, Uintah::Point>
intersectionLinePolyBinary(const Uintah::Point& P0,
                           const Uintah::Point& P1,
                           std::vector<Uintah::Point>::const_iterator poly_begin,
                           std::vector<Uintah::Point>::const_iterator poly_end,
                           const std::vector<Uintah::Point>& poly_orig)
{
  bool status_left, status_right;
  double side_p0_left, side_p0_right;
  double side_p1_left, side_p1_right;
  double t_p_left, t_p_right;
  double t_q_left, t_q_right;
  Uintah::Point intersection_left, intersection_right;

  auto num_pts = poly_end - poly_begin;
  auto mid_index = std::round(num_pts/2);

  auto index_left = poly_begin - poly_orig.begin();
  auto index_right = poly_begin + mid_index - poly_orig.begin();

  if (num_pts == 2) {
    std::tie(status_left, side_p0_left, side_p1_left, t_p_left, t_q_left, intersection_left) = 
      intersectionPointAndSide(P0, P1, *poly_begin, *(poly_end-1));

    return std::make_tuple(
      status_left, index_left, t_p_left, t_q_left, intersection_left,
      status_left, index_right, t_p_left, t_q_left, intersection_left);
  }

  std::tie(status_left, side_p0_left, side_p1_left, t_p_left, t_q_left, intersection_left) = 
      intersectionPointAndSide(P0, P1, *poly_begin, *(poly_begin+mid_index));

  std::tie(status_right, side_p0_right, side_p1_right, t_p_right, t_q_right, intersection_right) = 
      intersectionPointAndSide(P0, P1, *(poly_begin+mid_index), *(poly_end-1));

  if ((side_p0_left  > 0 && side_p1_left  > 0) || 
      (side_p0_right < 0 && side_p1_right < 0)) {
   
    std::tie(status_left,  index_left,  t_p_left,  t_q_left,  intersection_left,
             status_right, index_right, t_p_right, t_q_right, intersection_right) = 
      intersectionLinePolyBinary(P0, P1, poly_begin, poly_begin+mid_index+1, poly_orig);
  } else if ((side_p0_left  < 0 && side_p1_left  < 0) || 
             (side_p0_right > 0 && side_p1_right > 0)) {
    std::tie(status_left,  index_left,  t_p_left,  t_q_left,  intersection_left,
             status_right, index_right, t_p_right, t_q_right, intersection_right) = 
      intersectionLinePolyBinary(P0, P1, poly_begin+mid_index, poly_end, poly_orig);
  } else {
    if (t_q_left < 1.0) {
      std::tie(status_left,  index_left,  t_p_left,  t_q_left,  intersection_left,
               status_right, index_right, t_p_right, t_q_right, intersection_right) = 
        intersectionLinePolyBinary(P0, P1, poly_begin, poly_begin+mid_index+1, poly_orig);
    } else {
      std::tie(status_left,  index_left, t_p_left,  t_q_left,  intersection_left,
               status_right, index_right, t_p_right, t_q_right, intersection_right) = 
        intersectionLinePolyBinary(P0, P1, poly_begin+mid_index, poly_end, poly_orig);
    }
  }

  return std::make_tuple(
      status_left,  index_left,  t_p_left,  t_q_left,  intersection_left,
      status_right, index_right, t_p_right, t_q_right, intersection_right);
}


/* 
 * Find point of intersection between a quadratic B-spline and a line segment 
 * using Newton's method
 * Returns:
 *  bool   = true  if the point of intersection is inside the 
 *                 B-spline and the line segment
 *         = false otherwise
 *  double = value of segment interpolation parameter t at the point of intersection
 *  Point  = point of intersection of the B-spline and the line segment
 */
std::tuple<bool, Uintah::Vector, Uintah::Point> 
intersectionPointBSpline(const Uintah::Point& bezier_p0,
                         const Uintah::Point& bezier_p1, 
                         const Uintah::Point& bezier_p2, 
                         const Uintah::Point& seg_p0,
                         const Uintah::Point& seg_p1)
{
  Uintah::Vector t_old(0.5, 0.5, 0);
  Uintah::Vector t_new(-1, -1, 0);
  double dist = 1.0;
  int k = 0;
  constexpr double KMAX = 20;
  constexpr double TOLERANCE = 1.0e-10;
  while (dist > TOLERANCE && k < KMAX) {

    auto F_Jinv = evalFunctionJacobianInverse(bezier_p0, bezier_p1, bezier_p2, 
                                              seg_p0, seg_p1, t_old);
    auto F = F_Jinv.first;
    auto Jinv = F_Jinv.second;

    t_new = t_old - Jinv * F;
    //std::cout << "t_old = " << t_old << " t_new = " << t_new << "\n";
    dist = (t_new - t_old).length();
    t_old = t_new;
    ++k;
  }

  auto intersection = seg_p0 * (1 - t_new.y()) + seg_p1 * t_new.y();

  if (!(isInBounds<double>(t_new.x(), 0, 1) && isInBounds<double>(t_new.y(), 0, 1))) {
    //std::cout << "Bezier segment failed " << t_new << " " << intersection << "\n";
    return std::make_tuple(false, t_new, intersection.asPoint());
  }

  return std::make_tuple(true, t_new, intersection.asPoint());
}

/* 
 * Compute function and jacobian inverse for Newton method
 * to find intersection between a quadratic B-spline and a line segment 
 * Returns:
 *  Vector = value of function
 *  Matrix3 = value of jacobian inverse
 */
std::pair<Uintah::Vector, Uintah::Matrix3> 
evalFunctionJacobianInverse(const Uintah::Point& bezier_p0,
                            const Uintah::Point& bezier_p1, 
                            const Uintah::Point& bezier_p2, 
                            const Uintah::Point& seg_p0,
                            const Uintah::Point& seg_p1,
                            const Uintah::Vector& t)
{
  Uintah::Point B(0, 0, 0);
  Uintah::Vector T(0, 0, 0);
  Uintah::Vector N(0, 0, 0);
  computeOpenUniformQuadraticBSpline(t.x(), quadBSpline, bezier_p0, bezier_p1,
                                     bezier_p2, B, T, N);

  auto L = seg_p0 * (1 - t.y()) + seg_p1 * t.y();
  auto F = B.asVector() - L;
  //std::cout << "B = " << B << " T = " << T << " L = " << L << " F = " << F ;

  Uintah::Matrix3 J(0.0);
  J(0,0) = T.x();
  J(1,0) = T.y();
  J(0,1) = (seg_p0 - seg_p1).x();
  J(1,1) = (seg_p0 - seg_p1).y();
  J(2,2) = 1.0;
  //std::cout << " J = " << J << "\n";

  Uintah::Matrix3 Jinv;
  if (J.Determinant() == 0) {
    Jinv = Uintah::Matrix3(1.0e100);
  } else {
    Jinv = J.Inverse();
  }
  return std::make_pair(F, Jinv);
}

/* 
 * Find point of intersection between a quadratic B-spline fit to a polyline
 * and a line segment using Newton's method
 * Returns:
 *  bool   = true  if the point of intersection is inside the 
 *                 polyline and the line segment
 *         = false otherwise
 *  double = value of segment interpolation parameter t at the point of intersection
 *  Point  = point of intersection of the B-spline and the line segment
 */
std::tuple<bool, double, Uintah::Point> 
intersectionPointBSpline(const std::vector<Uintah::Point>& polyline,
                         const Uintah::Point& seg_p0,
                         const Uintah::Point& seg_p1)
{
  bool status;
  std::size_t index;
  double t_p, t_q;
  Uintah::Point intersection;

  std::tie(status, index, t_p, t_q, intersection) = intersectionPointBinarySearch(polyline, seg_p0, seg_p1);

  // The segment does not intersect the polyline
  if (!status) {
    //std::cout << "Does not intersect polyline " << t_p << " " << intersection << "\n";
    return std::make_tuple(false, t_p, intersection);
  }

  // Identify the three Bezier control points
  Uintah::Point bezier_p0, bezier_p1, bezier_p2;
  if (index == 0) {
    bezier_p0 = polyline[index];
    bezier_p1 = polyline[index+1];
    bezier_p2 = polyline[index+2];
  } else if (index == polyline.size() - 1) {
    bezier_p0 = polyline[index-2];
    bezier_p1 = polyline[index-1];
    bezier_p2 = polyline[index];
  } else {
    bezier_p0 = polyline[index-1];
    bezier_p1 = polyline[index];
    bezier_p2 = polyline[index+1];
  }

  // If the three control points are collinear just return the
  // intersection point already computed
  if (isCollinear(bezier_p0, bezier_p1, bezier_p2)) {
    std::cout << "Collinear " << t_p << " " << intersection << "\n";
    return std::make_tuple(true, t_p, intersection);
  }
  
  // Compute the intersection of the Bezier with segment 
  Uintah::Vector t1t2;
  std::tie(status, t1t2, intersection) = 
    intersectionPointBSpline(bezier_p0, bezier_p1, bezier_p2, seg_p0, seg_p1);

  // The segment does not intersect the Bezier curve in t \in [0, 1]
  // Look for the next three control points
  if (!status) {
    if (!(isInBounds<double>(t1t2.x(), 0, 1))) {
      if (index > 0 && index < polyline.size() - 2) {
        bezier_p0 = polyline[index];
        bezier_p1 = polyline[index+1];
        bezier_p2 = polyline[index+2];
        std::tie(status, t1t2, intersection) = 
          intersectionPointBSpline(bezier_p0, bezier_p1, bezier_p2, 
                                   seg_p0, seg_p1);
        if (!status) {
          std::cout << "Second segment failed " << t1t2 << " " << intersection << "\n";
          return std::make_tuple(false, t1t2.y(), intersection);
        }
      }
    } else {
      std::cout << "First segment failed " << t1t2 << " " << intersection << "\n";
      return std::make_tuple(false, t1t2.y(), intersection);
    }
  }

  return std::make_tuple(true, t1t2.y(), intersection);
}

/* 
 * Check if three points are collinear
 */
bool isCollinear(const Uintah::Point& p0, const Uintah::Point& p1,
                 const Uintah::Point& p2)
{
  Uintah::Vector vec_01 = (p1 - p0);
  Uintah::Vector vec_02 = (p2 - p0);
  vec_01.normalize();
  vec_02.normalize();
  double dot_product_minus_one = Uintah::Dot(vec_01, vec_02) - 1.0;  
  if (std::abs(dot_product_minus_one) < 1.0e-4) {
    return true;
  }
  return false; 
}

/*
 * Integrate the normalized deviatoric stress rate
 *
 *  sigma_s \dot(s_hat) = 2 G [eta - s_hat (s_hat : eta)]
 *
 * where
 *
 *  sigma_s = closest point on yield surface
 *  G = shear modulus
 *  dt = time increment
 *  sigma_elastic = stress at the start of the elastic-plastic step
 *  D = rate of deformation 
 *  sigma_init = initial value for Newton iterations
 *
 * Returns:
 *  s_hat
 */
Uintah::Matrix3 
  integrateNormalizedDeviatoricStressRate(double sigma_s,
                                          double G,
                                          double dt,
                                          const Uintah::Matrix3& sigma_elastic,
                                          const Uintah::Matrix3& D,
                                          const Uintah::Matrix3& sigma_init)
{
  double sig_e_m, sig_e_s, sig_i_m, sig_i_s;
  Tensor::Vector6Mandel I_hat, S_e_hat, S_i_hat;
  std::tie(sig_e_m, sig_e_s, I_hat, S_e_hat) = 
    Tensor::computeIsomorphicDecomposition(sigma_elastic);
  std::tie(sig_i_m, sig_i_s, I_hat, S_i_hat) = 
    Tensor::computeIsomorphicDecomposition(sigma_init);

  double d_m, d_s;
  Tensor::Vector6Mandel I, delta_eta;
  std::tie(d_m, d_s, I, delta_eta) = Tensor::computeVolDevDecomposition(D*dt);

  auto I6x6 = Tensor::IdentityMatrix6Mandel();
  auto S_old = S_i_hat;
  auto S_new = S_old;
  double dist = 1.0;
  int k = 0;
  constexpr double KMAX = 20;
  constexpr double TOLERANCE = 1.0e-10;
  while (dist > TOLERANCE && k < KMAX) {

    auto S_eta_inner = S_old.transpose() * delta_eta;
    auto S_eta_outer = Tensor::constructMatrix6Mandel(S_old, delta_eta);
    auto f = sigma_s * (S_old - S_e_hat) - 
                              2.0 * G * (delta_eta - S_old * S_eta_inner);
    auto J = (sigma_s + 2.0 * G * S_eta_inner(0,0)) * I6x6 + 2.0 * G * S_eta_outer;
    auto J_inv = J.inverse();
    auto S_new = S_old - J_inv * f;
    //std::cout << "S_new = " << S_new.transpose() << "\n"
    //          << " S_old = " << S_old.transpose() << "\n"
    //          << " J_inv = \n" << J_inv << "\n"
    //          << " f = " << f.transpose() << "\n";
    dist = (S_new - S_old).norm();
    S_old = S_new;
    ++k;
  }

  S_new /= S_new.norm();
  auto s_hat = Tensor::constructMatrix3(S_new);
  return s_hat; 
}
                                       

} // End namespace Util

} // end namespace Vaango
