#include <CCA/Components/MPM/ConstitutiveModel/Models/YieldCondUtils.h>
#include <Core/Geometry/Vector.h>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iterator>

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

/* Find two yield surface segment that is closest to input point */
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
intersectionPoint(const Uintah::Point& start_pt_1,
                  const Uintah::Point& end_pt_1,
                  const Uintah::Point& start_pt_2,
                  const Uintah::Point& end_pt_2)
{
  auto x11 = start_pt_1.x();
  auto y11 = start_pt_1.y();
  auto x12 = end_pt_1.x();
  auto y12 = end_pt_1.y();

  auto x21 = start_pt_2.x();
  auto y21 = start_pt_2.y();
  auto x22 = end_pt_2.x();
  auto y22 = end_pt_2.y();

  auto a11 = x12 - x11;
  auto a12 = -(x22 - x21);
  auto a21 = y12 - y11;
  auto a22 = -(y22 - y21);

  auto b1 = x21 - x11;
  auto b2 = y21 - y11;

  auto detA = a11 * a22 - a12 * a21;

  // If the segments are parallel or collinear
  if (detA == 0.0) {
    return std::make_tuple(false, 0.0, 0.0, start_pt_1);
  }

  auto t1 = (a22*b1 - a12*b2)/detA;
  auto t2 = (-a21*b1 + a11*b2)/detA;

  auto x1 = (1 - t1) * x11 + t1 * x12;
  auto y1 = (1 - t1) * y11 + t1 * y12;

  // Point of intersection is outside segments
  if (t1 < 0 || t1 > 1 || t2 < 0 || t2 > 1) {
    return std::make_tuple(false, t1, t2, Uintah::Point(x1, y1, 0));
  }

  return std::make_tuple(true, t1, t2, Uintah::Point(x1, y1, 0));
}

/* 
 * Find point of intersection between a polyline and a line segment 
 * Returns:
 *  bool   = true  if the point of intersection is inside the 
 *                 polyline and the line segment
 *         = false otherwise
 *  size_t = index of first point on polyline segment
 *  double = value of interpolation parameter t at the point of intersection
 *  Point  = point of intersection of the polyline and the line segment
 */
std::tuple<bool, std::size_t, double, Uintah::Point> 
intersectionPoint(const std::vector<Uintah::Point>& polyline,
                  const Uintah::Point& start_pt,
                  const Uintah::Point& end_pt)
{
  bool status = false;
  double t1, t2;
  Uintah::Point intersect;

  std::size_t index = 0;
  for (auto iter = polyline.begin(); iter < polyline.end()-1; ++iter) {
    auto start_seg = *iter;
    auto end_seg = *(iter+1);

    //std::cout << start_seg << " " << end_seg << " ";
    std::tie(status, t1, t2, intersect) = 
      intersectionPoint(start_seg, end_seg, start_pt, end_pt);
    //std::cout << status << " " << t1 << " " << t2 << " " << intersect << "\n";
    if (status) {
      return std::make_tuple(true, index, t2, intersect);
    }
    ++index;
  }

  return std::make_tuple(false, index, t2, intersect);
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
  double t;
  Uintah::Point intersection;

  std::tie(status, index, t, intersection) = intersectionPoint(polyline, seg_p0, seg_p1);

  // The segment does not intersect the polyline
  if (!status) {
    return std::make_tuple(false, t, intersection);
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
    return std::make_tuple(true, t, intersection);
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
          return std::make_tuple(false, t1t2.y(), intersection);
        }
      }
    } else {
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

} // End namespace Util

} // end namespace Vaango
