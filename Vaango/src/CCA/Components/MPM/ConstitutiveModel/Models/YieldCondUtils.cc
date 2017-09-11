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

  //std::cout << "Normals:";
  //std::copy(normals.begin(), normals.end(),
  //          std::ostream_iterator<Uintah::Vector>(std::cout, " "));
  //std::cout << std::endl;

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
void
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
  return;
}

/* Find two yield surface segments that are closest to input point */
void
getClosestSegments(const Uintah::Point& pt,
                   const std::vector<Uintah::Point>& poly,
                   std::vector<Uintah::Point>& segments)
{
  // Set up the first segment to start from the end of the polygon
  // **TODO** Make sure that the second to last point is being chosen because
  // the
  //          polygon has been closed
  Uintah::Point p_prev = *(poly.rbegin() + 1);

  // Set up the second segment to start from the beginning of the polygon
  auto iterNext = poly.begin();
  ++iterNext;
  Uintah::Point p_next = *iterNext;
  Uintah::Point min_p_prev, min_p, min_p_next;

  double min_dSq = std::numeric_limits<double>::max();

  // Loop through the polygon
  Uintah::Point closest;
  for (const auto& poly_pt : poly) {

    // std::cout << "Pt = " << pt << std::endl
    //          << " Poly_pt = " << poly_pt << std::endl
    //          << " Prev = " << p_prev << std::endl
    //          << " Next = " << p_next << std::endl;

    std::vector<Uintah::Point> segment = { poly_pt, p_next };
    Vaango::Util::findClosestPoint(pt, segment, closest);

    // Compute distance sq
    double dSq = (pt - closest).length2();
    // std::cout << " distance = " << dSq << std::endl;
    // std::cout << " min_distance = " << min_dSq << std::endl;

    if (dSq - min_dSq < 1.0e-16) {
      min_dSq = dSq;
      min_p = closest;
      min_p_prev = p_prev;
      min_p_next = p_next;
    }

    ++iterNext;

    // Update prev and next
    p_prev = poly_pt;
    p_next = *iterNext;
  }

  // Return the three points
  segments.push_back(min_p_prev);
  segments.push_back(min_p);
  segments.push_back(min_p_next);
  // std::cout << "Closest_segments = " << min_p_prev << std::endl
  //          << min_p << std::endl
  //          << min_p_next << std::endl;

  return;
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
  linspaced.clear();
  if (num > 0) {
    double delta = (end - start) / (double)num;
    for (int i = 0; i < num + 1; ++i) {
      linspaced.push_back(start + delta * (double)i);
    }
  }
  return;
}

} // End namespace Util

} // end namespace Vaango