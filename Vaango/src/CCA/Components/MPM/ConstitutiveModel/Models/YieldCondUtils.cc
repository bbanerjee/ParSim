#include <CCA/Components/MPM/ConstitutiveModel/Models/YieldCondUtils.h>
#include <Core/Geometry/Vector.h>
#include <algorithm>
#include <sstream>

namespace Vaango {

namespace Util {

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
  double delta = (end - start) / (double)num;

  std::vector<double> linspaced;
  for (int i = 0; i < num + 1; ++i) {
    linspaced.push_back(start + delta * (double)i);
  }
  return linspaced;
}

/* linspace function */
void
linspace(const double& start, const double& end, const int& num,
         std::vector<double>& linspaced)
{
  double delta = (end - start) / (double)num;

  for (int i = 0; i < num + 1; ++i) {
    linspaced.push_back(start + delta * (double)i);
  }
  return;
}

} // End namespace Util

} // end namespace Vaango