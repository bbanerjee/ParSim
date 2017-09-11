#ifndef VAANGO_MPM_CONSTITUTIVE_MODEL_YIELDCOND_UTIL_H
#define VAANGO_MPM_CONSTITUTIVE_MODEL_YIELDCOND_UTIL_H

#include <Core/Geometry/Point.h>
#include <vector>

namespace Vaango {

namespace Util {

/* Type trait for static asserts in class methods that should not be called */
template <typename T>
struct DoNotUse : std::false_type
{
};

/* A reverse range iterator 
   std::vector<int> v = {1, 2, 3, 4, 5};
   for (auto x : reverse(v)) {}
   From : https://gist.github.com/arvidsson/7231973 */
template <typename T>
class reverse_range
{
  T &x;

public:
  reverse_range(T &x) : x(x) {}
  auto begin() const -> decltype(this->x.rbegin())
  {
    return x.rbegin();
  }
  auto end() const -> decltype(this->x.rend())
  {
    return x.rend();
  }
};
template <typename T>
reverse_range<T> reverse(T &x)
{
  return reverse_range<T>(x);
}

/* Compute outward normal at each point of a polyline */
std::vector<Uintah::Vector>
computeNormals(const std::vector<Uintah::Point>& polyline);

/* Check whether three points are clockwise or anticloskwise order */
double
checkOrder(const Uintah::Point& p1, const Uintah::Point& p2,
           const Uintah::Point& p3);

/* Compute the convex hull of a set of points in 2D */
std::vector<Uintah::Point> 
convexHull2D(const std::vector<Uintah::Point>& points);

/* Get the closest point on a polyline from a given point */
void 
findClosestPoint(const Uintah::Point& p,
                 const std::vector<Uintah::Point>& poly,
                 Uintah::Point& min_p);

/* Get closest segments */
void 
getClosestSegments(const Uintah::Point& pt,
                   const std::vector<Uintah::Point>& poly,
                   std::vector<Uintah::Point>& segments);

/* linspace functions */
void 
linspace(const double& start, const double& end, const int& num,
         std::vector<double>& linspaced);
std::vector<double> 
linspace(double start, double end, int num);
}
}

#endif // VAANGO_MPM_CONSTITUTIVE_MODEL_YIELDCOND_UTIL_H
