#ifndef VAANGO_MPM_CONSTITUTIVE_MODEL_YIELDCOND_UTIL_H
#define VAANGO_MPM_CONSTITUTIVE_MODEL_YIELDCOND_UTIL_H

#include <Core/Geometry/Point.h>
#include <Core/Math/Matrix3.h>
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
double 
findClosestPoint(const Uintah::Point& p,
                 const std::vector<Uintah::Point>& poly,
                 Uintah::Point& min_p);

/* Get closest segments */
std::size_t 
getClosestSegments(const Uintah::Point& pt,
                   const std::vector<Uintah::Point>& poly,
                   std::vector<Uintah::Point>& segments);

/* linspace functions */
void 
linspace(const double& start, const double& end, const int& num,
         std::vector<double>& linspaced);
std::vector<double> 
linspace(double start, double end, int num);

/* Quadratic B-spline matrices */
static Uintah::Matrix3 quadBSplineLo(2, 0, 0, -4, 4, 0, 2, -3, 1);
static Uintah::Matrix3 quadBSplineHi(1, 1, 0, -2, 2, 0, 1, -3, 2);
static Uintah::Matrix3 quadBSpline(1, 1, 0, -2, 2, 0, 1, -2, 1);

/* Create open quadratic uniform B-spline approximating a polyline */
void
computeOpenUniformQuadraticBSpline(const std::vector<Uintah::Point>& polyline,
                                   size_t ptsPerSegment,
                                   std::vector<Uintah::Point>& spline);

/* Create open quadratic uniform B-spline approximating a segment of a polyline */
void
computeOpenUniformQuadraticBSpline(const std::vector<Uintah::Point>& polyline,
                                   size_t segmentStartIndex,
                                   size_t segmentEndIndex,
                                   size_t ptsPerSegment,
                                   std::vector<Uintah::Point>& spline);

/* Get a single point on open quadratic uniform B-spline between three points */
Uintah::Point
computeOpenUniformQuadraticBSpline(const double& t,
                                   const Uintah::Matrix3& splineMatrix,
                                   const Uintah::Point& point_k,
                                   const Uintah::Point& point_k1,
                                   const Uintah::Point& point_k2);

/* Find closest point to open quadratic uniform B-spline approximating a 
   segment of a polyline */
std::tuple<Uintah::Point, Uintah::Vector>
computeClosestPointQuadraticBSpline(const Uintah::Point pt,
                                    const std::vector<Uintah::Point>& polyline,
                                    size_t segmentStartIndex,
                                    size_t segmentEndIndex);

/* Find closest point on a quadratic B-spline using Newton iterations */
std::tuple<double, Uintah::Point, Uintah::Vector, Uintah::Vector>
findClosestPointToQuadraticBSplineNewton(const Uintah::Point pt, 
                                         const Uintah::Matrix3& splineMatrix,
                                         const Uintah::Point& point_k,
                                         const Uintah::Point& point_k1,
                                         const Uintah::Point& point_k2);

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
                                   Uintah::Vector& N);
}
}

#endif // VAANGO_MPM_CONSTITUTIVE_MODEL_YIELDCOND_UTIL_H
