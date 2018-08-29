#ifndef VAANGO_MPM_CONSTITUTIVE_MODEL_YIELDCOND_UTIL_H
#define VAANGO_MPM_CONSTITUTIVE_MODEL_YIELDCOND_UTIL_H

#include <Core/Geometry/Point.h>
#include <Core/Math/Matrix3.h>
#include <vector>

namespace Vaango {

namespace Util {

/* Constants */
static const double sqrt_two = std::sqrt(2.0);
static const double one_sqrt_two = 1.0 / sqrt_two;
static const double sqrt_three = std::sqrt(3.0);
static const double one_sqrt_three = 1.0 / sqrt_three;
static const double sqrt_two_third = std::sqrt(2.0/3.0);
static const double one_third(1.0 / 3.0);
static const double two_third(2.0 / 3.0);
static const double four_third = 4.0 / 3.0;
static const double one_sixth = 1.0 / 6.0;
static const double one_ninth = 1.0 / 9.0;
static const double pi = M_PI;
static const double pi_fourth = 0.25 * pi;
static const double pi_half = 0.5 * pi;
static const double large_number = 1.0e100;
static const Uintah::Matrix3 Identity(1,0,0,0,1,0,0,0,1);
static const Uintah::Matrix3 Zero(0.0);

/* Quadratic B-spline matrices */
static Uintah::Matrix3 quadBSplineLo(2, 0, 0, -4, 4, 0, 2, -3, 1);
static Uintah::Matrix3 quadBSplineHi(1, 1, 0, -2, 2, 0, 1, -3, 2);
static Uintah::Matrix3 quadBSpline(1, 1, 0, -2, 2, 0, 1, -2, 1);

/* Return codes */
enum class YieldStatus 
{
  IS_ELASTIC,
  HAS_YIELDED
};

/* Type trait for static asserts in class methods that should not be called */
template <typename T>
struct DoNotUse : std::false_type
{
};

/* Check whether a value is within bounds : value in [low, high] */
template <typename T>
bool isInBounds(const T& value, const T& low, const T& high, double epsilon=1.0e-4) {
  //constexpr double epsilon = 1.0e-4;
  return !(value < low - epsilon) && !(high + epsilon < value);
}

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
std::tuple<Uintah::Point, Uintah::Vector, Uintah::Vector>
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

/* Convert a pbar-sqrtJ2 point to z-rprime coordinates and vice-versa */
inline
void convertToZRprime(const double& sqrtKG, 
                      const double& p_bar, const double& sqrt_J2,
                      double& z, double& r_prime)
{
  z = -sqrt_three * p_bar;
  r_prime = sqrt_two * sqrt_J2 * sqrtKG;
}

inline
void revertFromZRprime(const double& sqrtKG, 
                       const double& z, const double& r_prime,
                       double& p_bar, double& sqrt_J2)
{
  p_bar = -z / sqrt_three;
  sqrt_J2 = r_prime /( sqrt_two * sqrtKG);
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
                    const Uintah::Point& end_pt_2);

/* 
 * Find point of intersection between a polyline and a line segment 
 * Returns:
 *  bool   = true  if the point of intersection is inside the 
 *                 polyline and the line segment
 *         = false otherwise
 *  size_t = index of first point of intersecting segment in polyline
 *  double = value of interpolation parameter t at the point of intersection
 *  Point  = point of intersection of the polyline and the line segment
 */
std::tuple<bool, std::size_t, double, Uintah::Point> 
  intersectionPoint(const std::vector<Uintah::Point>& poly,
                    const Uintah::Point& start_pt,
                    const Uintah::Point& end_pt);

/* 
 * Find point of intersection between a quadratic B-spline and a line segment 
 * Returns:
 *  bool   = true  if the point of intersection is inside the 
 *                 polyline and the line segment
 *         = false otherwise
 *  double = value of segment interpolation parameter t at the point of intersection
 *  Point  = point of intersection of the B-spline and the line segment
 */
std::tuple<bool, Uintah::Vector, Uintah::Point> 
  intersectionPointBSpline(const Uintah::Point& bezier_p0,
                           const Uintah::Point& bezier_p1, 
                           const Uintah::Point& bezier_p2, 
                           const Uintah::Point& seg_p0,
                           const Uintah::Point& seg_p1);

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
                              const Uintah::Vector& t);

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
                           const Uintah::Point& seg_p1);

/* 
 * Check if three points are collinear
 */
bool isCollinear(const Uintah::Point& p0, const Uintah::Point& p1,
                 const Uintah::Point& p2);

} // End namespace Util
} // End namespace Vaango

#endif // VAANGO_MPM_CONSTITUTIVE_MODEL_YIELDCOND_UTIL_H
