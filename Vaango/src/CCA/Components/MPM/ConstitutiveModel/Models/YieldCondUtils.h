#ifndef VAANGO_MPM_CONSTITUTIVE_MODEL_YIELDCOND_UTIL_H
#define VAANGO_MPM_CONSTITUTIVE_MODEL_YIELDCOND_UTIL_H

#include <Core/Geometry/Point.h>
#include <vector>

namespace Vaango {

namespace Util {

/* Get the closest point on a polyline from a given point */
void findClosestPoint(const Uintah::Point& p,
                      const std::vector<Uintah::Point>& poly,
                      Uintah::Point& min_p);

/* Get closest segments */
void getClosestSegments(const Uintah::Point& pt,
                        const std::vector<Uintah::Point>& poly,
                        std::vector<Uintah::Point>& segments);

/* linspace functions */
void linspace(const double& start, const double& end, const int& num,
              std::vector<double>& linspaced);
std::vector<double> linspace(double start, double end, int num);
}
}

#endif // VAANGO_MPM_CONSTITUTIVE_MODEL_YIELDCOND_UTIL_H
