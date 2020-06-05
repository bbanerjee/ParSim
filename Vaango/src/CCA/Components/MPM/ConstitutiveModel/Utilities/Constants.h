#ifndef VAANGO_MPM_CONSTITUTIVE_MODEL_CONSTANTS_H
#define VAANGO_MPM_CONSTITUTIVE_MODEL_CONSTANTS_H

#include <Core/Geometry/Point.h>
#include <Core/Math/Matrix3.h>

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

} // end namespace Util

} // end namespace Vaango

#endif //VAANGO_MPM_CONSTITUTIVE_MODEL_CONSTANTS_H
