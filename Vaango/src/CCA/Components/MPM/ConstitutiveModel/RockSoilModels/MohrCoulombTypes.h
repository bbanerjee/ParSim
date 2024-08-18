/*
 * The MIT License
 *
 * Copyright (c) 2015-2024 Biswajit Banerjee
 *
 * License for the specific language governing rights and limitations under
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef __MPM_CONSTITUTIVEMODEL_MODELS_TYPES_MOHRCOULOMB__
#define __MPM_CONSTITUTIVEMODEL_MODELS_TYPES_MOHRCOULOMB__

#include <Eigen/Dense>
#include <Core/Math/Matrix3.h>

namespace Uintah {

using Vector3 = Eigen::Matrix<double, 3, 1>;
using Vector6 = Eigen::Matrix<double, 6, 1>;
using Vector7 = Eigen::Matrix<double, 7, 1>;
using Vector8 = Eigen::Matrix<double, 8, 1>;

using Matrix33 = Eigen::Matrix<double, 3, 3>;
using Matrix66 = Eigen::Matrix<double, 6, 6>;
using Matrix67 = Eigen::Matrix<double, 6, 7>;
using Matrix77 = Eigen::Matrix<double, 7, 7>;
using Matrix88 = Eigen::Matrix<double, 8, 8>;

enum class RegionType
{
  APEX_REGION_1 = 1,
  APEX_REGION_2 = 2,
  APEX_REGION_3 = 3,
  APEX_REGION_4 = 4
};

std::ostream& operator<<(std::ostream& out, const RegionType& dc);

enum class DriftCorrection
{
  NO_CORRECTION = 1,
  CORRECTION_AT_BEGIN = 2,
  CORRECTION_AT_END = 3
};

std::ostream& operator<<(std::ostream& out, const DriftCorrection& dc);

enum class ToleranceMethod
{
  EPUS_RELATIVE_ERROR = 0,
  SLOAN = 1
};

enum class SolutionAlgorithm
{
  RUNGE_KUTTA_SECOND_ORDER_MODIFIED_EULER = 1,
  RUNGE_KUTTA_THIRD_ORDER_NYSTROM = 2,
  RUNGE_KUTTA_THIRD_ORDER_BOGACKI = 3,
  RUNGE_KUTTA_FOURTH_ORDER = 4,
  RUNGE_KUTTA_FIFTH_ORDER_ENGLAND = 5,
  RUNGE_KUTTA_FIFTH_ORDER_CASH = 6,
  RUNGE_KUTTA_FIFTH_ORDER_DORMAND = 7,
  RUNGE_KUTTA_FIFTH_ORDER_BOGACKI = 8,
  EXTRAPOLATION_BULIRSCH = 9
};

std::ostream& operator<<(std::ostream& out, const SolutionAlgorithm& sa);

enum class RetentionModel
{
  STATE_SURFACE = 1,
  VAN_GENUCHTEN = 2,
  GALLIPOLI = 3
};

Matrix33 toMatrix33(const Vector6& vec);
Vector6 toVector6(const Matrix33& mat);

Uintah::Matrix3 toMatrix3(const Vector6& vec);
Vector6 toVector6(const Uintah::Matrix3& mat);

Matrix33 toMatrix33Strain(const Vector6& vec);
Vector6 toVector6Strain(const Matrix33& mat);

Uintah::Matrix3 toMatrix3Strain(const Vector6& vec);
Vector6 toVector6Strain(const Uintah::Matrix3& mat);

} // end namespace Uintah

#endif //__MPM_CONSTITUTIVEMODEL_MODELS_TYPES_MOHRCOULOMB__
