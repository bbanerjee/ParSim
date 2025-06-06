#ifndef VAANGO_MPM_CONSTITUTIVE_MODEL_TENSOR_UTIL_H
#define VAANGO_MPM_CONSTITUTIVE_MODEL_TENSOR_UTIL_H

#include <Core/Geometry/Point.h>
#include <Core/Math/Matrix3.h>
#include <vector>
#include <Eigen/Core>

namespace Vaango {

namespace Tensor {

/* Types */
using Vector6Mandel = Eigen::Matrix<double, 6, 1>;
using Matrix6Mandel = Eigen::Matrix<double, 6, 6>;
using Vector9Mandel = Eigen::Matrix<double, 9, 1>;
using Matrix9Mandel = Eigen::Matrix<double, 9, 9>;

/* Constants */
static const double sqrt_two = std::sqrt(2.0);
static const double one_sqrt_two = 1.0 / sqrt_two;
static const double sqrt_three = std::sqrt(3.0);
static const double one_sqrt_three = 1.0 / sqrt_three;
static constexpr double half = 0.5;
static constexpr double third = 1.0/3.0;
static constexpr double two_third = 2.0/3.0;
static constexpr double pi = M_PI;
static constexpr double large_number = 1.0e100;
static const Uintah::Matrix3 Identity(1,0,0,0,1,0,0,0,1);
static const Uintah::Matrix3 Zero(0.0);

/**
 * Convert a symmetric rank-2 tensor (3 x 3 symmetric matrix) into a 6x1 vector 
 * using Mandel notation (11, 22, 33, 23, 31, 12)
 */
inline Vector6Mandel
constructVector6Mandel(const Uintah::Matrix3& A) {
  Vector6Mandel vec;
  vec << A(0,0), A(1,1), A(2,2), 
         sqrt_two * A(1,2), sqrt_two * A(2,0), sqrt_two * A(0,1);
  return vec;
}

/**
 * Convert a 6x1 vector in Mandel notation to a 3x3 symmetric matrix
 * using Mandel notation (11, 22, 33, 23, 31, 12)
 */
inline Uintah::Matrix3
constructMatrix3(const Vector6Mandel& a) {
  Uintah::Matrix3 mat(a(0), one_sqrt_two * a(5), one_sqrt_two * a(4),
       one_sqrt_two * a(5),                a(1), one_sqrt_two * a(3),
       one_sqrt_two * a(4), one_sqrt_two * a(3),                a(2));
  return mat;
}

/**
 * Convert a rank-2 tensor (3 x 3 matrix) into a 9x1 vector 
 * using Mandel notation (11, 22, 33, 23, 31, 12, 32, 13, 21)
 */
inline Vector9Mandel
constructVector9Mandel(const Uintah::Matrix3& A) {
  Vector9Mandel vec;
  vec << A(0,0), A(1,1), A(2,2), 
         one_sqrt_two * (A(1,2) + A(2,1)), 
         one_sqrt_two * (A(2,0) + A(0,2)), 
         one_sqrt_two * (A(0,1) + A(1,0)),
         one_sqrt_two * (A(2,1) - A(1,2)), 
         one_sqrt_two * (A(0,2) - A(2,0)), 
         one_sqrt_two * (A(1,0) - A(0,1));
  return vec;
}

/**
 * Construct a 9x1 vector from the outer product of 2 3x1 vectors
 * using Mandel notation (11, 22, 33, 23, 31, 12)
 */
inline Vector9Mandel
constructVector9Mandel(const Uintah::Vector& a, const Uintah::Vector& b) {
  Vector9Mandel vec;
  vec << a[0]*b[0], a[1]*b[1], a[2]*b[2], 
         one_sqrt_two * (a[1]*b[2] + a[2]*b[1]), 
         one_sqrt_two * (a[2]*b[0] + a[0]*b[2]), 
         one_sqrt_two * (a[0]*b[1] + a[1]*b[0]),
         one_sqrt_two * (a[2]*b[1] - a[1]*b[2]), 
         one_sqrt_two * (a[0]*b[2] - a[2]*b[0]), 
         one_sqrt_two * (a[1]*b[0] - a[0]*b[1]);
  return vec;
}

/**
 * Convert a rank-4 tensor with major and minor symmetry (outer product of two 
 * symmetric 3 x 3 matrices) into a 6x6 matrix
 * using Mandel notation (11, 22, 33, 23, 31, 12)
 */
inline Matrix6Mandel
constructMatrix6Mandel(const Uintah::Matrix3& A, const Uintah::Matrix3& B) {
  Matrix6Mandel mat;
  mat << A(0,0) * B(0,0), A(0,0) * B(1,1), A(0,0) * B(2,2), 
         sqrt_two * A(0,0) * B(1,2), sqrt_two * A(0,0) * B(2,0), sqrt_two * A(0,0) * B(0,1),
         A(1,1) * B(0,0), A(1,1) * B(1,1), A(1,1) * B(2,2), 
         sqrt_two * A(1,1) * B(1,2), sqrt_two * A(1,1) * B(2,0), sqrt_two * A(1,1) * B(0,1),
         A(2,2) * B(0,0), A(2,2) * B(1,1), A(2,2) * B(2,2), 
         sqrt_two * A(2,2) * B(1,2), sqrt_two * A(2,2) * B(2,0), sqrt_two * A(2,2) * B(0,1),
         sqrt_two * A(1,2) * B(0,0), sqrt_two * A(1,2) * B(1,1), sqrt_two * A(1,2) * B(2,2), 
         2 * A(1,2) * B(1,2), 2 * A(1,2) * B(2,0), 2 * A(1,2) * B(0,1),
         sqrt_two * A(2,0) * B(0,0), sqrt_two * A(2,0) * B(1,1), sqrt_two * A(2,0) * B(2,2), 
         2 * A(2,0) * B(1,2), 2 * A(2,0) * B(2,0), 2 * A(2,0) * B(0,1),
         sqrt_two * A(0,1) * B(0,0), sqrt_two * A(0,1) * B(1,1), sqrt_two * A(0,1) * B(2,2), 
         2 * A(0,1) * B(1,2), 2 * A(0,1) * B(2,0), 2 * A(0,1) * B(0,1);
  return mat;
}

/**
 * Convert a rank-4 tensor with major and minor symmetry (outer product of two 
 * symmetric 6 x 1 vectors) into a 6x6 matrix
 * using Mandel notation (11, 22, 33, 23, 31, 12) -> (0, 1, 2, 3, 4, 5)
 */
inline Matrix6Mandel
constructMatrix6Mandel(const Vector6Mandel& a, const Vector6Mandel& b) {
  Matrix6Mandel mat;
  mat << a(0) * b(0), a(0) * b(1), a(0) * b(2), a(0) * b(3), a(0) * b(4), a(0) * b(5),
         a(1) * b(0), a(1) * b(1), a(1) * b(2), a(1) * b(3), a(1) * b(4), a(1) * b(5),
         a(2) * b(0), a(2) * b(1), a(2) * b(2), a(2) * b(3), a(2) * b(4), a(2) * b(5),
         a(3) * b(0), a(3) * b(1), a(3) * b(2), a(3) * b(3), a(3) * b(4), a(3) * b(5),
         a(4) * b(0), a(4) * b(1), a(4) * b(2), a(4) * b(3), a(4) * b(4), a(4) * b(5),
         a(5) * b(0), a(5) * b(1), a(5) * b(2), a(5) * b(3), a(5) * b(4), a(5) * b(5);
  return mat;
}

/**
 * Contruct a 9x9 matrix from the outer product of two 3x3 matrices
 * using Mandel notation (11, 22, 33, 23, 31, 12)
 */
inline Matrix9Mandel
constructMatrix9Mandel(const Uintah::Matrix3& A, const Uintah::Matrix3& B) {
  Matrix9Mandel mat;
  mat << A(0,0) * B(0,0), A(0,0) * B(1,1), A(0,0) * B(2,2), 
         one_sqrt_two * A(0,0) * (B(1,2) + B(2,1)), 
         one_sqrt_two * A(0,0) * (B(2,0) + B(0,2)), 
         one_sqrt_two * A(0,0) * (B(0,1) + B(1,0)),
         -one_sqrt_two * A(0,0) * (B(1,2) - B(2,1)), 
         -one_sqrt_two * A(0,0) * (B(2,0) - B(0,2)), 
         -one_sqrt_two * A(0,0) * (B(0,1) - B(1,0)),

         A(1,1) * B(0,0), A(1,1) * B(1,1), A(1,1) * B(2,2), 
         one_sqrt_two * A(1,1) * (B(1,2) + B(2,1)), 
         one_sqrt_two * A(1,1) * (B(2,0) + B(0,2)), 
         one_sqrt_two * A(1,1) * (B(0,1) + B(1,0)),
         -one_sqrt_two * A(1,1) * (B(1,2) - B(2,1)), 
         -one_sqrt_two * A(1,1) * (B(2,0) - B(0,2)), 
         -one_sqrt_two * A(1,1) * (B(0,1) - B(1,0)),

         A(2,2) * B(0,0), A(2,2) * B(1,1), A(2,2) * B(2,2), 
         one_sqrt_two * A(2,2) * (B(1,2) + B(2,1)), 
         one_sqrt_two * A(2,2) * (B(2,0) + B(0,2)), 
         one_sqrt_two * A(2,2) * (B(0,1) + B(1,0)),
         -one_sqrt_two * A(2,2) * (B(1,2) - B(2,1)), 
         -one_sqrt_two * A(2,2) * (B(2,0) - B(0,2)), 
         -one_sqrt_two * A(2,2) * (B(0,1) - B(1,0)),

         one_sqrt_two * (A(1,2) + A(2,1)) * B(0,0), 
         one_sqrt_two * (A(1,2) + A(2,1)) * B(1,1), 
         one_sqrt_two * (A(1,2) + A(2,1)) * B(2,2), 
         half * (A(1,2) + A(2,1)) * (B(1,2) + B(2,1)), 
         half * (A(1,2) + A(2,1)) * (B(2,0) + B(0,2)), 
         half * (A(1,2) + A(2,1)) * (B(0,1) + B(1,0)),
         -half * (A(1,2) + A(2,1)) * (B(1,2) - B(2,1)), 
         -half * (A(1,2) + A(2,1)) * (B(2,0) - B(0,2)), 
         -half * (A(1,2) + A(2,1)) * (B(0,1) - B(1,0)),

         one_sqrt_two * (A(2,0) + A(0,2)) * B(0,0), 
         one_sqrt_two * (A(2,0) + A(0,2)) * B(1,1), 
         one_sqrt_two * (A(2,0) + A(0,2)) * B(2,2), 
         half * (A(2,0) + A(0,2)) * (B(1,2) + B(2,1)), 
         half * (A(2,0) + A(0,2)) * (B(2,0) + B(0,2)), 
         half * (A(2,0) + A(0,2)) * (B(0,1) + B(1,0)),
         -half * (A(2,0) + A(0,2)) * (B(1,2) - B(2,1)), 
         -half * (A(2,0) + A(0,2)) * (B(2,0) - B(0,2)), 
         -half * (A(2,0) + A(0,2)) * (B(0,1) - B(1,0)),

         one_sqrt_two * (A(0,1) + A(1,0)) * B(0,0), 
         one_sqrt_two * (A(0,1) + A(1,0)) * B(1,1), 
         one_sqrt_two * (A(0,1) + A(1,0)) * B(2,2), 
         half * (A(0,1) + A(1,0)) * (B(1,2) + B(2,1)), 
         half * (A(0,1) + A(1,0)) * (B(2,0) + B(0,2)), 
         half * (A(0,1) + A(1,0)) * (B(0,1) + B(1,0)),
         -half * (A(0,1) + A(1,0)) * (B(1,2) - B(2,1)), 
         -half * (A(0,1) + A(1,0)) * (B(2,0) - B(0,2)), 
         -half * (A(0,1) + A(1,0)) * (B(0,1) - B(1,0)),

         one_sqrt_two * (A(2,1) - A(1,2)) * B(0,0), 
         one_sqrt_two * (A(2,1) - A(1,2)) * B(1,1), 
         one_sqrt_two * (A(2,1) - A(1,2)) * B(2,2), 
         half * (A(2,1) - A(1,2)) * (B(1,2) + B(2,1)), 
         half * (A(2,1) - A(1,2)) * (B(2,0) + B(0,2)), 
         half * (A(2,1) - A(1,2)) * (B(0,1) + B(1,0)),
         -half * (A(2,1) - A(1,2)) * (B(1,2) - B(2,1)), 
         -half * (A(2,1) - A(1,2)) * (B(2,0) - B(0,2)), 
         -half * (A(2,1) - A(1,2)) * (B(0,1) - B(1,0)),

         one_sqrt_two * (A(0,2) - A(2,0)) * B(0,0), 
         one_sqrt_two * (A(0,2) - A(2,0)) * B(1,1), 
         one_sqrt_two * (A(0,2) - A(2,0)) * B(2,2), 
         half * (A(0,2) - A(2,0)) * (B(1,2) + B(2,1)), 
         half * (A(0,2) - A(2,0)) * (B(2,0) + B(0,2)), 
         half * (A(0,2) - A(2,0)) * (B(0,1) + B(1,0)),
         -half * (A(0,2) - A(2,0)) * (B(1,2) - B(2,1)), 
         -half * (A(0,2) - A(2,0)) * (B(2,0) - B(0,2)), 
         -half * (A(0,2) - A(2,0)) * (B(0,1) - B(1,0)),

         one_sqrt_two * (A(1,0) - A(0,1)) * B(0,0), 
         one_sqrt_two * (A(1,0) - A(0,1)) * B(1,1), 
         one_sqrt_two * (A(1,0) - A(0,1)) * B(2,2), 
         half * (A(1,0) - A(0,1)) * (B(1,2) + B(2,1)), 
         half * (A(1,0) - A(0,1)) * (B(2,0) + B(0,2)), 
         half * (A(1,0) - A(0,1)) * (B(0,1) + B(1,0)),
         -half * (A(1,0) - A(0,1)) * (B(1,2) - B(2,1)), 
         -half * (A(1,0) - A(0,1)) * (B(2,0) - B(0,2)), 
         -half * (A(1,0) - A(0,1)) * (B(0,1) - B(1,0));
  return mat;
}

/**
 * Contruct a 9x9 matrix from the outer product of two 3x1 vectors and a 3x3 matrix 
 * using Mandel notation (11, 22, 33, 23, 31, 12)
 */
inline Matrix9Mandel
constructMatrix9Mandel(const Uintah::Vector& a, const Uintah::Vector& b,
                       const Uintah::Matrix3& B) {
  Uintah::Matrix3 A(a, b);
  return constructMatrix9Mandel(A, B);
}

/**
 * Contruct a 9x9 matrix from the outer product of two 6x1 vectors
 * using Mandel notation (11, 22, 33, 23, 31, 12)
 */
inline Matrix9Mandel
constructMatrix9Mandel(const Vector9Mandel& a, const Vector9Mandel& b) {
  Matrix9Mandel mat;
  mat << a(0) * b(0), a(0) * b(1), a(0) * b(2), 
         a(0) * b(3), a(0) * b(4), a(0) * b(5),
         a(0) * b(6), a(0) * b(7), a(0) * b(8),

         a(1) * b(0), a(1) * b(1), a(1) * b(2), 
         a(1) * b(3), a(1) * b(4), a(1) * b(5),
         a(1) * b(6), a(1) * b(7), a(1) * b(8),

         a(2) * b(0), a(2) * b(1), a(2) * b(2), 
         a(2) * b(3), a(2) * b(4), a(2) * b(5),
         a(2) * b(6), a(2) * b(7), a(2) * b(8),

         a(3) * b(0), a(3) * b(1), a(3) * b(2), 
         a(3) * b(3), a(3) * b(4), a(3) * b(5),
         a(3) * b(6), a(3) * b(7), a(3) * b(8),

         a(4) * b(0), a(4) * b(1), a(4) * b(2), 
         a(4) * b(3), a(4) * b(4), a(4) * b(5),
         a(4) * b(6), a(4) * b(7), a(4) * b(8),

         a(5) * b(0), a(5) * b(1), a(5) * b(2), 
         a(5) * b(3), a(5) * b(4), a(5) * b(5),
         a(5) * b(6), a(5) * b(7), a(5) * b(8),

         a(6) * b(0), a(6) * b(1), a(6) * b(2), 
         a(6) * b(3), a(6) * b(4), a(6) * b(5),
         a(6) * b(6), a(6) * b(7), a(6) * b(8),

         a(7) * b(0), a(7) * b(1), a(7) * b(2), 
         a(7) * b(3), a(7) * b(4), a(7) * b(5),
         a(7) * b(6), a(7) * b(7), a(7) * b(8),

         a(8) * b(0), a(8) * b(1), a(8) * b(2), 
         a(8) * b(3), a(8) * b(4), a(8) * b(5),
         a(8) * b(6), a(8) * b(7), a(8) * b(8);
  return mat;
}

/**
 * Contruct a 6x6 matrix corresponding to symmetric rank-4 identity tensor
 * using Mandel notation (11, 22, 33, 23, 31, 12) -> (0, 1, 2, 3, 4, 5)
 */
inline Matrix6Mandel IdentityMatrix6Mandel() {
  return Matrix6Mandel::Identity();
}

/**
 * Compute volumetric-deviatoric decomposition of rank-2 matrix in Mandel form
 */
inline std::tuple<double, double, Vector6Mandel, Vector6Mandel>
  computeVolDevDecomposition(const Uintah::Matrix3& tensor2)
{
  Vector6Mandel tensor2_mandel = constructVector6Mandel(tensor2);
  double I1 = tensor2_mandel(0,0) + tensor2_mandel(1,0) + tensor2_mandel(2,0);
  Vector6Mandel I_mandel = constructVector6Mandel(Identity);
  Vector6Mandel dev_mandel = tensor2_mandel - I_mandel * (third * I1);
  double sqrt_J2 = std::sqrt(0.5 * dev_mandel.transpose() * dev_mandel);
  return std::make_tuple(I1, sqrt_J2, I_mandel, dev_mandel);
}

/**
 * Compute isomorphic decomposition of stress matrix in Mandel form
 */
inline std::tuple<double, double, Vector6Mandel, Vector6Mandel>
  computeIsomorphicDecomposition(const Uintah::Matrix3& stress)
{
  Vector6Mandel sigma_mandel = constructVector6Mandel(stress);
  double I1 = sigma_mandel(0,0) + sigma_mandel(1,0) + sigma_mandel(2,0);
  Vector6Mandel I_mandel = constructVector6Mandel(Identity);
  Vector6Mandel s_mandel = sigma_mandel - I_mandel * (third * I1);

  double sigma_m = one_sqrt_three * I1;
  double sigma_s = std::sqrt(s_mandel.transpose() * s_mandel);
  I_mandel *= one_sqrt_three;
  s_mandel /= (sigma_s + 1.0e-16);

  return std::make_tuple(sigma_m, sigma_s, I_mandel, s_mandel);
}

} // End namespace Tensor
} // End namespace Vaango

#endif // VAANGO_MPM_CONSTITUTIVE_MODEL_TENSOR_UTIL_H
