#ifndef __MATRIX3_ROTATION_H__
#define __MATRIX3_ROTATION_H__

#include <Core/Geometry/Vector.h>
#include <Core/Math/Matrix3.h>
#include <Eigen/Dense>

namespace Vaango {

  using SCIRun::Vector;
  using Uintah::Matrix3;

  // These are for rotation of matrices given an angle of rotation 
  // and an axis of rotation
  typedef Eigen::Matrix<double, 9, 1, Eigen::DontAlign> Vector9d;
  typedef Eigen::Matrix<double, 9, 9, Eigen::DontAlign> Matrix9d;

  /*!  Functions that rotate a second order tensor in Matrix3 form by an angle
       around an axis vector. 
       ** For Further information **  
          See Norris, Andrew N. "Euler-Rodrigues and Cayley Formulae for Rotation of Elasticity 
              Tensors." Mathematics and Mechanics of Solids 13.6 (2008): 465-498.
  */
  void rotateMatrix(const double& angle, const Vector& axis, const Matrix3& mat, Matrix3& mat_rot);
  void rotateMatrix(const Matrix9d& QQ, const Matrix3& mat, Matrix3& mat_rot);
  void formRotationMatrix(const double& angle, const Vector& axis, Matrix9d& QQ);
  void formSkewSymmetricMatrices(const Vector& vec, Matrix9d& PP, Matrix9d& PP1, Matrix9d& PP2);

}

#endif
