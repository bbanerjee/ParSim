/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

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
