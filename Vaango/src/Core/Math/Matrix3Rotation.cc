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

#include <Core/Math/Matrix3Rotation.h>

#include <cmath>

namespace Vaango {

 
  // Utility function:
  //   This function rotates a 3x3 matrix given an angle of rotation and an axis of rotation
  //   (See Norris, Andrew N. "Euler-Rodrigues and Cayley Formulae for Rotation of Elasticity 
  //        Tensors." Mathematics and Mechanics of Solids 13.6 (2008): 465-498.)
  void 
  rotateMatrix(const double& angle, const Vector& axis, const Matrix3& mat, Matrix3& mat_rot)
  {
    // Compute the rotation matrix
    Matrix9d QQ;
    formRotationMatrix(angle, axis, QQ);

    // Rotate the matrix 
    rotateMatrix(QQ, mat, mat_rot);
  }

  // Utility function:
  //   This function rotates a 3x3 matrix given a 9x9 rotation matrix Q
  //   (See Norris, Andrew N. "Euler-Rodrigues and Cayley Formulae for Rotation of Elasticity 
  //        Tensors." Mathematics and Mechanics of Solids 13.6 (2008): 465-498.)
  void 
  rotateMatrix(const Matrix9d& QQ, const Matrix3& mat, Matrix3& mat_rot)
  {
    // Convert the input matrix into a 9-vector
    Vector9d mat_vec;
    mat_vec << mat(0,0), mat(1,1), mat(2,2), mat(1,2), mat(2,2), mat(0,1), mat(2,1), mat(0,2), mat(1,0);

    // Compute the rotated vector
    Vector9d mat_vec_rot = QQ*mat_vec;

    // Write out the rotated matrix in Matrix3 format
    mat_rot(0,0) = mat_vec_rot(0);
    mat_rot(1,1) = mat_vec_rot(1);
    mat_rot(2,2) = mat_vec_rot(2);
    mat_rot(1,2) = mat_vec_rot(3);
    mat_rot(2,2) = mat_vec_rot(4);
    mat_rot(0,1) = mat_vec_rot(5);
    mat_rot(2,1) = mat_vec_rot(6);
    mat_rot(0,2) = mat_vec_rot(7);
    mat_rot(1,0) = mat_vec_rot(8);
  }

  // Utility function:
  //   This function creates a 9x9 rotation matrix Q given an angle of rotation and 
  //   and axis of rotation
  //   (See Norris, Andrew N. "Euler-Rodrigues and Cayley Formulae for Rotation of Elasticity 
  //        Tensors." Mathematics and Mechanics of Solids 13.6 (2008): 465-498.)
  void 
  formRotationMatrix(const double& angle, const Vector& axis, Matrix9d& QQ)
  {
    Matrix9d PP, P1, P2;
    formSkewSymmetricMatrices(axis, PP, P1, P2);

    QQ = Matrix9d::Identity();

    // First term 
    double s1 = std::sin(angle);
    double c1 = std::cos(angle);
    QQ += (P1*s1 + P1*P1*(1.0 - c1));

    // Second term 
    double s2 = std::sin(2.0*angle);
    double c2 = std::cos(2.0*angle);
    QQ += (P2*s2 + P2*P2*(1.0 - c2));
  }

  // Utility function:
  //   This function creates three 9x9 skew symmetric matrices (P, P1, P2) given 
  //   an input axial vector (vec).
  //   (See Norris, Andrew N. "Euler-Rodrigues and Cayley Formulae for Rotation of Elasticity 
  //        Tensors." Mathematics and Mechanics of Solids 13.6 (2008): 465-498.)
  void
  formSkewSymmetricMatrices(const Vector& vec, Matrix9d& PP, Matrix9d& PP1, Matrix9d& PP2) 
  {
    double p1 = vec[0];
    double p2 = vec[1];
    double p3 = vec[2];

    Eigen::Matrix3d Zero = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d XX, YY, ZZ;
    XX << 0, p3, 0, 0, 0, p1, p2, 0, 0;
    YY << 0, p2, 0, 0, 0, p3, p1, 0, 0;
    ZZ << 0, p1, 0, 0, 0, p2, p3, 0, 0;

    Matrix9d RR;
    RR << Zero , YY , YY , ZZ , Zero , XX , ZZ , XX , Zero;
  
    PP = RR - RR.transpose();

    Matrix9d One = Matrix9d::Identity();
    PP1 = PP*(PP*PP+One*4.0)*(1.0/3.0);

    PP2 = PP*(PP*PP+One)*(-1.0/6.0);
  }

} // end namespace Vaango

