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

#include <Core/Math/SymmMatrix6.h>

using namespace Uintah;


/*! Matrix-vector product:
    1) The input matrix is converted to a 6x1 vector in Voigt notation.
    2) The vector is multiplied with the 6x6 matrix
    3) The resulting vector is converted back into a 3x3 matrix
 */
Matrix3 
SymmMatrix6::operator*(const Matrix3& mat) const
{
  Vector6d vec;
  vec << mat(0,0), mat(1,1), mat(2,2), 2.0*mat(1,2), 2.0*mat(2,0), 2.0*mat(0,1);
  Vector6d prod = d_mat6*vec;
  return Matrix3(prod(0), prod(5), prod(4), prod(5), prod(1), prod(3), prod(4), prod(3), prod(2));
}

/*! Compute the inverse of a SymmMatrix6 */
void 
SymmMatrix6::inverse(SymmMatrix6& inv)
{
  Matrix6d invMat = d_mat6.inverse();
  inv = invMat;
}

/*! Compute the eigenvalues and eigenvectors of a SymmMatrix6 */
void
SymmMatrix6::eigen(Vector6d& eval, Matrix6d& evec)
{
  Eigen::SelfAdjointEigenSolver<Matrix6d> es(d_mat6);

  // The eigenvalues (unsorted)
  eval = es.eigenvalues();

  // The columns of evec are the eigenvectors
  evec = es.eigenvectors();
  
}
