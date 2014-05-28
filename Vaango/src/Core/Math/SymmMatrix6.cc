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
