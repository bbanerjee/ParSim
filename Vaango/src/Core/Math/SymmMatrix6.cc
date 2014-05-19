#include <Core/Math/SymmMatrix6.h>

using namespace Uintah;

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
