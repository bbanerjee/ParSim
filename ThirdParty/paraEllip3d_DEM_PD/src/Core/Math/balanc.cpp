#include <Core/Types/realtypes.h>

namespace dem {
  
#define RADIX 2.0
  
  void balanc(REAL **a, int n) {
    int last,j,i;
    REAL s,r,g,f,c,sqrdx;
  
    sqrdx=RADIX*RADIX;
    last=0;
    while (last == 0) {
      last=1;
      for (i=1;i<=n;i++) {
	r=c=0.0;
	for (j=1;j<=n;j++)
	  if (j != i) {
	    c += fabs(a[j][i]);
	    r += fabs(a[i][j]);
	  }
	if (c && r) {
	  g=r/RADIX;
	  f=1.0;
	  s=c+r;
	  while (c<g) {
	    f *= RADIX;
	    c *= sqrdx;
	  }
	  g=r*RADIX;
	  while (c>g) {
	    f /= RADIX;
	    c /= sqrdx;
	  }
	  if ((c+r)/f < 0.95*s) {
	    last=0;
	    g=1.0/f;
	    for (j=1;j<=n;j++) a[i][j] *= g;
	    for (j=1;j<=n;j++) a[j][i] *= f;
	  }
	}
      }
    }
  }
  
#undef RADIX
  
} // namespace dem ends

/*
  11.6.1 Balancing
  The sensitivity of eigenvalues to rounding errors during the execution of some
  algorithms can be reduced by the procedure of balancing. The errors in the eigensystem
  found by a numerical procedure are generally proportional to the Euclidean norm
  of the matrix, that is, to the square root of the sum of the squares of the elements.
  The idea of balancing is to use similarity transformations to make corresponding
  rows and columns of the matrix have comparable norms, thus reducing the overall
  norm of the matrix while leaving the eigenvalues unchanged. A symmetric matrix is
  already balanced.

  Balancing is a procedure with of order N2 operations. Thus, the time taken
  by the procedure balance, given below, should never be significant compared to
  the total time required to find the eigenvalues. It is therefore recommended that
  you always balance nonsymmetric matrices. It never hurts, and it can substantially
  improve the accuracy of the eigenvalues computed for a badly balanced matrix.
  The actual algorithm used is due to Osborne, as discussed in [1]. It consists of a
  sequence of similarity transformations by diagonal matrices D. To avoid introducing
  rounding errors during the balancing process, the elements of D are restricted to be
  exact powers of the radix base employed for floating-point arithmetic (i.e., 2 for all
  modern machines, but 16 for some historical mainframe architectures). The output
  is a matrix that is balanced in the norm given by summing the absolute magnitudes
  of the matrix elements. This is more efficient than using the Euclidean norm, and
  equally effective: A large reduction in one norm implies a large reduction in the
  other.

  Note that if the off-diagonal elements of any row or column of a matrix are
  all zero, then the diagonal element is an eigenvalue. If the eigenvalue happens to
  be ill-conditioned (sensitive to small changes in the matrix elements), it will have
  relatively large errors when determined by the routine hqr (11.7). Had we merely
  inspected the matrix beforehand, we could have determined the isolated eigenvalue
  exactly and then deleted the corresponding row and column from the matrix. You
  should consider whether such a pre-inspection might be useful in your application.
  (For symmetric matrices, the routines we gave will determine isolated eigenvalues
  accurately in all cases.)

  The routine balance keeps track of the scale factors used in the balancing. If
  you are computing eigenvectors as well as eigenvalues, then the accumulated similarity
  transformation of the original matrix is undone by applying these scale factors
  in the routine balbak
*/
