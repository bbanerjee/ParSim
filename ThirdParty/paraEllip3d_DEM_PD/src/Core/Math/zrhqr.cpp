#include <Core/Math/nrutil.h>

namespace dem {

#define NRANSI
#define MAXM 50

bool
zrhqr(REAL a[], int m, REAL rtr[], REAL rti[])
{
  void balanc(REAL * *a, int n);
  bool hqr(REAL * *a, int n, REAL wr[], REAL wi[]);
  int j, k;
  REAL **hess, xr, xi;

  if (m > MAXM || a[m] == 0.0) // nrerror("bad args in zrhqr");
    return false;

  hess = matrix(1, MAXM, 1, MAXM);
  for (k = 1; k <= m; k++) {
    hess[1][k] = -a[m - k] / a[m];
    for (j = 2; j <= m; j++)
      hess[j][k] = 0.0;
    if (k != m)
      hess[k + 1][k] = 1.0;
  }
  balanc(hess, m);
  if (!hqr(hess, m, rtr, rti)) {
    free_matrix(hess, 1, MAXM, 1, MAXM);
    return false;
  }

  for (j = 2; j <= m; j++) {
    xr = rtr[j];
    xi = rti[j];
    for (k = j - 1; k >= 1; k--) {
      if (rtr[k] <= xr)
        break;
      rtr[k + 1] = rtr[k];
      rti[k + 1] = rti[k];
    }
    rtr[k + 1] = xr;
    rti[k + 1] = xi;
  }
  free_matrix(hess, 1, MAXM, 1, MAXM);
  return true;
}

#undef MAXM
#undef NRANS

} // namespace dem ends
