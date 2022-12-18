/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
   utility file nrutil.h.  Do not confuse this file with the same-named
   file nrutil.h that is supplied in the 'misc' subdirectory.
   *That* file is the one from the book, and contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only ANSI C.               */

#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

#include <Core/Types/RealTypes.h>

namespace dem {

[[maybe_unused]] static REAL sqrarg;
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)

[[maybe_unused]] static REAL dsqrarg;
#define DSQR(a) ((dsqrarg = (a)) == 0.0 ? 0.0 : dsqrarg * dsqrarg)

[[maybe_unused]] static REAL dmaxarg1, dmaxarg2;
#define DMAX(a, b)                                                             \
  (dmaxarg1 = (a), dmaxarg2 = (b),                                             \
   (dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))

[[maybe_unused]] static REAL dminarg1, dminarg2;
#define DMIN(a, b)                                                             \
  (dminarg1 = (a), dminarg2 = (b),                                             \
   (dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))

[[maybe_unused]] static REAL maxarg1, maxarg2;
#define FMAX(a, b)                                                             \
  (maxarg1 = (a), maxarg2 = (b), (maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

[[maybe_unused]] static REAL minarg1, minarg2;
#define FMIN(a, b)                                                             \
  (minarg1 = (a), minarg2 = (b), (minarg1) < (minarg2) ? (minarg1) : (minarg2))

[[maybe_unused]] static long lmaxarg1, lmaxarg2;
#define LMAX(a, b)                                                             \
  (lmaxarg1 = (a), lmaxarg2 = (b),                                             \
   (lmaxarg1) > (lmaxarg2) ? (lmaxarg1) : (lmaxarg2))

[[maybe_unused]] static long lminarg1, lminarg2;
#define LMIN(a, b)                                                             \
  (lminarg1 = (a), lminarg2 = (b),                                             \
   (lminarg1) < (lminarg2) ? (lminarg1) : (lminarg2))

[[maybe_unused]] static int imaxarg1, imaxarg2;
#define IMAX(a, b)                                                             \
  (imaxarg1 = (a), imaxarg2 = (b),                                             \
   (imaxarg1) > (imaxarg2) ? (imaxarg1) : (imaxarg2))

[[maybe_unused]] static int iminarg1, iminarg2;
#define IMIN(a, b)                                                             \
  (iminarg1 = (a), iminarg2 = (b),                                             \
   (iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))

#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void nrerror(char error_text[]);
REAL* vector(long nl, long nh);
int* ivector(long nl, long nh);
unsigned char* cvector(long nl, long nh);
unsigned long* lvector(long nl, long nh);
REAL* dvector(long nl, long nh);
REAL** matrix(long nrl, long nrh, long ncl, long nch);
REAL** dmatrix(long nrl, long nrh, long ncl, long nch);
int** imatrix(long nrl, long nrh, long ncl, long nch);
REAL** submatrix(REAL** a, long oldrl, long oldrh, long oldcl, long oldch,
                 long newrl, long newcl);
REAL** convert_matrix(REAL* a, long nrl, long nrh, long ncl, long nch);
REAL*** f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_vector(REAL* v, long nl, long nh);
void free_ivector(int* v, long nl, long nh);
void free_cvector(unsigned char* v, long nl, long nh);
void free_lvector(unsigned long* v, long nl, long nh);
void free_dvector(REAL* v, long nl, long nh);
void free_matrix(REAL** m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(REAL** m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int** m, long nrl, long nrh, long ncl, long nch);
void free_submatrix(REAL** b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(REAL** b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(REAL*** t, long nrl, long nrh, long ncl, long nch, long ndl,
                   long ndh);

} // namespace dem ends

#endif /* _NR_UTILS_H_ */
