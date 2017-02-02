#ifndef REALTYPES_H
#define REALTYPES_H

#ifndef QUADMATH
#include <cmath>
typedef double REAL;

#else
extern "C" {
#include <quadmath.h>
}
#include <iostream>

typedef __float128 REAL;

extern std::ostream& operator << (std::ostream& os, const __float128& number);
extern std::istream& operator >> (std::istream& is, __float128& number);

extern REAL fabs(REAL x);
extern REAL sqrt(REAL x);
extern REAL fmax(REAL x, REAL y);
extern REAL fmin(REAL x, REAL y);
extern REAL pow(REAL x, REAL y);
extern REAL nearbyint(REAL x);

extern REAL sin(REAL x);
extern REAL asin(REAL x);
extern REAL cos(REAL x);
extern REAL acos(REAL x);
extern REAL tan(REAL x);
extern REAL atan(REAL x);
#endif

#endif
