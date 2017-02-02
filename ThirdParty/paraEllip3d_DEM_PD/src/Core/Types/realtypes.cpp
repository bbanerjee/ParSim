#include "realtypes.h"

#ifdef QUADMATH

std::ostream& operator << (std::ostream& os, const __float128& number) {
  os << (double) number;
  return os;
}

std::istream& operator >> (std::istream& is, __float128& number) {
  double val;
  is >> val;
  number = (__float128) val;
  return is;
}

REAL fabs(REAL x) {
  return fabsq(x);
}

REAL sqrt(REAL x) {
  return sqrtq(x);
}

REAL fmax(REAL x, REAL y) {
  return fmaxq(x, y);
}

REAL fmin(REAL x, REAL y) {
  return fminq(x, y);
}

REAL pow(REAL x, REAL y) {
  return powq(x,y);
}

REAL nearbyint(REAL x) {
  return nearbyintq(x);
}

REAL sin(REAL x) {
  return sinq(x);
}

REAL asin(REAL x) {
  return asinq(x);
}

REAL cos(REAL x) {
  return cosq(x);
}

REAL acos(REAL x) {
  return acosq(x);
}

REAL tan(REAL x) {
  return tanq(x);
}

REAL atan(REAL x) {
  return atanq(x);
}

#endif
