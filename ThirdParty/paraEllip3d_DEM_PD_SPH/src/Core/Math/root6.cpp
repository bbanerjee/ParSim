//    Numerical Recipes in C
//    void zrhqr(float a[], int m, float rtr[], float rti[])
//    Find all the roots of a polynomial with real coefficients, i=0~m, a(i)x^i,
//    given the degree m
//    and the coefficients a[0..m]. The method is to construct an upper
//    Hessenberg matrix whose
//    eigenvalues are the desired roots, and then use the routines balanc and
//    hqr. The real and
//    imaginary parts of the roots are returned in rtr[1..m] and rti[1..m],
//    respectively

//    The purpose of this function is to find the point on surface of particle
//    2, which penetrates
//    deepest into surface of particle 1. However, the algorithm in the thesis
//    can't guarantee that
//    the point found is inside particle 1, which means, even if the two
//    particles are completely
//    separated (not-in-touch), the algorithm could still find a point!
//
//    Input: coef1[0] and coef2[0] are always 1.0, which has been guaranteed by
//    the caller.
//
//    Return values:
//      false - non-overlapped
//      true  - overlapped and vector point returns the deepest penetrated
//      point.

#include <Core/Const/Constants.h>
#include <Core/Math/nr.h>
#include <Core/Math/root6.h>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace dem {

bool
root6(REAL coef1[], REAL coef2[], Vec& point)
{
  return root6(coef1, coef2, point, 1.0, 1, 2);
}

bool
root6(REAL coef1[], REAL coef2[], Vec& point, REAL radius, std::size_t partID1,
      std::size_t partID2)
{

  if ((coef1[0] == 1.0 && coef1[1] == 1.0 && coef1[2] == 1.0 &&
       coef1[3] == 0.0 && coef1[4] == 0.0 && coef1[5] == 0.0) &&
      (coef2[0] == 1.0 && coef2[1] == 1.0 && coef2[2] == 1.0 &&
       coef2[3] == 0.0 && coef2[4] == 0.0 && coef2[5] == 0.0)) {
    REAL x1 = -coef1[6] / 2;
    REAL y1 = -coef1[7] / 2;
    REAL z1 = -coef1[8] / 2;
    REAL R1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1 - coef1[9]);
    REAL x2 = -coef2[6] / 2;
    REAL y2 = -coef2[7] / 2;
    REAL z2 = -coef2[8] / 2;
    REAL R2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2 - coef2[9]);
    Vec dirc = Vec(x1 - x2, y1 - y2, z1 - z2);
    REAL dist = vnormL2(dirc);
    dirc = normalize(dirc);
    point = Vec(x2, y2, z2) + R2 * dirc;
    if (R1 + R2 - dist > EPS)
      return true; // overlapped
    else
      return false; // non-overlapped
  }

  REAL a1 = coef1[0];
  REAL b1 = coef1[1];
  REAL c1 = coef1[2];
  REAL d1 = coef1[3];
  REAL e1 = coef1[4];
  REAL f1 = coef1[5];
  REAL g1 = coef1[6];
  REAL h1 = coef1[7];
  REAL i1 = coef1[8];
  REAL j1 = coef1[9];

  REAL a2 = coef2[0];
  REAL b2 = coef2[1];
  REAL c2 = coef2[2];
  REAL d2 = coef2[3];
  REAL e2 = coef2[4];
  REAL f2 = coef2[5];
  REAL g2 = coef2[6];
  REAL h2 = coef2[7];
  REAL i2 = coef2[8];
  REAL j2 = coef2[9];
  REAL rtc[7], rti[7], rtr[7];

  REAL b1c1 = b1 * c1;
  REAL a1b1c1 = a1 * b1c1;
  REAL a2b1c1 = a2 * b1c1;
  REAL d1e1 = d1 * e1;
  REAL f1g1 = f1 * g1;
  REAL f1g1h1 = f1g1 * h1;
  REAL a1c1 = a1 * c1;
  REAL e2g1 = e2 * g1;
  REAL a1b2 = a1 * b2;
  REAL c1e1 = c1 * e1;
  REAL a1b1 = a1 * b1;
  REAL c2e1 = c2 * e1;
  REAL c1d1 = c1 * d1;
  REAL d2e1 = d2 * e1;
  REAL e1f1 = e1 * f1;
  REAL e1f2 = e1 * f2;
  REAL d1f1 = d1 * f1;
  REAL f2g1 = f2 * g1;
  REAL g1h1 = g1 * h1;
  REAL g2h1 = g2 * h1;
  REAL f1h1 = f1 * h1;
  REAL g1i1 = g1 * i1;
  REAL d1d2 = d1 * d2;
  REAL e1e2 = e1 * e2;
  REAL h1i1 = h1 * i1;
  REAL g2i1 = g2 * i1;
  REAL d1f2 = d1 * f2;
  REAL h2i1 = h2 * i1;
  REAL g1i2 = g1 * i2;
  REAL f1i1 = f1 * i1;
  REAL h1i2 = h1 * i2;
  REAL f1j2 = f1 * j2;
  REAL a1a1 = a1 * a1;
  REAL b1b1 = b1 * b1;
  REAL c1c1 = c1 * c1;
  REAL d1d1 = d1 * d1;
  REAL d1d1d1 = d1d1 * d1;
  REAL d1d1d1d1 = d1d1d1 * d1;
  REAL e1e1 = e1 * e1;
  REAL e1e1e1 = e1e1 * e1;
  REAL e1e1e1e1 = e1e1e1 * e1;
  REAL f1f1 = f1 * f1;
  REAL f1f1f1 = f1f1 * f1;
  REAL f1f1f1f1 = f1f1f1 * f1;
  REAL g1g1 = g1 * g1;
  REAL h1h1 = h1 * h1;
  REAL i1i1 = i1 * i1;
  REAL a2b1 = a2 * b1;
  REAL a2c1 = a2 * c1;
  REAL a1c2 = a1 * c2;
  REAL b1c2 = b1 * c2;
  REAL b2c1 = b2 * c1;
  REAL c1d2 = c1 * d2;
  REAL g1g2 = g1 * g2;
  REAL h1h2 = h1 * h2;
  REAL i1i2 = i1 * i2;
  REAL g1h2 = g1 * h2;
  REAL a1d1 = a1 * d1;
  REAL a1d2 = a1 * d2;
  REAL f1f2 = f1 * f2;
  REAL b1c1d1e1 = b1c1 * d1e1;
  REAL f1g1g2 = f1g1 * g2;
  REAL a1c1d1e1 = a1c1 * d1e1;
  REAL e2g1h1 = e2g1 * h1;
  REAL c2e1f1g1h1 = c2e1 * f1g1h1;
  REAL e2f1g1h1 = e2 * f1g1h1;
  REAL e1f2g1h1 = e1f2 * g1h1;
  REAL f2g1h1 = f2g1 * h1;
  REAL e1f1g2h1 = e1f1 * g2h1;
  REAL e1f1g1 = e1 * f1g1;
  REAL e1f1g1h2 = e1f1g1 * h2;
  REAL d1e1g1i1 = d1e1 * g1i1;
  REAL a1b2c1 = a1b2 * c1;
  REAL a1b1c2 = a1b1 * c2;
  REAL d2e1g1i1 = d2e1 * g1i1;
  REAL e2g1i1 = e2g1 * i1;
  REAL d1e2g1i1 = d1 * e2g1i1;
  REAL f1g1i1 = f1g1 * i1;
  REAL f2g1i1 = f2g1 * i1;
  REAL b1d1e1 = b1 * d1e1;
  REAL d1e1g2i1 = d1e1 * g2i1;
  REAL a1c1d1 = a1c1 * d1;
  REAL d1f1h1i1 = d1f1 * h1i1;
  REAL f1h1i1 = f1h1 * i1;
  REAL d1f2h1i1 = d1f2 * h1i1;
  REAL f2h1i1 = f2 * h1i1;
  REAL d1f1h2i1 = d1f1 * h2i1;
  REAL d1e1g1i2 = d1e1 * g1i2;
  REAL d1f1h1i2 = d1f1 * h1i2;
  REAL a1b1d1e1 = a1b1 * d1e1;
  REAL d1e1f1j2 = d1e1 * f1j2;
  REAL h1i1a1a1 = h1i1 * a1a1;
  REAL h2i1a1a1 = h2i1 * a1a1;
  REAL h1i2a1a1 = h1i2 * a1a1;
  REAL b1c1e1 = b1c1 * e1;
  REAL f1g1i1b1b1 = f1g1i1 * b1b1;
  REAL f2g1i1b1b1 = f2g1i1 * b1b1;
  REAL f1g2i1 = f1 * g2i1;
  REAL f1g1i2 = f1g1 * i2;
  REAL a1d2e1 = a1d2 * e1;

  rtc[0] =
    -8 * b1c1d1e1 * f1g1g2 +
    8 * (a1c1d1e1 * e2g1h1 + a2b1c1 * e1 * f1g1h1 + a1b2 * c1e1 * f1g1h1 +
         a1b1 * c2e1f1g1h1 - 0.5 * c1d1 * d2e1 * f1g1h1 - a1b1c1 * e2f1g1h1 -
         a1b1c1 * e1f2g1h1 + b1c1 * d1f1 * f2g1h1 - a1b1c1 * e1f1g2h1 -
         a1b1c1 * e1f1g1h2 - a1c1d1e1 * f1h1 * h2) +
    8 * (a2b1c1 * d1e1g1i1 + a1b2c1 * d1e1g1i1 + a1b1c2 * d1e1g1i1 -
         a1b1c1 * d2e1g1i1 - a1b1c1 * d1e2g1i1 + b1c1 * d1d2 * f1g1i1 +
         a1b1 * e1e2 * f1g1i1 - 0.5 * b1d1e1 * f1 * f2g1i1 - a1b1c1 * d1e1g2i1 +
         a1c1d1 * d2e1 * h1i1 + a2b1c1 * d1f1h1i1 + a1b2 * c1d1 * f1h1i1 +
         a1b1c2 * d1f1h1i1 - a1b1c1 * d2 * f1h1i1 -
         0.5 * a1 * d1e1 * e2 * f1h1i1 - a1b1c1 * d1f2h1i1 +
         a1b1 * e1f1 * f2h1i1 - a1b1c1 * d1f1h2i1) -
    8 * a1b1c1 * d1e1g1i2 - 8 * a1b1c1 * d1f1h1i2 - 8 * a1b1d1e1 * f1i1 * i2 +
    32 * a1b1c1 * d1e1f1j2 - 16 * b2 * c1e1 * h1i1a1a1 -
    16 * b1 * c2e1 * h1i1a1a1 + 16 * b1c1 * e2 * h1i1a1a1 +
    16 * b1c1e1 * h2i1a1a1 + 16 * b1c1e1 * h1i2a1a1 - 16 * a2c1 * f1g1i1b1b1 -
    16 * a1c2 * f1g1i1b1b1 + 16 * a1c1 * f2g1i1b1b1 +
    16 * a1c1 * f1g2i1 * b1b1 + 16 * a1c1 * f1g1i2 * b1b1 -
    32 * c1 * i1i2 * a1a1 * b1b1 +
    16 * (-a2b1 * d1 * g1h1 - a1b2 * d1 * g1h1 + a1b1 * d2 * g1h1 +
          a1b1 * d1 * g2h1 + a1b1 * d1 * g1h2 - 2 * b1 * h1h2 * a1a1 -
          2 * a1 * g1g2 * b1b1 + 4 * j2 * a1a1 * b1b1) *
      c1c1 +
    2 * (c2e1f1g1h1 - c1 * e2f1g1h1 - c1e1 * f2g1h1 + 3 * c1e1 * f1 * g2h1 +
         3 * c1e1 * f1g1 * h2 - c1 * d2e1g1i1 - 2 * b2c1 * f1g1i1 +
         2 * b1c2 * f1g1i1 - 2 * b1c1 * f2g1i1 - 2 * b1c1 * f1g2i1 -
         2 * a2 * c1e1 * h1i1 + 2 * a1 * c2e1 * h1i1 - 2 * a1c1 * e2 * h1i1 -
         c1d2 * f1h1i1 - 2 * a1c1 * e1 * h2i1 - 2 * b1c1 * f1g1i2 -
         2 * a1c1 * e1 * h1i2 + 8 * a1b1c1 * i1i2 + 4 * b1 * g1g2 * c1c1 +
         2 * d2 * g1h1 * c1c1 + 4 * a1 * h1h2 * c1c1 - 16 * a1b1 * j2 * c1c1) *
      d1d1 +
    2 * (-c2e1 * g1i1 + c1 * e2g1i1 + c1e1 * g2i1 - c2 * f1h1i1 + c1 * f2h1i1 +
         c1 * f1 * h2i1 + c1e1 * g1i2 + c1 * f1h1 * i2 + e1f1 * i1i2 -
         4 * c1e1 * f1j2 - 2 * g2h1 * c1c1 - 2 * g1h2 * c1c1) *
      d1d1d1 -
    2 * c1 * i1i2 * d1d1d1d1 + 4 * j2 * c1c1 * d1d1d1d1 +
    (16 * a1b1c1 * g1g2 + 4 * a2 * c1d1 * g1h1 - 4 * a1c2 * d1 * g1h1 -
     4 * a1c1 * d2 * g1h1 - 2 * a1 * e2f1g1h1 - 4 * a1c1d1 * g2h1 -
     4 * a1c1d1 * g1h2 - 2 * a1d1 * e2g1i1 + 4 * a2b1 * f1g1i1 -
     4 * a1b2 * f1g1i1 - 4 * a1b1 * f2g1i1 - 4 * a1b1 * f1g2i1 +
     2 * a2 * d1f1h1i1 - 2 * a1d2 * f1h1i1 - 2 * a1 * d1f2h1i1 +
     6 * a1 * d1f1h2i1 - 4 * a1b1 * f1g1i2 + 6 * a1 * d1f1h1i2 +
     8 * c1 * h1h2 * a1a1 + 4 * e2 * h1i1a1a1 + 8 * b1 * i1i2 * a1a1 -
     32 * b1c1 * j2 * a1a1 - 2 * c1 * g1g2 * d1d1 + 2 * f2g1i1 * d1d1 -
     2 * f1g2i1 * d1d1 - 2 * f1g1i2 * d1d1 - 2 * a1 * i1i2 * d1d1 +
     8 * a1c1 * j2 * d1d1) *
      e1e1 +
    (2 * d1 * f1g1g2 - 2 * a2 * f1g1h1 + 2 * a1 * f2g1h1 + 2 * a1 * f1 * g2h1 +
     2 * a1 * f1g1 * h2 - 2 * a2 * d1 * g1i1 + 2 * a1d2 * g1i1 +
     2 * a1d1 * g2i1 + 2 * a1d1 * g1i2 - 8 * a1 * d1f1 * j2 - 4 * h2i1a1a1 -
     4 * h1i2a1a1) *
      e1e1e1 -
    2 * a1 * g1g2 * e1e1e1e1 + 4 * j2 * a1a1 * e1e1e1e1 +
    2 * (2 * b2 * c1d1 * g1h1 - 2 * b1c2 * d1 * g1h1 - 2 * b1c1 * d2 * g1h1 -
         b1 * e1f2g1h1 - 2 * b1c1 * d1 * g2h1 - 2 * b1c1 * d1 * g1h2 +
         8 * a1b1c1 * h1h2 + b2 * d1e1g1i1 - b1 * d2e1g1i1 - b1 * d1e2g1i1 +
         3 * b1d1e1 * g2i1 - 2 * a2b1 * e1 * h1i1 + 2 * a1b2 * e1 * h1i1 -
         2 * a1b1 * e2 * h1i1 - b1 * d1f2h1i1 - 2 * a1b1 * e1 * h2i1 +
         3 * b1d1e1 * g1i2 - 2 * a1b1 * e1 * h1i2 + 4 * c1 * g1g2 * b1b1 +
         2 * f2g1i1b1b1 + 4 * a1 * i1i2 * b1b1 - 16 * a1c1 * j2 * b1b1 -
         c1 * h1h2 * d1d1 + e2 * h1i1 * d1d1 - e1 * h2i1 * d1d1 -
         e1 * h1i2 * d1d1 - b1 * i1i2 * d1d1 + 4 * b1c1 * j2 * d1d1 -
         b1 * g1g2 * e1e1 + d2 * g1h1 * e1e1 - d1 * g2h1 * e1e1 -
         d1 * g1h2 * e1e1 - a1 * h1h2 * e1e1 + 4 * a1b1 * j2 * e1e1 +
         2 * j2 * d1d1 * e1e1) *
      f1f1 +
    2 * (-b2 * e1 * g1h1 + b1 * e2g1h1 + b1 * e1 * g2h1 + b1 * e1 * g1h2 +
         d1e1 * h1h2 - b2 * d1 * h1i1 + b1 * d2 * h1i1 + b1 * d1 * h2i1 +
         b1 * d1 * h1i2 - 4 * b1d1e1 * j2 - 2 * g2i1 * b1b1 - 2 * g1i2 * b1b1) *
      f1f1f1 -
    2 * b1 * h1h2 * f1f1f1f1 + 4 * j2 * b1b1 * f1f1f1f1 +
    (-4 * b2c1 * d1e1 * f1 - 4 * b1c2 * d1e1 * f1 + 4 * b1c1 * d2e1 * f1 +
     4 * b1c1 * d1 * e2 * f1 + 4 * b1c1d1e1 * f2 - 8 * c1 * f1f2 * b1b1 -
     8 * b1 * d1d2 * c1c1 + 16 * a2 * b1b1 * c1c1 - 2 * c1e1 * e2 * d1d1 +
     4 * b2 * c1c1 * d1d1 - 8 * a2b1c1 * e1e1 + 2 * c1d1 * d2 * e1e1 +
     d1 * e2 * f1 * e1e1 + 2 * b1 * f1f2 * e1e1 + c2 * d1d1 * e1e1 -
     d2 * f1 * e1e1e1 - d1f2 * e1e1e1 + a2 * e1e1e1e1 - 2 * b1 * e1e2 * f1f1 +
     4 * c2 * b1b1 * f1f1 + b2 * e1e1 * f1f1) *
      g1g1 +
    (-4 * a2c1 * d1e1 * f1 - 4 * a1c2 * d1e1 * f1 + 4 * a1c1 * d2e1 * f1 +
     4 * a1c1d1 * e2 * f1 + 4 * a1c1d1e1 * f2 - 8 * c1e1 * e2 * a1a1 -
     8 * a1 * d1d2 * c1c1 + 16 * b2 * a1a1 * c1c1 - 2 * c1 * f1f2 * d1d1 +
     4 * a2 * c1c1 * d1d1 - 2 * a1 * f1f2 * e1e1 + 4 * c2 * a1a1 * e1e1 -
     8 * a1b2c1 * f1f1 + 2 * c1d1 * d2 * f1f1 + 2 * a1 * e1e2 * f1f1 +
     d1e1 * f2 * f1f1 + c2 * d1d1 * f1f1 + a2 * e1e1 * f1f1 - d2e1 * f1f1f1 -
     d1 * e2 * f1f1f1 + b2 * f1f1f1f1) *
      h1h1 +
    (-4 * a2b1 * d1e1 * f1 - 4 * a1b2 * d1e1 * f1 + 4 * a1b1 * d2e1 * f1 +
     4 * a1b1 * d1 * e2 * f1 + 4 * a1b1d1e1 * f2 - 8 * b1 * e1e2 * a1a1 -
     8 * a1 * f1f2 * b1b1 + 16 * c2 * a1a1 * b1b1 - 8 * a1b1c2 * d1d1 +
     2 * a1 * e1e2 * d1d1 + d2e1 * f1 * d1d1 + 2 * b1 * f1f2 * d1d1 -
     e2 * f1 * d1d1d1 - e1f2 * d1d1d1 + c2 * d1d1d1d1 - 2 * a1 * d1d2 * e1e1 +
     4 * b2 * a1a1 * e1e1 + a2 * d1d1 * e1e1 - 2 * b1 * d1d2 * f1f1 +
     4 * a2 * b1b1 * f1f1 + b2 * d1d1 * f1f1) *
      i1i1;

  REAL term11 = -(d2 * e1 * g1) - d1 * e2g1 + 2 * b2 * f1g1 + 2 * b1 * f2g1;
  REAL term12 =
    -d1e1 * g2 + 2 * b1 * f1 * g2 + 2 * a2 * e1 * h1 + 2 * a1 * e2 * h1;
  REAL term13 = -d2 * f1h1 - d1f2 * h1 + 2 * a1 * e1 * h2 - d1f1 * h2;
  REAL term14 = -4 * a2b1 * i1 - 4 * a1b2 * i1 + 2 * d1d2 * i1;
  REAL term15 = -4 * a1b1 * i2 + i2 * d1d1;
  REAL term10 = term11 + term12 + term13 + term14 + term15;
  REAL term21 = -4 * b2c1 * g1 - 4 * b1c2 * g1 + 2 * e1 * e2g1;
  REAL term22 = -4 * b1c1 * g2 + 2 * c2 * d1 * h1 + 2 * c1d2 * h1 - e2 * f1h1;
  REAL term23 = -e1f2 * h1 + 2 * c1d1 * h2 - e1f1 * h2 - d2e1 * i1;
  REAL term24 = -d1 * e2 * i1 + 2 * b2 * f1i1 + 2 * b1 * f2 * i1 - d1e1 * i2;
  REAL term25 = 2 * b1 * f1 * i2 + g2 * e1e1;
  REAL term20 = term21 + term22 + term23 + term24 + term25;
  REAL term31 = 8 * (a2b1c1 + a1b2 * c1 + a1b1 * c2);
  REAL term32 = -4 * c1d1 * d2 - 4 * a1 * e1e2 + 2 * d2e1 * f1;
  REAL term33 = 2 * d1 * e2 * f1 + 2 * d1e1 * f2 - 4 * b1 * f1f2;
  REAL term34 = -2 * c2 * d1d1 - 2 * a2 * e1e1 - 2 * b2 * f1f1;
  REAL term30 = term31 + term32 + term33 + term34;
  REAL term41 = 2 * c2 * d1 * g1 + 2 * c1d2 * g1 - e2 * f1g1 - e1f2 * g1;
  REAL term42 = 2 * c1d1 * g2 - e1f1 * g2 - 4 * a2c1 * h1 - 4 * a1c2 * h1;
  REAL term43 = 2 * f1f2 * h1 - 4 * a1c1 * h2 + 2 * a2 * e1 * i1;
  REAL term44 = 2 * a1 * e2 * i1 - d2 * f1i1 - d1f2 * i1 + 2 * a1 * e1 * i2;
  REAL term45 = -d1f1 * i2 + h2 * f1f1;
  REAL term40 = term41 + term42 + term43 + term44 + term45;
  rtc[1] =
    (-2 * c2 * d1e1 * g1 + 2 * c1d1 * e2g1 + 4 * b1c2 * f1g1 - e1e2 * f1g1 -
     4 * b1c1 * f2g1 + 4 * a1 * c2e1 * h1 - 4 * a1c1 * e2 * h1 -
     2 * c2 * d1f1 * h1 + 2 * c1d1 * f2 * h1 - e1f1 * f2 * h1 -
     8 * a1b1c2 * i1 + 2 * a1 * e1e2 * i1 - d1 * e2 * f1i1 - d1e1 * f2 * i1 +
     2 * b1 * f1f2 * i1 + 8 * a1b1c1 * i2 + 2 * d1e1 * f1 * i2 +
     2 * c2 * i1 * d1d1 - 2 * c1 * i2 * d1d1 + f2g1 * e1e1 -
     2 * a1 * i2 * e1e1 + e2 * h1 * f1f1 - 2 * b1 * i2 * f1f1) *
      term10 +
    (-8 * a2b1c1 * g1 + 2 * c1d1 * d2 * g1 - d2e1 * f1g1 - d1e1 * f2g1 +
     2 * b1 * f1 * f2g1 + 8 * a1b1c1 * g2 + 2 * d1e1 * f1 * g2 +
     4 * a2 * c1d1 * h1 - 4 * a1c1 * d2 * h1 - 2 * a2 * e1f1 * h1 +
     2 * a1 * e1f2 * h1 - d1f1 * f2 * h1 - 2 * a2 * d1e1 * i1 +
     2 * a1d2e1 * i1 + 4 * a2b1 * f1i1 - d1d2 * f1i1 - 4 * a1b1 * f2 * i1 -
     2 * c1 * g2 * d1d1 + f2 * i1 * d1d1 + 2 * a2 * g1 * e1e1 -
     2 * a1 * g2 * e1e1 + -2 * b1 * g2 * f1f1 + d2 * h1 * f1f1) *
      term20 +
    (-4 * b1c1 * g1g2 + 2 * c1d1 * g2h1 - e1f1g2h1 + 2 * c1d1 * g1h2 -
     e1f1g1h2 - 4 * a1c1 * h1h2 - d1e1g2i1 + 2 * b1 * f1g2i1 +
     2 * a1 * e1 * h2i1 - d1f1h2i1 - d1e1g1i2 + 2 * b1 * f1g1i2 +
     2 * a1 * e1 * h1i2 - d1f1h1i2 - 4 * a1b1 * i1i2 + 16 * a1b1c1 * j2 +
     4 * d1e1f1j2 + i1i2 * d1d1 - 4 * c1 * j2 * d1d1 + g1g2 * e1e1 -
     4 * a1 * j2 * e1e1 + h1h2 * f1f1 - 4 * b1 * j2 * f1f1) *
      term30 +
    (4 * b2 * c1d1 * g1 - 4 * b1c1 * d2 * g1 - d1e1 * e2g1 - 2 * b2 * e1f1g1 +
     2 * b1 * e2 * f1g1 - 8 * a1b2c1 * h1 + 2 * c1d1 * d2 * h1 +
     2 * a1 * e1e2 * h1 - d2e1 * f1h1 - d1 * e2 * f1h1 + 8 * a1b1c1 * h2 +
     2 * d1e1 * f1 * h2 + 4 * a1b2 * e1 * i1 - d1 * d2e1 * i1 -
     4 * a1b1 * e2 * i1 - 2 * b2 * d1f1 * i1 + 2 * b1 * d2 * f1i1 -
     2 * c1 * h2 * d1d1 + e2 * i1 * d1d1 + d2 * g1 * e1e1 - 2 * a1 * h2 * e1e1 +
     2 * b2 * h1 * f1f1 - 2 * b1 * h2 * f1f1) *
      term40;

  REAL d2d2 = d2 * d2;
  REAL e2e2 = e2 * e2;
  REAL f2f2 = f2 * f2;

  REAL term111 = -(d2 * e2g1) + 2 * b2 * f2g1 - d2e1 * g2 - d1 * e2 * g2;
  REAL term112 =
    2 * b2 * f1 * g2 + 2 * b1 * f2 * g2 + 2 * a2 * e2 * h1 - d2 * f2 * h1;
  REAL term113 = 2 * a2 * e1 * h2 + 2 * a1 * e2 * h2 - d2 * f1 * h2 - d1f2 * h2;
  REAL term114 = -4 * a2 * b2 * i1 - 4 * a2b1 * i2 - 4 * a1b2 * i2;
  REAL term115 = 2 * d1d2 * i2 + i1 * d2d2;
  REAL term110 = term111 + term112 + term113 + term114 + term115;
  REAL term121 = -d1e1 * g2 + 2 * b1 * f1 * g2 + 2 * a2 * e1 * h1;
  REAL term122 = 2 * a1 * e2 * h1 - d2 * f1h1 - d1f2 * h1 + 2 * a1 * e1 * h2;
  REAL term123 = -d1f1 * h2 - 4 * a2b1 * i1 - 4 * a1b2 * i1;
  REAL term124 = 2 * d1d2 * i1 + term15;
  REAL term150 = term11 + term121 + term122 + term123 + term124;
  REAL term131 = -4 * b2 * c2 * g1 - 4 * b2c1 * g2 - 4 * b1c2 * g2;
  REAL term132 =
    2 * e1e2 * g2 + 2 * c2 * d2 * h1 - e2 * f2 * h1 + 2 * c2 * d1 * h2;
  REAL term133 = 2 * c1d2 * h2 - e2 * f1 * h2 - e1f2 * h2 - d2 * e2 * i1;
  REAL term134 = 2 * b2 * f2 * i1 - d2e1 * i2 - d1 * e2 * i2 + 2 * b2 * f1 * i2;
  REAL term135 = 2 * b1 * f2 * i2 + g1 * e2e2;
  REAL term130 = term131 + term132 + term133 + term134 + term135;
  REAL term141 = -4 * b1c1 * g2 + 2 * c2 * d1 * h1 + 2 * c1d2 * h1;
  REAL term142 = -e2 * f1h1 - e1f2 * h1 + 2 * c1d1 * h2 - e1f1 * h2;
  REAL term143 = -d2e1 * i1 - d1 * e2 * i1 + 2 * b2 * f1i1 + 2 * b1 * f2 * i1;
  REAL term144 = -d1e1 * i2 + term25;
  REAL term160 = term21 + term141 + term142 + term143 + term144;
  REAL term171 = 8 * a2 * b2c1 + 8 * a2b1 * c2 + 8 * a1b2 * c2;
  REAL term172 = -4 * c2 * d1d2 - 4 * a2 * e1e2 + 2 * d2 * e2 * f1;
  REAL term173 = 2 * d2e1 * f2 + 2 * d1 * e2 * f2 - 4 * b2 * f1f2;
  REAL term174 = -2 * c1 * d2d2 - 2 * a1 * e2e2 - 2 * b1 * f2f2;
  REAL term170 = term171 + term172 + term173 + term174;
  REAL term181 =
    2 * c2 * d2 * g1 - e2 * f2g1 + 2 * c2 * d1 * g2 + 2 * c1d2 * g2;
  REAL term182 = -e2 * f1 * g2 - e1f2 * g2 - 4 * a2 * c2 * h1 - 4 * a2c1 * h2;
  REAL term183 =
    -4 * a1c2 * h2 + 2 * f1f2 * h2 + 2 * a2 * e2 * i1 - d2 * f2 * i1;
  REAL term184 = 2 * a2 * e1 * i2 + 2 * a1 * e2 * i2 - d2 * f1 * i2 - d1f2 * i2;
  REAL h1f2f2 = h1 * f2f2;
  REAL term180 = term181 + term182 + term183 + term184 + h1f2f2;
  rtc[2] =
    (-2 * c2 * d1 * e1 * g1 + 2 * c1 * d1 * e2 * g1 + 4 * b1 * c2 * f1 * g1 -
     e1 * e2 * f1 * g1 - 4 * b1 * c1 * f2 * g1 + 4 * a1 * c2 * e1 * h1 -
     4 * a1 * c1 * e2 * h1 - 2 * c2 * d1 * f1 * h1 + 2 * c1 * d1 * f2 * h1 -
     e1 * f1 * f2 * h1 - 8 * a1 * b1 * c2 * i1 + 2 * a1 * e1 * e2 * i1 -
     d1 * e2 * f1 * i1 - d1e1 * f2 * i1 + 2 * b1 * f1 * f2 * i1 +
     8 * a1b1c1 * i2 + 2 * d1e1 * f1 * i2 + 2 * c2 * i1 * d1d1 -
     2 * c1 * i2 * d1d1 + f2g1 * e1e1 - 2 * a1 * i2 * e1e1 + e2 * h1 * f1f1 -
     2 * b1 * i2 * f1f1) *
      term110 +
    f2 * term150 * term20 +
    (-8 * a2b1c1 * g1 + 2 * c1d1 * d2 * g1 - d2e1 * f1g1 - d1e1 * f2g1 +
     2 * b1 * f1 * f2g1 + 8 * a1b1c1 * g2 + 2 * d1e1 * f1 * g2 +
     4 * a2 * c1d1 * h1 - 4 * a1c1 * d2 * h1 - 2 * a2 * e1f1 * h1 +
     2 * a1 * e1f2 * h1 - d1f1 * f2 * h1 - 2 * a2 * d1e1 * i1 +
     2 * a1d2e1 * i1 + 4 * a2b1 * f1i1 - d1d2 * f1i1 - 4 * a1b1 * f2 * i1 -
     2 * c1 * g2 * d1d1 + f2 * i1 * d1d1 + 2 * a2 * g1 * e1e1 -
     2 * a1 * g2 * e1e1 - 2 * b1 * g2 * f1f1 + d2 * h1 * f1f1) *
      term130 +
    (i2 * term150 + g2 * term160) * term30 +
    (e2 * term150 + d2 * term160 + h2 * term30) * term40 +
    (-4 * b1c1 * g1g2 + 2 * c1d1 * g2h1 - e1f1g2h1 + 2 * c1d1 * g1h2 -
     e1f1g1h2 - 4 * a1c1 * h1h2 - d1e1g2i1 + 2 * b1 * f1g2i1 +
     2 * a1 * e1 * h2i1 - d1f1h2i1 - d1e1g1i2 + 2 * b1 * f1g1i2 +
     2 * a1 * e1 * h1i2 - d1f1h1i2 - 4 * a1b1 * i1i2 + 16 * a1b1c1 * j2 +
     4 * d1e1f1j2 + i1i2 * d1d1 - 4 * c1 * j2 * d1d1 + g1g2 * e1e1 -
     4 * a1 * j2 * e1e1 + h1h2 * f1f1 - 4 * b1 * j2 * f1f1) *
      term170 +
    (4 * b2 * c1d1 * g1 - 4 * b1c1 * d2 * g1 - d1e1 * e2g1 - 2 * b2 * e1f1g1 +
     2 * b1 * e2 * f1g1 - 8 * a1b2c1 * h1 + 2 * c1d1 * d2 * h1 +
     2 * a1 * e1e2 * h1 - d2e1 * f1h1 - d1 * e2 * f1h1 + 8 * a1b1c1 * h2 +
     2 * d1e1 * f1 * h2 + 4 * a1b2 * e1 * i1 - d1 * d2e1 * i1 -
     4 * a1b1 * e2 * i1 - 2 * b2 * d1f1 * i1 + 2 * b1 * d2 * f1i1 -
     2 * c1 * h2 * d1d1 + e2 * i1 * d1d1 + d2 * g1 * e1e1 - 2 * a1 * h2 * e1e1 +
     2 * b2 * h1 * f1f1 - 2 * b1 * h2 * f1f1) *
      term180 +
    c2 * term150 * term150 + a2 * term160 * term160 + j2 * term30 * term30 +
    b2 * term40 * term40;

  REAL term211 =
    -(d2 * e2 * g2) + 2 * b2 * f2 * g2 + 2 * a2 * e2 * h2 - d2 * f2 * h2;
  REAL term212 = -4 * a2 * b2 * i2 + i2 * d2d2;
  REAL term210 = term211 + term212;
  REAL term231 =
    -4 * b2 * c2 * g2 + 2 * c2 * d2 * h2 - e2 * f2 * h2 - d2 * e2 * i2;
  REAL term232 = 2 * b2 * f2 * i2 + g2 * e2e2;
  REAL term230 = term231 + term232;
  REAL term241 = 8 * a2 * b2 * c2 + 2 * d2 * e2 * f2 - 2 * c2 * d2d2;
  REAL term242 = -2 * a2 * e2e2 - 2 * b2 * f2f2;
  REAL term240 = term241 + term242;
  REAL term251 = 2 * c2 * d2 * g2 - e2 * f2 * g2 - 4 * a2 * c2 * h2 +
                 2 * a2 * e2 * i2 - d2 * f2 * i2 + h2 * f2f2;

  rtc[3] =
    2 * c2 * term10 * term110 +
    (-2 * c2 * d1e1 * g1 + 2 * c1d1 * e2g1 + 4 * b1c2 * f1g1 - e1e2 * f1g1 -
     4 * b1c1 * f2g1 + 4 * a1 * c2e1 * h1 - 4 * a1c1 * e2 * h1 -
     2 * c2 * d1f1 * h1 + 2 * c1d1 * f2 * h1 - e1f1 * f2 * h1 -
     8 * a1b1c2 * i1 + 2 * a1 * e1e2 * i1 - d1 * e2 * f1i1 - d1e1 * f2 * i1 +
     2 * b1 * f1f2 * i1 + 8 * a1b1c1 * i2 + 2 * d1e1 * f1 * i2 +
     2 * c2 * i1 * d1d1 - 2 * c1 * i2 * d1d1 + f2g1 * e1e1 -
     2 * a1 * i2 * e1e1 + e2 * h1 * f1f1 - 2 * b1 * i2 * f1f1) *
      term210 +
    f2 * term110 * term20 + (f2 * term150 + 2 * a2 * term20) * term130 +
    (-8 * a2b1c1 * g1 + 2 * c1d1 * d2 * g1 - d2e1 * f1g1 - d1e1 * f2g1 +
     2 * b1 * f1 * f2g1 + 8 * a1b1c1 * g2 + 2 * d1e1 * f1 * g2 +
     4 * a2 * c1d1 * h1 - 4 * a1c1 * d2 * h1 - 2 * a2 * e1f1 * h1 +
     2 * a1 * e1f2 * h1 - d1f1 * f2 * h1 - 2 * a2 * d1e1 * i1 +
     2 * a1d2e1 * i1 + 4 * a2b1 * f1i1 - d1d2 * f1i1 - 4 * a1b1 * f2 * i1 -
     2 * c1 * g2 * d1d1 + f2 * i1 * d1d1 + 2 * a2 * g1 * e1e1 -
     2 * a1 * g2 * e1e1 - 2 * b1 * g2 * f1f1 + d2 * h1 * f1f1) *
      term230 +
    (i2 * term110 + g2 * term130) * term30 +
    (e2 * term110 + d2 * term130) * term40 +
    (i2 * term150 + g2 * term160 + 2 * j2 * term30 + h2 * term40) * term170 +
    (-4 * b1c1 * g1g2 + 2 * c1d1 * g2h1 - e1f1g2h1 + 2 * c1d1 * g1h2 -
     e1f1g1h2 - 4 * a1c1 * h1h2 - d1e1g2i1 + 2 * b1 * f1g2i1 +
     2 * a1 * e1 * h2i1 - d1f1h2i1 - d1e1g1i2 + 2 * b1 * f1g1i2 +
     2 * a1 * e1 * h1i2 - d1f1h1i2 - 4 * a1b1 * i1i2 + 16 * a1b1c1 * j2 +
     4 * d1e1f1j2 + i1i2 * d1d1 - 4 * c1 * j2 * d1d1 + g1g2 * e1e1 -
     4 * a1 * j2 * e1e1 + h1h2 * f1f1 - 4 * b1 * j2 * f1f1) *
      term240 +
    (e2 * term150 + d2 * term160 + h2 * term30 + 2 * b2 * term40) * term180 +
    (4 * b2 * c1d1 * g1 - 4 * b1c1 * d2 * g1 - d1e1 * e2g1 - 2 * b2 * e1f1g1 +
     2 * b1 * e2 * f1g1 - 8 * a1b2c1 * h1 + 2 * c1d1 * d2 * h1 +
     2 * a1 * e1e2 * h1 - d2e1 * f1h1 - d1 * e2 * f1h1 + 8 * a1b1c1 * h2 +
     2 * d1e1 * f1 * h2 + 4 * a1b2 * e1 * i1 - d1 * d2e1 * i1 -
     4 * a1b1 * e2 * i1 - 2 * b2 * d1f1 * i1 + 2 * b1 * d2 * f1i1 -
     2 * c1 * h2 * d1d1 + e2 * i1 * d1d1 + d2 * g1 * e1e1 - 2 * a1 * h2 * e1e1 +
     2 * b2 * h1 * f1f1 - 2 * b1 * h2 * f1f1) *
      term251;

  rtc[4] =
    2 * c2 * term10 * term210 + f2 * term210 * term20 + f2 * term110 * term130 +
    f2 * term150 * term230 + 2 * a2 * term20 * term230 + i2 * term210 * term30 +
    g2 * term230 * term30 + e2 * term210 * term40 + d2 * term230 * term40 +
    i2 * term110 * term170 + g2 * term130 * term170 + i2 * term150 * term240 +
    g2 * term160 * term240 + 2 * j2 * term30 * term240 + h2 * term40 * term240 +
    e2 * term110 * term180 + d2 * term130 * term180 + h2 * term170 * term180 +
    e2 * term150 * term251 + d2 * term160 * term251 + h2 * term30 * term251 +
    2 * b2 * term40 * term251 + c2 * term110 * term110 +
    a2 * term130 * term130 + j2 * term170 * term170 + b2 * term180 * term180;

  rtc[5] = 2 * c2 * term110 * term210 + f2 * term210 * term130 +
           f2 * term110 * term230 + 2 * a2 * term130 * term230 +
           i2 * term210 * term170 + g2 * term230 * term170 +
           i2 * term110 * term240 + g2 * term130 * term240 +
           2 * j2 * term170 * term240 + e2 * term210 * term180 +
           d2 * term230 * term180 + h2 * term240 * term180 +
           e2 * term110 * term251 + d2 * term130 * term251 +
           h2 * term170 * term251 + 2 * b2 * term180 * term251;

  rtc[6] = f2 * term210 * term230 + i2 * term210 * term240 +
           g2 * term230 * term240 + e2 * term210 * term251 +
           d2 * term230 * term251 + h2 * term240 * term251 +
           c2 * term210 * term210 + a2 * term230 * term230 +
           j2 * term240 * term240 + b2 * term251 * term251;

  int order = 0;
  for (int i = 0; i <= 6; i++) {
    if (fabs(rtc[i]) > EPS) {
      order = i;
    }
  }

  for (int i = 0; i <= order; i++) {
    rtc[i] /= rtc[order];
  }

#ifndef NDEBUG // tested: order == 6 whether or not in contact
  debugInf << "root6.cpp: iter=" << iteration << " order=" << order
           << std::endl;
#endif

  // find roots for a polynomial using eigenvalues method
  if (!zrhqr(rtc, order, rtr, rti)) {
    return false;
  }

/*
if ((partID1 == 2 && partID2 == 95) ||
    (partID1 == 95 && partID2 == 2)) {
  //std::cout << "rtc = [";
  for (int ii = 0; ii < 7; ii++) {
    //std::cout << std::setprecision(15) << rtc[ii] << " ";
  }
  //std::cout << "]\n";
}
*/

#ifndef NDEBUG
  for (int i = 1; i <= order; i++)
    debugInf << " rtr=" << std::setw(16) << rtr[i] << " rti=" << std::setw(16)
             << rti[i] << std::endl;
#endif

  REAL lamda[6];
  int nlamda = 0;
  for (int i = 1; i <= order; i++)
    if (fabs(rti[i]) < EPS)
      lamda[nlamda++] = rtr[i];

#ifndef NDEBUG
  debugInf << " nlamda=" << nlamda << std::endl;
#endif

  REAL x, y, z, det, within, itself;
  bool found = false;
  REAL deepest = 0;
  REAL lambda = 0.0, lambdaSq = 0.0, lambdaCu = 0.0;
  /*
  if ((partID1 == 2 && partID2 == 95) ||
      (partID1 == 95 && partID2 == 2)) {
    //std::cout << std::setprecision(15) << "det = [";
  }
  */
  for (int k = 0; k < nlamda; k++) {
    x = 0;
    y = 0;
    z = 0;
    within = 0;
    lambda = lamda[k];
    lambdaSq = lambda * lambda;
    lambdaCu = lambdaSq * lambda;
    det = 8 * a1b1c1 + 2 * d1e1 * f1 - 2 * c1 * d1d1 - 2 * a1 * e1e1 -
          2 * b1 * f1f1 + lambda * term30 + term170 * lambdaSq +
          term240 * lambdaCu;

    // if determinant is zero, there are infinite solutions
    // it is necessary to skip this case, otherwise computation will fail.
    if (fabs(det) > EPS) {
      x = ((-4 * b1c1 * g1 + 2 * c1d1 * h1 - e1f1 * h1 - d1e1 * i1 +
            2 * b1 * f1i1 + g1 * e1e1 + lambda * term20 + term130 * lambdaSq +
            term230 * lambdaCu)) /
          det;
      y = ((2 * c1d1 * g1 - e1f1g1 - 4 * a1c1 * h1 + 2 * a1 * e1 * i1 -
            d1f1 * i1 + h1 * f1f1 + lambda * term40 + term180 * lambdaSq +
            term251 * lambdaCu)) /
          det;
      z = ((-(d1e1 * g1) + 2 * b1 * f1g1 + 2 * a1 * e1 * h1 - d1f1 * h1 -
            4 * a1b1 * i1 + i1 * d1d1 + lambda * term10 + term110 * lambdaSq +
            term210 * lambdaCu)) /
          det;

      within = a1 * x * x + b1 * y * y + c1 * z * z + d1 * x * y + e1 * y * z +
               f1 * z * x + g1 * x + h1 * y + i1 * z + j1;
      itself = a2 * x * x + b2 * y * y + c2 * z * z + d2 * x * y + e2 * y * z +
               f2 * z * x + g2 * x + h2 * y + i2 * z + j2;

/*
if ((partID1 == 2 && partID2 == 95) ||
    (partID1 == 95 && partID2 == 2)) {
  //std::cout << " within = " << within
            << " deepest = " << deepest
            << " itself = " << itself << "\n\t";
}
*/
#ifndef NDEBUG
      debugInf << " k=" << k << " det=" << std::setw(16) << det
               << " x=" << std::setw(16) << x << " y=" << std::setw(16) << y
               << " z=" << std::setw(16) << z << " within=" << std::setw(16)
               << within << " itself=" << std::setw(16) << itself << std::endl;
#endif

      // first, the found point must be on the surface of particle 2, i.e.,
      // itself < EPS
      // second, the point must penetrate deepest into particle 1, i.e.,
      // smallest negative within
      if (fabs(itself) < radius * EPS && within < deepest) {
        deepest = within;
        point = Vec(x, y, z);
        found = true;
#ifndef NDEBUG
        debugInf << " within selected=" << std::setw(16) << deepest
                 << std::endl;
#endif
      }
    } else {
#ifndef NDEBUG
      debugInf << " k=" << k << " det=" << det << " determinant is 0!"
               << std::endl;
#endif
    }
    /*
    if ((partID1 == 2 && partID2 == 95) ||
        (partID1 == 95 && partID2 == 2)) {
      //std::cout << "(" << lambda << ":" << det << ":"
                << x << ":" << y << ":" << z << ")\n\t ";
    }
    */
  }
  /*
  if ((partID1 == 2 && partID2 == 95) ||
      (partID1 == 95 && partID2 == 2)) {
    //std::cout << "]: Found? " << std::boolalpha << found << " \n";
  }
  */
  return found;
}

} // namespace dem ends
