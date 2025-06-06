/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

/*
 *  PiecewiseInterp.h: base class for local family of interpolation techniques
 *
 *  Written by:
 *   Alexei Samsonov
 *   Department of Computer Science
 *   University of Utah
 *   July 2000
 *
 */

#ifndef __CORE_MATH_PIECEWISE_INTERP_H__
#define __CORE_MATH_PIECEWISE_INTERP_H__

#include <Core/Containers/Array1.h>
#include <Core/Math/MiscMath.h>
#include <iostream>

namespace Uintah {

template<class T>
class PiecewiseInterp
{
protected:
  bool data_valid;         // set to true in an interpolation specific cases
  int curr_intrv;          // number of last interval used
  double min_bnd, max_bnd; // and its boundaries (for camera path optimizations)
  Array1<double> points;   // sorted parameter points
  inline bool
  fill_data(const Array1<double>&);

public:
  virtual ~PiecewiseInterp(){};
  inline int
  get_interval(double);
  virtual bool
  get_value(double, T&) = 0;
  virtual bool
  set_data(const Array1<double>& pts, const Array1<T>& vals) = 0;
  void
  reset()
  {
    data_valid = false;
    curr_intrv = -1;
    min_bnd = max_bnd = 0;
  };
};

// returns -1 if data is not valid or the value is out of boundaries
template<class T>
inline int
PiecewiseInterp<T>::get_interval(double w)
{
  if (data_valid) {
    if (w < min_bnd ||
        w >= max_bnd) { // taking advantage of smooth parameter changing
      int lbnd = 0, rbnd = points.size() - 1, delta;

      if (w >= points[lbnd] && w <= points[rbnd]) {

        if (curr_intrv >= 0) {
          if (w < min_bnd) {
            // the series of optimizations that will take advantage
            // on monotonic parameter changing (camera path interp.)
            if (w >= points[curr_intrv - 1]) {
              lbnd = curr_intrv - 1;
            }

            rbnd = curr_intrv;

          } else {
            if (w <= points[curr_intrv + 1]) {
              rbnd = curr_intrv + 1;
            }
            lbnd = curr_intrv;
          }
        }

        while ((delta = rbnd - lbnd) != 1) {
          if (w < points[lbnd + delta / 2]) {
            rbnd = lbnd + delta / 2;
          } else {
            lbnd += delta / 2;
          }
        }

        curr_intrv = lbnd;
        min_bnd    = points[curr_intrv];
        max_bnd    = points[curr_intrv + 1];
      } else {
        curr_intrv = (w < points[lbnd]) ? lbnd : rbnd;
        min_bnd    = 0;
        max_bnd    = 0;
      }
    }
  } else {
    reset();
  }

  return curr_intrv;
}

template<class T>
inline bool
PiecewiseInterp<T>::fill_data(const Array1<double>& pts)
{
  for (int i = 1; i < pts.size(); i++) {
    if ((pts[i] - pts[i - 1]) < 1e-7) {
      return false;
    }
  }

  points.resize(pts.size());
  points = pts;

  return true;
}

} // End namespace Uintah

#endif // __CORE_MATH_PIECEWISE_INTERP_H__
