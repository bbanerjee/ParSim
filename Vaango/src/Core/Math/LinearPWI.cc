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
 *  LinearPWI.cc: linear piecewise interpolation
 *
 *  Written by:
 *   Alexei Samsonov
 *   Department of Computer Science
 *   University of Utah
 *   July 2000
 *
 */

#include <Core/Containers/Array1.h>
#include <Core/Math/LinearPWI.h>

namespace Uintah {

LinearPWI::LinearPWI() {}

LinearPWI::LinearPWI(const Array1<double>& pts, const Array1<double>& vals)
{
  set_data(pts, vals);
}

// takes sorted array of points
bool
LinearPWI::set_data(const Array1<double>& pts, const Array1<double>& vals)
{
  reset();
  if (fill_data(pts) && points.size() > 1) {
    p.resize(points.size());
    for (int i = 0; i < points.size() - 1; i++) {
      p[i].a = (vals[i] * points[i + 1] - vals[i + 1] * points[i]) /
               (points[i + 1] - points[i]);
      p[i].b = (vals[i + 1] - vals[i]) / (points[i + 1] - points[i]);
    }
    return data_valid = true;
  } else {
    return data_valid = false;
  }
}

} // End namespace Uintah
