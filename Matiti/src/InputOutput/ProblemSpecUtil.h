/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef __MATITI_PROBLEMSPEC_UTIL_H__
#define __MATITI_PROBLEMSPEC_UTIL_H__

#include <Geometry/Polygon3D.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Geometry/Vector.h>
#include <string>

namespace Matiti_ProblemSpecUtil
{
  // A routine to read in the boundary of a two-dimensional region
  void readVector(Uintah::ProblemSpecP& ps, std::vector<double>& coeffVector);
  void readBoundary(Uintah::ProblemSpecP& ps, Matiti::Polygon3D& boundary);
  void parseVector(const std::string& stringValue, Uintah::Vector& value);
  void checkForInputError(const std::string & stringValue);
}

#endif