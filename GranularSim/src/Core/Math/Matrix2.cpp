/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
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

#include <Core/Math/Matrix2.h>

#include <cstdlib>

#include <iostream>
#include <fstream>

using namespace dem;

void Matrix2::set(int i, int j, double value)
{
  // Assign the Matrix2 the value components
  mat2[i][j] = value;
}


inline void swap(double& v1, double& v2)
{
  double tmp = v1;
  v1 = v2;
  v2 = tmp;
}


void Matrix2::prettyPrint(std::ostream &out_file) const
{
    Matrix2 m2 = *this;
    out_file <<  m2(0,0) << ' ' << m2(0,1) << std::endl;
    out_file <<  m2(1,0) << ' ' << m2(1,1) << std::endl;
}

namespace dem {
  std::ostream &
  operator << (std::ostream &out_file, const Matrix2 &m2)
  {
    // Overload the output stream << operator
    
    out_file <<  m2(0,0) << ' ' << m2(0,1) << ' ';
    out_file <<  m2(1,0) << ' ' << m2(1,1) << ' ';
    
    return out_file;
  }
}
