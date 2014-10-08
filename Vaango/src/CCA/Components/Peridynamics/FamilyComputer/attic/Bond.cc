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

#include <CCA/Components/Peridynamics/FamilyComputer/Bond.h>

#include <iostream>
#define _USE_MATH_DEFINE
#include <cmath>

using namespace Vaango;

Bond::Bond()
  :d_start(0),d_end(0),d_broken(false)
{
}

Bond::Bond(const Uintah::ParticleID& start, const Uintah::ParticleID& end)
  :d_start(start),d_end(end),d_broken(false)
{
}

Bond::~Bond()
{
}

bool 
Bond::operator==(const Bond& bond) const
{
  return ((d_start == bond.d_start && d_end == bond.d_end) ||
          (d_start == bond.d_end && d_end == bond.d_start));
}

namespace Vaango {

  std::ostream& operator<<(std::ostream& out, const Bond& bond)
  {
    out.setf(std::ios::floatfield);
    out.precision(3);
    out << "Bond: [" << bond.d_start << " - " << bond.d_end 
        << "], broken = " << std::boolalpha << bond.d_broken 
        << std::endl;
    return out;
  }
}
