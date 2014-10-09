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

#ifndef __MATITI_DENSITY_H__
#define __MATITI_DENSITY_H__

#include <Pointers/NodeP.h>
#include <Pointers/DensitySP.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <vector>

namespace Matiti {

  
  class Density {

  public:

    Density();
    Density(const Density& den);
    virtual ~Density();

    void clone(const DensitySP den);

    void initialize(Uintah::ProblemSpecP& ps);

    void nodeDensity (const NodeP node, double& node_density);

    double remind (double lengthPeriod, double nodePos);

    inline const double& ringWidth() const { return d_ring_width; }
//    inline void ringWidth(const double& width) { d_ring_width = width; }

  protected:

    

    double density (const std::vector<double>& polyCoeff, double reminder);

  private:

    double d_ring_width;
    std::vector<double> d_poly_coeffs;

 }; // end class

}; // end namespace

#endif

