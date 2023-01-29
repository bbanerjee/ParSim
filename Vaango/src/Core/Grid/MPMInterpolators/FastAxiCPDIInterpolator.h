/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
 * Copyright (c) 2018-2023 Parresia Research Limited, NZ
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
#ifndef FAST_AXICPDI_INTERPOLATOR_H
#define FAST_AXICPDI_INTERPOLATOR_H

#include <Core/Grid/MPMInterpolators/ParticleInterpolator.h>
#include <Core/Grid/MPMInterpolators/FastCPDIInterpolator.h>
namespace Uintah {

  class Patch;

  class FastAxiCPDIInterpolator : public FastCPDIInterpolator {
    
  public:
    
    FastAxiCPDIInterpolator();
    FastAxiCPDIInterpolator(const Patch* patch);
    virtual ~FastAxiCPDIInterpolator();
    
    virtual std::unique_ptr<ParticleInterpolator> clone(const Patch*);
    
    virtual void findCellAndWeights(const Point& p,vector<IntVector>& ni, 
                                    std::vector<double>& S, const Matrix3& size,
                                    const Matrix3& defgrad);

    virtual void findCellAndShapeDerivatives(const Point& pos,
                                             std::vector<IntVector>& ni,
                                             std::vector<Vector>& d_S,
                                             const Matrix3& size,
                                             const Matrix3& defgrad);

    virtual void findCellAndWeightsAndShapeDerivatives(const Point& pos,
                                                       std::vector<IntVector>& ni,
                                                       std::vector<double>& S,
                                                       std::vector<Vector>& d_S,
                                                       const Matrix3& size,
                                                       const Matrix3& defgrad);
    virtual int size();
    
  private:
    const Patch* d_patch;
    int d_size;
    
  };
}

#endif
