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
#include <Core/Grid/MPMInterpolators/FastAxiCPDIInterpolator.h>
#include <Core/Grid/MPMInterpolators/CPDIInterpolator.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Level.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/MiscMath.h>

using namespace Uintah;
using namespace Uintah;

FastAxiCPDIInterpolator::FastAxiCPDIInterpolator()
{
  d_size = 27;
  d_patch = 0;
}

FastAxiCPDIInterpolator::FastAxiCPDIInterpolator(const Patch* patch):FastCPDIInterpolator(patch)
{
  d_size = 27;
  d_patch = patch;
}

FastAxiCPDIInterpolator::~FastAxiCPDIInterpolator()
{
}

std::unique_ptr<ParticleInterpolator> 
FastAxiCPDIInterpolator::clone(const Patch* patch)
{
  return std::make_unique<FastAxiCPDIInterpolator>(patch);
}

void FastAxiCPDIInterpolator::findCellAndWeights(const Point& pos,
					    std::vector<IntVector>& ni,
					    std::vector<double>& S,
					    const Matrix3& size,
                                            const Matrix3& defgrad)
{
  Matrix3 defgrad1=Matrix3(defgrad(0,0),defgrad(0,1),defgrad(0,2),
                           defgrad(1,0),defgrad(1,1),defgrad(1,2),
                           defgrad(2,0),defgrad(2,1),1);
  FastCPDIInterpolator::findCellAndWeights(pos,ni,S,size,defgrad1);
}

void FastAxiCPDIInterpolator::findCellAndShapeDerivatives(const Point& pos,
						     std::vector<IntVector>& ni,
						     std::vector<Vector>& d_S,
						     const Matrix3& size,
                                                     const Matrix3& defgrad)
{
  Matrix3 defgrad1=Matrix3(defgrad(0,0),defgrad(0,1),defgrad(0,2),
                           defgrad(1,0),defgrad(1,1),defgrad(1,2),
                           defgrad(2,0),defgrad(2,1),1);
  FastCPDIInterpolator::findCellAndShapeDerivatives(pos,ni,d_S,size,defgrad1);
}

void FastAxiCPDIInterpolator::findCellAndWeightsAndShapeDerivatives(const Point& pos,
							  std::vector<IntVector>& ni,
							  std::vector<double>& S,
							  std::vector<Vector>& d_S,
							  const Matrix3& size,
                                                          const Matrix3& defgrad)
{
  Matrix3 defgrad1=Matrix3(defgrad(0,0),defgrad(0,1),defgrad(0,2),
                           defgrad(1,0),defgrad(1,1),defgrad(1,2),
                           defgrad(2,0),defgrad(2,1),1);
  FastCPDIInterpolator::findCellAndWeightsAndShapeDerivatives(pos,ni,S,d_S,
                                                          size,defgrad1);
}

int FastAxiCPDIInterpolator::size()
{
  return d_size;
}
