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

#ifndef FAST_CPDI_INTERPOLATOR_H
#define FAST_CPDI_INTERPOLATOR_H

#include <Core/Grid/MPMInterpolators/ParticleInterpolator.h>

namespace Uintah {

class Patch;

class FastCPDIInterpolator : public ParticleInterpolator
{

public:
  FastCPDIInterpolator();
  FastCPDIInterpolator(const Patch* patch);
  virtual ~FastCPDIInterpolator();

  virtual std::unique_ptr<ParticleInterpolator>
  clone(const Patch*);

  virtual void
  findCellAndWeights(const Point& p,
                     vector<IntVector>& ni,
                     std::vector<double>& S,
                     const Matrix3& size,
                     const Matrix3& defgrad);
  virtual void
  findCellAndShapeDerivatives(const Point& pos,
                              std::vector<IntVector>& ni,
                              std::vector<Vector>& d_S,
                              const Matrix3& size,
                              const Matrix3& defgrad);
  virtual void
  findCellAndWeightsAndShapeDerivatives(const Point& pos,
                                        std::vector<IntVector>& ni,
                                        std::vector<double>& S,
                                        std::vector<Vector>& d_S,
                                        const Matrix3& size,
                                        const Matrix3& defgrad);
  //__________________________________
  //  Needed for AMRMPM
  virtual void
  findCellAndWeights([[maybe_unused]] const Point& pos,
                     [[maybe_unused]] std::vector<IntVector>& ni,
                     [[maybe_unused]] std::vector<double>& S,
                     [[maybe_unused]] constNCVariable<Stencil7>& zoi,
                     [[maybe_unused]] constNCVariable<Stencil7>& zoi_fine,
                     [[maybe_unused]] const bool& getFiner,
                     [[maybe_unused]] int& num_cur,
                     [[maybe_unused]] int& num_fine,
                     [[maybe_unused]] int& num_coarse,
                     [[maybe_unused]] const Vector& size,
                     [[maybe_unused]] bool coarse_part,
                     [[maybe_unused]] const Patch* patch)
  {
  }

  virtual void
  findCellAndWeights_CFI([[maybe_unused]] const Point& pos,
                         [[maybe_unused]] std::vector<IntVector>& ni,
                         [[maybe_unused]] std::vector<double>& S,
                         [[maybe_unused]] constNCVariable<Stencil7>& zoi)
  {
  }

  virtual void
  findCellAndWeightsAndShapeDerivatives_CFI(
    [[maybe_unused]] const Point& pos,
    [[maybe_unused]] std::vector<IntVector>& CFI_ni,
    [[maybe_unused]] std::vector<double>& S,
    [[maybe_unused]] std::vector<Vector>& d_S,
    [[maybe_unused]] constNCVariable<Stencil7>& zoi)
  {
  }
  virtual int
  size();

private:
  const Patch* d_patch;
  int d_size;
};
} // namespace Uintah

#endif
