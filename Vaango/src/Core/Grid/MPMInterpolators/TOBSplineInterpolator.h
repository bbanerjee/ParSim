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

#ifndef TO_BSPLINE_INTERPOLATOR_H
#define TO_BSPLINE_INTERPOLATOR_H

#include <Core/Grid/MPMInterpolators/ParticleInterpolator.h>

namespace Uintah {

class Patch;

class TOBSplineInterpolator : public ParticleInterpolator
{

  // TO = ThirdOrder B Splines

public:
  TOBSplineInterpolator();
  TOBSplineInterpolator(const Patch* patch);
  virtual ~TOBSplineInterpolator();

  virtual std::unique_ptr<ParticleInterpolator>
  clone(const Patch*);

  virtual void
  findCellAndWeights(const Point& p,
                     std::vector<IntVector>& ni,
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

  void
  findNodeComponents(const int& idx,
                     int* xn,
                     int& count,
                     const int& low,
                     const int& hi,
                     const double& cellpos);

  void
  getBSplineWeights(double* Sd,
                    const int* xn,
                    const int& low,
                    const int& hi,
                    const int& count,
                    const double& cellpos);

  void
  getBSplineGrads(double* dSd,
                  const int* xn,
                  const int& low,
                  const int& hi,
                  const int& count,
                  const double& cellpos);

  double
  evalType1BSpline(const double& cp);
  double
  evalType2BSpline(const double& cp);
  double
  evalType3BSpline(const double& cp);

  double
  evalType1BSplineGrad(const double& cp);
  double
  evalType2BSplineGrad(const double& cp);
  double
  evalType3BSplineGrad(const double& cp);

private:
  const Patch* d_patch;
  int d_size;
};
} // namespace Uintah

#endif
