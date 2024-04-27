/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#ifndef __VELOCITY_GRADIENT_COMPUTER_H__
#define __VELOCITY_GRADIENT_COMPUTER_H__

#include <CCA/Components/MPM/GradientComputer/GradientComputer.h>

namespace Uintah {

//////////////////////////////////////////////////////////////////////////
/*!
  \class VelocityGradientComputer
  \brief Class for computing velocity gradients
*/
//////////////////////////////////////////////////////////////////////////

class VelocityGradientComputer : public GradientComputer
{

public:
  VelocityGradientComputer(const MPMFlags* MFlag);
  VelocityGradientComputer(const VelocityGradientComputer* gc);
  virtual ~VelocityGradientComputer();

  // Make a clone of the gradient computer
  VelocityGradientComputer*
  clone();

  // Actually compute velocity gradient
  void
  computeVelGrad(ParticleInterpolator* interpolator,
                 const double* oodx,
                 const short pgFld[],
                 const Point& px,
                 const Matrix3& pSize,
                 const Matrix3& pDefGrad_old,
                 constNCVariable<Vector>& gVelocity,
                 constNCVariable<Vector>& GVelocity,
                 Matrix3& velGrad_new);

protected:
  void
  computeAxiSymVelocityGradient(Matrix3& velGrad,
                                std::vector<IntVector>& ni,
                                std::vector<Vector>& d_S,
                                std::vector<double>& S,
                                const double* oodx,
                                constNCVariable<Vector>& gVelocity,
                                const Point& px);

  void
  computeVelocityGradient(Matrix3& velGrad,
                          std::vector<IntVector>& ni,
                          std::vector<Vector>& d_S,
                          const double* oodx,
                          const short pgFld[],
                          constNCVariable<Vector>& gVelocity,
                          constNCVariable<Vector>& GVelocity);
};
} // End namespace Uintah

#endif // __VELOCITY_GRADIENT_COMPUTER_H__
