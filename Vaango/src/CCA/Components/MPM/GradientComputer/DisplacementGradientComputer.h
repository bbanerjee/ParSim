/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2022 Parresia Research Limited, New Zealand
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

#ifndef __DISPLACEMENT_GRADIENT_COMPUTER_H__
#define __DISPLACEMENT_GRADIENT_COMPUTER_H__

#include <CCA/Components/MPM/GradientComputer/GradientComputer.h>


namespace Uintah {

  //////////////////////////////////////////////////////////////////////////
  /*!
    \class DisplacementGradientComputer
    \brief Class for computing displacement gradients
  */
  //////////////////////////////////////////////////////////////////////////

  class DisplacementGradientComputer : public GradientComputer {

  public:
         
    DisplacementGradientComputer(MPMFlags* MFlag);
    DisplacementGradientComputer(const DisplacementGradientComputer* gc);
    virtual ~DisplacementGradientComputer();

    // Make a clone of the gradient computer
    DisplacementGradientComputer* clone();

    // Actually compute displacement gradient
    void computeDispGrad(ParticleInterpolator* interp,
                         const double* oodx,
                         const Point& px,
                         const Matrix3& psize,
                         const Matrix3& pDefGrad_old,
                         constNCVariable<Vector> gDisp,
                         Matrix3& dispGrad_new);

  };
} // End namespace Uintah
      


#endif  // __DISPLACEMENT_GRADIENT_COMPUTER_H__

