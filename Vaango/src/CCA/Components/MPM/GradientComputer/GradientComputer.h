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

#ifndef __GRADIENT_COMPUTER_H__
#define __GRADIENT_COMPUTER_H__

#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>
#include <vector>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/MPMInterpolators/LinearInterpolator.h>
#include <CCA/Components/MPM/MPMFlags.h>
#include <vector>


namespace Uintah {

  class MPMFlags;
  class ParticleVariableBase;

  //////////////////////////////////////////////////////////////////////////
  /*!
    \class GradientComputer
    \brief Base class for computing gradients
  */
  //////////////////////////////////////////////////////////////////////////

  class GradientComputer {
  public:
         
    GradientComputer(MPMFlags* MFlag);
    GradientComputer(const GradientComputer* gc);
    virtual ~GradientComputer();

    // Make a clone of the gradient computer
    virtual GradientComputer* clone() = 0;

  public:

    /*! Calculate gradient of vector field for 8 noded interpolation, B matrix
        for Kmat and B matrix for Kgeo */
    void computeGradAndBmats(Matrix3& grad,
                             std::vector<IntVector>& ni,
                             std::vector<Vector>& d_S,
                             const double* oodx, 
                             constNCVariable<Vector>& gVec,
                             const Array3<int>& l2g,
                             double B[6][24],
                             double Bnl[3][24],
                             int* dof);

    void computeBmats(vector<IntVector>& ni,
                      std::vector<Vector>& d_S,
                      const double* oodx, 
                      const Array3<int>& l2g,
                      double B[6][24],
                      double Bnl[3][24],
                      int* dof);

  protected:
    
    /*! Calculate gradient of a vector field for 8 noded interpolation */
    void computeGrad(Matrix3& grad,
                     std::vector<IntVector>& ni,
                     std::vector<Vector>& d_S,
                     const double* oodx, 
                     constNCVariable<Vector>& gVec);

    MPMFlags* flag;

  };
} // End namespace Uintah
      

#endif  // __GRADIENT_COMPUTER_H__

