/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
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

#ifndef __STABILITY_CHECK_H__
#define __STABILITY_CHECK_H__

#include <CCA/Components/MPM/ConstitutiveModel/Utilities/TensorUtils.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/TangentModulusTensor.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {

/*! \class StabilityCheck
  \brief  A generic wrapper for various methods of checking stability.

  \author  Biswajit Banerjee, \n
  C-SAFE and Department of Mechanical Engineering,\n
  University of Utah.\n

  Examples: loss of hyperbolicity/ellipticity, Drucker
  stability criterion, Hill condition etc.
  Provides an abstract base class for various methods of checking
  the stability of motion/bifurcation points
*/
class StabilityCheck
{

public:
  //! Construct an object that can be used to check stability
  StabilityCheck() = default;

  //! Destructor of stability check
  virtual ~StabilityCheck() = default;

  virtual void outputProblemSpec(ProblemSpecP& ps) = 0;

  // Determine if we do the stability.  Instead of checking for the
  // existence of d_stable in ElasticPlastic.cc, instead check
  // do() which is true except for NoneCheck.cc.
  virtual bool doIt() { return true; }

  /*! Check the stability and return the direction of instability
    if any */
  virtual bool checkStability(const Matrix3& cauchyStress,
                              const Matrix3& deformRate,
                              const TangentModulusTensor& tangentModulus,
                              Vector& direction) = 0;

  virtual bool checkStability(const Matrix3& cauchyStress,
                              const Matrix3& deformRate,
                              const Vaango::Tensor::Matrix6Mandel& C_e,
                              const Vaango::Tensor::Vector6Mandel& P_vec,
                              const Vaango::Tensor::Vector6Mandel& N_vec,
                              double H,
                              Vector& direction) = 0;

protected:

  /* Compute A_e = n . C_e . n */
  Matrix3 
  elasticAcousticTensor(const Vaango::Tensor::Matrix6Mandel& C_e,
                        const Vector& n) const;

  /* Compute A = n . C_e . n - n . C_p . n
     C_p = (P otimes (N:C_e))/ (N:P + H) */
  Matrix3 
  elasticPlasticAcousticTensor(const Vaango::Tensor::Matrix6Mandel& C_e,
                               const Vaango::Tensor::Vector6Mandel& P_vec,
                               const Vaango::Tensor::Vector6Mandel& N_vec,
                               double H,
                               const Vector& n) const;

  /* Compute A = n . C_e . n - n . C_p . n, and
             J = det(A) _{mkln} A^{-1}_{lk} */
  std::tuple<Matrix3, double, Matrix3>
  computeAandJ(const Vaango::Tensor::Matrix6Mandel& C_e,
               const Vaango::Tensor::Vector6Mandel& P_vec,
               const Vaango::Tensor::Vector6Mandel& N_vec,
               double H,
               const Vector& n) const;

  /* Compute 6x6 Q matrix for coordinate transformations of second- and fourth-order
     tensors */
  Vaango::Tensor::Matrix6Mandel
  coordTransformMatrix(const Eigen::Matrix3d& eigenvecs) const;

};
} // End namespace Uintah

#endif // __STABILITY_CHECK_H__
