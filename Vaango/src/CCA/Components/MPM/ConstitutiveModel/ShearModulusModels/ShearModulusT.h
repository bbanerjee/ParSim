/*
 * The MIT License
 *
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#ifndef __SHEAR_MODULUS_TEMPLATED_MODEL_H__
#define __SHEAR_MODULUS_TEMPLATED_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateT.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

using ParameterDict = std::map<std::string, double>;

/*! \class ShearModulusT
 *  \brief A generic wrapper for various shear modulus models
 *  \author Biswajit Banerjee,
 *  \author C-SAFE and Department of Mechanical Engineering,
 *  \author University of Utah.
 *
 * Provides an abstract base class for various shear modulus models
*/
template <typename DerivedT, typename StateT, typename EquationOfStateT>
class ShearModulusT
{
public:
  ~ShearModulusT() = default;

  void
  outputProblemSpec(Uintah::ProblemSpecP& ps)
  {
    derived()->l_outputProblemSpec(ps);
  }

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Get the pressure model
   */
  /////////////////////////////////////////////////////////////////////////
  EquationOfStateT*
  getPressureModel() const
  {
    return d_eos;
  }

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Get the model parameters
   */
  /////////////////////////////////////////////////////////////////////////
  std::map<std::string, double>
  getParameters() const
  {
    return derived()->l_getParameters();
  }

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Compute the shear modulus
  */
  /////////////////////////////////////////////////////////////////////////
  double
  computeInitialShearModulus()
  {
    return derived()->l_computeInitialShearModulus();
  }

  double
  computeShearModulus(const StateT* state)
  {
    return derived()->l_computeShearModulus(state);
  }

  double
  computeShearModulus(const StateT* state) const
  {
    return derived()->l_computeShearModulus(state);
  }

  /*! Compute the shear strain energy */
  double
  computeStrainEnergy(const StateT* state)
  {
    return derived()->l_computeStrainEnergy(state);
  }

  /////////////////////////////////////////////////////////////////////////
  /*
    Compute q = 3 mu epse_s
       where mu = shear modulus
             epse_s = sqrt{2/3} ||ee||
             ee = deviatoric part of elastic strain = epse - 1/3 epse_v I
             epse = total elastic strain
             epse_v = tr(epse)
  */
  /////////////////////////////////////////////////////////////////////////
  double
  computeQ(const StateT* state) const
  {
    return derived()->l_computeQ(state);
  }

  /////////////////////////////////////////////////////////////////////////
  /*
    Compute dq/depse_s
  */
  /////////////////////////////////////////////////////////////////////////
  double
  computeDqDepse_s(const StateT* state) const
  {
    return derived()->l_computeDqDepse_s(state);
  }

  /////////////////////////////////////////////////////////////////////////
  /*
    Compute dq/depse_v
  */
  /////////////////////////////////////////////////////////////////////////
  double
  computeDqDepse_v(const StateT* state) const
  {
    return derived()->l_computeDqDepse_v(state);
  }

protected:
  double d_shearModulus;   // the initial shear modulus
  EquationOfStateT* d_eos; // the associated Pressure EOS model

private:
  ShearModulusT() { d_shearModulus = 0.0; }

  DerivedT*
  derived()
  {
    return static_cast<DerivedT*>(this);
  }

  const DerivedT*
  derived() const
  {
    return static_cast<const DerivedT*>(this);
  }

  friend DerivedT;
};
} // End namespace Uintah

#endif // __SHEAR_MODULUS_TEMPLATED_MODEL_H__
