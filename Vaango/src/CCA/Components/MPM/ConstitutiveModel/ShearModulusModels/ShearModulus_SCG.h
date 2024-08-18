/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef __SCG_SHEAR_MODEL_H__
#define __SCG_SHEAR_MODEL_H__

#include "ShearModulusModel.h"
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

/*! \class ShearModulus_SCG
 *  \brief The shear modulus model used by Steinberg,Cochran,Guinan in
 *         the SCG plasticity model.
 *  \author Biswajit Banerjee,
 *  \author C-SAFE and Department of Mechanical Engineering,
 *  \author University of Utah.
 *
*/
class ShearModulus_SCG : public ShearModulusModel
{

private:
  double d_mu0; // Material constant (also in SCG model)
  double d_A;   // Material constant (also in SCG model)
  double d_B;   // Material constant (also in SCG model)

  ShearModulus_SCG&
  operator=(const ShearModulus_SCG& smm);

public:
  /*! Construct a constant shear modulus model. */
  ShearModulus_SCG(Uintah::ProblemSpecP& ps);

  /*! Construct a copy of constant shear modulus model. */
  ShearModulus_SCG(const ShearModulus_SCG* smm);

  /*! Destructor of constant shear modulus model.   */
  ~ShearModulus_SCG() override;

  void
  outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  /*! Get parameters */
  std::map<std::string, double>
  getParameters() const override
  {
    std::map<std::string, double> params;
    params["mu_0"] = d_mu0;
    params["A"]    = d_A;
    params["B"]    = d_B;
    return params;
  }

  /*! Compute the shear modulus */
  double
  computeInitialShearModulus() override
  {
    return d_mu0;
  }
  double
  computeShearModulus(const ModelStateBase* state) override;
  double
  computeShearModulus(const ModelStateBase* state) const override;

  /*! Compute the shear strain energy */
  double
  computeStrainEnergy(const ModelStateBase* state) override
  {
    return 0.0;
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
  computeQ(const ModelStateBase* state) const override
  {
    return 0.0;
  }

  /////////////////////////////////////////////////////////////////////////
  /*
    Compute dq/depse_s
  */
  /////////////////////////////////////////////////////////////////////////
  double
  computeDqDepse_s(const ModelStateBase* state) const override
  {
    return 0.0;
  }

  /////////////////////////////////////////////////////////////////////////
  /*
    Compute dq/depse_v
  */
  /////////////////////////////////////////////////////////////////////////
  double
  computeDqDepse_v(const ModelStateBase* state) const override
  {
    return 0.0;
  }

private:
  double
  evalShearModulus(double temperature,
                   double density,
                   double initialDensity,
                   double pressure) const;
};
} // End namespace Vaango

#endif // __SCG_SHEAR_MODEL_H__
