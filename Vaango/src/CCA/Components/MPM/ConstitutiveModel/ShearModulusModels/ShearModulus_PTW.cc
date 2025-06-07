/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#include "ShearModulus_PTW.h"
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <cmath>
#include <iostream>
#include <sstream>

using namespace Vaango;

// Construct a shear modulus model.
ShearModulus_PTW::ShearModulus_PTW(Uintah::ProblemSpecP& ps)
{
  ps->require("mu_0", d_mu0);
  ps->require("alpha", d_alpha);
  ps->require("alphap", d_alphap);
  ps->require("slope_mu_p_over_mu0", d_slope_mu_p_over_mu0);
}

// Construct a copy of a shear modulus model.
ShearModulus_PTW::ShearModulus_PTW(const ShearModulus_PTW* smm)
{
  d_mu0                 = smm->d_mu0;
  d_alpha               = smm->d_alpha;
  d_alphap              = smm->d_alphap;
  d_slope_mu_p_over_mu0 = smm->d_slope_mu_p_over_mu0;
}

// Destructor of shear modulus model.
ShearModulus_PTW::~ShearModulus_PTW() = default;

void
ShearModulus_PTW::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP shear_ps = ps->appendChild("shear_modulus_model");
  shear_ps->setAttribute("type", "ptw_shear");

  shear_ps->appendElement("mu_0", d_mu0);
  shear_ps->appendElement("alpha", d_alpha);
  shear_ps->appendElement("alphap", d_alphap);
  shear_ps->appendElement("slope_mu_p_over_mu0", d_slope_mu_p_over_mu0);
}

// Compute the shear modulus
double
ShearModulus_PTW::computeShearModulus(const ModelStateBase* state)
{
  return evalShearModulus(state->temperature,
                          state->meltingTemp,
                          state->density,
                          state->initialDensity,
                          state->pressure);
}

double
ShearModulus_PTW::computeShearModulus(const ModelStateBase* state) const
{
  return evalShearModulus(state->temperature,
                          state->meltingTemp,
                          state->density,
                          state->initialDensity,
                          state->pressure);
}

double
ShearModulus_PTW::evalShearModulus(double temperature,
                                   double meltingTemp,
                                   double density,
                                   double initialDensity,
                                   double pressure) const
{
  double eta = density / initialDensity;
  ASSERT(eta > 0.0);
  eta         = pow(eta, 1.0 / 3.0);
  double That = temperature / meltingTemp;
  double P    = -pressure;
  double mu0P = d_mu0 * (1.0 + d_slope_mu_p_over_mu0 * P / eta);
  double mu   = mu0P * (1.0 - d_alphap * That);
  return mu;
}
