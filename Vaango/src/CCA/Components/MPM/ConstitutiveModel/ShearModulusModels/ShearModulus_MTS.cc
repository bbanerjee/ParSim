/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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

#include "ShearModulus_MTS.h"
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Default.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <cmath>
#include <iostream>
#include <sstream>

using namespace Vaango;

// Construct a shear modulus model.
ShearModulus_MTS::ShearModulus_MTS(Uintah::ProblemSpecP& ps)
{
  ps->require("mu_0", d_mu0);
  ps->require("D", d_D);
  ps->require("T_0", d_T0);
}

// Construct a copy of a shear modulus model.
ShearModulus_MTS::ShearModulus_MTS(const ShearModulus_MTS* smm)
{
  d_mu0 = smm->d_mu0;
  d_D = smm->d_D;
  d_T0 = smm->d_T0;
}

// Destructor of shear modulus model.
ShearModulus_MTS::~ShearModulus_MTS() = default;

void
ShearModulus_MTS::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP shear_ps = ps->appendChild("shear_modulus_model");
  shear_ps->setAttribute("type", "mts_shear");

  shear_ps->appendElement("mu_0", d_mu0);
  shear_ps->appendElement("D", d_D);
  shear_ps->appendElement("T_0", d_T0);
}

// Compute the shear modulus
double
ShearModulus_MTS::computeShearModulus(const ModelStateBase* state_in)
{
  auto state = static_cast<const ModelState_Default*>(state_in);
  double T = state->temperature;
  ASSERT(T > 0.0);
  return evalShearModulus(T);
}

// Compute the shear modulus
double
ShearModulus_MTS::computeShearModulus(const ModelStateBase* state_in) const
{
  auto state = static_cast<const ModelState_Default*>(state_in);
  double T = state->temperature;
  ASSERT(T > 0.0);
  return evalShearModulus(T);
}

double
ShearModulus_MTS::evalShearModulus(double T) const
{
  double expT0_T = exp(d_T0 / T) - 1.0;
  ASSERT(expT0_T != 0);
  double mu = d_mu0 - d_D / expT0_T;
  if (!(mu > 0.0)) {
    std::ostringstream desc;
    desc << "**Compute MTS Shear Modulus ERROR** Shear modulus <= 0." << "\n";
    desc << "T = " << T << " mu0 = " << d_mu0 << " T0 = " << d_T0
         << " exp(To/T) = " << expT0_T << " D = " << d_D << "\n";
    throw Uintah::InvalidValue(desc.str(), __FILE__, __LINE__);
  }
  return mu;
}
