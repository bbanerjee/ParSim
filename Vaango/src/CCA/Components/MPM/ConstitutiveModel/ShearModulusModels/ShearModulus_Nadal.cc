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

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <CCA/Components/MPM/ConstitutiveModel/ShearModulusModels/ShearModulus_Nadal.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <cmath>
#include <iostream>
#include <sstream>

using namespace Uintah;
using namespace Vaango;

// Construct a shear modulus model.
ShearModulus_Nadal::ShearModulus_Nadal(ProblemSpecP& ps,
                                       MPMEquationOfState* eos)
{
  d_eos = eos;

  ps->require("mu_0", d_mu0);
  ps->require("zeta", d_zeta);
  ps->require("slope_mu_p_over_mu0", d_slope_mu_p_over_mu0);
  ps->require("C", d_C);
  ps->require("m", d_m);
}

// Construct a copy of a shear modulus model.
ShearModulus_Nadal::ShearModulus_Nadal(const ShearModulus_Nadal* smm)
{
  d_eos = smm->d_eos;

  d_mu0                 = smm->d_mu0;
  d_zeta                = smm->d_zeta;
  d_slope_mu_p_over_mu0 = smm->d_slope_mu_p_over_mu0;
  d_C                   = smm->d_C;
  d_m                   = smm->d_m;
}

// Destructor of shear modulus model.
ShearModulus_Nadal::~ShearModulus_Nadal() = default;

void
ShearModulus_Nadal::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP shear_ps = ps->appendChild("shear_modulus_model");
  shear_ps->setAttribute("type", "np_shear");

  shear_ps->appendElement("mu_0", d_mu0);
  shear_ps->appendElement("zeta", d_zeta);
  shear_ps->appendElement("slope_mu_p_over_mu0", d_slope_mu_p_over_mu0);
  shear_ps->appendElement("C", d_C);
  shear_ps->appendElement("m", d_m);
}

// Compute the shear modulus
double
ShearModulus_Nadal::computeInitialShearModulus()
{
  return d_mu0;
}

double
ShearModulus_Nadal::computeShearModulus(const ModelStateBase* state)
{
  return evalShearModulus(state->temperature,
                          state->meltingTemp,
                          state->density,
                          state->initialDensity,
                          state->pressure);
}

double
ShearModulus_Nadal::computeShearModulus(const ModelStateBase* state) const
{
  return evalShearModulus(state->temperature,
                          state->meltingTemp,
                          state->density,
                          state->initialDensity,
                          state->pressure);
}

double
ShearModulus_Nadal::evalShearModulus(double temperature,
                                     double meltingTemp,
                                     double density,
                                     double initialDensity,
                                     double pressure) const
{
  double That = temperature / meltingTemp;
  if (That <= 0)
    return d_mu0;

  double mu = 1.0e-8; // Small value to keep the code from crashing
  if (That > 1.0 + d_zeta)
    return mu;

  double j_denom = d_zeta * (1.0 - That / (1.0 + d_zeta));
  double J       = 1.0 + exp((That - 1.0) / j_denom);
  if (!finite(J))
    return mu;

  double eta = density / initialDensity;
  ASSERT(eta > 0.0);
  eta = pow(eta, 1.0 / 3.0);

  // Pressure is +ve in this calculation
  double P     = -pressure;
  double t1    = d_mu0 * (1.0 + d_slope_mu_p_over_mu0 * P / eta);
  double t2    = 1.0 - That;
  double k_amu = 1.3806503e4 / 1.6605402;
  double t3    = density * k_amu * temperature / (d_C * d_m);
  mu           = 1.0 / J * (t1 * t2 + t3);

  if (mu < 1.0e-8) {
    std::cout << "mu = " << mu << " T = " << temperature
              << " Tm = " << meltingTemp << " T/Tm = " << That << " J = " << J
              << " rho/rho_0 = " << eta << " p = " << P << " t1 = " << t1
              << " t2 = " << t2 << " t3 = " << t3 << "\n";
    mu = 1.0e-8;
  }
  return mu;
}
