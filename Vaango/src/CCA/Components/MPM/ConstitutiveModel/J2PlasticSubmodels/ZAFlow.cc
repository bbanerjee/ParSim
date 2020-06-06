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

#include <CCA/Components/MPM/ConstitutiveModel/J2PlasticSubmodels/ZAFlow.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_Default.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <cmath>

using namespace Uintah;
using Vaango::ModelState_Default;

////////////////////////////////////////////////////////////////////////////////

ZAFlow::ZAFlow(ProblemSpecP& ps)
{
  d_CM.c_0 = 0.0;
  ps->get("c_0", d_CM.c_0);
  if (d_CM.c_0 == 0.0) {
    ps->require("sigma_g", d_CM.sigma_g);
    ps->require("k_H", d_CM.k_H);
    ps->require("sqrt_l_inv", d_CM.sqrt_l);
    d_CM.c_0 = d_CM.sigma_g + d_CM.k_H * d_CM.sqrt_l;
  }
  ps->require("B", d_CM.B);
  ps->require("beta_0", d_CM.beta_0);
  ps->require("beta_1", d_CM.beta_1);
  ps->require("B_0", d_CM.B_0);
  ps->require("alpha_0", d_CM.alpha_0);
  ps->require("alpha_1", d_CM.alpha_1);
  ps->require("K", d_CM.K);
  ps->require("n", d_CM.n);
}

ZAFlow::ZAFlow(const ZAFlow* cm)
{
  d_CM.c_0 = cm->d_CM.c_0;
  d_CM.sigma_g = cm->d_CM.sigma_g;
  d_CM.k_H = cm->d_CM.k_H;
  d_CM.sqrt_l = cm->d_CM.sqrt_l;
  d_CM.B = cm->d_CM.B;
  d_CM.beta_0 = cm->d_CM.beta_0;
  d_CM.beta_1 = cm->d_CM.beta_1;
  d_CM.B_0 = cm->d_CM.B_0;
  d_CM.alpha_0 = cm->d_CM.alpha_0;
  d_CM.alpha_1 = cm->d_CM.alpha_1;
  d_CM.K = cm->d_CM.K;
  d_CM.n = cm->d_CM.n;
}

ZAFlow::~ZAFlow() = default;

void
ZAFlow::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP flow_ps = ps->appendChild("flow_model");
  flow_ps->setAttribute("type", "zerilli_armstrong");

  flow_ps->appendElement("c_0", d_CM.c_0);
  if (d_CM.c_0 == 0.0) {
    flow_ps->appendElement("sigma_g", d_CM.sigma_g);
    flow_ps->appendElement("k_H", d_CM.k_H);
    flow_ps->appendElement("sqrt_l_inv", d_CM.sqrt_l);
  }
  flow_ps->appendElement("B", d_CM.B);
  flow_ps->appendElement("beta_0", d_CM.beta_0);
  flow_ps->appendElement("beta_1", d_CM.beta_1);
  flow_ps->appendElement("B_0", d_CM.B_0);
  flow_ps->appendElement("alpha_0", d_CM.alpha_0);
  flow_ps->appendElement("alpha_1", d_CM.alpha_1);
  flow_ps->appendElement("K", d_CM.K);
  flow_ps->appendElement("n", d_CM.n);
}

double
ZAFlow::computeFlowStress(const ModelStateBase* state_in, const double&,
                          const double&, const MPMMaterial*,
                          const particleIndex idx)
{
  auto state = static_cast<const ModelState_Default*>(state_in);
  double epdot = state->plasticStrainRate;
  double ep = state->plasticStrain;
  double T = state->temperature;
  epdot = (epdot == 0.0) ? 1.0e-8 : epdot;
  ep = (ep < 0.0) ? 0.0 : ep;
  T = (T < 0.0) ? 0.0 : T;

  double sigma_a = d_CM.c_0 + d_CM.K * pow(ep, d_CM.n);
  double alpha = d_CM.alpha_0 - d_CM.alpha_1 * log(epdot);
  double beta = d_CM.beta_0 - d_CM.beta_1 * log(epdot);
  double sigma_y =
    sigma_a + d_CM.B * exp(-beta * T) + d_CM.B_0 * sqrt(ep) * exp(-alpha * T);
  if (std::isnan(sigma_y)) {
    std::cout << "ZA_Flow_Stress:: idx = " << idx << " epdot = " << epdot
         << " ep = " << ep << " T = " << T << "\n";
    std::cout << " idx = " << idx << " sigma_a = " << sigma_a << " alpha = " << alpha
         << " beta = " << beta << " sigma_y = " << sigma_y << "\n";
  }

  return sigma_y;
}

//////////
/*! \brief Calculate the plastic strain rate [epdot(tau,ep,T)] */
//////////
double
ZAFlow::computeEpdot(const ModelStateBase* state_in, const double&,
                     const double& tolerance, const MPMMaterial*,
                     const particleIndex)
{
  auto state = static_cast<const ModelState_Default*>(state_in);
  double tau = state->yieldStress;
  double ep = state->plasticStrain;
  double T = state->temperature;
  ep = (ep < 0.0) ? 0.0 : ep;
  T = (T < 0.0) ? 0.0 : T;

  double sigma_a = d_CM.c_0 + d_CM.K * pow(ep, d_CM.n);
  double B0_sqrtEp = d_CM.B_0 * sqrt(ep);

  // Do Newton iteration
  double epdot = 1.0;
  double f = 0.0;
  double fPrime = 0.0;
  do {
    ASSERT(epdot >= 0.0);
    double alpha = d_CM.alpha_0 - d_CM.alpha_1 * log(epdot);
    double beta = d_CM.beta_0 - d_CM.beta_1 * log(epdot);
    double tAlpha = B0_sqrtEp * exp(-alpha * T);
    double tBeta = d_CM.B * exp(-beta * T);
    double sigma_y = sigma_a + tBeta + tAlpha;
    double term1 = d_CM.beta_1 * T * tBeta;
    double term2 = d_CM.alpha_1 * T * tAlpha;
    f = tau - sigma_y;
    fPrime = -(term1 + term2) / epdot;
    epdot -= f / fPrime;
  } while (std::abs(f) > tolerance);

  return epdot;
}

void
ZAFlow::computeTangentModulus(const Matrix3& stress, const ModelStateBase*,
                              const double&, const MPMMaterial*,
                              const particleIndex, TangentModulusTensor&,
                              TangentModulusTensor&)
{
  throw InternalError("Empty Function: ZAFlow::computeTangentModulus", __FILE__,
                      __LINE__);
}

void
ZAFlow::evalDerivativeWRTScalarVars(const ModelStateBase* state_in,
                                    const particleIndex idx, Vector& derivs)
{
  auto state = static_cast<const ModelState_Default*>(state_in);
  derivs[0] = evalDerivativeWRTStrainRate(state, idx);
  derivs[1] = evalDerivativeWRTTemperature(state, idx);
  derivs[2] = evalDerivativeWRTPlasticStrain(state, idx);
}

double
ZAFlow::evalDerivativeWRTPlasticStrain(const ModelStateBase* state_in,
                                       const particleIndex)
{
  auto state = static_cast<const ModelState_Default*>(state_in);
  // Get the state data
  double ep = state->plasticStrain;
  double epdot = state->plasticStrainRate;
  double T = state->temperature;
  epdot = (epdot == 0.0) ? 1.0e-8 : epdot;
  ep = (ep <= 0.0) ? 1.0e-8 : ep;
  T = (T < 0.0) ? 0.0 : T;

  double alpha = d_CM.alpha_0 - d_CM.alpha_1 * log(epdot);
  double term1 = d_CM.K * d_CM.n * pow(ep, d_CM.n - 1.0);
  double term2 = 0.5 * d_CM.B_0 * exp(-alpha * T) / sqrt(ep);
  double deriv = term1 + term2;
  ;

  return deriv;
}

///////////////////////////////////////////////////////////////////////////
/*  Compute the shear modulus. */
///////////////////////////////////////////////////////////////////////////
double
ZAFlow::computeShearModulus(const ModelStateBase* state_in)
{
  auto state = static_cast<const ModelState_Default*>(state_in);
  return state->shearModulus;
}

///////////////////////////////////////////////////////////////////////////
/* Compute the melting temperature */
///////////////////////////////////////////////////////////////////////////
double
ZAFlow::computeMeltingTemp(const ModelStateBase* state_in)
{
  auto state = static_cast<const ModelState_Default*>(state_in);
  return state->meltingTemp;
}

double
ZAFlow::evalDerivativeWRTTemperature(const ModelStateBase* state_in,
                                     const particleIndex)
{
  auto state = static_cast<const ModelState_Default*>(state_in);
  // Get the state data
  double ep = state->plasticStrain;
  double epdot = state->plasticStrainRate;
  double T = state->temperature;
  epdot = (epdot == 0.0) ? 1.0e-8 : epdot;
  ep = (ep < 0.0) ? 0.0 : ep;
  T = (T < 0.0) ? 0.0 : T;

  double alpha = d_CM.alpha_0 - d_CM.alpha_1 * log(epdot);
  double beta = d_CM.beta_0 - d_CM.beta_1 * log(epdot);
  double term1 = -d_CM.B_0 * alpha * sqrt(ep) * exp(-alpha * T);
  double term2 = -d_CM.B * beta * exp(-beta * T);
  double deriv = term1 + term2;

  return deriv;
}

double
ZAFlow::evalDerivativeWRTStrainRate(const ModelStateBase* state_in,
                                    const particleIndex)
{
  auto state = static_cast<const ModelState_Default*>(state_in);
  // Get the state data
  double ep = state->plasticStrain;
  double epdot = state->plasticStrainRate;
  double T = state->temperature;
  epdot = (epdot == 0.0) ? 1.0e-8 : epdot;
  ep = (ep < 0.0) ? 0.0 : ep;
  T = (T < 0.0) ? 0.0 : T;

  double alpha = d_CM.alpha_0 - d_CM.alpha_1 * log(epdot);
  double beta = d_CM.beta_0 - d_CM.beta_1 * log(epdot);
  double term1 = d_CM.B * d_CM.beta_1 * T * exp(-beta * T);
  double term2 = d_CM.B_0 * sqrt(ep) * d_CM.alpha_1 * T * exp(-alpha * T);
  double deriv = (term1 + term2) / epdot;

  return deriv;
}
