/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015 Parresia Research Limited, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/Models/PressureModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ShearModulus_Borja.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <cmath>
#include <iostream>
#include <sstream>

using namespace Uintah;
using namespace Vaango;
using namespace std;

// Construct a shear modulus model.
ShearModulus_Borja::ShearModulus_Borja(ProblemSpecP& ps, PressureModel* eos)
{
  d_eos = eos;

  if (!eos) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Cannot initialize shear mdoulus model unless pressure model has"
        << " been initialized properly";
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  ParameterDict eosParams = d_eos->getParameters();
  d_alpha = eosParams["alpha"];
  d_p0 = eosParams["p0"];
  d_kappatilde = eosParams["kappatilde"];
  d_epse_v0 = eosParams["epse_v0"];

  ps->require("mu0", d_mu0);
}

// Construct a copy of a shear modulus model.
ShearModulus_Borja::ShearModulus_Borja(const ShearModulus_Borja* smm)
{
  d_eos = smm->d_eos;

  d_mu0 = smm->d_mu0;
  d_alpha = smm->d_alpha;
  d_p0 = smm->d_p0;
  d_epse_v0 = smm->d_epse_v0;
  d_kappatilde = smm->d_kappatilde;
}

// Destructor of shear modulus model.
ShearModulus_Borja::~ShearModulus_Borja() = default;

void
ShearModulus_Borja::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP shear_ps = ps->appendChild("elastic_shear_modulus_model");
  shear_ps->setAttribute("type", "borja_shear_modulus");

  shear_ps->appendElement("mu0", d_mu0);
}

// Compute the shear modulus
double
ShearModulus_Borja::computeInitialShearModulus()
{
  double mu_vol = evalShearModulus(0.0);
  return (d_mu0 - mu_vol);
}

double
ShearModulus_Borja::computeShearModulus(const ModelStateBase* state_input)
{
  const ModelState_CamClay* state =
    dynamic_cast<const ModelState_CamClay*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_CamClay.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  double mu_vol = evalShearModulus(state->epse_v);
  return (d_mu0 - mu_vol);
}

double
ShearModulus_Borja::computeShearModulus(const ModelStateBase* state_input) const
{
  const ModelState_CamClay* state =
    dynamic_cast<const ModelState_CamClay*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_CamClay.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  double mu_vol = evalShearModulus(state->epse_v);
  return (d_mu0 - mu_vol);
}

// Compute the shear strain energy
// W = 3/2 mu epse_s^2
double
ShearModulus_Borja::computeStrainEnergy(const ModelStateBase* state_input)
{
  const ModelState_CamClay* state =
    dynamic_cast<const ModelState_CamClay*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_CamClay.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  double mu_vol = evalShearModulus(state->epse_v);
  double W = 1.5 * (d_mu0 - mu_vol) * (state->epse_s * state->epse_s);
  return W;
}

/* Compute q = 3 mu epse_s
         where mu = shear modulus
               epse_s = sqrt{2/3} ||ee||
               ee = deviatoric part of elastic strain = epse - 1/3 epse_v I
               epse = total elastic strain
               epse_v = tr(epse) */
double
ShearModulus_Borja::computeQ(const ModelStateBase* state_input) const
{
  const ModelState_CamClay* state =
    dynamic_cast<const ModelState_CamClay*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_CamClay.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  return evalQ(state->epse_v, state->epse_s);
}

/* Compute dq/depse_s */
double
ShearModulus_Borja::computeDqDepse_s(const ModelStateBase* state_input) const
{
  const ModelState_CamClay* state =
    dynamic_cast<const ModelState_CamClay*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_CamClay.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  return evalDqDepse_s(state->epse_v, state->epse_s);
}

/* Compute dq/depse_v */
double
ShearModulus_Borja::computeDqDepse_v(const ModelStateBase* state_input) const
{
  const ModelState_CamClay* state =
    dynamic_cast<const ModelState_CamClay*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_CamClay.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  return evalDqDepse_v(state->epse_v, state->epse_s);
}

// Private methods below:

//  Shear modulus computation (only pressure contribution)
double
ShearModulus_Borja::evalShearModulus(const double& epse_v) const
{
  double mu_vol = d_alpha * d_p0 * exp(-(epse_v - d_epse_v0) / d_kappatilde);
  return mu_vol;
}

//  Shear stress magnitude computation
double
ShearModulus_Borja::evalQ(const double& epse_v, const double& epse_s) const
{
  double mu_vol = evalShearModulus(epse_v);
  double q = 3.0 * (d_mu0 - mu_vol) * epse_s;

  return q;
}

//  volumetric derivative computation
double
ShearModulus_Borja::evalDqDepse_v(const double& epse_v,
                                  const double& epse_s) const
{
  double mu_vol = evalShearModulus(epse_v);
  double dmu_depse_v = mu_vol / d_kappatilde;
  double dq_depse_v = 3.0 * dmu_depse_v * epse_s;
  return dq_depse_v;
}

//  deviatoric derivative computation
double
ShearModulus_Borja::evalDqDepse_s(const double& epse_v,
                                  const double& epse_s) const
{
  double mu_vol = evalShearModulus(epse_v);
  double dq_depse_s = 3.0 * (d_mu0 - mu_vol);
  return dq_depse_s;
}
