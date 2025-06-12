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

#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCond_CamClay.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <cmath>

using namespace Uintah;
using namespace Vaango;

YieldCond_CamClay::YieldCond_CamClay(Uintah::ProblemSpecP& ps,
                                     IntVar_BorjaPressure* intvar)
{
  d_intvar = intvar;
  ps->require("M", d_M);
}

YieldCond_CamClay::YieldCond_CamClay(const YieldCond_CamClay* yc)
{
  d_intvar = yc->d_intvar;
  d_M = yc->d_M;
}

YieldCond_CamClay::~YieldCond_CamClay() = default;

void
YieldCond_CamClay::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP yield_ps = ps->appendChild("yield_condition");
  yield_ps->setAttribute("type", "camclay");
  yield_ps->appendElement("M", d_M);
}

//--------------------------------------------------------------
// Evaluate yield condition (q = state->q
//                           p = state->p
//                           p_c = state->p_c)
//--------------------------------------------------------------
std::pair<double, Util::YieldStatus>
YieldCond_CamClay::evalYieldCondition(const ModelStateBase* state)
{
  double f = computeYieldFunction(state);
  if (f > 0.0) {
    return std::make_pair(f, Util::YieldStatus::HAS_YIELDED);
  }
  return std::make_pair(f, Util::YieldStatus::IS_ELASTIC);
}

double
YieldCond_CamClay::computeYieldFunction(const ModelStateBase* state_input) const
{
  const ModelState_CamClay* state =
    static_cast<const ModelState_CamClay*>(state_input);

  double p   = state->p;
  double q   = state->q;
  double p_c = state->p_c;
  double f   = q * q / (d_M * d_M) + p * (p - p_c);
  return f;
}

// Evaluate yield condition (s = deviatoric stress
//                           p = state->p
//                           p_c = state->p_c)
double
YieldCond_CamClay::evalYieldCondition(const Uintah::Matrix3&,
                                      const ModelStateBase* state)
{
  return computeYieldFunction(state);
}

//--------------------------------------------------------------
// Evaluate yield condition max (q = state->q
//                               p = state->p
//                               p_c = state->p_c)
//--------------------------------------------------------------
double
YieldCond_CamClay::evalYieldConditionMax(const ModelStateBase* state_input)
{
  const ModelState_CamClay* state =
    static_cast<const ModelState_CamClay*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_CamClay.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  double p_c  = state->p_c;
  double qmax = std::abs(0.5 * d_M * p_c);
  return qmax * qmax / (d_M * d_M);
}

//--------------------------------------------------------------
// Derivatives needed by return algorithms and Newton iterations

//--------------------------------------------------------------
// Compute df/dp  where p = volumetric stress = 1/3 Tr(sigma)
//   df/dp = 2p - p_c
//--------------------------------------------------------------
double
YieldCond_CamClay::df_dp(const ModelStateBase* state_input)
{
  const ModelState_CamClay* state =
    static_cast<const ModelState_CamClay*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_CamClay.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  // std::std::cout << " p = " << state->p << " pc = " << state->p_c << " dfdp =
  // " <<
  // 2*state->p-state->p_c << "\n";
  return (2.0 * state->p - state->p_c);
}

//--------------------------------------------------------------
// Compute df/dq
//   df/dq = 2q/M^2
//--------------------------------------------------------------
double
YieldCond_CamClay::df_dq(const ModelStateBase* state_input)
{
  const ModelState_CamClay* state =
    static_cast<const ModelState_CamClay*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_CamClay.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  return 2.0 * state->q / (d_M * d_M);
}

//--------------------------------------------------------------
// Compute d/depse_v(df/dp)
//   df/dp = 2p(epse_v, epse_s) - p_c(epse_v)
//   d/depse_v(df/dp) = 2dp/depse_v - dp_c/depse_v
//
// Requires:  Equation of state and internal variable
//--------------------------------------------------------------
double
YieldCond_CamClay::d2f_dp_depsVol(const ModelStateBase* state_input,
                                  const MPMEquationOfState* eos,
                                  const ShearModulusModel*)
{
  double dpdepsev = eos->computeDpDepse_v(state_input);
  double dpcdepsev =
    d_intvar->computeVolStrainDerivOfInternalVariable("borja_pressure", state_input);
  return 2.0 * dpdepsev - dpcdepsev;
}

//--------------------------------------------------------------
// Compute d/depse_s(df/dp)
//   df/dp = 2p(epse_v, epse_s) - p_c(epse_v)
//   d/depse_s(df/dp) = 2dp/depse_s
//
// Requires:  Equation of state
//--------------------------------------------------------------
double
YieldCond_CamClay::d2f_dp_depsDev(const ModelStateBase* state_input,
                                  const MPMEquationOfState* eos,
                                  const ShearModulusModel*)
{
  double dpdepses = eos->computeDpDepse_s(state_input);
  return 2.0 * dpdepses;
}

//--------------------------------------------------------------
// Compute d/depse_v(df/dq)
//   df/dq = 2q(epse_v, epse_s)/M^2
//   d/depse_v(df/dq) = 2/M^2 dq/depse_v
//
// Requires:  Shear modulus model
//--------------------------------------------------------------
double
YieldCond_CamClay::d2f_dq_depsVol(const ModelStateBase* state_input,
                                  const MPMEquationOfState*,
                                  const ShearModulusModel* shear)
{
  double dqdepsev = shear->computeDqDepse_v(state_input);
  return (2.0 * dqdepsev) / (d_M * d_M);
}

//--------------------------------------------------------------
// Compute d/depse_s(df/dq)
//   df/dq = 2q(epse_v, epse_s)/M^2
//   d/depse_s(df/dq) = 2/M^2 dq/depse_s
//
// Requires:  Shear modulus model
//--------------------------------------------------------------
double
YieldCond_CamClay::d2f_dq_depsDev(const ModelStateBase* state_input,
                                  const MPMEquationOfState*,
                                  const ShearModulusModel* shear)
{
  double dqdepses = shear->computeDqDepse_s(state_input);
  return (2.0 * dqdepses) / (d_M * d_M);
}

//--------------------------------------------------------------
// Compute df/depse_v
//   df/depse_v = df/dq dq/depse_v + df/dp dp/depse_v - p dp_c/depse_v
//
// Requires:  Equation of state, shear modulus model, internal variable model
//--------------------------------------------------------------
double
YieldCond_CamClay::df_depsVol(const ModelStateBase* state_input,
                              const MPMEquationOfState* eos,
                              const ShearModulusModel* shear)
{
  const ModelState_CamClay* state =
    static_cast<const ModelState_CamClay*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_CamClay.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  double dfdq     = df_dq(state_input);
  double dfdp     = df_dp(state_input);
  double dqdepsev = shear->computeDqDepse_v(state_input);
  double dpdepsev = eos->computeDpDepse_v(state_input);
  double dpcdepsev =
    d_intvar->computeVolStrainDerivOfInternalVariable("borja_pressure", state_input);
  double dfdepsev = dfdq * dqdepsev + dfdp * dpdepsev - state->p * dpcdepsev;

  return dfdepsev;
}

//--------------------------------------------------------------
// Compute df/depse_s
//   df/depse_s = df/dq dq/depse_s + df/dp dp/depse_s
//
// Requires:  Equation of state, shear modulus model
//--------------------------------------------------------------
double
YieldCond_CamClay::df_depsDev(const ModelStateBase* state_input,
                              const MPMEquationOfState* eos,
                              const ShearModulusModel* shear)
{
  double dfdq     = df_dq(state_input);
  double dfdp     = df_dp(state_input);
  double dqdepses = shear->computeDqDepse_s(state_input);
  double dpdepses = eos->computeDpDepse_s(state_input);
  double dfdepses = dfdq * dqdepses + dfdp * dpdepses;

  return dfdepses;
}

//--------------------------------------------------------------
// Other yield condition functions

// Compute df/dsigma
//    df/dsigma = (2p - p_c)/3 I + sqrt(3/2) 2q/M^2 s/||s||
//              = 1/3 df/dp I + sqrt(3/2) df/dq s/||s||
//              = 1/3 df/dp I + df/ds
// where
//    s = sigma - 1/3 tr(sigma) I
Matrix3
YieldCond_CamClay::df_dsigma(const ModelStateBase* state) 
{
  return df_dsigma(Vaango::Util::Identity, state);
}

/*! Derivative with respect to the Cauchy stress (\f$\sigma \f$) */
//   p_c = state->p_c
Uintah::Matrix3
YieldCond_CamClay::df_dsigma(const Matrix3& stress,
                             const ModelStateBase* state_input)
{
  const ModelState_CamClay* state =
    static_cast<const ModelState_CamClay*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_CamClay.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  double p       = stress.Trace() / 3.0;
  double df_dp   = 2.0 * p - state->p_c;
  Matrix3 df_ds  = df_dxi(stress, nullptr);
  Matrix3 df_dsigma = Vaango::Util::Identity * (df_dp / 3.0) + df_ds;
  return df_dsigma;
}

/*! Derivative with respect to the \f$xi\f$ where \f$\xi = s \f$
    where \f$s\f$ is deviatoric part of Cauchy stress */
Uintah::Matrix3
YieldCond_CamClay::df_dxi(const Matrix3& stress,
                          const ModelStateBase*)
{
  Matrix3 sigDev = stress.Deviator();
  double sigDevNorm = sigDev.Norm();
  Matrix3 n         = sigDev / sigDevNorm;
  double q_scaled   = 3.0 * sigDevNorm;
  Matrix3 df_ds     = n * (q_scaled / d_M * d_M);
  return df_ds;
}

/* Derivative with respect to \f$ s \f$ and \f$ \beta \f$ */
std::pair<Uintah::Matrix3, Uintah::Matrix3>
YieldCond_CamClay::df_dsigmaDev_dbeta(const Matrix3& stress,
                                      const ModelStateBase* state)
{
  Matrix3 df_ds = df_dxi(stress, state);
  Matrix3 df_dbeta(0.0);
  return std::make_pair(df_ds, df_dbeta);
}

//--------------------------------------------------------------
// Tangent moduli
void
YieldCond_CamClay::computeElasPlasTangentModulus([[maybe_unused]] const TangentModulusTensor& Ce,
                                                 [[maybe_unused]] const Matrix3& sigma,
                                                 [[maybe_unused]] double sigY,
                                                 [[maybe_unused]] double dsigYdep,
                                                 [[maybe_unused]] double porosity,
                                                 double,
                                                 [[maybe_unused]] TangentModulusTensor& Cep)
{
  std::cout
    << "YieldCond_CamClay: computeElasPlasTangentModulus not implemented yet "
    << "\n";
  return;
}

void
YieldCond_CamClay::computeTangentModulus([[maybe_unused]] const TangentModulusTensor& Ce,
                                         [[maybe_unused]] const Matrix3& f_sigma,
                                         [[maybe_unused]] double f_q1,
                                         [[maybe_unused]] double h_q1,
                                         [[maybe_unused]] TangentModulusTensor& Cep)
{
  std::cout << "YieldCond_CamClay: computeTangentModulus not implemented yet "
            << "\n";
  return;
}
