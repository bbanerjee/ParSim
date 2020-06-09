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

#include "YieldCond_Rousselier.h"
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <cmath>
#include <sstream>

using namespace Vaango;

YieldCond_Rousselier::YieldCond_Rousselier(Uintah::ProblemSpecP& ps)
{
  ps->require("D", d_params.D);
  ps->require("sigma_1", d_params.sigma_1);
}

YieldCond_Rousselier::YieldCond_Rousselier(const YieldCond_Rousselier* cm)
{
  d_params.D     = cm->d_params.D;
  d_params.sigma_1 = cm->d_params.sigma_1;
}

void
YieldCond_Rousselier::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP yield_ps = ps->appendChild("plastic_yield_condition");
  yield_ps->setAttribute("type", "rousselier");
  yield_ps->appendElement("D", d_params.D);
  yield_ps->appendElement("sigma_1", d_params.sigma_1);
}

double
YieldCond_Rousselier::evalYieldCondition(const double q,
                                         const double sigma_f,
                                         const double I1,
                                         const double phi,
                                         double& sig)
{
  double D    = d_params.D;
  double sig1 = d_params.sigma_1;
  double p_phi = I1 / (3.0 * (1.0 - phi)) ;
  double q_phi = q / (1.0 - phi);
  double f = q_phi + D * sig1 * phi * std::exp(p_phi / sig1) - sigma_f;
  sig = f + sigma_f;
  return f;
}

/**
 *  Assumes: xi = dev(sigma) - beta, beta = kinematic hardening backstress 
 *  sigma = dev(sigma) + Identity * state->pressure;
 */
double
YieldCond_Rousselier::evalYieldCondition(const Uintah::Matrix3& xi,
                                         const ModelStateBase* state)
{
  Uintah::Matrix3 s_dev = xi;
  if (state->backStress) {
    s_dev = xi + *(state->backStress);
  } 

  double phi = state->porosity;
  double sigma_f = state->yieldStress;
  double I1 = 3.0 * state->pressure; 
  double J2 = 0.5 * s_dev.Contract(s_dev);
  double q = std::sqrt(3.0 * J2);
  double sig = 0.0;

  return evalYieldCondition(q, sigma_f, I1, phi, sig);
}

std::pair<double, Util::YieldStatus>
YieldCond_Rousselier::evalYieldCondition(const ModelStateBase* state)
{
  double phi = state->porosity;
  double sigma_f = state->yieldStress;
  double I1 = state->I1;
  double q = std::sqrt(3.0 * state->J2);
  double sig = 0.0;
  double f = evalYieldCondition(q, sigma_f, I1, phi, sig);
  if (f > 0) {
    return std::make_pair(f, Util::YieldStatus::HAS_YIELDED);
  }
  return std::make_pair(f, Util::YieldStatus::IS_ELASTIC);
}

double
YieldCond_Rousselier::evalYieldConditionMax(const ModelStateBase* state)
{
  double phi = state->porosity;
  double sigma_f = state->yieldStress;
  double p = state->I1 / 3.0;

  double D    = d_params.D;
  double sig1 = d_params.sigma_1;

  double p_phi = p / (1.0 - phi) ;
  double q_max =  (1.0 - phi) * (sigma_f - D * sig1 * phi * std::exp(p_phi / sig1));
  return std::abs(q_max);
}

void
YieldCond_Rousselier::df_dsigma(const Uintah::Matrix3& sig,
                                const double /*sigFlow*/,
                                const double phi,
                                Uintah::Matrix3& dfdsigma)
{
  double p  = sig.Trace() / 3.0;
  Uintah::Matrix3 s_dev = sig - Vaango::Util::Identity * p;
  dfdsigma = df_dsigma(s_dev, p, phi);
}

void
YieldCond_Rousselier::df_dsigma(const Uintah::Matrix3& xi,
                                const ModelStateBase* state,
                                Uintah::Matrix3& dfdsigma)
{
  dfdsigma = df_dsigma(xi, state->pressure, state->porosity);
}

Uintah::Matrix3
YieldCond_Rousselier::df_dsigma(const Uintah::Matrix3& s_dev,
                                double p,
                                double phi) const
{
  double ss = s_dev.NormSquared();
  Uintah::Matrix3 dp_dsigma = Vaango::Util::Identity * (1.0 / 3.0);
  Uintah::Matrix3 dq_dsigma = s_dev / std::sqrt(ss * 2.0 / 3.0);

  double D    = d_params.D;
  double sig1 = d_params.sigma_1;
  double p_phi = p / (1.0 - phi) ;
  double dfdp = (D * phi) / (1.0 - phi) * std::exp(p_phi / sig1);
  double dfdq = 1.0 / (1.0 - phi);

  return (dp_dsigma * dfdp + dq_dsigma * dfdq);
}

/*! \warning Derivative is taken assuming sig_eq^2 - sig_Y(p)^2 = 0 form.
  This is needed for the HypoElasticPlastic algorithm.  Needs
  to be more generalized if possible. */
void
YieldCond_Rousselier::df_dsigmaDev(const Uintah::Matrix3& sig,
                                   const double,
                                   const double,
                                   Matrix3& derivative)
{
  double p  = sig.Trace() / 3.0;
  Uintah::Matrix3 s_dev = sig - Vaango::Util::Identity * p;
  derivative = s_dev * 3.0;
}

/**
 * Assumes that yield condition has q replaced with q_xi 
 * q_xi = sqrt(3/2 (sig_dev - beta) : (sig_dev - beta)) 
 */
void
YieldCond_Rousselier::df_dxi(const Uintah::Matrix3& xi,
                             const ModelStateBase* state,
                             Uintah::Matrix3& f_xi)
{
  double phi = state->porosity;
  double df_dq = 1.0 / (1 - phi);
  Uintah::Matrix3 dq_dxi = xi * std::sqrt(1.5 / xi.NormSquared()); 
  f_xi = dq_dxi * df_dq;
}

void
YieldCond_Rousselier::df_dsigmaDev_dbeta(const Uintah::Matrix3& xi,
                                         const ModelStateBase* state,
                                         Uintah::Matrix3& df_ds,
                                         Uintah::Matrix3& df_dbeta)
{
  df_dxi(xi, state, df_ds);
  df_dbeta = df_ds * (-1.0); 
}

double
YieldCond_Rousselier::df_dp(const ModelStateBase* state)
{
  double D    = d_params.D;
  double sig1 = d_params.sigma_1;
  double phi = state->porosity;
  double p_phi = state->pressure / (1.0 - phi) ;
  double dfdp = (D * phi) / (1.0 - phi) * std::exp(p_phi / sig1);
  return dfdp;
}

double
YieldCond_Rousselier::df_dq(const ModelStateBase* state)
{
  double phi = state->porosity;
  double dfdq = 1.0 / (1.0 - phi);
  return dfdq;
}

double
YieldCond_Rousselier::df_dplasticStrain(const Uintah::Matrix3& xi,
                                        const double& d_sigy_dep,
                                        const ModelStateBase* state)
{
  return -d_sigy_dep;
}

double
YieldCond_Rousselier::df_dporosity(const Uintah::Matrix3& xi,
                                   const ModelStateBase* state)
{
  double D    = d_params.D;
  double sig1 = d_params.sigma_1;

  double phi = state->porosity;
  double overphi = 1.0 / (1.0 - phi);

  double p = state->pressure;
  double exp_term = D * sig1 * std::exp(p * overphi / sig1);

  double ss = xi.NormSquared();
  double q  = std::sqrt(ss * 1.5);

  double df_dphi = - q * overphi * overphi  
                   + exp_term  
                   - phi * exp_term * overphi * overphi;
  return df_dphi;
}

double
YieldCond_Rousselier::d2f_dp_depsVol(const ModelStateBase* state,
                                     const PressureModel* eos,
                                     const ShearModulusModel* shear,
                                     const InternalVariableModel* intvar)
{
  return 0.0;
}

double
YieldCond_Rousselier::d2f_dp_depsDev(const ModelStateBase* state,
                                     const PressureModel* eos,
                                     const ShearModulusModel* shear,
                                     const InternalVariableModel* intvar)
{
  return 0.0;
}

double
YieldCond_Rousselier::d2f_dq_depsVol(const ModelStateBase* state,
                                     const PressureModel* eos,
                                     const ShearModulusModel* shear,
                                     const InternalVariableModel* intvar)
{
  return 0.0;
}

double
YieldCond_Rousselier::d2f_dq_depsDev(const ModelStateBase* state,
                                     const PressureModel* eos,
                                     const ShearModulusModel* shear,
                                     const InternalVariableModel* intvar)
{
  return 0.0;
}

double
YieldCond_Rousselier::df_depsVol(const ModelStateBase* state,
                                 const PressureModel* eos,
                                 const ShearModulusModel* shear,
                                 const InternalVariableModel* intvar)
{
  return 0.0;
}

double
YieldCond_Rousselier::df_depsDev(const ModelStateBase* state,
                                 const PressureModel* eos,
                                 const ShearModulusModel* shear,
                                 const InternalVariableModel* intvar)
{
  return 0.0;
}

/*! Compute h_alpha  where \f$d/dt(ep) = d/dt(gamma)~h_{\alpha}\f$ */
double
YieldCond_Rousselier::eval_h_alpha(const Uintah::Matrix3& xi,
                                   const ModelStateBase* state)
{
  Uintah::Matrix3 dfdsigma;
  df_dsigma(xi, state, dfdsigma);
  double h_alpha = std::sqrt(2.0 / 3.0 * dfdsigma.Contract(dfdsigma));
  return h_alpha; 
}

/*! Compute h_phi  where \f$d/dt(phi) = d/dt(gamma)~h_{\phi}\f$ */
double
YieldCond_Rousselier::eval_h_phi(const Uintah::Matrix3& xi,
                                 const double& voidNucFac,
                                 const ModelStateBase* state)
{
  Uintah::Matrix3 s_dev = xi;
  if (state->backStress) {
    s_dev = xi + *(state->backStress);
  } 
  Uintah::Matrix3 sigma = s_dev + Vaango::Util::Identity * state->pressure;

  Uintah::Matrix3 dfdsigma;
  df_dsigma(xi, state, dfdsigma);

  double sigma_f_sigma = sigma.Contract(dfdsigma);
  return (1.0 - state->porosity) * dfdsigma.Trace() +
         voidNucFac * sigma_f_sigma / 
         ((1.0 - state->porosity) * state->yieldStress);
}


void
YieldCond_Rousselier::computeTangentModulus(const Uintah::TangentModulusTensor& Ce,
                                            const Uintah::Matrix3& f_sigma,
                                            double f_q1,
                                            double f_q2,
                                            double h_q1,
                                            double h_q2,
                                            Uintah::TangentModulusTensor& Cep)
{
  double fqhq = f_q1 * h_q1 + f_q2 * h_q2;
  Uintah::Matrix3 Cr(0.0), rC(0.0);
  double rCr = 0.0;
  for (int ii = 0; ii < 3; ++ii) {
    for (int jj = 0; jj < 3; ++jj) {
      Cr(ii, jj) = 0.0;
      rC(ii, jj) = 0.0;
      for (int kk = 0; kk < 3; ++kk) {
        for (int ll = 0; ll < 3; ++ll) {
          Cr(ii, jj) += Ce(ii, jj, kk, ll) * f_sigma(kk, ll);
          rC(ii, jj) += f_sigma(kk, ll) * Ce(kk, ll, ii, jj);
        }
      }
      rCr += rC(ii, jj) * f_sigma(ii, jj);
    }
  }
  for (int ii = 0; ii < 3; ++ii) {
    for (int jj = 0; jj < 3; ++jj) {
      for (int kk = 0; kk < 3; ++kk) {
        for (int ll = 0; ll < 3; ++ll) {
          Cep(ii, jj, kk, ll) =
            Ce(ii, jj, kk, ll) - Cr(ii, jj) * rC(kk, ll) / (-fqhq + rCr);
        }
      }
    }
  }
}

void
YieldCond_Rousselier::computeElasPlasTangentModulus(
  const Uintah::TangentModulusTensor& Ce,
  const Uintah::Matrix3& sigma,
  double sigY,
  double dsigYdep,
  double porosity,
  double voidNuclFac,
  Uintah::TangentModulusTensor& Cep)
{
  double p = sigma.Trace() / 3.0;
  Uintah::Matrix3 s_dev = sigma - Vaango::Util::Identity * (sigma.Trace() / 3.0);

  ModelStateBase state;
  state.yieldStress = sigY;
  state.porosity = porosity;
  state.pressure = p;

  // Calculate the derivative of the yield function wrt sigma
  Uintah::Matrix3 f_sigma = df_dsigma(s_dev, p, porosity);

  // Calculate derivative wrt porosity
  double f_q1 = df_dporosity(s_dev, &state);

  // Calculate derivative wrt plastic strain
  double f_q2 = df_dplasticStrain(s_dev, dsigYdep, &state);

  // Compute h_q1
  double sigma_f_sigma = sigma.Contract(f_sigma);
  double h_q1 = (1.0 - porosity) * f_sigma.Trace() +
         voidNuclFac * sigma_f_sigma / ((1.0 - porosity) * sigY);

  // Compute h_q2
  double h_q2 = std::sqrt(2.0 / 3.0 * f_sigma.Contract(f_sigma));

  // Calculate elastic-plastic tangent modulus
  computeTangentModulus(Ce, f_sigma, f_q1, f_q2, h_q1, h_q2, Cep);
}
