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

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCond_Gurson.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <cmath>
#include <iostream>

using namespace Uintah;
using namespace Vaango;

YieldCond_Gurson::YieldCond_Gurson(Uintah::ProblemSpecP& ps,
                                   IntVar_Metal* intvar,
                                   const FlowStressModel* flow)
{
  d_intvar = intvar;
  d_flow = flow;

  ps->require("q1", d_CM.q1);
  ps->require("q2", d_CM.q2);
  ps->require("q3", d_CM.q3);
  ps->require("k", d_CM.k);
  ps->require("f_c", d_CM.f_c);
}

YieldCond_Gurson::YieldCond_Gurson(const YieldCond_Gurson* cm)
{
  d_intvar = cm->d_intvar;
  d_flow = cm->d_flow;

  d_CM.q1  = cm->d_CM.q1;
  d_CM.q2  = cm->d_CM.q2;
  d_CM.q3  = cm->d_CM.q3;
  d_CM.k   = cm->d_CM.k;
  d_CM.f_c = cm->d_CM.f_c;
}

YieldCond_Gurson::~YieldCond_Gurson() = default;

void
YieldCond_Gurson::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP yield_ps = ps->appendChild("yield_condition");
  yield_ps->setAttribute("type", "gurson");
  yield_ps->appendElement("q1", d_CM.q1);
  yield_ps->appendElement("q2", d_CM.q2);
  yield_ps->appendElement("q3", d_CM.q3);
  yield_ps->appendElement("k", d_CM.k);
  yield_ps->appendElement("f_c", d_CM.f_c);
}

std::pair<double, Util::YieldStatus>
YieldCond_Gurson::evalYieldCondition(const ModelStateBase* state)
{
  double f = computeYieldFunction(state);
  if (f > 0) {
    return std::make_pair(f, Util::YieldStatus::HAS_YIELDED);
  }
  return std::make_pair(f, Util::YieldStatus::IS_ELASTIC);
}

double
YieldCond_Gurson::computeYieldFunction(const ModelStateBase* state) const
{
  return computeYieldFunction(state->getStress(), state);
}

double
YieldCond_Gurson::evalYieldCondition(const Matrix3& stress,
                                     const ModelStateBase* state)
{
  return computeYieldFunction(stress, state);
}

double
YieldCond_Gurson::computeYieldFunction(const Matrix3& stress,
                                       const ModelStateBase* state) const
{
  Matrix3 xi   = stress.Deviator() - state->backStress.Deviator();
  double sig_m = (stress.Trace() - state->backStress.Trace()) / 3.0;

  // Get the state data
  double sigy     = state->yieldStress;
  double porosity = state->porosity;
  double xi_xi    = xi.NormSquared();
  double sig_eq   = std::sqrt(1.5 * xi_xi);

  ASSERT(sigy != 0);

  double q1    = d_CM.q1;
  double q2    = d_CM.q2;
  double q3    = d_CM.q3;
  double k     = d_CM.k;
  double phi_c = d_CM.f_c;

  double phi_star = porosity;
  if (porosity > phi_c) {
    phi_star = phi_c + k * (porosity - phi_c);
  }

  double maxincosh = 30.0;
  double incosh = 1.5 * q2 * sig_m / sigy;
  if (std::abs(incosh) > 30.0) {
    incosh = copysign(maxincosh, incosh);
  }
  double a = 2.0 * q1 * phi_star * cosh(incosh);
  double b = 1.0 + q3 * phi_star * phi_star;
  double f = (sig_eq * sig_eq) / (sigy * sigy) + a - b;
  return f;
}

double
YieldCond_Gurson::evalYieldConditionMax(const ModelStateBase* state)
{
  return state->yieldStress * state->yieldStress;
}

/*! Derivative with respect to the Cauchy stress (\f$\sigma \f$)*/
Matrix3
YieldCond_Gurson::df_dsigma(const ModelStateBase* state) 
{
  return df_dsigma(state->getStress(), state);
}

Uintah::Matrix3
YieldCond_Gurson::df_dsigma(const Matrix3& stress,
                            const ModelStateBase* state)
{
  Matrix3 xi   = stress.Deviator() - state->backStress.Deviator();
  double sig_m = (stress.Trace() - state->backStress.Trace()) / 3.0;

  double sigy = state->yieldStress;
  ASSERT(sigy != 0);

  double porosity = state->porosity;

  double q1    = d_CM.q1;
  double q2    = d_CM.q2;
  double k     = d_CM.k;
  double phi_c = d_CM.f_c;

  double phi_star = porosity;
  if (porosity > phi_c) {
    phi_star = phi_c + k * (porosity - phi_c);
  }

  double maxinsinh = 30.0;
  double insinh    = 1.5 * q2 * sig_m / sigy;
  if (std::abs(insinh) > 30.0) {
    insinh  = copysign(maxinsinh, insinh);
  }
  double a  = 3.0 / (sigy * sigy);
  double b  = (q1 * q2 * phi_star / sigy) * sinh(insinh);
  Matrix3 df_dsigma = xi * a + Vaango::Util::Identity * b;
  return df_dsigma;
}

/*! Derivative with respect to the \f$xi\f$ where \f$\xi = s - \beta \f$
    where \f$s\f$ is deviatoric part of Cauchy stress and
    \f$\beta\f$ is the backstress */
Uintah::Matrix3
YieldCond_Gurson::df_dxi(const Matrix3& stress,
                         const ModelStateBase* state)
{
  auto xi = stress.Deviator() - state->backStress.Deviator();

  double sigy = state->yieldStress;
  ASSERT(sigy != 0);
  double a = 3.0 / (sigy * sigy);
  Matrix3 df_dxi   = xi * a;
  return df_dxi;
}

/* Derivative with respect to \f$ s \f$ and \f$ \beta \f$ */
std::pair<Uintah::Matrix3, Uintah::Matrix3>
YieldCond_Gurson::df_dsigmaDev_dbeta(const Matrix3& stress,
                                     const ModelStateBase* state)
{
  Matrix3 df_ds = df_dxi(stress, state);
  double f_beta_p = df_dbeta_p(stress, state);
  Matrix3 df_dbeta = df_ds * (-1.0) + Vaango::Util::Identity * (-f_beta_p);
  return std::make_pair(df_ds, df_dbeta);
}

double
YieldCond_Gurson::df_dbeta_p(const Matrix3& stress,
                             const ModelStateBase* state) const
{
  double sig_m = (stress.Trace() - state->backStress.Trace()) / 3.0;

  double sigy = state->yieldStress;
  ASSERT(sigy != 0);

  double phi   = state->porosity;

  double q1    = d_CM.q1;
  double q2    = d_CM.q2;
  double k     = d_CM.k;
  double phi_c = d_CM.f_c;

  double phi_star = phi;
  if (phi > phi_c) {
    phi_star = phi_c + k * (phi - phi_c);
  }

  double maxinsinh = 30.0;
  double insinh    = 1.5 * q2 * sig_m / sigy;
  if (std::abs(insinh) > 30.0) {
    insinh  = copysign(maxinsinh, insinh);
  }
  double f_beta_p = (q1 * q2 * phi_star / sigy) * sinh(insinh);
  return f_beta_p;
}

/*! Derivative with respect to internal variables */
void
YieldCond_Gurson::df_dintvar(const ModelStateBase* state,
                             MetalIntVar& df_dintvar) const
{
  double f_plasticStrain = df_dplasticStrain(state);
  double f_porosity = df_dporosity(state);
  df_dintvar = {f_plasticStrain, f_porosity};
}

/*! Derivative with respect to the equivalent plastic strain (\f$\epsilon^p \f$)
    The calling routine has to send in the correct stress state so that
    q = sqrt(3/2 ||s - beta||^2) , s = dev. stress, beta = backstress */
double
YieldCond_Gurson::df_dplasticStrain(const ModelStateBase* state) const
{
  double sigy = state->yieldStress;
  ASSERT(sigy != 0);

  double dsigy_dep = d_flow->evalDerivativeWRTPlasticStrain(state, 0);

  Matrix3 xi = state->devStress - state->backStress.Deviator();

  double sig_m    = state->p - state->backStress.Trace()/3.0;
  double sig_eq   = std::sqrt(1.5 * xi.Contract(xi));
  double porosity = state->porosity;

  double q1    = d_CM.q1;
  double q2    = d_CM.q2;
  double k     = d_CM.k;
  double phi_c = d_CM.f_c;

  double phi_star = porosity;
  if (porosity > phi_c) {
    phi_star = phi_c + k * (porosity - phi_c);
  }

  double maxinsinh = 30.0;
  double insinh    = 1.5 * q2 * sig_m / sigy;
  if (std::abs(insinh) > 30.0) {
    insinh = std::copysign(maxinsinh, insinh);
  }
  double a = (2.0 * sig_eq * sig_eq) / (sigy * sigy * sigy);
  double b = 3.0 * q1 * q2 * sig_m * phi_star * sinh(insinh) / (sigy * sigy);

  return -(a + b) * dsigy_dep;
}

/*! Derivative with respect to the porosity (\f$\phi \f$)*/
double
YieldCond_Gurson::df_dporosity(const ModelStateBase* state) const
{
  double sigy = state->yieldStress;
  ASSERT(sigy != 0);

  double sig_m    = state->p - state->backStress.Trace()/3.0;
  double porosity = state->porosity;

  double q1    = d_CM.q1;
  double q2    = d_CM.q2;
  double q3    = d_CM.q3;
  double k     = d_CM.k;
  double phi_c = d_CM.f_c;

  double phi_star       = porosity;
  double dphiStar_dphi = 1.0;
  if (porosity > phi_c) {
    phi_star      = phi_c + k * (porosity - phi_c);
    dphiStar_dphi = k;
  }

  double maxincosh = 30.0;
  double incosh    = 1.5 * q2 * sig_m / sigy;
  if (std::abs(incosh) > maxincosh) {
    incosh = std::copysign(maxincosh, incosh);
  }
  double a = 2.0 * q1 * cosh(incosh);
  double b = 2.0 * q3 * phi_star;

  return (a - b) * dphiStar_dphi;
}

void
YieldCond_Gurson::computeTangentModulus(const TangentModulusTensor& Ce,
                                        const Matrix3& f_sigma,
                                        double f_q1,
                                        double f_q2,
                                        double h_q1,
                                        double h_q2,
                                        TangentModulusTensor& Cep)
{
  // std::cerr <<  getpid() << " Ce = " << Ce;
  // std::cerr <<  getpid() << " f_sigma = " << f_sigma << std::endl;
  // std::cerr <<  getpid() << " f_q1 = " << f_q1 << " f_q2 = " << f_q2 << std::endl;
  // std::cerr <<  getpid() << " h_q1 = " << h_q1 << " h_q2 = " << h_q2 << std::endl <<
  // endl;
  double fqhq = f_q1 * h_q1 + f_q2 * h_q2;
  Matrix3 Cr(0.0), rC(0.0);
  double rCr = 0.0;
  for (int ii = 0; ii < 3; ++ii) {
    for (int jj = 0; jj < 3; ++jj) {
      Cr(ii, jj) = 0.0;
      rC(ii, jj) = 0.0;
      for (int kk = 0; kk < 3; ++kk) {
        for (int ll = 0; ll < 3; ++ll) {
          double Ce1 = Ce(ii, jj, kk, ll);
          double Ce2 = Ce(kk, ll, ii, jj);
          double fs  = f_sigma(kk, ll);
          Cr(ii, jj) += Ce1 * fs;
          rC(ii, jj) += fs * Ce2;
        }
      }
      rCr += rC(ii, jj) * f_sigma(ii, jj);
    }
  }
  double rCr_fqhq = rCr - fqhq;
  for (int ii = 0; ii < 3; ++ii) {
    for (int jj = 0; jj < 3; ++jj) {
      for (int kk = 0; kk < 3; ++kk) {
        for (int ll = 0; ll < 3; ++ll) {
          Cep(ii, jj, kk, ll) =
            Ce(ii, jj, kk, ll) - Cr(ii, jj) * rC(kk, ll) / rCr_fqhq;
        }
      }
    }
  }
}

void
YieldCond_Gurson::computeElasPlasTangentModulus(const TangentModulusTensor& Ce,
                                                const Matrix3& sigma,
                                                double sigY,
                                                double dsigYdep,
                                                double porosity,
                                                double voidNuclFac,
                                                TangentModulusTensor& Cep)
{
  ModelStateBase state;
  state.setStress(sigma);
  state.yieldStress = sigY;
  state.porosity    = porosity;
  state.backStress  = Matrix3(0.0);

  Matrix3 xi = state.devStress - state.backStress.Deviator();

  // Calculate the derivative of the yield function wrt sigma
  Matrix3 f_sigma = df_dsigma(xi, &state);

  // Calculate derivative wrt internal variables
  MetalIntVar f_intvar; 
  df_dintvar(&state, f_intvar);
  double f_q1 = f_intvar.plasticPorosity;
  double f_q2 = f_intvar.eqPlasticStrain;

  // Compute hardening moduli
  MetalIntVar h_eta;
  d_intvar->computeHardeningModulus(&state, h_eta);
  double h_q1 = h_eta.eqPlasticStrain;
  double h_q2 = h_eta.plasticPorosity;

  // Calculate elastic-plastic tangent modulus
  computeTangentModulus(Ce, f_sigma, f_q1, f_q2, h_q1, h_q2, Cep);
}

