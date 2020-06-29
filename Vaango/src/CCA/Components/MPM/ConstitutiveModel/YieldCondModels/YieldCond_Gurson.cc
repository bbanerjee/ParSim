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
  Matrix3 xi = state->devStress - state->backStress.Deviator();
  double f = evalYieldCondition(xi, state);
  if (f > 0) {
    return std::make_pair(f, Util::YieldStatus::HAS_YIELDED);
  }
  return std::make_pair(f, Util::YieldStatus::IS_ELASTIC);
}

double
YieldCond_Gurson::evalYieldCondition(const Matrix3& xi,
                                     const ModelStateBase* state)
{
  // Get the state data
  double sigy     = state->yieldStress;
  double porosity = state->porosity;
  double sig_m    = state->p;
  double xi_xi    = xi.NormSquared();
  double sig_eq   = std::sqrt(1.5 * xi_xi);
  double sig      = 0.0;

  return evalYieldCondition(sig_eq, sigy, 3.0*sig_m, porosity, sig);
}

double
YieldCond_Gurson::evalYieldCondition(const double sig_eq,
                                     const double sigy,
                                     const double sig_tr,
                                     const double porosity,
                                     double& sig)
{
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
  double incosh = 0.5 * q2 * sig_tr / sigy;
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
  Matrix3 xi = state->devStress - state->backStress.Deviator();
  return df_dsigma(xi, state);
}

Uintah::Matrix3
YieldCond_Gurson::df_dsigma(const Matrix3& xi,
                            const ModelStateBase* state)
{
  double sigy = state->yieldStress;
  ASSERT(sigy != 0);

  double sig_m    = state->p;
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
  double b  = q1 * q2 * phi_star * sinh(insinh);
  Matrix3 df_dsigma = xi * a + Vaango::Util::Identity * b;
  df_dsigma /= df_dsigma.Norm();
  return df_dsigma;
}

/*! Derivative with respect to the \f$xi\f$ where \f$\xi = s - \beta \f$
    where \f$s\f$ is deviatoric part of Cauchy stress and
    \f$\beta\f$ is the backstress */
Uintah::Matrix3
YieldCond_Gurson::df_dxi(const Matrix3& xi,
                         const ModelStateBase* state)
{
  double sigy = state->yieldStress;
  ASSERT(sigy != 0);
  double a = 3.0 / (sigy * sigy);
  Matrix3 df_dxi   = xi * a;
  df_dxi /= df_dxi.Norm();
  return df_dxi;
}

/* Derivative with respect to \f$ s \f$ and \f$ \beta \f$ */
std::pair<Uintah::Matrix3, Uintah::Matrix3>
YieldCond_Gurson::df_dsigmaDev_dbeta(const Matrix3& xi,
                                     const ModelStateBase* state)
{
  Matrix3 df_ds = df_dxi(xi, state);
  Matrix3 df_dbeta = df_ds * (-1.0);
  return std::make_pair(df_ds, df_dbeta);
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

  double sig_m    = state->p;
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

  double sig_m    = state->p;
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

inline double
YieldCond_Gurson::computePorosityFactor_h1(double sigma_f_sigma,
                                           double tr_f_sigma,
                                           double porosity,
                                           double sigma_Y,
                                           double A)
{
  return (1.0 - porosity) * tr_f_sigma +
         A * sigma_f_sigma / ((1.0 - porosity) * sigma_Y);
}

inline double
YieldCond_Gurson::computePlasticStrainFactor_h2(double sigma_f_sigma,
                                                double porosity,
                                                double sigma_Y)
{
  return sigma_f_sigma / ((1.0 - porosity) * sigma_Y);
}

/*! Compute h_alpha  where \f$d/dt(ep) = d/dt(gamma)~h_{\alpha}\f$ */
double
YieldCond_Gurson::eval_h_alpha(const Matrix3& xi, const ModelStateBase* state)
{
  double sigy = state->yieldStress;
  ASSERT(sigy != 0);
  double phi = state->porosity;
  ASSERT(phi != 1);
  double p         = state->pressure;
  double trBetaHat = state->backStress.Trace();

  Matrix3 One;
  One.Identity();
  Matrix3 xi_hat = xi + One * (p - 1.0 / 3.0 * trBetaHat);

  Matrix3 dfdsigma = df_dsigma(xi, state);

  double numer = xi_hat.Contract(dfdsigma);
  double denom = (1.0 - phi) * sigy;
  return numer / denom;
}

/*! Compute h_phi  where \f$d/dt(phi) = d/dt(gamma)~h_{\phi}\f$ */
double
YieldCond_Gurson::eval_h_phi(const Matrix3& xi,
                             const double& factorA,
                             const ModelStateBase* state)
{
  double sigy = state->yieldStress;
  ASSERT(sigy != 0);
  double phi = state->porosity;
  ASSERT(phi != 1);
  double p = state->pressure;

  Matrix3 dfdsigma = df_dsigma(xi, state);

  double tr_df_dsigma = dfdsigma.Trace();
  double a            = (1.0 - phi) * tr_df_dsigma;

  Matrix3 One;
  One.Identity();
  double trBetaHat = state->backStress.Trace();
  Matrix3 xi_hat   = xi + One * (p - 1.0 / 3.0 * trBetaHat);

  double numer = xi_hat.Contract(dfdsigma);
  double denom = (1.0 - phi) * sigy;
  double b     = factorA * numer / denom;

  return (a + b);
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
  // cerr << getpid() << " Ce = " << Ce;
  // cerr << getpid() << " f_sigma = " << f_sigma << endl;
  // cerr << getpid() << " f_q1 = " << f_q1 << " f_q2 = " << f_q2 << endl;
  // cerr << getpid() << " h_q1 = " << h_q1 << " h_q2 = " << h_q2 << endl <<
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

  // Compute h_q1
  double sigma_f_sigma = sigma.Contract(f_sigma);
  double tr_f_sigma    = f_sigma.Trace();
  double h_q1          = computePorosityFactor_h1(
    sigma_f_sigma, tr_f_sigma, porosity, sigY, voidNuclFac);

  // Compute h_q2
  double h_q2 = computePlasticStrainFactor_h2(sigma_f_sigma, porosity, sigY);

  // Calculate elastic-plastic tangent modulus
  computeTangentModulus(Ce, f_sigma, f_q1, f_q2, h_q1, h_q2, Cep);
}
