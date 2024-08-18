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
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCond_vonMises.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <cmath>

using namespace Uintah;
using namespace Vaango;

YieldCond_vonMises::YieldCond_vonMises(Uintah::ProblemSpecP& ps,
                                       IntVar_Metal* intvar,
                                       const FlowStressModel* flow)
{
  d_intvar = intvar;
  d_flow = flow;
}

YieldCond_vonMises::YieldCond_vonMises(const YieldCond_vonMises* ym)
{
  d_intvar = ym->d_intvar;
  d_flow = ym->d_flow;
}

YieldCond_vonMises::~YieldCond_vonMises() = default;

void
YieldCond_vonMises::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP yield_ps = ps->appendChild("yield_condition");
  yield_ps->setAttribute("type", "von_mises");
}

//--------------------------------------------------------------
// Compute value of yield function
//--------------------------------------------------------------
std::pair<double, Util::YieldStatus>
YieldCond_vonMises::evalYieldCondition(const ModelStateBase* state)
{
  double f = computeYieldFunction(state);
  if (f <= 0.0) {
    return std::make_pair(f, Util::YieldStatus::IS_ELASTIC);
  }
  return std::make_pair(f, Util::YieldStatus::HAS_YIELDED);
}

double
YieldCond_vonMises::computeYieldFunction(const ModelStateBase* state) const
{
  return computeYieldFunction(state->getStress(), state);
}

double
YieldCond_vonMises::evalYieldCondition(const Matrix3& stress,
                                       const ModelStateBase* state)
{
  return computeYieldFunction(stress, state);
}

double
YieldCond_vonMises::computeYieldFunction(const Matrix3& stress,
                                         const ModelStateBase* state) const
{
  Matrix3 xi    = stress.Deviator() - state->backStress.Deviator();
  double sigy   = state->yieldStress;
  double xiNorm = xi.Norm();
  double f      = Vaango::Util::sqrt_three_half * xiNorm - sigy;
  //std::cout << __FILE__ << "\n"
  //          << "xi = " << xi << " ||xi|| = " << xiNorm << " sigy = " << sigy
  //          << " f = " << f << "\n";
  return f;
}

double
YieldCond_vonMises::evalYieldConditionMax(const ModelStateBase* state)
{
  return state->yieldStress;
}

//--------------------------------------------------------------
// Compute derivatives of yield function
// df/dsigma = df/ds = sqrt(3/2) s / ||s||
// ||df/dsigma|| = ||df/ds|| = sqrt(3/2)
// df/dsigma / ||df/dsigma|| = s // ||s||
//--------------------------------------------------------------
/*! Derivative with respect to the Cauchy stress (\f$\sigma \f$)*/
Matrix3
YieldCond_vonMises::df_dsigma(const ModelStateBase* state) 
{
  return df_dsigma(state->getStress(), state);
}

/*! Derivative with respect to the Cauchy stress (\f$\sigma \f$)
    Assume f = sqrt{3/2} ||xi|| - sigma_y , xi = s - beta
    df/dsigma = sqrt(3/2) (xi / ||xi|| + 1/3 tr(beta) / ||xi|| I) */
Matrix3
YieldCond_vonMises::df_dsigma(const Matrix3& stress, const ModelStateBase* state)
{
  Matrix3 df_dsigma = df_dxi(stress, state);
  return df_dsigma;
}

/*! Derivative with respect to the \f$xi\f$ where \f$\xi = s - \beta \f$
    where \f$s\f$ is deviatoric part of Cauchy stress and
    \f$\beta\f$ is the backstress
    Assume f = sqrt{3/2} ||xi|| - sigma_y
    df/dxi = sqrt(3/2) xi / ||xi|| */
Matrix3
YieldCond_vonMises::df_dxi(const Matrix3& stress, const ModelStateBase* state)
{
  Matrix3 xi     = stress.Deviator() - state->backStress.Deviator();
  Matrix3 df_dxi = xi * (Vaango::Util::sqrt_three_half / xi.Norm());
  return df_dxi;
}

/* Derivative with respect to \f$ s \f$ and \f$ \beta \f$ */
std::pair<Matrix3, Matrix3>
YieldCond_vonMises::df_dsigmaDev_dbeta(const Matrix3& stress,
                                       const ModelStateBase* state)
{
  Matrix3 df_ds    = df_dxi(stress, state);
  Matrix3 df_dbeta = df_ds * (-1.0);
  return std::make_pair(df_ds, df_dbeta);
}

/*! Derivative with respect to internal variables 
    Assume f = sqrt{3/2} ||xi|| - sigma_y */
void
YieldCond_vonMises::df_dintvar(const ModelStateBase* state,
                               MetalIntVar& df_dintvar) const
{
  double dsigy_dep = d_flow->evalDerivativeWRTPlasticStrain(state, 0);
  double df_dplasticStrain = -dsigy_dep;
  double df_dporosity = 0.0;
  df_dintvar = {df_dplasticStrain, df_dporosity};
}

void
YieldCond_vonMises::computeElasPlasTangentModulus(
  const TangentModulusTensor& Ce,
  const Matrix3& sigma,
  double sigY,
  double dsigYdep,
  double porosity,
  double,
  TangentModulusTensor& Cep)
{
  ModelStateBase state;
  state.setStress(sigma);
  state.yieldStress = sigY;
  state.porosity    = 0.0;
  state.backStress  = Matrix3(0.0);

  Matrix3 xi = state.devStress - state.backStress.Deviator();

  // Calculate the derivative of the yield function wrt sigma
  Matrix3 f_sigma = df_dsigma(xi, &state);

  // Calculate derivative wrt plastic strain
  double f_q1 = dsigYdep;

  // Compute h_q1
  double sigma_f_sigma = sigma.Contract(f_sigma);
  double h_q1          = computePlasticStrainFactor(sigma_f_sigma, sigY);

  // Calculate elastic-plastic tangent modulus
  computeTangentModulus(Ce, f_sigma, f_q1, h_q1, Cep);
}

double
YieldCond_vonMises::computePlasticStrainFactor(double sigma_f_sigma,
                                               double sigma_Y)
{
  return sigma_f_sigma / sigma_Y;
}

void
YieldCond_vonMises::computeTangentModulus(const TangentModulusTensor& Ce,
                                          const Matrix3& f_sigma,
                                          double f_q1,
                                          double h_q1,
                                          TangentModulusTensor& Cep)
{
  double fqhq = f_q1 * h_q1;
  Matrix3 Cr(0.0), rC(0.0);
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
