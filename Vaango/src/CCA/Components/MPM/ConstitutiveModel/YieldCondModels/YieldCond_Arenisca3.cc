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

#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCond_Arenisca3.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <cmath>

using namespace Vaango;

const double YieldCond_Arenisca3::sqrt_three     = std::sqrt(3.0);
const double YieldCond_Arenisca3::one_sqrt_three = 1.0 / sqrt_three;

YieldCond_Arenisca3::YieldCond_Arenisca3(Uintah::ProblemSpecP& ps)
{
  // Nonlinear Drucker-Prager parameters
  ps->require("PEAKI1", d_inputParam.PEAKI1); // Shear Limit Surface Parameter
  ps->require("FSLOPE", d_inputParam.FSLOPE); // Shear Limit Surface Parameter
  ps->require("STREN", d_inputParam.STREN);   // Shear Limit Surface Parameter
  ps->require("YSLOPE", d_inputParam.YSLOPE); // Shear Limit Surface Parameter
  ps->require("BETA_nonassociativity",
              d_inputParam.BETA_nonassociativity); // Nonassociativity Parameter

  // Cap parameters: CR = (peakI1-kappa)/(peakI1-X)
  ps->require("CR", d_capParam.CR); // Cap Shape Parameter

  // Check the input parameters
  checkInputParameters();

  // Compute the model parameters from the input parameters
  computeModelParameters();
}

YieldCond_Arenisca3::YieldCond_Arenisca3(const YieldCond_Arenisca3* yc)
{
  d_inputParam = yc->d_inputParam;
  d_modelParam = yc->d_modelParam;
}

YieldCond_Arenisca3::~YieldCond_Arenisca3() = default;

void
YieldCond_Arenisca3::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP yield_ps = ps->appendChild("yield_condition");
  yield_ps->setAttribute("type", "arenisca3");
  yield_ps->appendElement("FSLOPE", d_inputParam.FSLOPE);
  yield_ps->appendElement("PEAKI1", d_inputParam.PEAKI1);
  yield_ps->appendElement("STREN", d_inputParam.STREN);
  yield_ps->appendElement("YSLOPE", d_inputParam.YSLOPE);
  yield_ps->appendElement("BETA_nonassociativity",
                          d_inputParam.BETA_nonassociativity);

  yield_ps->appendElement("CR", d_capParam.CR);
}

//--------------------------------------------------------------
// Check that the input parameters are reasonable
//--------------------------------------------------------------
void
YieldCond_Arenisca3::checkInputParameters()
{
  std::ostringstream warn;
  if (d_inputParam.PEAKI1 < 0.0) {
    warn << "PEAKI1 must be nonnegative. PEAKI1 = " << d_inputParam.PEAKI1
         << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_inputParam.FSLOPE < 0.0) {
    warn << "FSLOPE must be nonnegative. FSLOPE = " << d_inputParam.FSLOPE
         << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_inputParam.FSLOPE < d_inputParam.YSLOPE) {
    warn << "FSLOPE must be greater than YSLOPE. FSLOPE = "
         << d_inputParam.FSLOPE << ", YSLOPE = " << d_inputParam.YSLOPE
         << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_inputParam.BETA_nonassociativity <= 0.0) {
    warn << "BETA_nonassociativity must be positive. BETA_nonassociativity = "
         << d_inputParam.BETA_nonassociativity << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_capParam.CR >= 1 || d_capParam.CR <= 0.0) {
    warn << "CR must be 0<CR<1. CR = " << d_capParam.CR << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
}

//--------------------------------------------------------------
// Compute the model parameters a1, a2, a3, a4, beta from the
// input parameters FSLOPE, PEAKI1, STREN, SLOPE, BETA_nonassoc
// The shear limit surface is defined in terms of the a1,a2,a3,a4 parameters,
// but
// the user inputs are the more intuitive set of FSLOPE. YSLOPE, STREN, and
// PEAKI1.
//
// Note: This routine computes the a_i parameters from the user inputs.  The
// code was
// originally written by R.M. Brannon, with modifications by M.S. Swan.
//--------------------------------------------------------------
void
YieldCond_Arenisca3::computeModelParameters(double factor)
{
  double FSLOPE = d_inputParam.FSLOPE, // Slope at I1=PEAKI1
    STREN       = d_inputParam.STREN,  // Value of rootJ2 at I1=0
    YSLOPE      = d_inputParam.YSLOPE, // High pressure slope
    PEAKI1      = d_inputParam.PEAKI1; // Value of I1 at strength=0

  if (FSLOPE > 0.0 && PEAKI1 >= 0.0 && STREN == 0.0 &&
      YSLOPE == 0.0) { // ----------------------------------------------Linear
                       // Drucker Prager
    d_modelParam.a1 = PEAKI1 * FSLOPE;
    d_modelParam.a2 = 0.0;
    d_modelParam.a3 = 0.0;
    d_modelParam.a4 = FSLOPE;
  } else if (FSLOPE == 0.0 && PEAKI1 == 0.0 && STREN > 0.0 &&
             YSLOPE ==
               0.0) { // -------------------------------------------------------
                      // Von Mises
    d_modelParam.a1 = STREN;
    d_modelParam.a2 = 0.0;
    d_modelParam.a3 = 0.0;
    d_modelParam.a4 = 0.0;
  } else if (FSLOPE > 0.0 && YSLOPE == 0.0 && STREN > 0.0 &&
             PEAKI1 ==
               0.0) { // -------------------------------------------------------
                      // 0 PEAKI1 to vonMises
    d_modelParam.a1 = STREN;
    d_modelParam.a2 = FSLOPE / STREN;
    d_modelParam.a3 = STREN;
    d_modelParam.a4 = 0.0;
  } else if (FSLOPE > YSLOPE && YSLOPE > 0.0 && STREN > YSLOPE * PEAKI1 &&
             PEAKI1 >=
               0.0) { // -------------------------------------------------------
                      // Nonlinear Drucker-Prager
    d_modelParam.a1 = STREN;
    d_modelParam.a2 = (FSLOPE - YSLOPE) / (STREN - YSLOPE * PEAKI1);
    d_modelParam.a3 =
      (STREN - YSLOPE * PEAKI1) * exp(-d_modelParam.a2 * PEAKI1);
    d_modelParam.a4 = YSLOPE;
  } else {
    // Bad inputs, throw exception:
    std::ostringstream warn;
    warn << "Bad input parameters for shear limit surface. "
         << "FSLOPE = " << FSLOPE << ", YSLOPE = " << YSLOPE
         << ", PEAKI1 = " << PEAKI1 << ", STREN = " << STREN << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
}

//--------------------------------------------------------------
// Evaluate yield condition (q = state->q
//                           p = state->p)
// f := J2 - Ff^2*Fc^2 = 0
// where
//     J2 = q^2/3
//     I1 = 3*(p - zeta)
//     Ff := a1 - a3*exp(a2*I1) - a4*I1
//     Fc^2 := 1 - (kappa - 3*p)^2/(kappa - X)^2
//     kappa = I1_0 - CR*(I1_0 - X)
//--------------------------------------------------------------
std::pair<double, Util::YieldStatus>
YieldCond_Arenisca3::evalYieldCondition(const ModelStateBase* state_input)
{
  const ModelState_Arenisca3* state =
    static_cast<const ModelState_Arenisca3*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Arenisca3.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  // Get the local vars from the model state
  double zeta  = state->zeta;
  double kappa = state->kappa;
  double capX  = state->capX;

  // Initialize hasYielded to 0
  double hasYielded = -1.0;

  // Cauchy stress invariants: I1 = 3*p, J2 = q^2/3
  double I1      = state->I1;
  double sqrt_J2 = state->sqrt_J2;

  // Shift stress to evaluate yield condition where zeta is the isotropic back
  // stress
  double I1_eff = (I1 - zeta);

  // --------------------------------------------------------------------
  // *** SHEAR LIMIT FUNCTION (Ff) ***
  // --------------------------------------------------------------------
  double Ff = d_modelParam.a1 -
              d_modelParam.a3 * exp(d_modelParam.a2 * I1_eff) -
              d_modelParam.a4 * I1_eff;

  // --------------------------------------------------------------------
  // *** Branch Point (Kappa) ***
  // --------------------------------------------------------------------
  kappa = d_inputParam.PEAKI1 -
          d_capParam.CR * (d_inputParam.PEAKI1 - capX); // Branch Point

  // --------------------------------------------------------------------
  // *** COMPOSITE YIELD FUNCTION ***
  // --------------------------------------------------------------------
  // Evaluate Composite Yield Function F(I1) = Ff(I1)*fc(I1) in each region.
  // The elseif statements have nested if statements, which is not equivalent
  // to them having a single elseif(A&&B&&C)
  //---------------------------------------------------(I1<X)
  if (I1_eff < capX) {
    hasYielded = 1.0;
    return std::make_pair(hasYielded, Util::YieldStatus::HAS_YIELDED);
  }

  // **Elliptical Cap Function: (fc)**
  // fc = sqrt(1.0 - Pow((Kappa-I1mZ)/(Kappa-X)),2.0);
  // faster version: fc2 = fc^2
  // **WARNING** p3 is the maximum achievable volumetric plastic strain in
  // compresson
  // so if a value of 0 has been specified this indicates the user
  // wishes to run without porosity, and no cap function is used, i.e. fc=1
  if ((capX < I1_eff) && (I1_eff < kappa)) { // ---------------(X<I1<kappa)

    double kappaRatio = (kappa - I1_eff) / (kappa - capX);
    double fc2        = 1.0 - kappaRatio * kappaRatio;
    if (sqrt_J2 * sqrt_J2 > Ff * Ff * fc2) {
      hasYielded = 1.0;
      return std::make_pair(hasYielded, Util::YieldStatus::HAS_YIELDED);
    }
  } else { // --------- X >= I1 or kappa <= I1

    if (I1_eff <= d_inputParam.PEAKI1) { // ----- (kappa <= I1 <= PEAKI1)
      if (sqrt_J2 > Ff) {
        hasYielded = 1.0;
        return std::make_pair(hasYielded, Util::YieldStatus::HAS_YIELDED);
      }
    } else { // I1 > PEAKI1
      hasYielded = 1.0;
      return std::make_pair(hasYielded, Util::YieldStatus::HAS_YIELDED);
    }
  }

  return std::make_pair(hasYielded, Util::YieldStatus::IS_ELASTIC);
}

//--------------------------------------------------------------
// Evaluate yield condition max (q = state->q
//                               p = state->p)
//--------------------------------------------------------------
double
YieldCond_Arenisca3::evalYieldConditionMax(const ModelStateBase*)
{
  std::ostringstream out;
  out << "**ERROR** evalYieldConditionMax should not be called by "
      << " models that use the Arenisca3 yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return 0.0;
}

//--------------------------------------------------------------
// Derivatives needed by return algorithms and Newton iterations

//--------------------------------------------------------------
// Compute df/dp  where p = volumetric stress = 1/3 Tr(sigma)
//   df/dp = derivative of the yield function wrt p
//
// f := J2 - Ff^2*Fc^2 = 0
// where
//     J2 = q^2/3
//     I1 = 3*(p - zeta)
//     Ff := a1 - a3*exp(a2*I1) - a4*I1
//     Fc^2 := 1 - (kappa - 3*p)^2/(kappa - X)^2
//     kappa = I1_0 - CR*(I1_0 - X)
//
// df/dp = 6*Ff*(a2*a3*exp(a2*I1) + a4)*Fc^2 -
//             6*Ff^2*(kappa - I1)/(kappa - X)^2
//--------------------------------------------------------------
double
YieldCond_Arenisca3::df_dp(const ModelStateBase* state_input)
{
  const ModelState_Arenisca3* state =
    static_cast<const ModelState_Arenisca3*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Arenisca3.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  // Get the local vars from the model state
  double zeta  = state->zeta;
  double kappa = state->kappa;
  double capX  = state->capX;

  // Cauchy stress invariants: I1 = 3*p, J2 = q^2/3
  double I1 = state->I1;

  // Shift stress to evaluate yield condition where zeta is the isotropic back
  // stress
  double I1_eff = (I1 - zeta);

  // --------------------------------------------------------------------
  // *** SHEAR LIMIT FUNCTION (Ff) ***
  // --------------------------------------------------------------------
  double Ff = d_modelParam.a1 -
              d_modelParam.a3 * exp(d_modelParam.a2 * I1_eff) -
              d_modelParam.a4 * I1_eff;

  // --------------------------------------------------------------------
  // *** Branch Point (Kappa) ***
  // --------------------------------------------------------------------
  kappa = d_inputParam.PEAKI1 -
          d_capParam.CR * (d_inputParam.PEAKI1 - capX); // Branch Point

  // --------------------------------------------------------------------
  // **Elliptical Cap Function: (fc)**
  // --------------------------------------------------------------------
  double kappaRatio = (kappa - I1_eff) / (kappa - capX);
  double Fc_sq      = 1.0 - kappaRatio * kappaRatio;

  // --------------------------------------------------------------------
  // Derivatives
  // --------------------------------------------------------------------
  // term1 = 6*Ff*(a2*a3*exp(a2*I1) + a4)*Fc^2
  double term1 = 6.0 * Ff * Fc_sq * (d_modelParam.a2 * d_modelParam.a3 *
                                       exp(d_modelParam.a2 * I1_eff) +
                                     d_modelParam.a4);
  // term2 = 6*Ff^2*(kappa - I1)/(kappa - X)^2
  double term2 = 6.0 * Ff * Ff * kappaRatio / (kappa - capX);

  return (term1 - term2);
}

//--------------------------------------------------------------
// Compute df/dq
// f := J2 - Ff^2*Fc^2 = 0
// where
//     J2 = q^2/3
// df/dq = 2q/3
//--------------------------------------------------------------
double
YieldCond_Arenisca3::df_dq(const ModelStateBase* state_input)
{
  const ModelState_Arenisca3* state =
    static_cast<const ModelState_Arenisca3*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Arenisca3.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  return (2.0 / std::sqrt(3.0)) * state->sqrt_J2;
}

//--------------------------------------------------------------
// Compute d/depse_v(df/dp)
//   df/dp = 6*Ff*(a2*a3*exp(3*a2*p) + a4)*Fc^2 -
//             6*Ff^2*(kappa - I1)/(kappa - X)^2
//   d/depse_v(df/dp) =
//
// Requires:  Equation of state and internal variable
//--------------------------------------------------------------
double
YieldCond_Arenisca3::d2f_dp_depsVol(const ModelStateBase* state_input,
                                    const MPMEquationOfState* eos,
                                    const ShearModulusModel*)
{
  std::ostringstream out;
  out << "**ERROR** d2f_dp_depsVol should not be called by "
      << " models that use the Arenisca3 yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return 0.0;
}

//--------------------------------------------------------------
// Compute d/depse_s(df/dp)
//   df/dp =
//   d/depse_s(df/dp) =
//
// Requires:  Equation of state
//--------------------------------------------------------------
double
YieldCond_Arenisca3::d2f_dp_depsDev(const ModelStateBase* state_input,
                                    const MPMEquationOfState* eos,
                                    const ShearModulusModel*)
{
  std::ostringstream out;
  out << "**ERROR** d2f_dp_depsDev should not be called by "
      << " models that use the Arenisca3 yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return 0.0;
}

//--------------------------------------------------------------
// Compute d/depse_v(df/dq)
//   df/dq =
//   d/depse_v(df/dq) =
//
// Requires:  Shear modulus model
//--------------------------------------------------------------
double
YieldCond_Arenisca3::d2f_dq_depsVol(const ModelStateBase* state_input,
                                    const MPMEquationOfState*,
                                    const ShearModulusModel* shear)
{
  std::ostringstream out;
  out << "**ERROR** d2f_dq_depsVol should not be called by "
      << " models that use the Arenisca3 yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return 0.0;
}

//--------------------------------------------------------------
// Compute d/depse_s(df/dq)
//   df/dq =
//   d/depse_s(df/dq) =
//
// Requires:  Shear modulus model
//--------------------------------------------------------------
double
YieldCond_Arenisca3::d2f_dq_depsDev(const ModelStateBase* state_input,
                                    const MPMEquationOfState*,
                                    const ShearModulusModel* shear)
{
  std::ostringstream out;
  out << "**ERROR** d2f_dq_depsDev should not be called by "
      << " models that use the Arenisca3 yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return 0.0;
}

//--------------------------------------------------------------
// Compute df/depse_v
//   df/depse_v =
//
// Requires:  Equation of state, shear modulus model, internal variable model
//--------------------------------------------------------------
double
YieldCond_Arenisca3::df_depsVol(const ModelStateBase* state_input,
                                const MPMEquationOfState* eos,
                                const ShearModulusModel* shear)
{
  std::ostringstream out;
  out << "**ERROR** df_depsVol should not be called by "
      << " models that use the Arenisca3 yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return 0.0;
}

//--------------------------------------------------------------
// Compute df/depse_s
//   df/depse_s =
//
// Requires:  Equation of state, shear modulus model
//--------------------------------------------------------------
double
YieldCond_Arenisca3::df_depsDev(const ModelStateBase* state_input,
                                const MPMEquationOfState* eos,
                                const ShearModulusModel* shear)
{
  std::ostringstream out;
  out << "**ERROR** df_depsVol should not be called by "
      << " models that use the Arenisca3 yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return 0.0;
}

//--------------------------------------------------------------
// Other yield condition functions

// Evaluate the yield function.
double
YieldCond_Arenisca3::evalYieldCondition(const double p,
                                        const double q,
                                        const double dummy0,
                                        const double dummy1,
                                        double& dummy2)
{
  std::ostringstream out;
  out << "**ERROR** Deprecated evalYieldCondition with double arguments. "
      << " Should not be called by models that use the Arenisca3 yield "
         "criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return 0.0;
}

// Evaluate yield condition (s = deviatoric stress
//                           p = state->p)
double
YieldCond_Arenisca3::evalYieldCondition(const Uintah::Matrix3&,
                                        const ModelStateBase* state_input)
{
  std::ostringstream out;
  out << "**ERROR** evalYieldCondition with a Matrix3 argument should not be "
         "called by "
      << " models that use the Arenisca3 yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return 0.0;
}

//--------------------------------------------------------------
// Other derivatives

// Compute df/dsigma
//    df/dsigma =
// where
//    s = sigma - 1/3 tr(sigma) I
Uintah::Matrix3
YieldCond_Arenisca3::df_dsigma(const Uintah::Matrix3& sig,
                               const double p_c,
                               const double)
{
  std::ostringstream out;
  out << "**ERROR** evalDerivOfYieldCondition with a Matrix3 argument should "
         "not be "
      << "called by models that use the Arenisca3 yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return Uintah::Matrix3(0.0);
}

// Compute df/ds  where s = deviatoric stress
//    df/ds =
Uintah::Matrix3
YieldCond_Arenisca3::df_dsigmaDev(const Uintah::Matrix3& sigDev,
                                  const double,
                                  const double)
{
  std::ostringstream out;
  out << "**ERROR** evalDerivOfYieldCondition with a Matrix3 argument should "
         "not be "
      << "called by models that use the Arenisca3 yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return Uintah::Matrix3(0.0);
}

/*! Derivative with respect to the Cauchy stress (\f$\sigma \f$) */
Uintah::Matrix3
YieldCond_Arenisca3::df_dsigma(const Matrix3& sig,
                               const ModelStateBase* state_input)
{
  std::ostringstream out;
  out << "**ERROR** df_dsigma with a Matrix3 argument should not be "
      << "called by models that use the Arenisca3 yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return Uintah::Matrix3(0.0);
}

/*! Derivative with respect to the \f$xi\f$ where \f$\xi = s \f$
    where \f$s\f$ is deviatoric part of Cauchy stress */
Uintah::Matrix3
YieldCond_Arenisca3::df_dxi(const Matrix3& sigDev,
                            const ModelStateBase*)

{
  std::ostringstream out;
  out << "**ERROR** df_dxi with a Matrix3 argument should not be "
      << "called by models that use the Arenisca3 yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return Uintah::Matrix3(0.0);
}

/* Derivative with respect to \f$ s \f$ and \f$ \beta \f$ */
std::pair<Uintah::Matrix3, Uintah::Matrix3>
YieldCond_Arenisca3::df_dsigmaDev_dbeta(const Matrix3& sigDev,
                                        const ModelStateBase*)
{
  std::ostringstream out;
  out << "**ERROR** df_dsigmaDev_dbeta with a Matrix3 argument should not be "
      << "called by models that use the Arenisca3 yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return std::make_pair(Uintah::Matrix3(0.0), Uintah::Matrix3(0.0));
}

/*! Compute h_alpha  where \f$d/dt(ep) = d/dt(gamma)~h_{\alpha}\f$ */
double
YieldCond_Arenisca3::eval_h_alpha(const Matrix3&, const ModelStateBase*)
{
  std::ostringstream out;
  out << "**ERROR** eval_h_alpha with a Matrix3 argument should not be "
      << "called by models that use the Arenisca3 yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return 1.0;
}

/*! Compute h_phi  where \f$d/dt(phi) = d/dt(gamma)~h_{\phi}\f$ */
double
YieldCond_Arenisca3::eval_h_phi(const Matrix3&,
                                const double&,
                                const ModelStateBase*)
{
  std::ostringstream out;
  out << "**ERROR** eval_h_phi with a Matrix3 argument should not be "
      << "called by models that use the Arenisca3 yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return 0.0;
}

//--------------------------------------------------------------
// Tangent moduli
void
YieldCond_Arenisca3::computeElasPlasTangentModulus(
  const TangentModulusTensor& Ce,
  const Matrix3& sigma,
  double sigY,
  double dsigYdep,
  double porosity,
  double,
  TangentModulusTensor& Cep)
{
  std::ostringstream out;
  out << "**ERROR** computeElasPlasTangentModulus with a Matrix3 argument "
         "should not be "
      << "called by models that use the Arenisca3 yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return;
}

void
YieldCond_Arenisca3::computeTangentModulus(const TangentModulusTensor& Ce,
                                           const Matrix3& f_sigma,
                                           double f_q1,
                                           double h_q1,
                                           TangentModulusTensor& Cep)
{
  std::ostringstream out;
  out << "**ERROR** coputeTangentModulus with a Matrix3 argument should not be "
      << "called by models that use the Arenisca3 yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return;
}
