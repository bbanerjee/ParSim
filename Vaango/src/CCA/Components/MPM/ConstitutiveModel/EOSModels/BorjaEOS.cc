/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/BorjaEOS.h>
#include <Core/Exceptions/ConvergenceFailure.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Math/DEIntegrator.h>
#include <cmath>
#include <iostream>

using namespace Uintah;
using namespace Vaango;

BorjaEOS::BorjaEOS(ProblemSpecP& ps)
{
  ps->require("p0", d_p0);
  ps->require("alpha", d_alpha);
  ps->require("kappatilde", d_kappatilde);
  ps->require("epse_v0", d_epse_v0);

  // Compute kappa_hat for large deformation
  d_kappahat = d_kappatilde / (1.0 - d_kappatilde);
}

BorjaEOS::BorjaEOS(const BorjaEOS* cm)
{
  d_p0         = cm->d_p0;
  d_alpha      = cm->d_alpha;
  d_kappatilde = cm->d_kappatilde;
  d_kappahat   = cm->d_kappahat;
  d_epse_v0    = cm->d_epse_v0;
}

BorjaEOS::~BorjaEOS() = default;

void
BorjaEOS::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP eos_ps = ps->appendChild("equation_of_state");
  eos_ps->setAttribute("type", "borja_pressure");

  eos_ps->appendElement("p0", d_p0);
  eos_ps->appendElement("alpha", d_alpha);
  eos_ps->appendElement("kappatilde", d_kappatilde);
  eos_ps->appendElement("epse_v0", d_epse_v0);
}

//////////
// Calculate the pressure using the Borja pressure model
//  (look at the header file for the equation)
double
BorjaEOS::computePressure(const MPMMaterial*,
                          const ModelStateBase* state_input,
                          const Matrix3&,
                          const Matrix3&,
                          const double&)
{
  const ModelState_CamClay* state =
    static_cast<const ModelState_CamClay*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_CamClay.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  double p = evalPressure(state->epse_v, state->epse_s);
  return p;
}

// Calculate the derivative of p with respect to epse_v
//      where epse_v = tr(epse)
//            epse = total elastic strain
//   dp/depse_v = p0 beta/kappahat exp[(epse_v - epse_v0)/kappahat]
//              = p/kappahat
double
BorjaEOS::computeDpDepse_v(const ModelStateBase* state_input) const
{
  const ModelState_CamClay* state =
    static_cast<const ModelState_CamClay*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_CamClay.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  double dp_depse_v = evalDpDepse_v(state->epse_v, state->epse_s);
  return dp_depse_v;
}

// Calculate the derivative of p with respect to epse_s
//      where epse_s = sqrt{2}{3} ||ee||
//            ee = epse - 1/3 tr(epse) I
//            epse = total elastic strain
double
BorjaEOS::computeDpDepse_s(const ModelStateBase* state_input) const
{
  const ModelState_CamClay* state =
    static_cast<const ModelState_CamClay*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_CamClay.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  return evalDpDepse_s(state->epse_v, state->epse_s);
}

// Compute the derivative of p with respect to J
/* Assume dp/dJ = dp/deps_v (true for inifintesimal strains)
   Then
   dp/dJ = p0 beta/kappahat exp[(epse_v - epse_v0)/kappahat]
         = p/kappahat
*/
double
BorjaEOS::eval_dp_dJ(const MPMMaterial*,
                     const double&,
                     const ModelStateBase* state_input)
{
  return computeDpDepse_v(state_input);
}

// Set the initial value of the bulk modulus
void
BorjaEOS::setInitialBulkModulus()
{
  d_bulkModulus = evalDpDepse_v(0.0, 0.0);
  return;
}

// Compute incremental bulk modulus
double
BorjaEOS::computeBulkModulus(const ModelStateBase* state_input) const
{
  const ModelState_CamClay* state =
    static_cast<const ModelState_CamClay*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_CamClay.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  double K = evalDpDepse_v(state->epse_v, 0.0);
  return K;
}

// Compute volumetric strain energy
//   The strain energy function for the Borja model has the form
//      U(epse_v) = p0 kappahat exp[(epse_v - epse_v0)/kappahat]
double
BorjaEOS::computeStrainEnergy(const ModelStateBase* state_input)
{
  const ModelState_CamClay* state =
    static_cast<const ModelState_CamClay*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_CamClay.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  double Wvol =
    -d_p0 * d_kappahat * exp(-(state->epse_v - d_epse_v0) / d_kappahat);
  return Wvol;
}

// No isentropic increase in temperature with increasing strain
double
BorjaEOS::computeIsentropicTemperatureRate(const double,
                                           const double,
                                           const double,
                                           const double)
{
  return 0.0;
}

//--------------------------------------------------------------------------
// The following are needed for MPMICE coupling
//--------------------------------------------------------------------------

// Compute pressure (option 1) for MPMICE coupling
//   Assume epse_s = 0 for coupling purposes until the interface can be made
//   more general.
double
BorjaEOS::computePressure(const double& rho_orig, const double& rho_cur) const
{
  // Calculate epse_v
  double epse_v = rho_orig / rho_cur - 1.0;
  double p      = evalPressure(epse_v, 0.0);
  return p;
}

// Compute pressure (option 2) - for MPMICE coupling
//   Assume epse_s = 0 for coupling purposes until the interface can be made
//   more general.
//   c^2 = K/rho
//   dp/drho = -(J/rho) dp/depse_v = -(J/rho) K = -J c^2
void
BorjaEOS::computePressure(const double& rho_orig,
                          const double& rho_cur,
                          double& pressure,
                          double& dp_drho,
                          double& csquared)
{
  // Calculate J and epse_v
  double J      = rho_orig / rho_cur;
  double epse_v = J - 1.0;

  pressure = evalPressure(epse_v, 0.0);
  // std::cout << "J = " << J << " epse_v = " << epse_v << " pressure = " <<
  // pressure << std::endl;
  double K = computeBulkModulus(rho_orig, rho_cur);
  csquared = K / rho_cur;
  dp_drho  = -J * csquared;
  return;
}

// Compute the incremental bulk modulus
// The bulk modulus is defined at the tangent to the p - epse_v curve
// keeping epse_s fixed
// i.e., K = dp/depse_v
// For the purposes of coupling to MPMICE we assume that epse_s = 0
// and epse_v = J - 1
double
BorjaEOS::computeInitialBulkModulus() const
{
  double K = evalDpDepse_v(0.0, 0.0);
  return K;
}

double
BorjaEOS::computeBulkModulus(const double& rho_orig,
                             const double& rho_cur) const
{
  // Calculate epse_v
  double epse_v = rho_orig / rho_cur - 1.0;
  double K      = evalDpDepse_v(epse_v, 0.0);
  return K;
}

// Compute density given pressure (tension +ve)
//  rho = rho0/[1 + epse_v0 - kappahat ln(p/p0 beta)]
//  Assume epse_s = 0, i.e., beta = 1
double
BorjaEOS::computeDensity(const double& rho_orig, const double& pressure)
{
  if (pressure >= 0.0) {
    return rho_orig;
  }
  double denom = 1.0 + d_epse_v0 - d_kappahat * log(pressure / d_p0);
  // std::cout << "rho_orig = " << rho_orig << " pressure = " << pressure << "
  // denom = " << denom << std::endl;
  double rho = rho_orig / denom;
  return rho;
}

// Compute volumetric strain energy
//   The strain energy function for the Borja model has the form
//      U(epse_v) = p0 kappahat exp[(epse_v - epse_v0)/kappahat]
double
BorjaEOS::computeStrainEnergy(const double& rho_orig, const double& rho_cur)
{
  // Calculate epse_v
  double epse_v = rho_orig / rho_cur - 1.0;
  double Wvol =
    -d_p0 * d_kappahat * exp(-(epse_v - d_epse_v0) / d_kappahat);
  return Wvol;
}

double
BorjaEOS::computeElasticVolumetricStrain(const double& pp, const double& p0)
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute volume strain of Borja material."
         " It should be provided as an input."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}

double
BorjaEOS::computeExpElasticVolumetricStrain(const double& pp, const double& p0)
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute exp(volume strain) of Borja material."
         " It should be provided as an input."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}
double
BorjaEOS::computeDerivExpElasticVolumetricStrain(const double& pp,
                                                 const double& p0,
                                                 double& exp_eps_e_v)
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute derivative of exp(volume strain) of "
         " Borja material. It should be provided as an input."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}

//-------------------------------------------------------------------------
// Private methods below:
//-------------------------------------------------------------------------

//  Pressure computation
double
BorjaEOS::evalPressure(const double& epse_v, const double& epse_s) const
{
  double beta = 1.0 + 1.5 * (d_alpha / d_kappahat) * (epse_s * epse_s);
  double p    = d_p0 * beta * exp(-(epse_v - d_epse_v0) / d_kappahat);

  return p;
}

//  Pressure derivative computation
double
BorjaEOS::evalDpDepse_v(const double& epse_v, const double& epse_s) const
{
  double p = evalPressure(epse_v, epse_s);
  // std::cout << "p = " << p << " d_kappahat= " << d_kappahat
  //           << " p/kappahat = " << p/d_kappahat
  //           << " epse_v = " << epse_v << " epse_s = " << epse_s << "\n";
  return -p / d_kappahat;
}

//  Shear derivative computation
double
BorjaEOS::evalDpDepse_s(const double& epse_v, const double& epse_s) const
{
  double dbetaDepse_s = 3.0 * (d_alpha / d_kappahat) * epse_s;
  return d_p0 * dbetaDepse_s * exp(-(epse_v - d_epse_v0) / d_kappahat);
}
