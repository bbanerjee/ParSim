/*
 * The MIT License
 *
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

#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/WaterEOS.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Arena.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <cmath>

using namespace Uintah;
using namespace Vaango;

WaterEOS::WaterEOS()
{
  d_p0 = 0.0001; // Hardcoded (SI units).  *TODO* Get as input with ProblemSpec
                 // later.
  d_K0          = 2.21e9;
  d_n           = 6.029;
  d_bulkModulus = d_K0;
}

WaterEOS::WaterEOS(ProblemSpecP&)
{
  d_p0 = 0.0001; // Hardcoded (SI units).  *TODO* Get as input with ProblemSpec
                 // later.
  d_K0          = 2.21e9;
  d_n           = 6.029;
  d_bulkModulus = d_K0;
}

WaterEOS::WaterEOS(const WaterEOS* cm)
{
  d_p0          = cm->d_p0;
  d_K0          = cm->d_K0;
  d_n           = cm->d_n;
  d_bulkModulus = cm->d_bulkModulus;
}

WaterEOS::~WaterEOS() = default;

void
WaterEOS::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP eos_ps = ps->appendChild("equation_of_state");
  eos_ps->setAttribute("type", "water");
}

//////////
// Calculate the pressure using the elastic constitutive equation
double
WaterEOS::computePressure(const MPMMaterial* matl,
                          const ModelStateBase* state_input,
                          const Matrix3&,
                          const Matrix3& rateOfDeformation,
                          const double& delT)
{
  const ModelState_Arena* state =
    static_cast<const ModelState_Arena*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Arena.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  double rho_0 = matl->getInitialDensity();
  double rho   = state->density;
  double p     = computePressure(rho_0, rho);
  return p;
}

// Compute pressure (option 1)
double
WaterEOS::computePressure(const double& rho_orig, const double& rho_cur) const
{
  double J = rho_orig / rho_cur;
  double p = (J > 1.0) ? d_p0 : d_p0 + d_K0 / d_n * (std::pow(J, -d_n) - 1);
  return p;
}

// Compute pressure (option 2)
void
WaterEOS::computePressure(const double& rho_orig,
                          const double& rho_cur,
                          double& pressure,
                          double& dp_drho,
                          double& csquared)
{
  double J = rho_orig / rho_cur;
  if (J > 1.0) {
    pressure     = d_p0;
    double dp_dJ = -d_K0;
    dp_drho      = d_K0 / rho_cur;
    csquared     = dp_dJ / rho_cur;
  } else {
    pressure     = d_p0 + d_K0 / d_n * (std::pow(J, -d_n) - 1);
    double dp_dJ = -d_K0 * std::pow(J, -(d_n + 1));
    dp_drho      = d_K0 / rho_cur * std::pow(J, -d_n);
    csquared     = dp_dJ / rho_cur;
  }
}

// Compute derivative of pressure
double
WaterEOS::eval_dp_dJ(const MPMMaterial* matl,
                     const double& detF,
                     const ModelStateBase*)
{
  /*
  const ModelState_Arena* state =
    static_cast<const ModelState_Arena*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Arena.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  double J    = detF;
  double dpdJ = (J > 1.0) ? -d_K0 : -d_K0 * std::pow(J, -(d_n + 1));
  return dpdJ;
}

// Compute bulk modulus
double
WaterEOS::computeInitialBulkModulus() const
{
  return d_K0;
}

double
WaterEOS::computeBulkModulus(const double& pressure) const
{
  double bulkModulus =
    (pressure < d_p0) ? d_K0 : d_K0 + d_n * (pressure - d_p0);
  return bulkModulus;
}

double
WaterEOS::computeBulkModulus(const double& rho_orig,
                             const double& rho_cur) const
{
  double p           = computePressure(rho_orig, rho_cur);
  double bulkModulus = computeBulkModulus(p);
  return bulkModulus;
}

double
WaterEOS::computeBulkModulus(const ModelStateBase* state_input) const
{
  const ModelState_Arena* state =
    static_cast<const ModelState_Arena*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Arena.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  double p           = state->pbar_w;
  double bulkModulus = computeBulkModulus(p);
  return bulkModulus;
}

// Compute strain energy
double
WaterEOS::computeStrainEnergy(const double& rho_orig, const double& rho_cur)
{
  throw InternalError(
    "ComputeStrainEnergy has not been implemented yet for Water.",
    __FILE__,
    __LINE__);
  return 0.0;
}

double
WaterEOS::computeStrainEnergy(const ModelStateBase* state)
{
  throw InternalError(
    "ComputeStrainEnergy has not been implemented yet for Water.",
    __FILE__,
    __LINE__);
  return 0.0;
}

// Compute density given pressure (tension +ve)
double
WaterEOS::computeDensity(const double& rho_orig, const double& pressure)
{
  throw InternalError("ComputeDensity has not been implemented yet for Water.",
                      __FILE__,
                      __LINE__);
  return 0.0;
}

//  Calculate the derivative of p with respect to epse_v
double
WaterEOS::computeDpDepse_v(const ModelStateBase* state_input) const
{
  const ModelState_Arena* state =
    static_cast<const ModelState_Arena*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Arena.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  double p          = state->pbar_w;
  double dp_depse_v = (p < d_p0) ? d_K0 : d_K0 + d_n * (p - d_p0);
  return dp_depse_v;
}

////////////////////////////////////////////////////////////////////////
/**
 * Function: computeElasticVolumetricStrain
 *
 * Purpose:
 *   Compute the volumetric strain given a pressure (p)
 *
 * Inputs:
 *   pp  = current pressure
 *   p0 = initial pressure
 *
 * Returns:
 *   eps_e_v = current elastic volume strain
 */
////////////////////////////////////////////////////////////////////////
double
WaterEOS::computeElasticVolumetricStrain(const double& pp, const double& p0)
{

  // ASSERT(!(pp < 0))

  // Compute bulk modulus of water
  double Kw0 = computeBulkModulus(p0);
  double Kw  = computeBulkModulus(pp);

  // Compute volume strain
  double eps_e_v0 = (p0 < d_p0) ? 0.0 : -(p0 - d_p0) / Kw0;
  double eps_e_v  = (pp < d_p0) ? 0.0 : -(pp - d_p0) / Kw;
  return (eps_e_v - eps_e_v0);
}

////////////////////////////////////////////////////////////////////////
/**
 * Function: computeExpElasticVolumetricStrain
 *
 * Purpose:
 *   Compute the exponential of volumetric strain given a pressure (p)
 *
 * Inputs:
 *   pp  = current pressure
 *   p0 = initial pressure
 *
 * Returns:
 *   exp(eps_e_v) = exponential of the current elastic volume strain
 */
////////////////////////////////////////////////////////////////////////
double
WaterEOS::computeExpElasticVolumetricStrain(const double& pp, const double& p0)
{
  // ASSERT(!(pp < 0))

  // Compute volume strain
  double eps_e_v = computeElasticVolumetricStrain(pp, p0);
  return std::exp(eps_e_v);
}

////////////////////////////////////////////////////////////////////////
/**
 * Function: computeDerivExpElasticVolumetricStrain
 *
 * Purpose:
 *   Compute the pressure drivative of the exponential of
 *   the volumetric strain at a given pressure (p)
 *
 * Inputs:
 *   pp  = current pressure
 *   p0 = initial pressure
 *
 * Outputs:
 *   exp_eps_e_v = exp(eps_e_v) = exponential of elastic volumeric strain
 *
 * Returns:
 *   deriv = d/dp[exp(eps_e_v)] = derivative of the exponential of
 *                                current elastic volume strain
 */
////////////////////////////////////////////////////////////////////////
double
WaterEOS::computeDerivExpElasticVolumetricStrain(const double& pp,
                                                 const double& p0,
                                                 double& exp_eps_e_v)
{

  // ASSERT(!(pp < 0))

  // Compute the exponential of volumetric strain at pressure (pp)
  exp_eps_e_v = computeExpElasticVolumetricStrain(pp, p0);

  // Compute bulk modulus of water
  double Kw = computeBulkModulus(pp);
  std::cout << "p = " << pp << " Kw = " << Kw << std::endl;

  return -exp_eps_e_v / Kw;
}
