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

#include "MieGruneisenEOSEnergy.h"
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <Core/Exceptions/ConvergenceFailure.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Math/DEIntegrator.h>
#include <cmath>
#include <iostream>

using namespace Uintah;
using namespace Vaango;

MieGruneisenEOSEnergy::MieGruneisenEOSEnergy(ProblemSpecP& ps)
{
  ps->require("C_0", d_const.C_0);
  ps->require("Gamma_0", d_const.Gamma_0);
  ps->require("S_alpha", d_const.S_1);
  ps->getWithDefault("S_2", d_const.S_2, 0.0);
  ps->getWithDefault("S_3", d_const.S_3, 0.0);
  ps->require("rho_0", d_const.rho_0);
  d_J_min = 1.0 +
            (d_const.S_2 + std::sqrt(d_const.S_2 * d_const.S_2 -
                                     3.0 * d_const.S_1 * d_const.S_3)) /
              (3.0 * d_const.S_3);
  // std::cout << "J_min = " << d_J_min << "\n";
}

MieGruneisenEOSEnergy::MieGruneisenEOSEnergy(const MieGruneisenEOSEnergy* cm)
{
  d_const.C_0     = cm->d_const.C_0;
  d_const.Gamma_0 = cm->d_const.Gamma_0;
  d_const.S_1     = cm->d_const.S_1;
  d_const.S_2     = cm->d_const.S_2;
  d_const.S_3     = cm->d_const.S_3;
  d_const.rho_0   = cm->d_const.rho_0;
  d_J_min         = cm->d_J_min;
}

MieGruneisenEOSEnergy::~MieGruneisenEOSEnergy() = default;

void
MieGruneisenEOSEnergy::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP eos_ps = ps->appendChild("equation_of_state");
  eos_ps->setAttribute("type", "mie_gruneisen_energy");

  eos_ps->appendElement("C_0", d_const.C_0);
  eos_ps->appendElement("Gamma_0", d_const.Gamma_0);
  eos_ps->appendElement("S_alpha", d_const.S_1);
  eos_ps->appendElement("S_2", d_const.S_2);
  eos_ps->appendElement("S_3", d_const.S_3);
  eos_ps->appendElement("rho_0", d_const.rho_0);
}

std::map<std::string, double>
MieGruneisenEOSEnergy::getParameters() const
{
  std::map<std::string, double> params;
  params["C_0"]     = d_const.C_0;
  params["Gamma_0"] = d_const.Gamma_0;
  params["S_alpha"] = d_const.S_1;
  params["S_2"]     = d_const.S_2;
  params["S_3"]     = d_const.S_3;
  params["rho_0"]   = d_const.rho_0;
  return params;
}
//////////
// Calculate the pressure (tension +ve) using the Mie-Gruneisen equation of
// state
double
MieGruneisenEOSEnergy::computePressure(const MPMMaterial* matl,
                                       const ModelStateBase* state,
                                       const Matrix3&,
                                       const Matrix3&,
                                       const double&)
{
  // Get original density
  double rho_0 = matl->getInitialDensity();

  // Get the current density
  double rho = state->density;

  // Retrieve specific internal energy e
  double e = state->energy;

  // Constants
  double rho0_C0_sq = rho_0 * d_const.C_0 * d_const.C_0;
  double S_1        = d_const.S_1;
  double S_2        = d_const.S_2;
  double S_3        = d_const.S_3;
  double Gamma_0    = d_const.Gamma_0;

  // Calc. eta
  double J   = rho_0 / rho;
  J          = (J > d_J_min) ? J : d_J_min;
  double eta = 1. - J;

  // Calculate the pressure
  double p = rho_0 * Gamma_0 * e;
  if (eta >= 0.0) {
    double eta2  = eta * eta;
    double eta3  = eta2 * eta;
    double alpha = (1. - S_1 * eta - S_2 * eta2 - S_3 * eta3);
    double denom = alpha * alpha;
    p += rho0_C0_sq * eta * (1. - .5 * Gamma_0 * eta) / denom;
  } else {
    p += rho0_C0_sq * eta;
  }

  return -p;
}

double
MieGruneisenEOSEnergy::computeIsentropicTemperatureRate(const double T,
                                                        const double rho_0,
                                                        const double rho_cur,
                                                        const double Dtrace)
{
  double dTdt = -T * d_const.Gamma_0 * rho_0 * Dtrace / rho_cur;

  return dTdt;
}

double
MieGruneisenEOSEnergy::eval_dp_dJ(const MPMMaterial* matl,
                                  const double& detF,
                                  const ModelStateBase* state)
{
  double rho_0   = matl->getInitialDensity();
  double rho_cur = rho_0 / detF;
  return eval_dp_dJ(rho_0, rho_cur);
}

double
MieGruneisenEOSEnergy::eval_dp_dJ(double rho_0, double rho) const
{
  // Constants
  double rho0_C0_sq = rho_0 * d_const.C_0 * d_const.C_0;
  double S_1        = d_const.S_1;
  double S_2        = d_const.S_2;
  double S_3        = d_const.S_3;
  double Gamma_0    = d_const.Gamma_0;

  // Calc. eta
  double J   = rho_0 / rho;
  J          = (J > d_J_min) ? J : d_J_min;
  double eta = 1. - J;

  // Calculate the pressure
  double dp_dJ = 0.0;
  if (eta >= 0.0) {
    double eta2        = eta * eta;
    double eta3        = eta2 * eta;
    double alpha       = 1. - S_1 * eta - S_2 * eta2 - S_3 * eta3;
    double dalpha_deta = -S_1 - 2.0 * S_2 * eta - 3.0 * S_3 * eta2;
    double alpha2      = alpha * alpha;
    dp_dJ              = -(rho0_C0_sq / alpha2) *
            (1 - Gamma_0 * eta -
             2.0 * eta * (1 - Gamma_0 / 2 * eta) / alpha * dalpha_deta);
  } else {
    dp_dJ = -rho0_C0_sq;
  }
  // std::cout << "J = " << 1 - eta << " dp_dJ = " << dp_dJ << "\n";
  return dp_dJ;
}

// Compute bulk modulus
double
MieGruneisenEOSEnergy::computeInitialBulkModulus() const
{
  return computeBulkModulus(d_const.rho_0, d_const.rho_0);
}

double
MieGruneisenEOSEnergy::computeBulkModulus(const double& rho_orig,
                                          const double& rho_cur) const
{
  return -1.0 * eval_dp_dJ(rho_orig, rho_cur);
}

double
MieGruneisenEOSEnergy::computeBulkModulus(const ModelStateBase* state) const
{
  return computeBulkModulus(state->initialDensity, state->density);
}

// Compute pressure (option 1) - no internal energy contribution
// Compression part:
//  p = -((C0^2 Eta (2 - Eta Gamma0) rho0)/(2 (1 - Eta S1 - Eta^2 S2 - Eta^3
//  S3)^2))
// Tension part:
//  p = -((C0^2 Eta rho0)/(1 - Eta))
double
MieGruneisenEOSEnergy::computePressure(const double& rho_orig,
                                       const double& rho_cur) const
{
  // Calculate J
  double J = rho_orig / rho_cur;
  J        = (J > d_J_min) ? J : d_J_min;

  // Calc. eta = 1 - J  (Note that J = 1 - eta)
  double eta = 1. - J;

  // Calculate the pressure
  double p = 0.0;
  if (eta >= 0.0) {
    double etaMax =
      1.0 - 1.0e-16; // Hardcoded to take care of machine precision issues
    eta = (eta > etaMax) ? etaMax : eta;
    p   = pCompression(rho_orig, eta);
  } else {
    p = pTension(rho_orig, eta);
  }

  return p;
}

// Compute pressure (option 2) - no internal energy contribution
// Compression part:
//  p = -((C0^2 Eta (2 - Eta Gamma0) rho0)/(2 (1 - Eta S1 - Eta^2 S2 - Eta^3
//  S3)^2))
//  dp_dJ = (C0^2 rho0 (1 + S1 - (1 - Eta) S1 + 3 S2 - 6 (1 - Eta) S2 +
//     3 (1 - Eta)^2 S2 + 5 Eta^3 S3 - Eta Gamma0 (1 + Eta^2 S2 + 2 Eta^3 S3)))/
//       (1 - Eta S1 - Eta^2 S2 - Eta^3 S3)^3
//  K = dp_dJ
// Tension part:
//  p = -((C0^2 Eta rho0)/(1 - Eta))
//  dp_dJ = (C0^2 rho0)/(1 - Eta)^2
void
MieGruneisenEOSEnergy::computePressure(const double& rho_orig,
                                       const double& rho_cur,
                                       double& pressure,
                                       double& dp_drho,
                                       double& csquared)
{
  // Calculate J and dJ_drho
  double J       = rho_orig / rho_cur;
  J              = (J > d_J_min) ? J : d_J_min;
  double dJ_drho = -J / rho_cur;
  double dp_dJ   = 0.0;

  // Calc. eta = 1 - J  (Note that J = 1 - eta)
  double eta = 1. - J;

  // Calculate the pressure
  if (eta >= 0.0) {
    double etaMax =
      1.0 - 1.0e-16; // Hardcoded to take care of machine precision issues
    eta      = (eta > etaMax) ? etaMax : eta;
    pressure = pCompression(rho_orig, eta);
    dp_dJ    = dpdJCompression(rho_orig, eta);
  } else {
    pressure = pTension(rho_orig, eta);
    dp_dJ    = dpdJTension(rho_orig, eta);
  }
  dp_drho  = dp_dJ * dJ_drho;
  csquared = dp_dJ / rho_cur;

  if (std::isnan(pressure) || std::abs(dp_dJ) < 1.0e-30) {
     std::ostringstream desc;
    desc << "pressure = " << -pressure << " rho_cur = " << rho_cur
         << " dp_drho = " << -dp_drho << " c^2 = " << csquared << "\n";
    throw InvalidValue(desc.str(), __FILE__, __LINE__);
  }
  return;
}

// Compute density given pressure (tension +ve)
double
MieGruneisenEOSEnergy::computeDensity(const double& rho_orig,
                                      const double& pressure)
{
  double eta  = 0.0;
  double C0sq = d_const.C_0 * d_const.C_0;
  double bulk = C0sq * rho_orig;
  if (std::abs(pressure) < 0.1 * bulk) {
    // Use Newton's method for small pressures (less than 1 GPa hardcoded)
    // (should use p_ref instead for this this work under non-SI units - TO DO)

    if (pressure < 0.0) {
      // Compressive deformations
      const double J0        = 0.8;
      const double tolerance = 1.0e-3;
      const int maxIter      = 10;
      auto pFunc             = &MieGruneisenEOSEnergy::pCompression;
      auto dpdJFunc          = &MieGruneisenEOSEnergy::dpdJCompression;
      eta                    = findEtaNewton(
        pFunc, dpdJFunc, rho_orig, pressure, J0, tolerance, maxIter);
    } else {
      // Tensile deformations
      const double J0        = 1.5;
      const double tolerance = 1.0e-3;
      const int maxIter      = 10;
      auto pFunc             = &MieGruneisenEOSEnergy::pTension;
      auto dpdJFunc          = &MieGruneisenEOSEnergy::dpdJTension;
      eta                    = findEtaNewton(
        pFunc, dpdJFunc, rho_orig, pressure, J0, tolerance, maxIter);
    }
  } else {
    // Use Ridder's method for other pressures
    if (pressure < 0.0) {
      double etamin = 0.0;
      double etamax = 1.0 - 1.0e-16; // Hardcoded for machine precision issues
                                     // Needs to be resolved (TO DO)
      const double tolerance = 1.0e-3;
      const int maxIter      = 100;
      auto pFunc             = &MieGruneisenEOSEnergy::pCompression;
      eta                    = findEtaRidder(
        pFunc, rho_orig, pressure, etamin, etamax, tolerance, maxIter);
    } else {
      double etamin          = -5.0; // Hardcoded: Needs to be resolved (TO DO)
      double etamax          = 0.0;
      const double tolerance = 1.0e-3;
      const int maxIter      = 100;
      auto pFunc             = &MieGruneisenEOSEnergy::pTension;
      eta                    = findEtaRidder(
        pFunc, rho_orig, pressure, etamin, etamax, tolerance, maxIter);
    }
  }
  double J   = 1.0 - eta;
  double rho = rho_orig / J; // **TO DO** Infinity check

  return rho;
}

// Private method: Find root of p(eta) - p0 = 0 using Ridder's method
double
MieGruneisenEOSEnergy::findEtaRidder(pFuncPtr pFunc,
                                     const double& rho_orig,
                                     const double& p0,
                                     double& etamin,
                                     double& etamax,
                                     const double& tolerance,
                                     const int& maxIter) const
{
  double eta =
    1.0 - 1.0e-16; // Hardcoded to take care of machine precision issues
  double pp = (this->*pFunc)(rho_orig, eta);
  // if (0.1*std::abs(pp) < std::abs(p0)) return eta;
  if (0.1 * std::abs(pp) < std::abs(p0)) {
    pp = copysign(0.1 * std::abs(pp), p0);
  } else {
    pp = p0;
  }

  // etamin = 1.0 - Jmax;
  double pmin = (this->*pFunc)(rho_orig, etamin);
  double fmin = pmin - pp;
  if (fmin == 0.0) {
    return etamin;
  }

  // etamax = 1.0 - Jmin;
  double pmax = (this->*pFunc)(rho_orig, etamax);
  double fmax = pmax - pp;
  if (fmax == 0.0) {
    return etamax;
  }

  int count   = 1;
  double fnew = 0.0;
  while (count < maxIter) {

    // compute mid point
    double etamid = 0.5 * (etamin + etamax);
    double pmid   = (this->*pFunc)(rho_orig, etamid);
    double fmid   = pmid - pp;

    double ss = sqrt(fmid * fmid - fmin * fmax);
    if (ss == 0.0) {
      return eta;
    }

    // compute new point
    double dx = (etamid - etamin) * fmid / ss;
    if ((fmin - fmax) < 0.0) {
      dx = -dx;
    }
    double etanew = etamid + dx;
    double pnew   = (this->*pFunc)(rho_orig, etanew);
    fnew          = pnew - pp;

    // Test for convergence
    if (count > 1) {
      // if abs(etanew - eta) < tolerance*max(abs(etanew),1.0)
      if (std::abs(fnew) < tolerance * std::abs(pp)) {
        return etanew;
      }
    }
    eta = etanew;

    // Re-bracket the root as tightly as possible
    if (fmid * fnew > 0.0) {
      if (fmin * fnew < 0.0) {
        etamax = etanew;
        fmax   = fnew;
      } else {
        etamin = etanew;
        fmin   = fnew;
      }
    } else {
      etamin = etamid;
      fmin   = fmid;
      etamax = etanew;
      fmax   = fnew;
    }

    count++;
  }
   std::ostringstream desc;
  desc << "**ERROR** Ridder algorithm did not converge"
       << " pressure = " << p0 << " pp = " << pp << " eta = " << eta << "\n";
  throw ConvergenceFailure(desc.str(),
                           maxIter,
                           std::abs(fnew),
                           tolerance * std::abs(pp),
                           __FILE__,
                           __LINE__);

  return -1;
}

// Private method: Find root of p(eta) - p0 = 0 using Newton's method
double
MieGruneisenEOSEnergy::findEtaNewton(pFuncPtr pFunc,
                                     dpdJFuncPtr dpdJFunc,
                                     const double& rho_orig,
                                     const double& p0,
                                     const double& J0,
                                     const double& tolerance,
                                     const int& maxIter) const
{
  double p     = 0.0;
  double dp_dJ = 0.0;
  double J     = J0;
  double eta   = 1.0 - J;

  double f      = 0.0;
  double fPrime = 0.0;
  int iter      = 0;

  do {

    // Calculate p
    p = (this->*pFunc)(rho_orig, eta);

    // Calculate dp/dJ
    dp_dJ = (this->*dpdJFunc)(rho_orig, eta);

    // f(J) and f'(J) calc
    f      = p - p0;
    fPrime = dp_dJ;
    J -= f / fPrime;

    // Update eta
    eta = 1.0 - J;

    ++iter;
  } while (std::abs(f) > tolerance && iter < maxIter);

  if (iter >= maxIter) {
     std::ostringstream desc;
    desc << "**ERROR** Newton algorithm did not converge"
         << " pressure = " << p0 << " eta = " << eta << "\n";
    throw ConvergenceFailure(
      desc.str(), maxIter, std::abs(f), tolerance, __FILE__, __LINE__);
  }

  return eta;
}

// Private method: Compute p for compressive volumetric deformations
double
MieGruneisenEOSEnergy::pCompression(const double& rho_orig,
                                    const double& eta) const
{
  // Calc eta^2 and eta^3
  double etaSq = eta * eta;
  double etaCb = eta * eta * eta;

  // Calculate p
  double numer = rho_orig * d_const.C_0 * d_const.C_0 * eta *
                 (1.0 - 0.5 * d_const.Gamma_0 * eta);
  double denom =
    1.0 - d_const.S_1 * eta - d_const.S_2 * etaSq - d_const.S_3 * etaCb;
  double p = -numer / (denom * denom);

  return p;
}

// Private method: Compute dp/dJ for compressive volumetric deformations
double
MieGruneisenEOSEnergy::dpdJCompression(const double& rho_orig,
                                       const double& eta) const
{
  // Calc eta^2 and eta^3
  double etaSq = eta * eta;
  double etaCb = eta * eta * eta;

  // Calculate dp/dJ
  double J     = 1 - eta;
  double numer = 1.0 + eta * d_const.S_1 +
                 3.0 * (1.0 - 2.0 * J + J * J) * d_const.S_2 +
                 5.0 * etaCb * d_const.S_3 -
                 eta * d_const.Gamma_0 *
                   (1.0 + etaSq * d_const.S_2 + 2.0 * etaCb * d_const.S_3);
  double denom =
    1.0 - d_const.S_1 * eta - d_const.S_2 * etaSq - d_const.S_3 * etaCb;
  double dp_dJ =
    d_const.C_0 * d_const.C_0 * rho_orig * numer / (denom * denom * denom);

  return dp_dJ;
}

// Private method: Compute p for tensile volumetric deformations
double
MieGruneisenEOSEnergy::pTension(const double& rho_orig, const double& eta) const
{
  // Calculate p
  double p = -rho_orig * d_const.C_0 * d_const.C_0 * eta / (1.0 - eta);

  return p;
}

// Private method: Compute dp/dJ for tensile volumetric deformations
double
MieGruneisenEOSEnergy::dpdJTension(const double& rho_orig,
                                   const double& eta) const
{
  // Calculate dp/dJ
  double J     = 1 - eta;
  double dp_dJ = (rho_orig * d_const.C_0 * d_const.C_0) / (J * J);

  return dp_dJ;
}

// Compute strain energy
//   An exact integral does not exist and numerical integration is needed
//   Use double exponential integration because the function is smooth
//   (even though we use it only in a limited region)
//   **WARNING** Requires well behaved EOS that does not blow up to
//               infinity in the middle of the domain
double
MieGruneisenEOSEnergy::computeStrainEnergy(const ModelStateBase* state)
{
  return computeStrainEnergy(state->initialDensity, state->density);
}

double
MieGruneisenEOSEnergy::computeStrainEnergy(const double& rho_orig,
                                           const double& rho_cur)
{
  // Calculate J
  double J = rho_orig / rho_cur;

  // Calc. eta = 1 - J  (Note that J = 1 - eta)
  double eta = 1. - J;

  // Calculate the pressure
  double U    = 0.0;
  double C0sq = d_const.C_0 * d_const.C_0;
  if (eta >= 0.0) {
    int evals;
    double error;
    U = DEIntegrator<MieGruneisenEOSEnergy>::Integrate(
      this, 0, eta, 1.0e-6, evals, error);
    U *= rho_orig * C0sq;
  } else {
    U = C0sq * rho_orig * (J - 1.0 - log(J));
  }
  return U;
}

// Special operator for computing energy
double
MieGruneisenEOSEnergy::operator()(double eta) const
{
  // Calculate the pressure
  double etaSq = eta * eta;
  double etaCb = eta * eta * eta;
  double numer = eta * (1.0 - 0.5 * d_const.Gamma_0 * eta);
  double denom =
    1.0 - d_const.S_1 * eta - d_const.S_2 * etaSq - d_const.S_3 * etaCb;

  return numer / (denom * denom);
}

double
MieGruneisenEOSEnergy::computeDpDepse_v(const Vaango::ModelStateBase*) const
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute dp/deps_v of 5 parameter Mie-Gruneisen "
         "material"
         " unless the elastic part of J is provided."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}
double
MieGruneisenEOSEnergy::computeDpDepse_s(const Vaango::ModelStateBase*) const
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute dp/deps_s of 5 parameter Mie-Gruneisen "
         "material"
         " unless the elastic part of J is provided."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}
double
MieGruneisenEOSEnergy::computeElasticVolumetricStrain(const double& pp,
                                                      const double& p0)
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute volume strain of 5 parameter Mie-Gruneisen "
         "material."
         " It should be provided as an input."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}

double
MieGruneisenEOSEnergy::computeExpElasticVolumetricStrain(const double& pp,
                                                         const double& p0)
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute exp(volume strain) of 5 parameter "
         "Mie-Gruneisen material."
         " It should be provided as an input."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}
double
MieGruneisenEOSEnergy::computeDerivExpElasticVolumetricStrain(
  const double& pp,
  const double& p0,
  double& exp_eps_e_v)
{
  std::ostringstream err;
  err
    << "**ERROR** Cannot compute derivative of exp(volume strain) of "
       " 5 parameter Mie-Gruneisen material. It should be provided as an input."
       " Please change the equation_of_state if you need this "
       " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}
