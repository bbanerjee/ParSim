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

#include "MieGruneisenEOS.h"
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <Core/Exceptions/InvalidValue.h>
#include <cmath>

using namespace Uintah;
using namespace Vaango;

MieGruneisenEOS::MieGruneisenEOS(ProblemSpecP& ps)
{
  ps->require("C_0", d_const.C_0);
  ps->require("Gamma_0", d_const.Gamma_0);
  ps->require("S_alpha", d_const.S_alpha);
  ps->require("rho_0", d_const.rho_0);
}

MieGruneisenEOS::MieGruneisenEOS(const MieGruneisenEOS* cm)
{
  d_const.C_0     = cm->d_const.C_0;
  d_const.Gamma_0 = cm->d_const.Gamma_0;
  d_const.S_alpha = cm->d_const.S_alpha;
  d_const.rho_0   = cm->d_const.rho_0;
}

MieGruneisenEOS::~MieGruneisenEOS() = default;

void
MieGruneisenEOS::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP eos_ps = ps->appendChild("equation_of_state");
  eos_ps->setAttribute("type", "mie_gruneisen");

  eos_ps->appendElement("C_0", d_const.C_0);
  eos_ps->appendElement("Gamma_0", d_const.Gamma_0);
  eos_ps->appendElement("S_alpha", d_const.S_alpha);
  eos_ps->appendElement("rho_0", d_const.rho_0);
}

std::map<std::string, double>
MieGruneisenEOS::getParameters() const
{
  std::map<std::string, double> params;
  params["C_0"]     = d_const.C_0;
  params["Gamma_0"] = d_const.Gamma_0;
  params["S_alpha"] = d_const.S_alpha;
  params["rho_0"]   = d_const.rho_0;
  return params;
}

//////////
// Calculate the pressure using the Mie-Gruneisen equation of state
double
MieGruneisenEOS::computePressure(const MPMMaterial* matl,
                                 const ModelStateBase* state,
                                 const Matrix3&,
                                 const Matrix3&,
                                 const double&)
{
  // Get the state data
  double rho = state->density;
  double T   = state->temperature;
  double T_0 = state->initialTemperature;

  // Get original density
  double rho_0 = matl->getInitialDensity();

  // Calc. zeta
  double zeta = (rho / rho_0 - 1.0);

  // Calculate internal energy E
  double E = (state->specificHeat) * (T - T_0) * rho_0;

  // Calculate the pressure
  double p = d_const.Gamma_0 * E;
  if (rho != rho_0) {
    double numer = rho_0 * (d_const.C_0 * d_const.C_0) *
                   (1.0 / zeta + (1.0 - 0.5 * d_const.Gamma_0));
    double denom = 1.0 / zeta - (d_const.S_alpha - 1.0);
    if (denom == 0.0) {
      std::cout << "rh0_0 = " << rho_0 << " zeta = " << zeta
                << " numer = " << numer << "\n";
      denom = 1.0e-5;
    }
    p += numer / (denom * denom);
  }
  return -p;
}

double
MieGruneisenEOS::eval_dp_dJ(const MPMMaterial* matl,
                            const double& detF,
                            [[maybe_unused]] const ModelStateBase* state)
{
  double rho_0   = matl->getInitialDensity();
  double C_0     = d_const.C_0;
  double S_alpha = d_const.S_alpha;
  double Gamma_0 = d_const.Gamma_0;

  double J      = detF;
  double numer  = rho_0 * C_0 * C_0 * (1.0 + (S_alpha - Gamma_0) * (1.0 - J));
  double denom  = (1.0 - S_alpha * (1.0 - J));
  double denom3 = (denom * denom * denom);
  if (denom3 == 0.0) {
    std::cout << "rh0_0 = " << rho_0 << " J = " << J << " numer = " << numer
              << "\n";
    denom3 = 1.0e-5;
  }

  return (numer / denom3);
}

// Compute pressure (option 1)  (tension is +ve)
// The expression for p has been written in terms of J = det(F) = rho0/rho
// instead of zeta
//  p = (C0^2 (2 + Gamma0 (-1 + J)) (-1 + J) rho0)/(2 (1 + (-1 + J) Salpha)^2)
//  dp/dJ = (C0^2 rho0 (1 + Gamma0 (-1 + J) + Salpha - J Salpha))/(1 + (-1 + J)
//  Salpha)^3
//  K = J dp/dJ  (based on Scott, N.H., (2007) Math. and Mech. Solids, 12, pp.
//  526)
//    = dp/dJ  (based on Ogden) is more realistic
//  c^2 = K/rho
double
MieGruneisenEOS::computePressure(const double& rho_orig,
                                 const double& rho_cur) const
{
  // Calc. J
  double J     = rho_orig / rho_cur;
  double J_min = 1.0 - 1.0 / d_const.S_alpha;
  J            = (J > J_min) ? J : J_min;

  /*
  if (J < 1.0 - 1.0 / d_const.S_alpha) {
    throw InvalidValue("**ERROR: EOS Model invalid for extreme compression",
                       __FILE__, __LINE__);
  }
  */

  // Calculate the pressure
  double J_one = J - 1.0;
  double numer = rho_orig * (d_const.C_0 * d_const.C_0) * J_one *
                 (2.0 + d_const.Gamma_0 * J_one);
  double denom = 1.0 + J_one * d_const.S_alpha;
  denom        = (denom == 0.0) ? 1.0e-3 : denom;
  double p     = numer / (2.0 * denom * denom);
  return p;
}

// Compute pressure (option 2)  (tension is +ve)
// The expression for p has been written in terms of J = det(F) = rho0/rho
// instead of zeta
//  p = (C0^2 rho0 (-1 + J) (2 + Gamma0 (-1 + J)))/(2 (1 + (-1 + J) Salpha)^2)
//  dp/dJ = (C0^2 rho0 (1 + Gamma0 (-1 + J) - (-1 + J) Salpha))/(1 + (-1 + J)
//  Salpha)^3
//  dp/drho = (C0^2 J^2 (-1 - (-1 + J) Gamma0 + (-1 + J) Salpha))/(1 + (-1 + J)
//  Salpha)^3
//  K = J dp/dJ  (based on Scott, N.H., (2007) Math. and Mech. Solids, 12, pp.
//  526)
//    = dp/dJ  (based on Ogden) is more realistic
//  c^2 = K/rho
void
MieGruneisenEOS::computePressure(const double& rho_orig,
                                 const double& rho_cur,
                                 double& pressure,
                                 double& dp_drho,
                                 double& csquared)
{
  // Calc. J
  double J     = rho_orig / rho_cur;
  double J_min = 1.0 - 1.0 / d_const.S_alpha;
  J            = (J > J_min) ? J : J_min;

  /*
  if (J < 1.0 - 1.0 / d_const.S_alpha) {
    throw InvalidValue("**ERROR: EOS Model invalid for extreme compression",
                       __FILE__, __LINE__);
  }
  */

  // Calculate the pressure
  double J_one = J - 1.0;
  double c0sq  = d_const.C_0 * d_const.C_0;
  double numer = rho_orig * c0sq * J_one * (2.0 + d_const.Gamma_0 * J_one);
  double denom = 1.0 + J_one * d_const.S_alpha;
  denom        = (denom == 0.0) ? 1.0e-3 : denom;
  pressure     = numer / (2.0 * denom * denom);

  // Calculate dp/dJ and csquared
  numer = c0sq * rho_orig * (1.0 + J_one * (d_const.Gamma_0 - d_const.S_alpha));
  double dp_dJ = numer / (denom * denom * denom);
  // double bulk = J*dp_dJ;
  double bulk = dp_dJ;
  csquared    = bulk / rho_cur;

  // Calculate dp/drho
  dp_drho = -J * J * dp_dJ / rho_orig;

  return;
}

// Compute bulk modulus
double
MieGruneisenEOS::computeInitialBulkModulus() const
{
  return computeBulkModulus(d_const.rho_0, d_const.rho_0);
}

double
MieGruneisenEOS::computeBulkModulus(const double& rho_orig,
                                    const double& rho_cur) const
{
  // Calc. J
  double J     = rho_orig / rho_cur;
  double J_min = 1.0 - 1.0 / d_const.S_alpha;
  J            = (J > J_min) ? J : J_min;

  /*
  if (J < 1.0 - 1.0 / d_const.S_alpha) {
    std::ostringstream err;
    err << "**ERROR** J = " << J << " is less than (1-1/S_alpha) = "
        <<  1.0 - 1.0 / d_const.S_alpha
        << " EOS Model invalid for extreme compression",
    throw InvalidValue(err.str(), __FILE__, __LINE__);
  }
  */

  // Calculate dp/dJ
  double J_one = J - 1.0;
  double c0sq  = d_const.C_0 * d_const.C_0;
  double numer =
    c0sq * rho_orig * (1.0 + J_one * (d_const.Gamma_0 - d_const.S_alpha));
  double denom = 1.0 + J_one * d_const.S_alpha;
  denom        = (denom == 0.0) ? 1.0e-3 : denom;
  double dp_dJ = numer / (denom * denom * denom);
  double bulk = J*dp_dJ;
  //double bulk = dp_dJ;

  return bulk;
}

double
MieGruneisenEOS::computeBulkModulus(const ModelStateBase* state) const
{
  return computeBulkModulus(state->initialDensity, state->density);
}

// Compute strain energy
//  (Integrate the expression for p(J))
//  U = (1/(2 S_alpha^3))C0^2 rho0 (((-1 + J) S_alpha (-2 S_alpha +
//       Gamma0 (2 + (-1 + J) Salpha)))/(1 + (-1 + J) S_alpha) +
//       2 (-Gamma0 + S_alpha) Log[1 + (-1 + J) S_alpha])
//  conditional upon: (1 - J) S_alpha <= 1
// Compute strain energy
double
MieGruneisenEOS::computeStrainEnergy(const ModelStateBase* state)
{
  return computeStrainEnergy(state->initialDensity, state->density);
}

double
MieGruneisenEOS::computeStrainEnergy(const double& rho_orig,
                                     const double& rho_cur)
{
  // Calc. J
  double J    = rho_orig / rho_cur;
  double Sa   = d_const.S_alpha;
  double G0   = d_const.Gamma_0;
  double C0sq = d_const.C_0 * d_const.C_0;

  // Check validity condition
  double J_min = 1.0 - 1.0 / d_const.S_alpha;
  J            = (J > J_min) ? J : J_min;

  /*
  if (J < 1.0 - 1.0 / Sa) {
    throw InvalidValue("**ERROR: EOS Model invalid for extreme compression",
                       __FILE__, __LINE__);
  }
  */

  double J_one = J - 1.0;
  double numer = J_one * Sa * (-2.0 * Sa + G0 * (2.0 + J_one * Sa));
  double denom = 1.0 + J_one * Sa;
  double U     = rho_orig * C0sq / (2.0 * Sa * Sa * Sa) *
             (numer / denom + 2 * (Sa - G0) * log(1.0 + J_one * Sa));

  return U;
}

// Compute density given pressure (use Newton iterations)
double
MieGruneisenEOS::computeDensity(const double& rho_orig, const double& pressure)
{
  double J     = 0.8;
  double J_one = J - 1.0;
  double numer = 0.0;
  double denom = 1.0;
  double c0sq  = d_const.C_0 * d_const.C_0;
  double p     = 0.0;
  double dp_dJ = 0.0;

  double f      = 0.0;
  double fPrime = 0.0;
  int iter      = 0;
  int max_iter  = 100;
  double tol    = 1.0e-6 * pressure;
  do {

    // Calculate p
    J_one = J - 1.0;
    numer = rho_orig * c0sq * J_one * (2.0 + d_const.Gamma_0 * J_one);
    denom = 1.0 + J_one * d_const.S_alpha;
    p     = numer / (2.0 * denom * denom);

    // Calculate dp/dJ
    numer =
      c0sq * rho_orig * (1.0 + J_one * (d_const.Gamma_0 - d_const.S_alpha));
    dp_dJ = numer / (denom * denom * denom);

    // f(J) and f'(J) calc
    f      = p - pressure;
    fPrime = dp_dJ;
    J -= f / fPrime;
    ++iter;
  } while (fabs(f) > tol && iter < max_iter);

  double J_min = 1.0 - 1.0 / d_const.S_alpha;
  J            = (J > J_min) ? J : J_min;
  double rho   = rho_orig / J;

  return rho;
}

double
MieGruneisenEOS::computeDpDepse_v(const Vaango::ModelStateBase*) const
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute dp/deps_v of Mie-Gruneisen material"
         " unless the elastic part of J is provided."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}

double
MieGruneisenEOS::computeDpDepse_s(const Vaango::ModelStateBase*) const
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute dp/deps_s of Mie-Gruneisen material"
         " unless the elastic part of J is provided."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}

double
MieGruneisenEOS::computeElasticVolumetricStrain([[maybe_unused]] const double& pp,
                                                [[maybe_unused]] const double& p0)
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute volume strain of Mie-Gruneisen material"
         " unless the elastic part of J is provided."
         " It should be provided as an input."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}

double
MieGruneisenEOS::computeExpElasticVolumetricStrain([[maybe_unused]] const double& pp,
                                                   [[maybe_unused]] const double& p0)
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute exp(volume strain) of Mie-Gruneisen material"
         " unless the elastic part of J is provided."
         " It should be provided as an input."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}

double
MieGruneisenEOS::computeDerivExpElasticVolumetricStrain([[maybe_unused]] const double& pp,
                                                        [[maybe_unused]] const double& p0,
                                                        [[maybe_unused]] double& exp_eps_e_v)
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute derivative of exp(volume strain) of "
         " Mie-Gruneisen material unless the elastic part of J is provided."
         " It should be provided as an input."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}
