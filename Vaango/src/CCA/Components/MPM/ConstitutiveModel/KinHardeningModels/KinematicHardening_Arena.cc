/*
 * The MIT License
 *
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

#include <CCA/Components/MPM/ConstitutiveModel/KinHardeningModels/KinematicHardening_Arena.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Arena.h>
#include <Core/Exceptions/InternalError.h>
#include <cmath>

//#define NUMERICALLY_INTEGRATE_BACKSTRESS

using namespace Uintah;
using namespace Vaango;

const Matrix3 KinematicHardening_Arena::Identity(1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                                                 0.0, 0.0, 1.0);

KinematicHardening_Arena::KinematicHardening_Arena(
  ProblemSpecP& ps, InternalVariableModel* intvar)
{
  // Copy of internal variables object
  d_intvar = intvar;

  // Initial fluid pressure
  ps->require("fluid_pressure_initial", d_cm.fluid_pressure_initial);

  // Local labels
  initializeLocalMPMLabels();
}

KinematicHardening_Arena::KinematicHardening_Arena(
  const KinematicHardening_Arena* cm)
{
  // Copy of internal variables object
  d_intvar = cm->d_intvar;

  // Initial fluid pressure
  d_cm.fluid_pressure_initial = cm->d_cm.fluid_pressure_initial;

  // Local labels
  initializeLocalMPMLabels();
}

KinematicHardening_Arena::~KinematicHardening_Arena()
{
  // VarLabel::destroy(pZetaLabel);
  // VarLabel::destroy(pZetaLabel_preReloc);
}

void
KinematicHardening_Arena::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP plastic_ps = ps->appendChild("kinematic_hardening_model");
  plastic_ps->setAttribute("type", "arena_pore_pressure");
  plastic_ps->appendElement("fluid_pressure_initial",
                            d_cm.fluid_pressure_initial);
}

void
KinematicHardening_Arena::computeBackStress(const ModelStateBase* state_input,
                                            Matrix3& backStress_new)
{
  const ModelState_Arena* state =
    static_cast<const ModelState_Arena*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Arena.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  // If the state is tensile the back stress does not change. Return old
  // backtress.
  if (state->I1_eff > 0.0) {
    backStress_new = Identity * (state->pbar_w);
    return;
  }

  // Get the variables of interest
  double p0 = d_cm.fluid_pressure_initial;
  double phi0 = state->phi0;
  double Sw0 = state->Sw0;
  double phat_old = -state->pbar_w;
  double ep_v_bar_new = -state->ep_v;

#ifndef NUMERICALLY_INTEGRATE_BACKSTRESS

  double phat_new = phat_old;
  if (phat_old < 1.0e11) {
    phat_new = computePressureUnloaded(ep_v_bar_new, Sw0, phi0, p0, phat_old);
  } else {
    phat_new = computePressureBisection(ep_v_bar_new, Sw0, phi0, p0, phat_old);
  }

// double ep_v_bar_old  = ep_v_bar_new - dep_v_bar;
// std::cout << "\t Backstress: "
//          << " ep_v_old = " << ep_v_bar_old
//          << " delta ep_v = " << dep_v_bar
//          << " ep_v_new = " << ep_v_bar_new
//          << " phat_old = " << phat_old
//          << " phat_new = " << phat_new << std::endl;

#else
  // Compute volumetric strains in air, water, and matrix material at p = pbar_w
  double exp_ev = 0.0;
  double dexp_ev_air =
    d_air.computeDerivExpElasticVolumetricStrain(phat_old, 0.0, exp_ev);
  double dexp_ev_water =
    d_water.computeDerivExpElasticVolumetricStrain(phat_old, p0, exp_ev);
  double dexp_ev_matrix =
    d_granite.computeDerivExpElasticVolumetricStrain(phat_old, 0.0, exp_ev);

  // Compute denominator of rate equation
  double BB = phi0 * ((1.0 - Sw0) * dexp_ev_air + Sw0 * dexp_ev_water) +
              (1 - phi0) * dexp_ev_matrix;

  // Compute derivative of pbar_w wrt plastic vol strain
  double dep_v_bar = -state->dep_v;
  double ep_v_bar_old = ep_v_bar_new - dep_v_bar;
  double dpbar_w_dep_v = std::exp(-ep_v_bar_old) / (3.0 * BB);

  // Compute new pbar_w using forward Euler
  double phat_new = phat_old + dpbar_w_dep_v * dep_v_bar;
// std::cout << "\t Backstress: "
//          << " air = " << dexp_ev_air
//          << " water = " << dexp_ev_water
//          << " matrix = " << dexp_ev_matrix << std::endl
//          << "\t\t BB = " << BB
//          << " dpbar_w/dep_v = " << dpbar_w_dep_v
//          << " delta ep_v = " << dep_v_bar
//          << " ep_v = " << ep_v_bar_old
//          << " phat_old = " << phat_old
//          << " delta pbar_w = " << dpbar_w_dep_v*dep_v_bar
//          << " phat_new = " << phat_new << std::endl;
#endif

  // Compute backstress tensor
  backStress_new = -phat_new * Identity;

  return;
}

//------------------------------------------------------
// Newton solve for pressure
//------------------------------------------------------
double
KinematicHardening_Arena::computePressureUnloaded(const double& ev_p,
                                                  const double& S0,
                                                  const double& phi0,
                                                  const double& p0,
                                                  const double& p_init)
{
  int maxIter = 100;
  double tol = 1.0e-6;
  // double p = p0;
  double p_old = p_init;
  double p = p_init;
  int k = 1;
  while (k < maxIter) {
    double g_Dg = computeGpByDgp(p, ev_p, S0, phi0, p0);
    std::cout << "p = " << p << " ev_p = " << ev_p << " S0 = " << S0
              << " phi0 = " << phi0 << " p0 = " << p0 << " g_Dg = " << g_Dg
              << std::endl;
    p_old = p;
    p -= g_Dg;
    k++;
    if ((std::abs((p_old - p) / p) < tol) || (std::abs(g_Dg) < tol)) {
      // if (std::abs(g_Dg/p) < tol) {
      break;
    }
  }

  if (k >= maxIter) {
    // std::ostringstream out;
    // out << "**ERROR** Too many iterations needed to find backstress.";
    // throw InternalError(out.str(), __FILE__, __LINE__);

    // Instead of throwing an exception try bisection first
    p = computePressureBisection(ev_p, S0, phi0, p0, p_init);
  }

  return p;
}

double
KinematicHardening_Arena::computeGpByDgp(const double& p, const double& ev_p,
                                         const double& S0, const double& phi0,
                                         const double& p0)
{

  // Compute volumetric strains in air, water, and matrix material at p = pbar_w
  double exp_ev_air = 0.0;
  double exp_ev_water = 0.0;
  double exp_ev_matrix = 0.0;
  double dexp_ev_air =
    d_air.computeDerivExpElasticVolumetricStrain(p, 0.0, exp_ev_air);
  double dexp_ev_water =
    d_water.computeDerivExpElasticVolumetricStrain(p, p0, exp_ev_water);
  double dexp_ev_matrix =
    d_granite.computeDerivExpElasticVolumetricStrain(p, 0.0, exp_ev_matrix);

  double air_term = (1.0 - S0) * phi0;
  double water_term = S0 * phi0;
  double solid_term = (1.0 - phi0);

  std::cout << "In computeGpByDgp:" << std::endl;
  std::cout << " exp_ev_air = " << exp_ev_air
            << " dexp_ev_air = " << dexp_ev_air << " air_term = " << air_term
            << " air Gp = " << air_term * exp_ev_air << std::endl;
  std::cout << " exp_ev_water = " << exp_ev_water
            << " dexp_ev_water = " << dexp_ev_water
            << " water_term = " << water_term
            << " water Gp = " << water_term * exp_ev_water << std::endl;
  std::cout << " exp_ev_matrix = " << exp_ev_matrix
            << " dexp_ev_matrix = " << dexp_ev_matrix
            << " solid_term = " << solid_term
            << " solid Gp = " << solid_term * exp_ev_matrix << std::endl;
  std::cout << " ev_p = " << ev_p << " rhs = "
            << std::log(air_term * exp_ev_air + water_term * exp_ev_water +
                        solid_term * exp_ev_matrix)
            << std::endl;

  double Gp = air_term * exp_ev_air + water_term * exp_ev_water +
              solid_term * exp_ev_matrix - std::exp(-ev_p);

  double dGp = air_term * dexp_ev_air + water_term * dexp_ev_water +
               solid_term * dexp_ev_matrix;
  std::cout << " Gp = " << Gp << " dGp = " << dGp << std::endl;

  double Gp_by_dGp = Gp / dGp;
  if (!std::isfinite(Gp_by_dGp)) {
    std::ostringstream out;
    out << "**ERROR** Increment in Newton iteration is not finite."
        << " p = " << p << " Gp = " << Gp << " dGp = " << dGp << std::endl;
    throw InternalError(out.str(), __FILE__, __LINE__);
  }

  return Gp_by_dGp;
}

//------------------------------------------------------
// Bisection solve for pressure
//------------------------------------------------------
double
KinematicHardening_Arena::computePressureBisection(const double& ev_p,
                                                   const double& S0,
                                                   const double& phi0,
                                                   const double& p0,
                                                   const double& p_init)
{
  // Set the initial range
  double low = 0.1 * p_init;
  double high = 10.0 * p_init;

  // Set tol and maxiter
  int maxIter = 100;
  double tol = 1.0e-6;

  // Start iterating
  int k = 1;
  while (k < maxIter) {

    // Get midpoint
    double mid = 0.5 * (low + high);

    // Compute function
    double Gp = computeGp(mid, ev_p, S0, phi0, p0);

    // Compute high - low
    double range = high - low;

    std::cout << "k = " << k << " Gp = " << Gp << " low = " << low
              << " mid = " << mid << " high = " << high << std::endl;

    // Solution found
    if (std::abs(Gp) < tol || std::abs(range / high) < tol) {
      return mid;
    }

    // Continue iterating
    k++;
    double Glow = computeGp(low, ev_p, S0, phi0, p0);
    if (std::signbit(Gp) == std::signbit(Glow)) {
      low = mid;
    } else {
      high = mid;
    }
  }

  // Failure
  std::ostringstream out;
  out << "**ERROR** Too many bisection iterations needed to find backstress.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return 0;
}

double
KinematicHardening_Arena::computeGp(const double& p, const double& ev_p,
                                    const double& S0, const double& phi0,
                                    const double& p0)
{
  // Compute volumetric strains in air, water, and matrix material at p = pbar_w
  double exp_ev_air = 0.0;
  double exp_ev_water = 0.0;
  double exp_ev_matrix = 0.0;
  double dexp_ev_air =
    d_air.computeDerivExpElasticVolumetricStrain(p, 0.0, exp_ev_air);
  double dexp_ev_water =
    d_water.computeDerivExpElasticVolumetricStrain(p, p0, exp_ev_water);
  double dexp_ev_matrix =
    d_granite.computeDerivExpElasticVolumetricStrain(p, 0.0, exp_ev_matrix);

  double air_term = (1.0 - S0) * phi0;
  double water_term = S0 * phi0;
  double solid_term = (1.0 - phi0);

  std::cout << "In computeGp:" << std::endl;
  std::cout << " exp_ev_air = " << exp_ev_air
            << " dexp_ev_air = " << dexp_ev_air << " air_term = " << air_term
            << " air Gp = " << air_term * exp_ev_air << std::endl;
  std::cout << " exp_ev_water = " << exp_ev_water
            << " dexp_ev_water = " << dexp_ev_water
            << " water_term = " << water_term
            << " water Gp = " << water_term * exp_ev_water << std::endl;
  std::cout << " exp_ev_matrix = " << exp_ev_matrix
            << " dexp_ev_matrix = " << dexp_ev_matrix
            << " solid_term = " << solid_term
            << " solid Gp = " << solid_term * exp_ev_matrix << std::endl;
  std::cout << " ev_p = " << ev_p << " rhs = "
            << std::log(air_term * exp_ev_air + water_term * exp_ev_water +
                        solid_term * exp_ev_matrix)
            << std::endl;

  double Gp = air_term * exp_ev_air + water_term * exp_ev_water +
              solid_term * exp_ev_matrix - std::exp(-ev_p);

  if (!std::isfinite(Gp)) {
    std::ostringstream out;
    out << "**ERROR** Function evaluation in bisection iteration is not finite."
        << " p = " << p << " Gp = " << Gp << std::endl;
    throw InternalError(out.str(), __FILE__, __LINE__);
  }

  return Gp;
}
