/*
 * The MIT License
 *
 * Copyright (c) 2015-2017 Parresia Research Limited, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/Models/YieldCond_Tabular.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <chrono>
#include <cmath>

#define USE_GEOMETRIC_BISECTION
//#define USE_ALGEBRAIC_BISECTION
//#define DEBUG_YIELD_BISECTION
//#define DEBUG_YIELD_BISECTION_I1_J2
//#define DEBUG_YIELD_BISECTION_R
//#define CHECK_FOR_NANS

using namespace Vaango;

const double YieldCond_Tabular::sqrt_two = std::sqrt(2.0);
const double YieldCond_Tabular::sqrt_three = std::sqrt(3.0);
const double YieldCond_Tabular::one_sqrt_three = 1.0 / sqrt_three;

YieldCond_Tabular::YieldCond_Tabular(Uintah::ProblemSpecP& ps)
  : d_yield(ps)
{
  // Check the input parameters
  checkInputParameters();
}

void
YieldCond_Tabular::checkInputParameters()
{
}

YieldCond_Tabular::YieldCond_Tabular(const YieldCond_Tabular* yc)
{
  d_yield = yc->d_yield;
}

YieldCond_Tabular::~YieldCond_Tabular()
{
}

void
YieldCond_Tabular::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP yield_ps = ps->appendChild("plastic_yield_condition");
  yield_ps->setAttribute("type", "tabular");

  d_yield.table.outputProblemSpec(yield_ps);
}

//--------------------------------------------------------------
// Evaluate yield condition
//
// f := J2 - g(p) = 0
// where
//     J2 = 1/2 s:s,  s = sigma - p I,  p = 1/3 Tr(sigma)
//     kappa = I1_peak - CR*(I1_peak - X_eff)
//     g(p) = table
//
// Returns:
//   hasYielded = -1.0 (if elastic)
//              =  1.0 (otherwise)
//--------------------------------------------------------------
double
YieldCond_Tabular::evalYieldCondition(const ModelStateBase* state_input)
{
  const ModelState_Tabular* state =
    dynamic_cast<const ModelState_Tabular*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  // Initialize hasYielded to -1
  double hasYielded = -1.0;

  // Cauchy stress invariants: I1 = 3*p, J2 = q^2/3
  double I1 = state->I1;
  double sqrt_J2 = state->sqrt_J2;

  double pp = -I1/3;
  DoubleVec1D gg = d_yield.table.interpolate<1>({{pp}}); 

  if (sqrt_2 > gg) {
    hasYielded = 1.0;
  }

  return hasYielded;
}

//--------------------------------------------------------------
// Derivatives needed by return algorithms and Newton iterations

//--------------------------------------------------------------
// Evaluate yield condition max  value of sqrtJ2
//--------------------------------------------------------------
double
YieldCond_Tabular::evalYieldConditionMax(const ModelStateBase* state_input)
{
  const ModelState_Tabular* state =
    dynamic_cast<const ModelState_Tabular*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  // Number of points
  int num_points = 10;

  // Set up I1 values
  double rad = 0.5 * (d_local.PEAKI1 - X_eff);
  double cen = 0.5 * (d_local.PEAKI1 + X_eff);
  double theta_min = 0.0;
  double theta_max = M_PI;
  std::vector<double> theta_vec;
  linspace(theta_min, theta_max, num_points, theta_vec);
  double J2_max = std::numeric_limits<double>::min();
  // for (auto I1_eff : I1_eff_vec) {
  for (auto theta : theta_vec) {

    double I1 = cen + rad * std::cos(theta);

    double pp = -I1/3;
    DoubleVec1D gg = d_yield.table.interpolate<1>({{pp}}); 

    // Compute J2
    J2_max = std::max(J2_max, gg);
  }

  return std::sqrt(J2_max);
}

//--------------------------------------------------------------
/*! Compute Derivative with respect to the Cauchy stress (\f$\sigma \f$)
 *  Compute df/dsigma
 *
 *  for the yield function
 *      f := J2 - Ff^2*Fc^2 = 0
 *  where
 *      J2 = 1/2 s:s,  s = sigma - p I,  p = 1/3 Tr(sigma)
 *      I1_eff = 3*(p + pbar_w)
 *      X_eff  = X + 3*pbar_w
 *      kappa = I1_peak - CR*(I1_peak - X_eff)
 *      Ff := a1 - a3*exp(a2*I1_eff) - a4*I1_eff
 *      Fc^2 := 1 - (kappa - I1_eff)^2/(kappa - X_eff)^2
 *
 *  The derivative is
 *      df/dsigma = df/dp dp/dsigma + df/ds : ds/dsigma
 *
 *  where
 *      df/dp = computeVolStressDerivOfYieldFunction
 *      dp/dsigma = 1/3 I
 *  and
 *      df/ds = df/dJ2 dJ2/ds
 *      df/dJ2 = computeDevStressDerivOfYieldFunction
 *      dJ2/ds = s
 *      ds/dsigma = I(4s) - 1/3 II
 *  which means
 *      df/dp dp/dsigma = 1/3 df/dp I
 *      df/ds : ds/dsigma = df/dJ2 s : [I(4s) - 1/3 II]
 *                        = df/dJ2 s
*/
void
YieldCond_Tabular::eval_df_dsigma(const Matrix3&,
                                  const ModelStateBase* state_input,
                                  Matrix3& df_dsigma)
{
  const ModelState_Tabular* state =
    dynamic_cast<const ModelState_Tabular*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  double df_dp = computeVolStressDerivOfYieldFunction(state_input);
  double df_dJ2 = computeDevStressDerivOfYieldFunction(state_input);

  Matrix3 One;
  One.Identity();
  Matrix3 p_term = One * (df_dp / 3.0);
  Matrix3 s_term = state->deviatoricStressTensor * (df_dJ2);

  df_dsigma = p_term + s_term;
  // df_dsigma /= df_dsigma.Norm();

  return;
}

//--------------------------------------------------------------
// Compute df/dp  where pI = volumetric stress = 1/3 Tr(sigma) I
//   df/dp = derivative of the yield function wrt p
//
// for the yield function
//     f := J2 - Ff^2*Fc^2 = 0
// where
//     J2 = 1/2 s:s,  s = sigma - p I,  p = 1/3 Tr(sigma)
//     I1_eff = 3*(p + pbar_w)
//     X_eff  = X + 3*pbar_w
//     kappa = I1_peak - CR*(I1_peak - X_eff)
//     Ff := a1 - a3*exp(a2*I1_eff) - a4*I1_eff
//     Fc^2 := 1 - (kappa - I1_eff)^2/(kappa - X_eff)^2
//
// the derivative is
//     df/dp = -2 Ff Fc^2 dFf/dp - Ff^2 dFc^2/dp
// where
//     dFf/dp = dFf/dI1_eff dI1_eff/dp
//            = -[a2 a3 exp(a2 I1_eff) + a4] dI1_eff/dp
//     dFc^2/dp = dFc^2/dI1_eff dI1_eff/dp
//            = 2 (kappa - I1_eff)/(kappa - X_eff)^2  dI1_eff/dp
// and
//    dI1_eff/dp = 1/3
//--------------------------------------------------------------
double
YieldCond_Tabular::computeVolStressDerivOfYieldFunction(
  const ModelStateBase* state_input)
{
  const ModelState_Tabular* state =
    dynamic_cast<const ModelState_Tabular*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  // Get the particle specific internal variables from the model state
  double PEAKI1 = state->yieldParams.at("PEAKI1");
  double FSLOPE = state->yieldParams.at("FSLOPE");
  double STREN = state->yieldParams.at("STREN");
  double YSLOPE = state->yieldParams.at("YSLOPE");
  double CR = state->yieldParams.at("CR");

  std::vector<double> limitParameters =
    computeModelParameters(PEAKI1, FSLOPE, STREN, YSLOPE);
  double a1 = limitParameters[0];
  double a2 = limitParameters[1];
  double a3 = limitParameters[2];
  double a4 = limitParameters[3];

  // Get the plastic internal variables from the model state
  double X_eff = state->capX + 3.0 * state->pbar_w;
  double kappa = state->kappa;

  // Cauchy stress invariants: I1 = 3*p, J2 = q^2/3
  double I1_eff = state->I1_eff;

  // --------------------------------------------------------------------
  // *** SHEAR LIMIT FUNCTION (Ff) ***
  // --------------------------------------------------------------------
  double Ff = a1 - a3 * exp(a2 * I1_eff) - a4 * I1_eff;

  // --------------------------------------------------------------------
  // *** Branch Point (Kappa) ***
  // --------------------------------------------------------------------
  kappa = PEAKI1 - CR * (PEAKI1 - X_eff); // Branch Point

  // --------------------------------------------------------------------
  // **Elliptical Cap Function: (fc)**
  // --------------------------------------------------------------------
  double kappa_I1_eff = kappa - I1_eff;
  double kappa_X_eff = kappa - X_eff;
  double kappaRatio = kappa_I1_eff / kappa_X_eff;
  double Fc_sq = 1.0 - kappaRatio * kappaRatio;

  // --------------------------------------------------------------------
  // Derivatives
  // --------------------------------------------------------------------
  // dI1_eff/dp = 1/3
  double dI1_eff_dp = 1.0 / 3.0;

  // dFf/dp = dFf/dI1_eff dI1_eff/dp
  //        = -[a2 a3 exp(a2 I1_eff) + a4] dI1_eff/dp
  double dFf_dp = -(a2 * a3 * std::exp(a2 * I1_eff) + a4) * dI1_eff_dp;

  // dFc^2/dp = dFc^2/dI1_eff dI1_eff/dp
  //        = 2 (kappa - I1_eff)/(kappa - X_eff)^2  dI1_eff/dp
  double dFc_sq_dp =
    (2.0 * kappa_I1_eff / (kappa_X_eff * kappa_X_eff)) * dI1_eff_dp;

  // df/dp = -2 Ff Fc^2 dFf/dp - 2 Ff^2 dFc^2/dp
  //       = -2 Ff (Fc^2 dFf/dp + Ff dFc^2/dp)
  double df_dp = -Ff * (2.0 * Fc_sq * dFf_dp + Ff * dFc_sq_dp);

  return df_dp;
}

//--------------------------------------------------------------
// Compute df/dJ2  where J2 = 1/2 s:s ,  s = sigma - p I,  p = 1/3 Tr(sigma)
//   s = derivatoric stress
//   df/dJ2 = derivative of the yield function wrt J2
//
// for the yield function
//     f := J2 - Ff^2*Fc^2 = 0
// where
//     J2 = 1/2 s:s,  s = sigma - p I,  p = 1/3 Tr(sigma)
//     I1_eff = 3*(p + pbar_w)
//     X_eff  = X + 3*pbar_w
//     kappa = I1_peak - CR*(I1_peak - X_eff)
//     Ff := a1 - a3*exp(a2*I1_eff) - a4*I1_eff
//     Fc^2 := 1 - (kappa - I1_eff)^2/(kappa - X_eff)^2
//
// the derivative is
//     df/dJ2 = 1
//--------------------------------------------------------------
double
YieldCond_Tabular::computeDevStressDerivOfYieldFunction(
  const ModelStateBase* state_input)
{
  const ModelState_Tabular* state =
    dynamic_cast<const ModelState_Tabular*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  return 1.0;
}

/**
 * Function: getInternalPoint
 *
 * Purpose: Get a point that is inside the yield surface
 *
 * Inputs:
 *  state = state at the current time
 *
 * Returns:
 *   I1 = value of tr(stress) at a point inside the yield surface
 */
double
YieldCond_Tabular::getInternalPoint(const ModelStateBase* state_old_input,
                                    const ModelStateBase* state_trial_input)
{
  const ModelState_Tabular* state_old =
    dynamic_cast<const ModelState_Tabular*>(state_old_input);
  const ModelState_Tabular* state_trial =
    dynamic_cast<const ModelState_Tabular*>(state_trial_input);
  if ((!state_old) || (!state_trial)) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  // Compute effective trial stress
  double I1_eff_trial =
    state_trial->I1_eff - state_trial->pbar_w + state_old->pbar_w;

  // Get the particle specific internal variables from the model state
  double PEAKI1 = state_old->yieldParams.at("PEAKI1");

  // It may be better to use an interior point at the center of the yield
  // surface, rather than at
  // pbar_w, in particular when PEAKI1=0.  Picking the midpoint between PEAKI1
  // and X would be
  // problematic when the user has specified some no porosity condition (e.g.
  // p0=-1e99)
  double I1_eff_interior = 0.0;
  double upperI1 = PEAKI1;
  if (I1_eff_trial < upperI1) {
    if (I1_eff_trial >
        state_old->capX +
          3.0 * state_old->pbar_w) { // Trial is above yield surface
      I1_eff_interior = state_trial->I1_eff;
    } else { // Trial is past X, use yield midpoint as interior point
      I1_eff_interior =
        -3.0 * state_old->pbar_w +
        0.5 * (PEAKI1 + state_old->capX + 3.0 * state_old->pbar_w);
    }
  } else { // I1_trial + pbar_w >= I1_peak => Trial is past vertex
    double lTrial = sqrt(I1_eff_trial * I1_eff_trial +
                         state_trial->sqrt_J2 * state_trial->sqrt_J2);
    double lYield = 0.5 * (PEAKI1 - state_old->capX - 3.0 * state_old->pbar_w);
    I1_eff_interior =
      -3.0 * state_old->pbar_w + upperI1 - std::min(lTrial, lYield);
  }

  return I1_eff_interior;
}

/**
 * Function: getClosestPoint
 *
 * Purpose: Get the point on the yield surface that is closest to a given point
 * (2D)
 *
 * Inputs:
 *  state = current state
 *  px = x-coordinate of point
 *  py = y-coordinate of point
 *
 * Outputs:
 *  cpx = x-coordinate of closest point on yield surface
 *  cpy = y-coordinate of closest point
 *
 */
bool
YieldCond_Tabular::getClosestPoint(const ModelStateBase* state_input,
                                   const double& px, const double& py,
                                   double& cpx, double& cpy)
{
  const ModelState_Tabular* state =
    dynamic_cast<const ModelState_Tabular*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

#ifdef USE_GEOMETRIC_BISECTION
  // std::chrono::time_point<std::chrono::system_clock> start, end;
  // start = std::chrono::system_clock::now();
  Point pt(px, py, 0.0);
  Point closest(0.0, 0.0, 0.0);
  getClosestPointGeometricBisect(state, pt, closest);
  cpx = closest.x();
  cpy = closest.y();
// end = std::chrono::system_clock::now();
// std::cout << "Geomeric Bisection : Time taken = " <<
//    std::chrono::duration<double>(end-start).count() << std::endl;
#else
  // std::chrono::time_point<std::chrono::system_clock> start, end;
  // start = std::chrono::system_clock::now();
  Point pt(px, py, 0.0);
  Point closest(0.0, 0.0, 0.0);
  getClosestPointAlgebraicBisect(state, pt, closest);
  cpx = closest.x();
  cpy = closest.y();
// end = std::chrono::system_clock::now();
// std::cout << "Algebraic Bisection : Time taken = " <<
//    std::chrono::duration<double>(end-start).count() << std::endl;
#endif

  return true;
}

void
YieldCond_Tabular::getClosestPointGeometricBisect(
  const ModelState_Tabular* state, const Uintah::Point& z_r_pt,
  Uintah::Point& z_r_closest)
{
  // Get the particle specific internal variables from the model state
  // Store in a local struct
  d_local.PEAKI1 = state->yieldParams.at("PEAKI1");
  d_local.FSLOPE = state->yieldParams.at("FSLOPE");
  d_local.STREN = state->yieldParams.at("STREN");
  d_local.YSLOPE = state->yieldParams.at("YSLOPE");
  d_local.BETA = state->yieldParams.at("BETA");
  d_local.CR = state->yieldParams.at("CR");

  std::vector<double> limitParameters = computeModelParameters(
    d_local.PEAKI1, d_local.FSLOPE, d_local.STREN, d_local.YSLOPE);
  d_local.a1 = limitParameters[0];
  d_local.a2 = limitParameters[1];
  d_local.a3 = limitParameters[2];
  d_local.a4 = limitParameters[3];

  // Get the plastic internal variables from the model state
  double pbar_w = state->pbar_w;
  double X_eff = state->capX + 3.0 * pbar_w;

  // Compute kappa
  double I1_diff = d_local.PEAKI1 - X_eff;
  double kappa = d_local.PEAKI1 - d_local.CR * I1_diff;

  // Get the bulk and shear moduli and compute sqrt(3/2 K/G)
  double sqrtKG = std::sqrt(1.5 * state->bulkModulus / state->shearModulus);

  // Compute diameter of yield surface in z-r space
  double sqrtJ2_diff = 2.0 * evalYieldConditionMax(state);
  double yield_surf_dia_zrprime =
    std::max(I1_diff * one_sqrt_three, sqrtJ2_diff * sqrt_two * sqrtKG);
  double dist_to_trial_zr =
    std::sqrt(z_r_pt.x() * z_r_pt.x() + z_r_pt.y() * z_r_pt.y());
  double dist_dia_ratio = dist_to_trial_zr / yield_surf_dia_zrprime;
  // int num_points = std::max(5, (int) std::ceil(std::log(dist_dia_ratio)));
  int num_points = std::max(5, (int)std::ceil(std::log(dist_dia_ratio)));

  // Set up I1 limits
  double I1eff_min = X_eff;
  double I1eff_max = d_local.PEAKI1;

  // Set up bisection
  double eta_lo = 0.0, eta_hi = 1.0;

  // Set up mid point
  double I1eff_mid = 0.5 * (I1eff_min + I1eff_max);
  double eta_mid = 0.5 * (eta_lo + eta_hi);

  // Do bisection
  int iters = 1;
  double TOLERANCE = 1.0e-10;
  std::vector<Uintah::Point> z_r_points;
  std::vector<Uintah::Point> z_r_segments;
  std::vector<Uintah::Point> z_r_segment_points;
  Uintah::Point z_r_closest_old;
  z_r_closest_old.x(std::numeric_limits<double>::max());
  z_r_closest_old.y(std::numeric_limits<double>::max());
  z_r_closest_old.z(0.0);
  while (std::abs(eta_hi - eta_lo) > TOLERANCE) {

    // Get the yield surface points
    z_r_points.clear();
    getYieldSurfacePointsAll_RprimeZ(X_eff, kappa, sqrtKG, I1eff_min, I1eff_max,
                                     num_points, z_r_points);

    // Find the closest point
    findClosestPoint(z_r_pt, z_r_points, z_r_closest);

#ifdef DEBUG_YIELD_BISECTION_R
    std::cout << "iteration = " << iters << std::endl;
    std::cout << "K = " << state->bulkModulus << std::endl;
    std::cout << "G = " << state->shearModulus << std::endl;
    std::cout << "X = " << state->capX << std::endl;
    std::cout << "pbar_w = " << state->pbar_w << std::endl;
    std::cout << "yieldParams = list(BETA = " << d_local.BETA << ", "
              << "CR = " << d_local.CR << ", "
              << "FSLOPE = " << d_local.FSLOPE << ", "
              << "PEAKI1 = " << d_local.PEAKI1 << ", "
              << "STREN = " << d_local.STREN << ", "
              << "YSLOPE = " << d_local.YSLOPE << ")" << std::endl;
    std::cout << "z_r_pt = c(" << z_r_pt.x() << "," << z_r_pt.y() << ")"
              << std::endl;
    std::cout << "z_r_closest = c(" << z_r_closest.x() << "," << z_r_closest.y()
              << ")" << std::endl;
    std::cout << "z_r_yield_z = c(";
    for (auto& pt : z_r_points) {
      if (pt == z_r_points.back()) {
        std::cout << pt.x();
      } else {
        std::cout << pt.x() << ",";
      }
    }
    std::cout << ")" << std::endl;
    std::cout << "z_r_yield_r = c(";
    for (auto& pt : z_r_points) {
      if (pt == z_r_points.back()) {
        std::cout << pt.y();
      } else {
        std::cout << pt.y() << ",";
      }
    }
    std::cout << ")" << std::endl;
    if (iters == 1) {
      std::cout << "zr_df = \n"
                << "  ComputeFullYieldSurface(yieldParams, X, pbar_w, K, G, "
                   "num_points,\n"
                << "                          z_r_pt, z_r_closest, "
                   "z_r_yield_z, z_r_yield_r,\n"
                << "                          iteration, consistency_iter)"
                << std::endl;
    } else {
      std::cout << "zr_df = rbind(zr_df,\n"
                << "  ComputeFullYieldSurface(yieldParams, X, pbar_w, K, G, "
                   "num_points,\n"
                << "                          z_r_pt, z_r_closest, "
                   "z_r_yield_z, z_r_yield_r,\n"
                << "                          iteration, consistency_iter))"
                << std::endl;
    }
#endif

#ifdef DEBUG_YIELD_BISECTION
    // if (state->particleID == 3377699720593411) {
    std::cout << "Iteration = " << iters << std::endl;
    std::cout << "State = " << *state << std::endl;
    std::cout << "z_r_pt = " << z_r_pt << ";" << std::endl;
    std::cout << "z_r_closest = " << z_r_closest << ";" << std::endl;
    std::cout << "z_r_yield_z = [";
    for (auto& pt : z_r_points) {
      std::cout << pt.x() << " ";
    }
    std::cout << "];" << std::endl;
    std::cout << "z_r_yield_r = [";
    for (auto& pt : z_r_points) {
      std::cout << pt.y() << " ";
    }
    std::cout << "];" << std::endl;
    std::cout << "plot(z_r_yield_z, z_r_yield_r); hold on;" << std::endl;
    std::cout << "plot(z_r_pt(1), z_r_pt(2));" << std::endl;
    std::cout << "plot(z_r_closest(1), z_r_closest(2));" << std::endl;
    std::cout
      << "plot([z_r_pt(1) z_r_closest(1)],[z_r_pt(2) z_r_closest(2)], '--');"
      << std::endl;
//}
#endif
#ifdef DEBUG_YIELD_BISECTION_I1_J2
    // if (state->particleID == 3377699720593411) {
    double fac_z = std::sqrt(3.0);
    double fac_r = d_local.BETA * sqrtKG * std::sqrt(2.0);
    std::cout << "Iteration = " << iters << std::endl;
    std::cout << "I1_J2_trial = [" << z_r_pt.x() * fac_z << " "
              << z_r_pt.y() / fac_r << "];" << std::endl;
    std::cout << "I1_J2_closest = [" << z_r_closest.x() * fac_z << " "
              << z_r_closest.y() / fac_r << "];" << std::endl;
    std::cout << "I1_J2_yield_I1 = [";
    for (auto& pt : z_r_points) {
      std::cout << pt.x() * fac_z << " ";
    }
    std::cout << "];" << std::endl;
    std::cout << "I1_J2_yield_J2 = [";
    for (auto& pt : z_r_points) {
      std::cout << pt.y() / fac_r << " ";
    }
    std::cout << "];" << std::endl;
    std::cout << "plot(I1_J2_yield_I1, I1_J2_yield_J2); hold on;" << std::endl;
    std::cout << "plot(I1_J2_trial(1), I1_J2_trial(2), 'ko');" << std::endl;
    std::cout << "plot(I1_J2_closest(1), I1_J2_closest(2));" << std::endl;
    std::cout << "plot([I1_J2_trial(1) I1_J2_closest(1)],[I1_J2_trial(2) "
                 "I1_J2_closest(2)], '--');"
              << std::endl;
//}
#endif

    // Compute I1 for the closest point
    double I1eff_closest = sqrt_three * z_r_closest.x();

    // If (I1_closest < I1_mid)
    if (I1eff_closest < I1eff_mid) {
      I1eff_max = I1eff_mid;
      eta_hi = eta_mid;
    } else {
      I1eff_min = I1eff_mid;
      eta_lo = eta_mid;
    }

    I1eff_mid = 0.5 * (I1eff_min + I1eff_max);
    eta_mid = 0.5 * (eta_lo + eta_hi);

    // Distance to old closest point
    if (iters > 10 && (z_r_closest - z_r_closest_old).length2() < 1.0e-16) {
      break;
    }
    z_r_closest_old = z_r_closest;

    ++iters;
  }

  return;
}

/* Get the points on the yield surface */
void
YieldCond_Tabular::getYieldSurfacePointsAll_RprimeZ(
  const double& X_eff, const double& kappa, const double& sqrtKG,
  const double& I1eff_min, const double& I1eff_max, const int& num_points,
  std::vector<Uintah::Point>& z_r_vec)
{
  // Compute z_eff and r'
  computeZeff_and_RPrime(X_eff, kappa, sqrtKG, I1eff_min, I1eff_max, num_points,
                         z_r_vec);

  return;
}

/* Get the points on two segments the yield surface */
void
YieldCond_Tabular::getYieldSurfacePointsSegment_RprimeZ(
  const double& X_eff, const double& kappa, const double& sqrtKG,
  const Uintah::Point& start_point, const Uintah::Point& end_point,
  const int& num_points, std::vector<Uintah::Point>& z_r_poly)
{

  // Find the start I1 and end I1 values of the segments
  // **TODO** make sure that the start and end points are differenet
  double z_effStart = start_point.x();
  double z_effEnd = end_point.x();
  double I1_effStart = sqrt_three * z_effStart;
  double I1_effEnd = sqrt_three * z_effEnd;

  // Compute z_eff and r'
  computeZeff_and_RPrime(X_eff, kappa, sqrtKG, I1_effStart, I1_effEnd,
                         num_points, z_r_poly);

  return;
}

/*! Compute a vector of z_eff, r' values given a range of I1_eff values */
void
YieldCond_Tabular::computeZeff_and_RPrime(
  const double& X_eff, const double& kappa, const double& sqrtKG,
  const double& I1eff_min, const double& I1eff_max, const int& num_points,
  std::vector<Uintah::Point>& z_r_vec)
{
  // Set up points
  double rad = 0.5 * (d_local.PEAKI1 - X_eff);
  double cen = 0.5 * (d_local.PEAKI1 + X_eff);
  double theta_max = std::acos(std::max((I1eff_min - cen) / rad, -1.0));
  double theta_min = std::acos(std::min((I1eff_max - cen) / rad, 1.0));
  std::vector<double> theta_vec;
  linspace(theta_min, theta_max, num_points, theta_vec);

  for (auto theta : theta_vec) {
    double I1_eff = std::max(cen + rad * std::cos(theta), X_eff);

    // Compute F_f
    double Ff = d_local.a1 - d_local.a3 * std::exp(d_local.a2 * I1_eff) -
                d_local.a4 * (I1_eff);
    double Ff_sq = Ff * Ff;

    // Compute Fc
    double Fc_sq = 1.0;
    if (I1_eff < kappa) {
      double ratio = (kappa - I1_eff) / (kappa - X_eff);
      Fc_sq = std::max(1.0 - ratio * ratio, 0.0);
    }

    // Compute J2
    double J2 = Ff_sq * Fc_sq;

// Check for nans
#ifdef CHECK_FOR_NANS
    if (std::isnan(I1_eff) || std::isnan(J2)) {
      double ratio = (kappa - I1_eff) / (kappa - X_eff);
      std::cout << "theta = " << theta << " kappa = " << kappa
                << " X_eff = " << X_eff << " I1_eff = " << I1_eff
                << " J2 = " << J2 << " Ff = " << Ff << " Fc_sq = " << Fc_sq
                << " ratio = " << ratio << std::endl;
      std::cout << "rad = " << rad << " cen = " << cen
                << " theta_max = " << theta_max << " theta_min = " << theta_min
                << " I1eff_max = " << I1eff_max << " I1eff_min = " << I1eff_min
                << std::endl;
    }
#endif

    z_r_vec.push_back(Uintah::Point(
      I1_eff / sqrt_three, d_local.BETA * std::sqrt(2.0 * J2) * sqrtKG, 0.0));
  }

  return;
}

//--------------------------------------------------------------
// Other yield condition functions

//--------------------------------------------------------------
// Compute d/depse_v(df/dp)
//   df/dp = 6*Ff*(a2*a3*exp(3*a2*p) + a4)*Fc^2 -
//             6*Ff^2*(kappa - I1)/(kappa - X)^2
//   d/depse_v(df/dp) =
//
// Requires:  Equation of state and internal variable
//--------------------------------------------------------------
double
YieldCond_Tabular::computeVolStrainDerivOfDfDp(
  const ModelStateBase* state_input, const PressureModel* eos,
  const ShearModulusModel*, const InternalVariableModel*)
{
  std::ostringstream out;
  out << "**ERROR** computeVolStrainDerivOfDfDp should not be called by "
      << " models that use the Tabular yield criterion.";
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
YieldCond_Tabular::computeDevStrainDerivOfDfDp(
  const ModelStateBase* state_input, const PressureModel* eos,
  const ShearModulusModel*, const InternalVariableModel*)
{
  std::ostringstream out;
  out << "**ERROR** computeDevStrainDerivOfDfDp should not be called by "
      << " models that use the Tabular yield criterion.";
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
YieldCond_Tabular::computeVolStrainDerivOfDfDq(
  const ModelStateBase* state_input, const PressureModel*,
  const ShearModulusModel* shear, const InternalVariableModel*)
{
  std::ostringstream out;
  out << "**ERROR** computeVolStrainDerivOfDfDq should not be called by "
      << " models that use the Tabular yield criterion.";
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
YieldCond_Tabular::computeDevStrainDerivOfDfDq(
  const ModelStateBase* state_input, const PressureModel*,
  const ShearModulusModel* shear, const InternalVariableModel*)
{
  std::ostringstream out;
  out << "**ERROR** computeDevStrainDerivOfDfDq should not be called by "
      << " models that use the Tabular yield criterion.";
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
YieldCond_Tabular::computeVolStrainDerivOfYieldFunction(
  const ModelStateBase* state_input, const PressureModel* eos,
  const ShearModulusModel* shear, const InternalVariableModel*)
{
  std::ostringstream out;
  out
    << "**ERROR** computeVolStrainDerivOfYieldFunction should not be called by "
    << " models that use the Tabular yield criterion.";
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
YieldCond_Tabular::computeDevStrainDerivOfYieldFunction(
  const ModelStateBase* state_input, const PressureModel* eos,
  const ShearModulusModel* shear, const InternalVariableModel*)
{
  std::ostringstream out;
  out
    << "**ERROR** computeVolStrainDerivOfYieldFunction should not be called by "
    << " models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return 0.0;
}

// Evaluate the yield function.
double
YieldCond_Tabular::evalYieldCondition(const double p, const double q,
                                      const double dummy0, const double dummy1,
                                      double& dummy2)
{
  std::ostringstream out;
  out
    << "**ERROR** Deprecated evalYieldCondition with double arguments. "
    << " Should not be called by models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return 0.0;
}

// Evaluate yield condition (s = deviatoric stress
//                           p = state->p)
double
YieldCond_Tabular::evalYieldCondition(const Uintah::Matrix3&,
                                      const ModelStateBase* state_input)
{
  std::ostringstream out;
  out << "**ERROR** evalYieldCondition with a Matrix3 argument should not be "
         "called by "
      << " models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return 0.0;
}

//--------------------------------------------------------------
// Other derivatives

// Compute df/dsigma
//    df/dsigma =
// where
//    s = sigma - 1/3 tr(sigma) I
void
YieldCond_Tabular::evalDerivOfYieldFunction(const Uintah::Matrix3& sig,
                                            const double p_c, const double,
                                            Uintah::Matrix3& derivative)
{
  std::ostringstream out;
  out << "**ERROR** evalDerivOfYieldCondition with a Matrix3 argument should "
         "not be "
      << "called by models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return;
}

// Compute df/ds  where s = deviatoric stress
//    df/ds =
void
YieldCond_Tabular::evalDevDerivOfYieldFunction(const Uintah::Matrix3& sigDev,
                                               const double, const double,
                                               Uintah::Matrix3& derivative)
{
  std::ostringstream out;
  out << "**ERROR** evalDerivOfYieldCondition with a Matrix3 argument should "
         "not be "
      << "called by models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return;
}

/*! Derivative with respect to the \f$xi\f$ where \f$\xi = s \f$
    where \f$s\f$ is deviatoric part of Cauchy stress */
void
YieldCond_Tabular::eval_df_dxi(const Matrix3& sigDev, const ModelStateBase*,
                               Matrix3& df_ds)

{
  std::ostringstream out;
  out << "**ERROR** eval_df_dxi with a Matrix3 argument should not be "
      << "called by models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return;
}

/* Derivative with respect to \f$ s \f$ and \f$ \beta \f$ */
void
YieldCond_Tabular::eval_df_ds_df_dbeta(const Matrix3& sigDev,
                                       const ModelStateBase*, Matrix3& df_ds,
                                       Matrix3& df_dbeta)
{
  std::ostringstream out;
  out << "**ERROR** eval_df_ds_df_dbeta with a Matrix3 argument should not be "
      << "called by models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return;
}

/*! Derivative with respect to the plastic strain (\f$\epsilon^p \f$) */
double
YieldCond_Tabular::eval_df_dep(const Matrix3&, const double& dsigy_dep,
                               const ModelStateBase*)
{
  std::ostringstream out;
  out << "**ERROR** eval_df_dep with a Matrix3 argument should not be "
      << "called by models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return 0.0;
}

/*! Derivative with respect to the porosity (\f$\epsilon^p \f$) */
double
YieldCond_Tabular::eval_df_dphi(const Matrix3&, const ModelStateBase*)
{
  std::ostringstream out;
  out << "**ERROR** eval_df_dphi with a Matrix3 argument should not be "
      << "called by models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return 0.0;
}

/*! Compute h_alpha  where \f$d/dt(ep) = d/dt(gamma)~h_{\alpha}\f$ */
double
YieldCond_Tabular::eval_h_alpha(const Matrix3&, const ModelStateBase*)
{
  std::ostringstream out;
  out << "**ERROR** eval_h_alpha with a Matrix3 argument should not be "
      << "called by models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return 1.0;
}

/*! Compute h_phi  where \f$d/dt(phi) = d/dt(gamma)~h_{\phi}\f$ */
double
YieldCond_Tabular::eval_h_phi(const Matrix3&, const double&,
                              const ModelStateBase*)
{
  std::ostringstream out;
  out << "**ERROR** eval_h_phi with a Matrix3 argument should not be "
      << "called by models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return 0.0;
}

//--------------------------------------------------------------
// Tangent moduli
void
YieldCond_Tabular::computeElasPlasTangentModulus(const TangentModulusTensor& Ce,
                                                 const Matrix3& sigma,
                                                 double sigY, double dsigYdep,
                                                 double porosity, double,
                                                 TangentModulusTensor& Cep)
{
  std::ostringstream out;
  out << "**ERROR** computeElasPlasTangentModulus with a Matrix3 argument "
         "should not be "
      << "called by models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return;
}

void
YieldCond_Tabular::computeTangentModulus(const TangentModulusTensor& Ce,
                                         const Matrix3& f_sigma, double f_q1,
                                         double h_q1, TangentModulusTensor& Cep)
{
  std::ostringstream out;
  out << "**ERROR** coputeTangentModulus with a Matrix3 argument should not be "
      << "called by models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return;
}
