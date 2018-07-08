/*
 * The MIT License
 *
 * Copyright (c) 2015-2018 Parresia Research Limited, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/Models/YieldCond_TabularCap.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/YieldCondUtils.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <chrono>
#include <cmath>

//#define DEBUG_CLOSEST_POINT

using namespace Vaango;
using Point = Uintah::Point;
using Vector = Uintah::Vector;
using Matrix3 = Uintah::Matrix3;

const double YieldCond_TabularCap::sqrt_two = std::sqrt(2.0);
const double YieldCond_TabularCap::sqrt_three = std::sqrt(3.0);
const double YieldCond_TabularCap::one_sqrt_three = 1.0 / sqrt_three;
const double YieldCond_TabularCap::large_number = 1.0e100;
const Matrix3 YieldCond_TabularCap::One(1,0,0,0,1,0,0,0,1);

YieldCond_TabularCap::YieldCond_TabularCap(Uintah::ProblemSpecP& ps)
  : d_yield(ps)
{
  // Check the input parameters
  checkInputParameters();
  setYieldConditionRange();

  // Save the data as a polyline
  saveAsPolyline();

  // Compute normals at each point of of the table
  computeNormals();

}

YieldCond_TabularCap::YieldCond_TabularCap(const YieldCond_TabularCap* yc)
{
  d_yield = yc->d_yield;
  d_I1bar_min = yc->d_I1bar_min;
  d_I1bar_max = yc->d_I1bar_max;
  d_sqrtJ2_max = yc->d_sqrtJ2_max;
  d_polyline = yc->d_polyline;
  d_normals = yc->d_normals;
}

void
YieldCond_TabularCap::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP yield_ps = ps->appendChild("plastic_yield_condition");
  yield_ps->setAttribute("type", "tabular");

  d_yield.table.outputProblemSpec(yield_ps);
}

void
YieldCond_TabularCap::checkInputParameters()
{
  if (d_yield.table.getNumIndependents() != 1) {
    std::ostringstream out;
    out << "**ERROR** The tabular yield data file contains more than one"
        << " independent variable.";
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  if (d_yield.table.getNumDependents() != 1) {
    std::ostringstream out;
    out << "**ERROR** The tabular yield data file contains more than one"
        << " dependent variable.";
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  // Check that the table contains the right variables
  for (const auto& var : d_yield.table.getIndependentVars()) {
    if (var->name != "Pressure") {
      std::ostringstream out;
      out << "**ERROR** The tabular yield data file does not contain the"
          << " independent variable \"Presssure\"";
      throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
    }
  }
  for (const auto& var : d_yield.table.getDependentVars()) {
    if (var->name != "SqrtJ2") {
      std::ostringstream out;
      out << "**ERROR** The tabular yield data file does not contain the"
          << " dependent variable \"SqrtJ2\"";
      throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
    }
  }

  // Check that independent and dependent var data contain the
  // same number of points
  DoubleVec1D xvals = 
    d_yield.table.getIndependentVarData("Pressure", IndexKey(0,0,0,0));
  DoubleVec1D yvals = 
    d_yield.table.getDependentVarData("SqrtJ2",IndexKey(0,0,0,0));
  if (xvals.size() != yvals.size()) {
    std::ostringstream out;
    out << "**ERROR** The tabular yield data file does not contain the"
          << " same number of \"Pressure\" and \"SqrtJ2\" data points";
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  // Check convexity; If not convex, make convex
  std::vector<Point> points;
  for (auto ii = 0u; ii < xvals.size(); ii++) {
    points.push_back(Point(xvals[ii], yvals[ii], 0));
  }
  std::vector<Point> hull = Vaango::Util::convexHull2D(points);
  //std::copy(hull.begin(), hull.end(),
  //          std::ostream_iterator<Point>(std::cout, " "));
  //std::cout << std::endl;

  if (points.size() != hull.size()) {
    std::sort(hull.begin(), hull.end(),
              [](const Point& a, const Point& b) {
                return a.x() < b.x();
              });
    //std::copy(hull.begin(), hull.end(),
    //          std::ostream_iterator<Point>(std::cout, " "));
    //std::cout << std::endl;
    xvals.clear();
    yvals.clear();
    for (const auto& point : hull) {
      xvals.push_back(point.x());
      yvals.push_back(point.y());
    }
    d_yield.table.setIndependentVarData("Pressure", IndexKey(0,0,0,0), xvals);
    d_yield.table.setDependentVarData("SqrtJ2", IndexKey(0,0,0,0), yvals);
  }

  //DoubleVec1D xvals_new = 
  //  d_yield.table.getIndependentVarData("Pressure", IndexKey(0,0,0,0));
  //DoubleVec1D yvals_new = 
  //  d_yield.table.getDependentVarData("SqrtJ2",IndexKey(0,0,0,0));
  //std::copy(xvals_new.begin(), xvals_new.end(),
  //          std::ostream_iterator<double>(std::cout, " "));
  //std::cout << std::endl;
  //std::copy(yvals_new.begin(), yvals_new.end(),
  //          std::ostream_iterator<double>(std::cout, " "));
  //std::cout << std::endl;

  if (d_yield.capEllipticityRatio >= 1 || d_yield.capEllipticityRatio <= 0.0) {
    std::ostringstream warn;
    warn << "capEllipticityRatio must be in [0, 1]. Input value = " 
         << d_yield.capEllipticityRatio << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

}

/* Yield surface does not change over time.  So we can caculate the
   max sqrt(J2) at the beginning */
void 
YieldCond_TabularCap::setYieldConditionRange() {

  d_I1bar_min = std::numeric_limits<double>::max();
  d_I1bar_max = std::numeric_limits<double>::min();
  d_sqrtJ2_max = std::numeric_limits<double>::min();

  DoubleVec1D pvals = 
    d_yield.table.getIndependentVarData("Pressure", IndexKey(0,0,0,0));
  for (auto p_bar : pvals) {
    DoubleVec1D gg = d_yield.table.interpolate<1>({{p_bar}}); 
    d_I1bar_min = std::min(d_I1bar_min, 3.0*p_bar);
    d_I1bar_max = std::max(d_I1bar_max, 3.0*p_bar);
    d_sqrtJ2_max = std::max(d_sqrtJ2_max, gg[0]);
  }
}

/* Save the yield surface points as a polyline for future computations */
void 
YieldCond_TabularCap::saveAsPolyline()
{
  DoubleVec1D xvals = 
    d_yield.table.getIndependentVarData("Pressure", IndexKey(0,0,0,0));
  DoubleVec1D yvals = 
    d_yield.table.getDependentVarData("SqrtJ2",IndexKey(0,0,0,0));

  if (xvals.size() < 3) {
    // Just extending the data so that normals can be computed
    d_polyline.push_back(Point(xvals[0], -10.0, 0));
    for (auto ii = 0u; ii < xvals.size(); ii++) {
      d_polyline.push_back(Point(xvals[ii], yvals[ii], 0));
    }
    double t = 1.1;
    Vector extra = d_polyline[1]*(1 - t) + d_polyline[2]*t;
    d_polyline.push_back(Point(extra));
  } else {
    Point first(xvals[0], 0, 0);
    Point second(xvals[1], -yvals[1], 0);
    double t = 0.01;
    Vector extra1 = first*(1 - t) + second*(t);
    d_polyline.push_back(second);
    d_polyline.push_back(Point(extra1));
    d_polyline.push_back(Point(xvals[0], yvals[0], 0));
    d_polyline.push_back(Point(extra1.x(), -extra1.y(), 0));
    for (auto ii = 1u; ii < xvals.size(); ii++) {
      d_polyline.push_back(Point(xvals[ii], yvals[ii], 0));
    }
    Point last = d_polyline[d_polyline.size()-1];
    Point secondlast = d_polyline[d_polyline.size()-2];
    t = 1.1;
    extra1 = secondlast*(1 - t) + last*t;
    Vector extra2 = secondlast*(1 - t)*t + last*(t*t + 1 - t);
    d_polyline.push_back(Point(extra1));
    d_polyline.push_back(Point(extra2));
  }
  //std::copy(d_polyline.begin(), d_polyline.end(),
  //          std::ostream_iterator<Point>(std::cout, " "));
  //std::cout << std::endl;
}

/* Compute normal at each point on the yield surface */
void
YieldCond_TabularCap::computeNormals()
{
  if (d_polyline.size() < 2) {
    std::ostringstream out;
    out << "**ERROR** The tabular yield data has not yet been saved"
        << " as a polyline";
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  d_normals = Vaango::Util::computeNormals(d_polyline);
  //std::cout << "Normals:";
  //std::copy(d_normals.begin(), d_normals.end(),
  //          std::ostream_iterator<Vector>(std::cout, " "));
  //std::cout << std::endl;
}

Point 
YieldCond_TabularCap::getClosestPoint(const double& p_bar, const double& sqrtJ2)
{
  Point curr(p_bar, sqrtJ2, 0.0);
  Point closest;
  Vaango::Util::findClosestPoint(curr, d_polyline, closest);
  //double distSq = Vaango::Util::findClosestPoint(curr, d_polyline, closest);
  //std::cout << "closest = " << closest << " distSq = " << distSq << "\n";
  return closest;
}

//--------------------------------------------------------------
// Evaluate yield condition
//
// f := sqrt(J2) - g(p)*Fc(p, X_p) = 0
// where
//     J2 = 1/2 s:s,  s = sigma - p I,  p = 1/3 Tr(sigma)
//     g(pbar) = table,  pbar = -p
//     Fc^2 = 1 - (kappa - p)^2/(kappa - X_p)^2
//     kappa = p_tension - R*(p_tension - X_p)
//     X_p = hydrostatic strength
//
// Returns:
//   hasYielded = -1.0 (if elastic)
//              =  1.0 (otherwise)
//--------------------------------------------------------------
double
YieldCond_TabularCap::evalYieldCondition(const ModelStateBase* state_input)
{
  const ModelState_Tabular* state =
    dynamic_cast<const ModelState_Tabular*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  // First check if the state is outside the tension cap
  double p_bar_min = d_I1bar_min/3.0;
  double p_bar = -state->I1/3;
  if (p_bar < p_bar_min) {
    return 1.0;
  }

  // Next check if the state is outside the compression cap
  double p_bar_max = -state->capX/3.0; 
  if (p_bar > p_bar_max) {
    return 1.0;
  }

  // Compute the value of the yield function for the yield function without cap
  DoubleVec1D gg;
  try {
    gg = d_yield.table.interpolate<1>({{p_bar}}); 
  } catch (Uintah::InvalidValue& e) {
    std::ostringstream out;
    out << "**ERROR** In evalYieldCondition:"
        << " p_bar = " << p_bar << "\n"
        << e.message() ;
    throw Uintah::InvalidValue(out.str(), __FILE__, __LINE__);
  }

  // Find location of the center of the cap ellipse
  double kappa_bar = p_bar_min + d_yield.capEllipticityRatio * (p_bar_max - p_bar_min);

  // For the cap region only
  if (p_bar > kappa_bar) {

    double ratio = (p_bar - kappa_bar)/(p_bar_max - kappa_bar);
    double Fc_sq = 1.0 - ratio*ratio;
    if (state->sqrt_J2 * state->sqrt_J2 > gg[0] * gg[0] * Fc_sq) {
      return 1.0;
    }
  } else {
    if (state->sqrt_J2 > gg[0]) {
      return 1.0;
    }
  }
  
  //std::cout << "p_bar = " << p_bar << " gg = " << gg[0] 
  //          << " sqrtJ2 = " << state->sqrt_J2 << std::endl;

  return -1.0;
}

//--------------------------------------------------------------
// Derivatives needed by return algorithms and Newton iterations

//--------------------------------------------------------------
// Evaluate yield condition max value of sqrtJ2
//--------------------------------------------------------------
double
YieldCond_TabularCap::evalYieldConditionMax(const ModelStateBase* )
{
  setYieldConditionRange();
  return d_sqrtJ2_max;
}

//--------------------------------------------------------------
/*! Compute Derivative with respect to the Cauchy stress (\f$\sigma \f$)
 *  Compute df/dsigma
 *
 *  for the yield function
 *      f := sqrt(J2(s)) - g(p)*Fc(p, X_p) = 0
 *  where
 *      J2 = 1/2 s:s,  s = sigma - p I,  p = 1/3 Tr(sigma)
 *      g(pbar) = table, pbar = -p
 *
 *  The derivative is
 *      df/dsigma = df/dp dp/dsigma + df/ds : ds/dsigma
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
YieldCond_TabularCap::eval_df_dsigma(const Matrix3&,
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

  //std::cout << "p = " << state->I1/3 << " sqrtJ2 = " << state->sqrt_J2 << "\n";

  double df_dp = computeVolStressDerivOfYieldFunction(state);
  double df_dJ2 = computeDevStressDerivOfYieldFunction(state);

  //std::cout << "df_dp = " << df_dp << " df_dJ2 = " << df_dJ2 << "\n";

  Matrix3 p_term = One * (df_dp / 3.0);
  Matrix3 s_term = state->deviatoricStressTensor * (df_dJ2);

  df_dsigma = p_term + s_term;

  df_dsigma /= df_dsigma.Norm();

  return;
}

//--------------------------------------------------------------
// Compute df/dp  where pI = volumetric stress = 1/3 Tr(sigma) I
//   df/dp = derivative of the yield function wrt p
//
// for the yield function
//     f := sqrt(J2) - g(p)*Fc(p, X_p) = 0
// where
//     J2 = 1/2 s:s,  s = sigma - p I,  p = 1/3 Tr(sigma)
//     g(pbar) := table,  pbar = -p
//     Fc^2 = 1 - (kappa - p)^2/(kappa - X_p)^2
//     kappa = p_tension - R*(p_tension - X_p)
//
// the derivative is
//     df/dp = -dg/dp * F_c(p, X_p) - g(p) * dF_c/dp 
//     dg/dp = -dg/dpbar dpbar/dp  = dg/dpbar
//     dFc/dp = (1/F_c)*(kappa - p)/(kappa - X_p)^2

//--------------------------------------------------------------
double
YieldCond_TabularCap::computeVolStressDerivOfYieldFunction(
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

  // Set up limits
  double p_bar_min = d_I1bar_min/3.0;
  double p_bar_max = -state->capX/3.0; 

  // Find location of the center of the cap ellipse
  double kappa_bar = p_bar_min + d_yield.capEllipticityRatio * (p_bar_max - p_bar_min);

  // Check that the closest point is not at the vertex
  double p_bar = -state->I1/3;
  double epsilon = 1.0e-6;
  Point closest = getClosestPoint(p_bar, state->sqrt_J2);
  if (closest.x() - epsilon < p_bar_min) {
    return large_number;
  }

  // Compute dg/dp
  DoubleVec1D gg, g_lo, g_hi;
  try {
    gg = d_yield.table.interpolate<1>({{closest.x()}});
    g_lo = d_yield.table.interpolate<1>({{closest.x()-epsilon}});
    g_hi = d_yield.table.interpolate<1>({{closest.x()+epsilon}});
  } catch (Uintah::InvalidValue& e) {
    std::ostringstream out;
    out << "**ERROR** In compute df/dp:"
        << " p_bar = " << closest.x() << "\n"
        << e.message() ;
    throw Uintah::InvalidValue(out.str(), __FILE__, __LINE__);
  }
  double dg_dpbar = (g_hi[0] - g_lo[0])/(2*epsilon);

  // Compute dFc/dp
  double ratio = (closest.x() - kappa_bar)/(p_bar_max - kappa_bar);
  double Fc = std::sqrt(1.0 - ratio * ratio);
  double dFc_dp = ratio/(Fc*(p_bar_max - kappa_bar));

  // Compute df_dp
  double df_dp = dg_dpbar * Fc - gg[0] * dFc_dp ;

  return df_dp;
}

//--------------------------------------------------------------
// Compute df/dJ2  where J2 = 1/2 s:s ,  s = sigma - p I,  p = 1/3 Tr(sigma)
//   s = derivatoric stress
//   df/dJ2 = derivative of the yield function wrt J2
//
// for the yield function
//     f := sqrt(J2) - g(p)*F_c(p,X_p) = 0
// where
//     J2 = 1/2 s:s,  s = sigma - p I,  p = 1/3 Tr(sigma)
//     g(pbar) := table, pbar = -p
//     Fc^2 = 1 - (kappa - p)^2/(kappa - X_p)^2
//     kappa = p_tension - R*(p_tension - X_p)
//
// the derivative is
//     df/dJ2 = 1/(2 sqrt(J2))
//--------------------------------------------------------------
double
YieldCond_TabularCap::computeDevStressDerivOfYieldFunction(
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

  double df_dJ2 = (state->sqrt_J2 == 0) ? large_number : 1/(2*state->sqrt_J2);

  return df_dJ2;
}

/**
 * Function: getClosestPoint
 * Purpose: Get the point on the yield surface that is closest to a given point
 * (2D)
 * Inputs:
 *  state = current state
 *  z = x-coordinate of point
 *  rprime = y-coordinate of point
 * Outputs:
 *  cz = x-coordinate of closest point on yield surface
 *  crprime = y-coordinate of closest point
 */
bool
YieldCond_TabularCap::getClosestPoint(const ModelStateBase* state_input,
                                   const double& z, const double& rprime,
                                   double& cz, double& crprime)
{
  const ModelState_Tabular* state =
    dynamic_cast<const ModelState_Tabular*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  //std::chrono::time_point<std::chrono::system_clock> start, end;
  //start = std::chrono::system_clock::now();
  Point pt(z, rprime, 0.0);
  Point closest(0.0, 0.0, 0.0);
  if (d_polyline.size() < 5) {
    closest = getClosestPointTable(state, pt);
  } else {
    closest = getClosestPointSpline(state, pt);
  }

  cz = closest.x();
  crprime = closest.y();
  //end = std::chrono::system_clock::now();
  //std::cout << "Geomeric Bisection : Time taken = " <<
  //std::chrono::duration<double>(end-start).count() << std::endl;

  return true;
}

Point
YieldCond_TabularCap::getClosestPointTable(const ModelState_Tabular* state, 
                                        const Point& z_r_pt)
{
  // Get the bulk and shear moduli and compute sqrt(3/2 K/G)
  double sqrtKG = std::sqrt(1.5 * state->bulkModulus / state->shearModulus);

  // Add cap points to tabular data
  Polyline p_q_table;
  computeCapPoints(-state->capX, p_q_table);

  // Convert tabular data to z-rprime coordinates
  Polyline z_r_table;
  convertToZRprime(sqrtKG, p_q_table, z_r_table);

  // Find the closest point
  Point z_r_closest;
  Vaango::Util::findClosestPoint(z_r_pt, z_r_table, z_r_closest);

  return z_r_closest;
}

Point
YieldCond_TabularCap::getClosestPointSpline(const ModelState_Tabular* state, 
                                            const Point& z_r_pt)
{
  // Get the bulk and shear moduli and compute sqrt(3/2 K/G)
  double sqrtKG = std::sqrt(1.5 * state->bulkModulus / state->shearModulus);

  // Add cap points to tabular data
  Polyline p_q_table;
  computeCapPoints(-state->capX, p_q_table);

  // Convert tabular data to z-rprime coordinates
  Polyline z_r_table;
  convertToZRprime(sqrtKG, p_q_table, z_r_table);

  // Find the closest segments
  Polyline z_r_segments;
  std::size_t closest_index = 
    Vaango::Util::getClosestSegments(z_r_pt, z_r_table, z_r_segments);

  // Get the yield surface points for the closest segments
  // (Fit quadratic B_spline)
  std::size_t numPts = z_r_table.size();
  std::size_t ptsPerSegment = 30;
  Polyline z_r_spline;
  if (closest_index == 0) {
    Vaango::Util::computeOpenUniformQuadraticBSpline(z_r_table, 
                                                     0, 1, ptsPerSegment,
                                                     z_r_spline);
  } else if (closest_index == numPts-2) {
    Vaango::Util::computeOpenUniformQuadraticBSpline(z_r_table, 
                                                     numPts-2, numPts-1,
                                                     ptsPerSegment,
                                                     z_r_spline);
  } else {
    Vaango::Util::computeOpenUniformQuadraticBSpline(z_r_table, 
                                                     closest_index-1,
                                                     closest_index+1,
                                                     ptsPerSegment,
                                                     z_r_spline);
  }

  // Find the closest point
  Point z_r_closest;
  Vaango::Util::findClosestPoint(z_r_pt, z_r_spline, z_r_closest);

  #ifdef DEBUG_CLOSEST_POINT
    std::cout << "ZRPoint = " << z_r_pt << std::endl;
    std::cout << "ZRTable = ";
    std::copy(z_r_table.begin(), z_r_table.end(),
              std::ostream_iterator<Point>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "ZRSpline = ";
    std::copy(z_r_spline.begin(), z_r_spline.end(),
              std::ostream_iterator<Point>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "ZRClose = " << z_r_closest << "\n";
  #endif

  return z_r_closest;
}

/* Add cap points to yield function table */
void
YieldCond_TabularCap::computeCapPoints(double X_bar, Polyline& p_q_all)
{
  // Set up limits
  double p_bar_min = d_I1bar_min/3.0;
  double p_bar_max = X_bar/3.0; 
  double kappa_bar = p_bar_min + d_yield.capEllipticityRatio * (p_bar_max - p_bar_min);

  // Find the location of the p_bar_start, p_bar_end on table polyline
  auto end_iter = std::find_if(d_polyline.begin(), d_polyline.end(),
                                [&kappa_bar](const auto& point) 
                                {
                                    return kappa_bar < point.x();
                                });
  // Copy the relevant points
  std::copy(d_polyline.begin(), end_iter, std::back_inserter(p_q_all));

  // Compute sqrtJ2 at that value of kappa
  auto start_iter = end_iter - 1;
  if (end_iter == d_polyline.end()) {
    start_iter--;
    end_iter--;
  }
  auto start = *start_iter;
  auto end = *end_iter;
  double t = (end.x() - kappa_bar)/(end.x() - start.x());
  double sqrtJ2_kappa = (1 - t)*start.y() + t*end.y();

  // Set up theta vector
  std::vector<double> theta_vec;
  Vaango::Util::linspace(0.0, M_PI/2, 100, theta_vec);

  // Set up ellipse axes
  double a = p_bar_max - kappa_bar;
  double b = sqrtJ2_kappa;

  // Compute ellipse points (**WARNING** may not return in order)
  Polyline p_q_cap;
  std::transform(theta_vec.begin(), theta_vec.end(), p_q_cap.begin(),
                 [&a, &b, &kappa_bar](const auto& theta) 
                 {
                   auto x = kappa_bar + a * cos(theta);
                   auto y = b * sin(theta);
                   return Uintah::Point(x, y, 0);
                 });
  
  // Sort cap points in increasing x
  std::sort(p_q_cap.begin(), p_q_cap.end(), 
            [](const Point& pt1, const Point& pt2) 
            {
              return pt1.x() < pt2.x();
            });

  // Concatenate the two vectors
  p_q_all.insert(p_q_all.end(), p_q_cap.begin(), p_q_cap.end());

}

/* Convert yield function data to z_rprime coordinates */
void 
YieldCond_TabularCap::convertToZRprime(const double& sqrtKG, 
                                       const Polyline& p_q_points, 
                                       Polyline& z_r_points) const
{
  // Compute z and r' for the yield surface points
  for (const auto& pt : p_q_points) {
    double p_bar = pt.x();
    double sqrt_J2 = pt.y();
    double z = -sqrt_three * p_bar;
    double rprime = sqrt_two * sqrt_J2 * sqrtKG;
    z_r_points.push_back(Point(z, rprime, 0));
  }
}
                   
//--------------------------------------------------------------
// Other yield condition functions (not used)

//--------------------------------------------------------------
// Compute d/depse_v(df/dp)
//   df/dp = 6*Ff*(a2*a3*exp(3*a2*p) + a4)*Fc^2 -
//             6*Ff^2*(kappa - I1)/(kappa - X)^2
//   d/depse_v(df/dp) =
//
// Requires:  Equation of state and internal variable
//--------------------------------------------------------------
double
YieldCond_TabularCap::computeVolStrainDerivOfDfDp(
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
YieldCond_TabularCap::computeDevStrainDerivOfDfDp(
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
YieldCond_TabularCap::computeVolStrainDerivOfDfDq(
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
YieldCond_TabularCap::computeDevStrainDerivOfDfDq(
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
YieldCond_TabularCap::computeVolStrainDerivOfYieldFunction(
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
YieldCond_TabularCap::computeDevStrainDerivOfYieldFunction(
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
YieldCond_TabularCap::evalYieldCondition(const double p, const double q,
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
YieldCond_TabularCap::evalYieldCondition(const Matrix3&,
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
YieldCond_TabularCap::evalDerivOfYieldFunction(const Matrix3& sig,
                                            const double p_c, const double,
                                            Matrix3& derivative)
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
YieldCond_TabularCap::evalDevDerivOfYieldFunction(const Matrix3& sigDev,
                                               const double, const double,
                                               Matrix3& derivative)
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
YieldCond_TabularCap::eval_df_dxi(const Matrix3& sigDev, const ModelStateBase*,
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
YieldCond_TabularCap::eval_df_ds_df_dbeta(const Matrix3& sigDev,
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
YieldCond_TabularCap::eval_df_dep(const Matrix3&, const double& dsigy_dep,
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
YieldCond_TabularCap::eval_df_dphi(const Matrix3&, const ModelStateBase*)
{
  std::ostringstream out;
  out << "**ERROR** eval_df_dphi with a Matrix3 argument should not be "
      << "called by models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return 0.0;
}

/*! Compute h_alpha  where \f$d/dt(ep) = d/dt(gamma)~h_{\alpha}\f$ */
double
YieldCond_TabularCap::eval_h_alpha(const Matrix3&, const ModelStateBase*)
{
  std::ostringstream out;
  out << "**ERROR** eval_h_alpha with a Matrix3 argument should not be "
      << "called by models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return 1.0;
}

/*! Compute h_phi  where \f$d/dt(phi) = d/dt(gamma)~h_{\phi}\f$ */
double
YieldCond_TabularCap::eval_h_phi(const Matrix3&, const double&,
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
YieldCond_TabularCap::computeElasPlasTangentModulus(const TangentModulusTensor& Ce,
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
YieldCond_TabularCap::computeTangentModulus(const TangentModulusTensor& Ce,
                                         const Matrix3& f_sigma, double f_q1,
                                         double h_q1, TangentModulusTensor& Cep)
{
  std::ostringstream out;
  out << "**ERROR** coputeTangentModulus with a Matrix3 argument should not be "
      << "called by models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return;
}

namespace Vaango {

std::ostream&
operator<<(std::ostream& out, const YieldCond_TabularCap& yc)
{
  DoubleVec1D pvals, qvals;
  try {
    pvals = yc.d_yield.table.getIndependentVarData("Pressure", 
                                                   IndexKey(0,0,0,0));
  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << std::endl;
  }

  try {
    qvals = yc.d_yield.table.getDependentVarData("SqrtJ2", 
                                                 IndexKey(0,0,0,0));
  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << std::endl;
  }

  out << "p:";
  std::copy(pvals.begin(), pvals.end(),
            std::ostream_iterator<double>(out, " "));
  out << std::endl;
  out << "sqrtJ2:";
  std::copy(qvals.begin(), qvals.end(),
            std::ostream_iterator<double>(out, " "));
  out << std::endl;

  out << "I1_bar_min = " << yc.d_I1bar_min
      << " I1_bar_max = " << yc.d_I1bar_max
      << " sqrtJ2_max = " << yc.d_sqrtJ2_max << std::endl;

  out << "Polyline:";
  std::copy(yc.d_polyline.begin(), yc.d_polyline.end(),
            std::ostream_iterator<Point>(out, " "));
  out << std::endl;

  out << "Normals:";
  std::copy(yc.d_normals.begin(), yc.d_normals.end(),
            std::ostream_iterator<Vector>(out, " "));
  out << std::endl;
  return out;
}

} // end namespace Vaango
