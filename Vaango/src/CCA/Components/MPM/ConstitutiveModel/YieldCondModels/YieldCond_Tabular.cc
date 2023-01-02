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

#include <CCA/Components/MPM/ConstitutiveModel/Utilities/YieldCondUtils.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCond_Tabular.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <chrono>
#include <cmath>

//#define DEBUG_CLOSEST_POINT
#define USE_NEWTON_CLOSEST_POINT

using namespace Vaango;
using Point   = Uintah::Point;
using Vector  = Uintah::Vector;
using Matrix3 = Uintah::Matrix3;

const double YieldCond_Tabular::sqrt_two       = std::sqrt(2.0);
const double YieldCond_Tabular::sqrt_three     = std::sqrt(3.0);
const double YieldCond_Tabular::one_sqrt_three = 1.0 / sqrt_three;
const double YieldCond_Tabular::large_number   = 1.0e100;
const Matrix3 YieldCond_Tabular::One(1, 0, 0, 0, 1, 0, 0, 0, 1);

YieldCond_Tabular::YieldCond_Tabular(Uintah::ProblemSpecP& ps,
                                     IntVar_TabularCap* intvar)
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

YieldCond_Tabular::YieldCond_Tabular(const YieldCond_Tabular* yc)
{
  d_yield      = yc->d_yield;
  d_I1bar_min  = yc->d_I1bar_min;
  d_I1bar_max  = yc->d_I1bar_max;
  d_sqrtJ2_max = yc->d_sqrtJ2_max;
  d_polyline   = yc->d_polyline;
  d_normals    = yc->d_normals;
}

void
YieldCond_Tabular::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP yield_ps = ps->appendChild("yield_condition");
  yield_ps->setAttribute("type", "tabular");

  d_yield.table.outputProblemSpec(yield_ps);
}

void
YieldCond_Tabular::checkInputParameters()
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
    d_yield.table.getIndependentVarData("Pressure", IndexKey(0, 0, 0, 0));
  DoubleVec1D yvals =
    d_yield.table.getDependentVarData("SqrtJ2", IndexKey(0, 0, 0, 0));
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
  // std::cout << "Orig = ";
  // std::copy(points.begin(), points.end(),
  //          std::ostream_iterator<Point>(std::cout, " "));
  // std::cout << std::endl;

  std::vector<Point> hull = Vaango::Util::convexHull2D(points);

  if (points.size() != hull.size()) {
    std::sort(hull.begin(), hull.end(), [](const Point& a, const Point& b) {
      return a.x() < b.x();
    });
    // std::cout << "Hull = ";
    // std::copy(hull.begin(), hull.end(),
    //          std::ostream_iterator<Point>(std::cout, " "));
    // std::cout << std::endl;

    // Interpolate points.y() to hull
    /*
    for (auto& point : points) {
      auto xval = point.x();
      auto end_iter = std::find_if(hull.begin(), hull.end(),
                                  [&xval](const auto& hull_pt)
                                  {
                                    return xval < hull_pt.x();
                                  });
      auto start = *(end_iter-1);
      auto end = *end_iter;
      auto t = (xval - start.x())/(end.x() - start.x());
      point.y((1- t)*start.y() + t*end.y());
    }
    */

    // Project points to hull
    for (auto& point : points) {
      auto xval = point.x();
      auto yval = point.y();
      point     = getClosestPoint(hull, xval, yval);
    }

    xvals.clear();
    yvals.clear();
    // for (const auto& point : hull) {
    for (const auto& point : points) {
      xvals.push_back(point.x());
      yvals.push_back(point.y());
    }
    d_yield.table.setIndependentVarData(
      "Pressure", IndexKey(0, 0, 0, 0), xvals);
    d_yield.table.setDependentVarData("SqrtJ2", IndexKey(0, 0, 0, 0), yvals);
  }

  // DoubleVec1D xvals_new =
  //  d_yield.table.getIndependentVarData("Pressure", IndexKey(0,0,0,0));
  // DoubleVec1D yvals_new =
  //  d_yield.table.getDependentVarData("SqrtJ2",IndexKey(0,0,0,0));
  // std::copy(xvals_new.begin(), xvals_new.end(),
  //          std::ostream_iterator<double>(std::cout, " "));
  // std::cout << std::endl;
  // std::copy(yvals_new.begin(), yvals_new.end(),
  //          std::ostream_iterator<double>(std::cout, " "));
  // std::cout << std::endl;
}

/* Yield surface does not change over time.  So we can caculate the
   max sqrt(J2) at the beginning */
void
YieldCond_Tabular::setYieldConditionRange()
{

  d_I1bar_min  = std::numeric_limits<double>::max();
  d_I1bar_max  = std::numeric_limits<double>::min();
  d_sqrtJ2_max = std::numeric_limits<double>::min();

  DoubleVec1D pvals =
    d_yield.table.getIndependentVarData("Pressure", IndexKey(0, 0, 0, 0));
  for (auto p_bar : pvals) {
    DoubleVec1D gg = d_yield.table.interpolate<1>({ { p_bar } });
    d_I1bar_min    = std::min(d_I1bar_min, 3.0 * p_bar);
    d_I1bar_max    = std::max(d_I1bar_max, 3.0 * p_bar);
    d_sqrtJ2_max   = std::max(d_sqrtJ2_max, gg[0]);
  }
}

/* Save the yield surface points as a polyline for future computations */
void
YieldCond_Tabular::saveAsPolyline()
{
  DoubleVec1D xvals =
    d_yield.table.getIndependentVarData("Pressure", IndexKey(0, 0, 0, 0));
  DoubleVec1D yvals =
    d_yield.table.getDependentVarData("SqrtJ2", IndexKey(0, 0, 0, 0));

  if (xvals.size() < 3) {
    // Just extending the data so that normals can be computed
    d_polyline.push_back(Point(xvals[0], -10.0, 0));
    for (auto ii = 0u; ii < xvals.size(); ii++) {
      d_polyline.push_back(Point(xvals[ii], yvals[ii], 0));
    }
    double t     = 1.1;
    Vector extra = d_polyline[1] * (1 - t) + d_polyline[2] * t;
    d_polyline.push_back(Point(extra));
  } else {
    Point first(xvals[0], 0, 0);
    Point second(xvals[1], -yvals[1], 0);
    double t      = 0.01;
    Vector extra1 = first * (1 - t) + second * (t);
    d_polyline.push_back(second);
    d_polyline.push_back(Point(extra1));
    d_polyline.push_back(Point(xvals[0], yvals[0], 0));
    d_polyline.push_back(Point(extra1.x(), -extra1.y(), 0));
    for (auto ii = 1u; ii < xvals.size(); ii++) {
      d_polyline.push_back(Point(xvals[ii], yvals[ii], 0));
    }
    Point last       = d_polyline[d_polyline.size() - 1];
    Point secondlast = d_polyline[d_polyline.size() - 2];
    t                = 1.1;
    extra1           = secondlast * (1 - t) + last * t;
    Vector extra2    = secondlast * (1 - t) * t + last * (t * t + 1 - t);
    d_polyline.push_back(Point(extra1));
    d_polyline.push_back(Point(extra2));
  }
  // std::copy(d_polyline.begin(), d_polyline.end(),
  //          std::ostream_iterator<Point>(std::cout, " "));
  // std::cout << std::endl;
}

/* Compute normal at each point on the yield surface */
void
YieldCond_Tabular::computeNormals()
{
  if (d_polyline.size() < 2) {
    std::ostringstream out;
    out << "**ERROR** The tabular yield data has not yet been saved"
        << " as a polyline";
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  d_normals = Vaango::Util::computeNormals(d_polyline);
  // std::cout << "Normals:";
  // std::copy(d_normals.begin(), d_normals.end(),
  //          std::ostream_iterator<Vector>(std::cout, " "));
  // std::cout << std::endl;
}

Point
YieldCond_Tabular::getClosestPoint(const double& p_bar, const double& sqrtJ2)
{
  Point curr(p_bar, sqrtJ2, 0.0);
  Point closest;
  Vaango::Util::findClosestPoint(curr, d_polyline, closest);
  // double distSq = Vaango::Util::findClosestPoint(curr, d_polyline, closest);
  // std::cout << "closest = " << closest << " distSq = " << distSq << "\n";
  return closest;
}

Point
YieldCond_Tabular::getClosestPoint(const Polyline& polyline,
                                   const double& p_bar,
                                   const double& sqrtJ2)
{
  Point curr(p_bar, sqrtJ2, 0.0);
  Point closest;
  Vaango::Util::findClosestPoint(curr, polyline, closest);
  // double distSq = Vaango::Util::findClosestPoint(curr, d_polyline, closest);
  // std::cout << "closest = " << closest << " distSq = " << distSq << "\n";
  return closest;
}

//--------------------------------------------------------------
// Evaluate yield condition
//
// f := sqrt(J2) - g(p) = 0
// where
//     J2 = 1/2 s:s,  s = sigma - p I,  p = 1/3 Tr(sigma)
//     g(pbar) = table,  pbar = -p
//
// Returns:
//   hasYielded = -1.0 (if elastic)
//              =  1.0 (otherwise)
//--------------------------------------------------------------
std::pair<double, Util::YieldStatus>
YieldCond_Tabular::evalYieldCondition(const ModelStateBase* state_input)
{
  const ModelState_Tabular* state =
    static_cast<const ModelState_Tabular*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  double p_bar = -state->I1 / 3;
  if (p_bar < d_I1bar_min / 3) {
    return std::make_pair(1.0, Util::YieldStatus::HAS_YIELDED);
  }

  DoubleVec1D gg;
  try {
    gg = d_yield.table.interpolate<1>({ { p_bar } });
  } catch (Uintah::InvalidValue& e) {
    std::ostringstream out;
    out << "**ERROR** In evalYieldCondition:"
        << " p_bar = " << p_bar << "\n"
        << e.message();
    throw Uintah::InvalidValue(out.str(), __FILE__, __LINE__);
  }
  // std::cout << "p_bar = " << p_bar << " gg = " << gg[0]
  //          << " sqrtJ2 = " << state->sqrt_J2 << std::endl;
  if (state->sqrt_J2 > gg[0]) {
    return std::make_pair(1.0, Util::YieldStatus::HAS_YIELDED);
  }

  return std::make_pair(-1.0, Util::YieldStatus::IS_ELASTIC);
}

double
YieldCond_Tabular::computeYieldFunction(const ModelStateBase* state) const
{
  std::ostringstream out;
  out << "**ERROR** The yield function for the tabular plasticity models"
      << " cannot be evaluated for a given stress state.\n";
  throw Uintah::InternalError(out.str(), __FILE__, __LINE__);

  return 0.0;
}

//--------------------------------------------------------------
// Evaluate yield condition max value of sqrtJ2
//--------------------------------------------------------------
double
YieldCond_Tabular::evalYieldConditionMax(const ModelStateBase*)
{
  setYieldConditionRange();
  return d_sqrtJ2_max;
}

//--------------------------------------------------------------
/*! Compute Derivative with respect to the Cauchy stress (\f$\sigma \f$)
 *  Compute df/dsigma
 *
 *  for the yield function
 *      f := sqrt(J2(s)) - g(p) = 0
 *  where
 *      J2 = 1/2 s:s,  s = sigma - p I,  p = 1/3 Tr(sigma)
 *      g(pbar) = table, pbar = -p
 *
 *  The derivative is
 *      df/dsigma = df/dp dp/dsigma + df/ds : ds/dsigma
 *  where
 *      df/dp = df_dp
 *      dp/dsigma = 1/3 I
 *  and
 *      df/ds = df/dJ2 dJ2/ds
 *      df/dJ2 = df_dq
 *      dJ2/ds = s
 *      ds/dsigma = I(4s) - 1/3 II
 *  which means
 *      df/dp dp/dsigma = 1/3 df/dp I
 *      df/ds : ds/dsigma = df/dJ2 s : [I(4s) - 1/3 II]
 *                        = df/dJ2 s
*/
/*! Derivative with respect to the Cauchy stress (\f$\sigma \f$)*/
Matrix3
YieldCond_Tabular::df_dsigma(const ModelStateBase* state) 
{
  return df_dsigma(Vaango::Util::Identity, state);
}

Uintah::Matrix3
YieldCond_Tabular::df_dsigma(const Matrix3&,
                             const ModelStateBase* state_input)
{
  const ModelState_Tabular* state =
    static_cast<const ModelState_Tabular*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  // std::cout << "p = " << state->I1/3 << " sqrtJ2 = " << state->sqrt_J2 <<
  // "\n";

  double dfdp  = df_dp(state);
  double dfdJ2 = df_dq(state);

  // std::cout << "df_dp = " << df_dp << " df_dJ2 = " << df_dJ2 << "\n";

  Matrix3 p_term = One * (dfdp / 3.0);
  Matrix3 s_term = state->deviatoricStressTensor * (dfdJ2);

  Matrix3 df_dsigma = p_term + s_term;

  return df_dsigma;
}

//--------------------------------------------------------------
// Compute df/dp  where pI = volumetric stress = 1/3 Tr(sigma) I
//   df/dp = derivative of the yield function wrt p
//
// for the yield function
//     f := sqrt(J2) - g(p) = 0
// where
//     J2 = 1/2 s:s,  s = sigma - p I,  p = 1/3 Tr(sigma)
//     g(pbar) := table,  pbar = -p
//
// the derivative is
//     df/dp = -dg/dp = -dg/dpbar dpbar/dp = dg/dpbar
//--------------------------------------------------------------
double
YieldCond_Tabular::df_dp(const ModelStateBase* state_input)
{
  const ModelState_Tabular* state =
    static_cast<const ModelState_Tabular*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  double epsilon = 1.0e-6;
  Point closest  = getClosestPoint(-state->I1 / 3, state->sqrt_J2);
  if (closest.x() - epsilon < d_I1bar_min / 3.0) {
    return large_number;
  }

  DoubleVec1D g_lo, g_hi;
  try {
    g_lo = d_yield.table.interpolate<1>({ { closest.x() - epsilon } });
    g_hi = d_yield.table.interpolate<1>({ { closest.x() + epsilon } });
  } catch (Uintah::InvalidValue& e) {
    std::ostringstream out;
    out << "**ERROR** In compute df/dp:"
        << " p_bar = " << closest.x() << "\n"
        << e.message();
    throw Uintah::InvalidValue(out.str(), __FILE__, __LINE__);
  }
  double dg_dpbar = (g_hi[0] - g_lo[0]) / (2 * epsilon);

  return dg_dpbar;
}

//--------------------------------------------------------------
// Compute df/dJ2  where J2 = 1/2 s:s ,  s = sigma - p I,  p = 1/3 Tr(sigma)
//   s = derivatoric stress
//   df/dJ2 = derivative of the yield function wrt J2
//
// for the yield function
//     f := sqrt(J2) - g(p) = 0
// where
//     J2 = 1/2 s:s,  s = sigma - p I,  p = 1/3 Tr(sigma)
//     g(pbar) := table, pbar = -p
//
// the derivative is
//     df/dJ2 = 1/(2 sqrt(J2))
//--------------------------------------------------------------
double
YieldCond_Tabular::df_dq(const ModelStateBase* state_input)
{
  const ModelState_Tabular* state =
    static_cast<const ModelState_Tabular*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  double df_dJ2 = (state->sqrt_J2 == 0) ? 0.0 : 1 / (2 * state->sqrt_J2);

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
YieldCond_Tabular::getClosestPoint(const ModelStateBase* state_input,
                                   const double& z,
                                   const double& rprime,
                                   double& cz,
                                   double& crprime)
{
  const ModelState_Tabular* state =
    static_cast<const ModelState_Tabular*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  // std::chrono::time_point<std::chrono::system_clock> start, end;
  // start = std::chrono::system_clock::now();
  Point pt(z, rprime, 0.0);
  Point closest(0.0, 0.0, 0.0);
  if (d_polyline.size() < 5) {
    closest = getClosestPointTable(state, pt);
  } else {
#ifdef USE_NEWTON_CLOSEST_POINT
    closest = getClosestPointSplineNewton(state, pt);
#else
    closest = getClosestPointSpline(state, pt);
#endif
  }

  cz      = closest.x();
  crprime = closest.y();
  // end = std::chrono::system_clock::now();
  // std::cout << "Geomeric Bisection : Time taken = " <<
  // std::chrono::duration<double>(end-start).count() << std::endl;

  return true;
}

Point
YieldCond_Tabular::getClosestPointTable(const ModelState_Tabular* state,
                                        const Point& z_r_pt)
{
  // Get the bulk and shear moduli and compute sqrt(3/2 K/G)
  double sqrtKG = std::sqrt(1.5 * state->bulkModulus / state->shearModulus);

  // Convert tabular data to z-rprime coordinates
  Polyline z_r_table;
  convertToZRprime(sqrtKG, z_r_table);

  // Find the closest point
  Point z_r_closest;
  Vaango::Util::findClosestPoint(z_r_pt, z_r_table, z_r_closest);

  return z_r_closest;
}

Point
YieldCond_Tabular::getClosestPointSpline(const ModelState_Tabular* state,
                                         const Point& z_r_pt)
{
  // Get the bulk and shear moduli and compute sqrt(3/2 K/G)
  double sqrtKG = std::sqrt(1.5 * state->bulkModulus / state->shearModulus);

  // Convert tabular data to z-rprime coordinates
  Polyline z_r_table;
  convertToZRprime(sqrtKG, z_r_table);

  // Find the closest segments
  Polyline z_r_segments;
  std::size_t closest_index =
    Vaango::Util::getClosestSegments(z_r_pt, z_r_table, z_r_segments);

  // Get the yield surface points for the closest segments
  // (Fit quadratic B_spline)
  std::size_t numPts        = z_r_table.size();
  std::size_t ptsPerSegment = 30;
  Polyline z_r_spline;
  if (closest_index == 0) {
    Vaango::Util::computeOpenUniformQuadraticBSpline(
      z_r_table, 0, 1, ptsPerSegment, z_r_spline);
  } else if (closest_index == numPts - 2) {
    Vaango::Util::computeOpenUniformQuadraticBSpline(
      z_r_table, numPts - 2, numPts - 1, ptsPerSegment, z_r_spline);
  } else {
    Vaango::Util::computeOpenUniformQuadraticBSpline(z_r_table,
                                                     closest_index - 1,
                                                     closest_index + 1,
                                                     ptsPerSegment,
                                                     z_r_spline);
  }

  // Find the closest point
  Point z_r_closest;
  Vaango::Util::findClosestPoint(z_r_pt, z_r_spline, z_r_closest);

#ifdef DEBUG_CLOSEST_POINT
  std::cout << "ZRPoint = " << z_r_pt << std::endl;
  std::cout << "ZRTable = ";
  std::copy(z_r_table.begin(),
            z_r_table.end(),
            std::ostream_iterator<Point>(std::cout, " "));
  std::cout << std::endl;
  std::cout << "ZRSpline = ";
  std::copy(z_r_spline.begin(),
            z_r_spline.end(),
            std::ostream_iterator<Point>(std::cout, " "));
  std::cout << std::endl;
  std::cout << "ZRClose = " << z_r_closest << "\n";
#endif

  return z_r_closest;
}

Point
YieldCond_Tabular::getClosestPointSplineNewton(const ModelState_Tabular* state,
                                               const Point& z_r_pt)
{
  // Get the bulk and shear moduli and compute sqrt(3/2 K/G)
  double sqrtKG = std::sqrt(1.5 * state->bulkModulus / state->shearModulus);

  // Convert tabular data to z-rprime coordinates
  Polyline z_r_table;
  convertToZRprime(sqrtKG, z_r_table);

  // Find the closest segments
  Polyline z_r_segments;
  std::size_t closest_index =
    Vaango::Util::getClosestSegments(z_r_pt, z_r_table, z_r_segments);

  // Get the yield surface points for the closest segments
  // (Fit quadratic B_spline)
  std::size_t numPts = z_r_table.size();
  auto seg_start     = closest_index - 1;
  auto seg_end       = closest_index + 1;
  if (closest_index < 2) {
    seg_start = 0;
    seg_end   = 2;
  } else if (closest_index > numPts - 3) {
    seg_start = numPts - 3;
    seg_end   = numPts - 1;
  }

#ifdef DEBUG_CLOSEST_POINT
  std::cout << "closest_index = " << closest_index << "\n";
  std::cout << "indices are (" << seg_start << "," << seg_end << ") from"
            << "(0," << numPts - 1 << ")\n";
#endif

  Point z_r_closest(0, 0, 0);
  Vector z_r_tangent(0, 0, 0);
  Vector z_r_dtangent(0, 0, 0);
  std::tie(z_r_closest, z_r_tangent, z_r_dtangent) =
    Vaango::Util::computeClosestPointQuadraticBSpline(
      z_r_pt, z_r_table, seg_start, seg_end);

#ifdef DEBUG_CLOSEST_POINT
  std::cout << "ZRClose = " << std::setprecision(16) << z_r_closest << "\n";
#endif

  return z_r_closest;
}

/* Convert yield function data to z_rprime coordinates */
void
YieldCond_Tabular::convertToZRprime(const double& sqrtKG,
                                    Polyline& z_r_points) const
{
  // Compute z and r' for the yield surface
  for (const auto& pt : d_polyline) {
    double p_bar   = pt.x();
    double sqrt_J2 = pt.y();
    double z       = -sqrt_three * p_bar;
    double rprime  = sqrt_two * sqrt_J2 * sqrtKG;
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
YieldCond_Tabular::d2f_dp_depsVol(const ModelStateBase* state_input,
                                  const MPMEquationOfState* eos,
                                  const ShearModulusModel*)
{
  std::ostringstream out;
  out << "**ERROR** d2f_dp_depsVol should not be called by "
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
YieldCond_Tabular::d2f_dp_depsDev(const ModelStateBase* state_input,
                                  const MPMEquationOfState* eos,
                                  const ShearModulusModel*)
{
  std::ostringstream out;
  out << "**ERROR** d2f_dp_depsDev should not be called by "
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
YieldCond_Tabular::d2f_dq_depsVol(const ModelStateBase* state_input,
                                  const MPMEquationOfState*,
                                  const ShearModulusModel* shear)
{
  std::ostringstream out;
  out << "**ERROR** d2f_dq_depsVol should not be called by "
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
YieldCond_Tabular::d2f_dq_depsDev(const ModelStateBase* state_input,
                                  const MPMEquationOfState*,
                                  const ShearModulusModel* shear)
{
  std::ostringstream out;
  out << "**ERROR** d2f_dq_depsDev should not be called by "
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
YieldCond_Tabular::df_depsVol(const ModelStateBase* state_input,
                              const MPMEquationOfState* eos,
                              const ShearModulusModel* shear)
{
  std::ostringstream out;
  out << "**ERROR** df_depsVol should not be called by "
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
YieldCond_Tabular::df_depsDev(const ModelStateBase* state_input,
                              const MPMEquationOfState* eos,
                              const ShearModulusModel* shear)
{
  std::ostringstream out;
  out << "**ERROR** df_depsVol should not be called by "
      << " models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return 0.0;
}

// Evaluate yield condition (s = deviatoric stress
//                           p = state->p)
double
YieldCond_Tabular::evalYieldCondition(const Matrix3&,
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

/*! Derivative with respect to the \f$xi\f$ where \f$\xi = s \f$
    where \f$s\f$ is deviatoric part of Cauchy stress */
Uintah::Matrix3
YieldCond_Tabular::df_dxi(const Matrix3& stress,
                          const ModelStateBase*)

{
  std::ostringstream out;
  out << "**ERROR** df_dxi with a Matrix3 argument should not be "
      << "called by models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  
  return Uintah::Matrix3(0.0);
}

/* Derivative with respect to \f$ s \f$ and \f$ \beta \f$ */
std::pair<Uintah::Matrix3, Uintah::Matrix3>
YieldCond_Tabular::df_dsigmaDev_dbeta(const Matrix3& stress,
                                      const ModelStateBase*)
{
  std::ostringstream out;
  out << "**ERROR** df_dsigmaDev_dbeta with a Matrix3 argument should not be "
      << "called by models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  
  return std::make_pair(Uintah::Matrix3(0.0), Uintah::Matrix3(0.0));
}

//--------------------------------------------------------------
// Tangent moduli
void
YieldCond_Tabular::computeElasPlasTangentModulus(const TangentModulusTensor& Ce,
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
      << "called by models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return;
}

void
YieldCond_Tabular::computeTangentModulus(const TangentModulusTensor& Ce,
                                         const Matrix3& f_sigma,
                                         double f_q1,
                                         double h_q1,
                                         TangentModulusTensor& Cep)
{
  std::ostringstream out;
  out << "**ERROR** coputeTangentModulus with a Matrix3 argument should not be "
      << "called by models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return;
}

namespace Vaango {

std::ostream&
operator<<(std::ostream& out, const YieldCond_Tabular& yc)
{
  DoubleVec1D pvals, qvals;
  try {
    pvals =
      yc.d_yield.table.getIndependentVarData("Pressure", IndexKey(0, 0, 0, 0));
  } catch (const Uintah::InvalidValue& e) {
    std::cout << e.message() << std::endl;
  }

  try {
    qvals =
      yc.d_yield.table.getDependentVarData("SqrtJ2", IndexKey(0, 0, 0, 0));
  } catch (const Uintah::InvalidValue& e) {
    std::cout << e.message() << std::endl;
  }

  out << "p:";
  std::copy(
    pvals.begin(), pvals.end(), std::ostream_iterator<double>(out, " "));
  out << std::endl;
  out << "sqrtJ2:";
  std::copy(
    qvals.begin(), qvals.end(), std::ostream_iterator<double>(out, " "));
  out << std::endl;

  out << "I1_bar_min = " << yc.d_I1bar_min << " I1_bar_max = " << yc.d_I1bar_max
      << " sqrtJ2_max = " << yc.d_sqrtJ2_max << std::endl;

  out << "Polyline:";
  std::copy(yc.d_polyline.begin(),
            yc.d_polyline.end(),
            std::ostream_iterator<Point>(out, " "));
  out << std::endl;

  out << "Normals:";
  std::copy(yc.d_normals.begin(),
            yc.d_normals.end(),
            std::ostream_iterator<Vector>(out, " "));
  out << std::endl;
  return out;
}

} // end namespace Vaango
