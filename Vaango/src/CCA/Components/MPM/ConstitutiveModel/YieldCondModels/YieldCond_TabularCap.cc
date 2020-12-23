/*
 * The MIT License
 *
 * Copyright (c) 2018-2020 Parresia Research Limited, New Zealand
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
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCond_TabularCap.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <chrono>
#include <cmath>

//#define DO_BINARY_CLOSEST_SEGMENT
#define DO_CLOSEST_POINT_ELASTIC_CHECK
//#define DEBUG_CLOSEST_POINT_ELASTIC
#define USE_NEWTON_CLOSEST_POINT
//#define DEBUG_CLOSEST_POINT
//#define DEBUG_EVAL_YIELD
//#define DEBUG_DF_DSIGMA
//#define DEBUG_POINT_ON_TANGENT
//#define TIME_POLY_SEARCH

using namespace Vaango;
using Point   = Uintah::Point;
using Vector  = Uintah::Vector;
using Matrix3 = Uintah::Matrix3;

YieldCond_TabularCap::YieldCond_TabularCap(Uintah::ProblemSpecP& ps,
                                           IntVar_TabularCap* intvar)
  : d_yield(ps)
{
  d_intvar = intvar;

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
  d_yield      = yc->d_yield;
  d_intvar     = yc->d_intvar;
  d_I1bar_min  = yc->d_I1bar_min;
  d_I1bar_max  = yc->d_I1bar_max;
  d_sqrtJ2_max = yc->d_sqrtJ2_max;
  d_polyline   = yc->d_polyline;
  d_normals    = yc->d_normals;
}

void
YieldCond_TabularCap::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP yield_ps = ps->appendChild("yield_condition");
  yield_ps->setAttribute("type", "tabular_cap");

  d_yield.table.outputProblemSpec(yield_ps);
  yield_ps->appendElement("cap_ellipticity_ratio", d_yield.capEllipticityRatio);
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
    d_yield.table.getIndependentVarData("Pressure", IndexKey(0, 0, 0, 0));
  DoubleVec1D yvals =
    d_yield.table.getDependentVarData("SqrtJ2", IndexKey(0, 0, 0, 0));
  if (xvals.size() != yvals.size()) {
    std::ostringstream out;
    out << "**ERROR** The tabular yield data file does not contain the"
        << " same number of \"Pressure\" and \"SqrtJ2\" data points";
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  // Copy the data and increase sampling if necessary
  std::vector<Point> points;
  auto num_pts = xvals.size();
  if (num_pts < 6) {
    for (auto ii = 0u; ii < num_pts - 1; ii++) {
      auto x_start = xvals[ii];
      auto x_end   = xvals[ii + 1];
      auto y_start = yvals[ii];
      auto y_end   = yvals[ii + 1];
      std::vector<double> x_vec;
      Vaango::Util::linspace(x_start, x_end, 100, x_vec);
      for (const auto& x : x_vec) {
        auto t = (x - x_start) / (x_end - x_start);
        auto y = (1 - t) * y_start + t * y_end;
        points.push_back(Point(x, y, 0));
      }
    }
  } else {
    for (auto ii = 0u; ii < xvals.size(); ii++) {
      points.push_back(Point(xvals[ii], yvals[ii], 0));
    }
  }
  // std::cout << "Orig = ";
  // std::copy(points.begin(), points.end(),
  //          std::ostream_iterator<Point>(std::cout, " "));
  // std::cout << std::endl;

  // Check convexity; If not convex, make convex
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

  if (d_yield.capEllipticityRatio >= 1 || d_yield.capEllipticityRatio <= 0.0) {
    std::ostringstream warn;
    warn << "capEllipticityRatio must be in [0, 1]. Input value = "
         << d_yield.capEllipticityRatio << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
}

/* If the yield surface does not change over time we can calculate the
   max sqrt(J2) at the beginning */
void
YieldCond_TabularCap::setYieldConditionRange()
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
YieldCond_TabularCap::saveAsPolyline()
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
  //std::copy(d_polyline.begin(),
  //          d_polyline.end(),
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
  // std::cout << "Normals:";
  // std::copy(d_normals.begin(), d_normals.end(),
  //          std::ostream_iterator<Vector>(std::cout, " "));
  // std::cout << std::endl;
}

Point
YieldCond_TabularCap::getClosestPoint(const Polyline& polyline,
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
// f := sqrt(J2) - g(p)*Fc(p, X_p) = 0
// where
//     J2 = 1/2 s:s,  s = sigma - p I,  p = 1/3 Tr(sigma)
//     g(pbar) = table,  pbar = -p
//     Fc^2 = 1 - (kappa - p)^2/(kappa - X_p)^2
//     kappa = p_tension - R*(p_tension - X_p)
//     X_p = hydrostatic strength = X/3
//
// Returns:
//   hasYielded = -1.0 (if elastic)
//              =  1.0 (otherwise)
//--------------------------------------------------------------
std::pair<double, Util::YieldStatus>
YieldCond_TabularCap::evalYieldCondition(const ModelStateBase* state_input)
{
  const ModelState_TabularCap* state =
    static_cast<const ModelState_TabularCap*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_TabularCap.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  // First check if the state is outside the tension cap
  double p_bar_min = d_I1bar_min / 3.0;
  double p_bar     = -state->I1 / 3;
  if (p_bar < p_bar_min) {
#ifdef DEBUG_EVAL_YIELD
    std::cout << "Tensile state: [" << p_bar << " < " << p_bar_min << "]\n";
#endif
    return std::make_pair(1.0, Util::YieldStatus::HAS_YIELDED);
  }

  // Next check if the state is outside the compression cap
  double p_bar_max = -state->capX / 3.0;
  if (p_bar > p_bar_max) {
#ifdef DEBUG_EVAL_YIELD
    std::cout << "Compressive state: [" << p_bar << " > " << p_bar_max << "]\n";
#endif
    return std::make_pair(1.0, Util::YieldStatus::HAS_YIELDED);
  }

  // Compute the value of the yield function for the yield function without cap
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

  // Find location of the center of the cap ellipse
  double kappa_bar =
    p_bar_min + d_yield.capEllipticityRatio * (p_bar_max - p_bar_min);

  // For the cap region only
  if (p_bar > kappa_bar) {

    double ratio = (p_bar - kappa_bar) / (p_bar_max - kappa_bar);
    double Fc_sq = 1.0 - ratio * ratio;
    if ((state->sqrt_J2 * state->sqrt_J2 - gg[0] * gg[0] * Fc_sq) > 1.0e-10) {
#ifdef DEBUG_EVAL_YIELD
      double J2     = state->sqrt_J2 * state->sqrt_J2;
      double Fc2_g2 = gg[0] * gg[0] * Fc_sq;
      std::cout << "Yielded: [" << p_bar << " > " << kappa_bar << "] "
                << " J2 = " << J2 << " Fc^2 g^2 = " << Fc2_g2
                << " diff = " << J2 - Fc2_g2 << "\n";
#endif
      return std::make_pair(1.0, Util::YieldStatus::HAS_YIELDED);
    }
  } else {
    if (state->sqrt_J2 > gg[0]) {
#ifdef DEBUG_EVAL_YIELD
      std::cout << "Yielded: [" << p_bar << " <= " << kappa_bar << "] "
                << " sqrt_J2 = " << state->sqrt_J2 << " g = " << gg[0] << "\n";
#endif
      return std::make_pair(1.0, Util::YieldStatus::HAS_YIELDED);
    }
  }

#ifdef DEBUG_EVAL_YIELD
  std::cout << "p_bar = " << p_bar << " gg = " << gg[0]
            << " sqrtJ2 = " << state->sqrt_J2 << std::endl;
#endif

#ifdef DO_CLOSEST_POINT_ELASTIC_CHECK
  auto status = checkClosestPointDistance(state);
  if (status == Util::YieldStatus::HAS_YIELDED) {
    return std::make_pair(1.0, Util::YieldStatus::HAS_YIELDED);
  }
#endif

  return std::make_pair(-1.0, Util::YieldStatus::IS_ELASTIC);
}

double
YieldCond_TabularCap::computeYieldFunction(const ModelStateBase* state) const
{
  std::ostringstream out;
  out << "**ERROR** The yield function for the tabular plasticity models"
      << " cannot be evaluated for a given stress state.\n";
  throw Uintah::InternalError(out.str(), __FILE__, __LINE__);

  return 0.0;
}

/* Add cap points to yield function table */
void
YieldCond_TabularCap::computeCapPoints(double X_bar, Polyline& p_q_all)
{
  // Set up limits
  double p_bar_min = d_I1bar_min / 3.0;
  double p_bar_max = X_bar / 3.0;
  double kappa_bar =
    p_bar_min + d_yield.capEllipticityRatio * (p_bar_max - p_bar_min);

  // Create a working copy of the polyline
  Polyline polyline_copy;
  std::copy(
    d_polyline.begin(), d_polyline.end(), std::back_inserter(polyline_copy));

  // Add further points beyond available data using linear interpolation if
  // needed
  auto end         = polyline_copy.end() - 1;
  double p_bar_end = (*end).x();
  if (p_bar_max > p_bar_end) {
    auto start          = end - 1;
    double p_bar_start  = (*start).x();
    double sqrtJ2_start = (*start).y();
    double sqrtJ2_end   = (*end).y();
    double dp_bar       = p_bar_end - p_bar_start;
    double curr_p_bar   = p_bar_end + dp_bar;
    double t            = (curr_p_bar - p_bar_start) / dp_bar;
    double curr_sqrt_J2 = (1 - t) * sqrtJ2_start + t * sqrtJ2_end;
    polyline_copy.emplace_back(Point(curr_p_bar, curr_sqrt_J2, 0.0));
    while (curr_p_bar < p_bar_max) {
      curr_p_bar += dp_bar;
      t            = (curr_p_bar - p_bar_start) / dp_bar;
      curr_sqrt_J2 = (1 - t) * sqrtJ2_start + t * sqrtJ2_end;
      polyline_copy.emplace_back(Point(curr_p_bar, curr_sqrt_J2, 0.0));
    }
  }

  // Find the location of the p_bar_start, p_bar_end on table polyline
  // (ignore first two points)
  auto end_iter = std::find_if(
    polyline_copy.begin() + 2,
    polyline_copy.end(),
    [&kappa_bar](const auto& point) { return kappa_bar < point.x(); });

  // Copy the relevant points
  std::copy(polyline_copy.begin(), end_iter, std::back_inserter(p_q_all));
  // std::cout << "p_q_cone = ";
  // std::copy(p_q_all.begin(), p_q_all.end(),
  //            std::ostream_iterator<Point>(std::cout, " "));
  // std::cout << std::endl;

  // Set up default theta incremenets
  int num_theta = 180;
  double theta_inc = M_PI / (2.0 * num_theta);

  // Set up ellipse axes
  double a = p_bar_max - kappa_bar;

  /*
  // Compute distance incremenet of polyline
  auto last = p_q_all.rbegin();
  //auto last = d_polyline.rbegin();
  auto last_but_one = last - 1;
  auto dist_inc_poly = (*last - *last_but_one).length();

  // Compute length of theta_inc arc
  auto x_inc = kappa_bar + a * cos(theta_inc);
  auto b_inc = computeEllipseHeight(d_polyline, x_inc);
  auto y_inc = b_inc * sin(theta_inc);

  auto dist_inc_theta = std::sqrt((x_inc - kappa_bar)*(x_inc - kappa_bar) + y_inc * y_inc);
  //std::cout << "dist_inc polyline = " << dist_inc_poly << " dist_inc theta = " << dist_inc_theta << "\n";

  // Adjust theta increments
  if (dist_inc_theta > dist_inc_poly) {
    num_theta = std::max(num_theta, num_theta*static_cast<int>(std::ceil(dist_inc_theta/dist_inc_poly)));
  }
  theta_inc = M_PI / (2.0 * num_theta);
  */

  /*
  std::cout << "dist_inc_theta = " << dist_inc_theta << " dist_inc_poly = " << dist_inc_poly
            << " num_theta = " << num_theta << " theta_inc = " << theta_inc << "\n";
  */

  // Set up theta vector
  std::vector<double> theta_vec;
  Vaango::Util::linspace(0, M_PI / 2, num_theta, theta_vec);
  std::reverse(std::begin(theta_vec), std::end(theta_vec));
  theta_vec.emplace_back(-theta_inc);
  theta_vec.emplace_back(-2.0*theta_inc);

  // Compute ellipse points
  Polyline p_q_cap;
  for (auto theta : theta_vec) {
    auto x = kappa_bar + a * cos(theta);
    auto b = computeEllipseHeight(d_polyline, x);
    auto y = b * sin(theta);
    p_q_cap.emplace_back(Uintah::Point(x, y, 0));
  }

  // Concatenate the two vectors
  p_q_all.insert(p_q_all.end(), p_q_cap.begin(), p_q_cap.end());
  // std::cout << "p_q_all = ";
  // std::copy(p_q_all.begin(), p_q_all.end(),
  //            std::ostream_iterator<Point>(std::cout, " "));
  // std::cout << std::endl;
}

double
YieldCond_TabularCap::computeEllipseHeight(const Polyline& p_q_points,
                                           double p_cap)
{
  // Find the location of the p_bar_start, p_bar_end on table polyline
  // (ignore first two points)
  auto end_iter =
    std::find_if(d_polyline.begin() + 2,
                 d_polyline.end(),
                 [&p_cap](const auto& point) { return p_cap < point.x(); });

  // Compute sqrtJ2 at that value of p_cap
  auto start_iter = end_iter - 1;
  if (end_iter == d_polyline.end()) {
    start_iter--;
    end_iter--;
  }
  auto start    = *start_iter;
  auto end      = *end_iter;
  double t      = (p_cap - start.x()) / (end.x() - start.x());
  double sqrtJ2 = (1 - t) * start.y() + t * end.y();

  return sqrtJ2;
}

Util::YieldStatus
YieldCond_TabularCap::checkClosestPointDistance(
  const ModelState_TabularCap* state)
{
  double p       = state->I1 / 3;
  double sqrt_J2 = state->sqrt_J2;
  Point trial_pt(p, sqrt_J2, 0.0);

  // Convert the yield surface points (p, sqrtJ2) into tension +ve form
  Polyline yield_f_pts = state->yield_f_pts;
  for (auto& pt : yield_f_pts) {
    pt.x(-pt.x());
  }

  // Find the closest segments
  Polyline p_J2_segments;
  std::size_t closest_index = 0;
#ifdef DO_BINARY_CLOSEST_SEGMENT
  closest_index = Vaango::Util::getClosestSegmentsBinarySearch(trial_pt, yield_f_pts, 
                                                               p_J2_segments);
#else
  if (yield_f_pts.size() < Vaango::Util::NUM_PTS_KDTREE_SWITCH) {
    closest_index = Vaango::Util::getClosestSegments(trial_pt, yield_f_pts, 
                                                     p_J2_segments);
  } else {
    closest_index = Vaango::Util::getClosestSegmentsKDTree(trial_pt, yield_f_pts, 
                                                           p_J2_segments);
  }
#endif

  // Find the closest point on the quadratic B-spline through closest segments
  std::size_t numPts = yield_f_pts.size();
  auto seg_start     = closest_index - 1;
  auto seg_end       = closest_index + 1;
  if (closest_index < 2) {
    seg_start = 0;
    seg_end   = 2;
  } else if (closest_index > numPts - 3) {
    seg_start = numPts - 3;
    seg_end   = numPts - 1;
  }

  Point closest_to_spline(0, 0, 0);
  Vector tangent(0, 0, 0);
  Vector dtangent(0, 0, 0);
  std::tie(closest_to_spline, tangent, dtangent) =
    Vaango::Util::computeClosestPointQuadraticBSpline(
      trial_pt, yield_f_pts, seg_start, seg_end);

  // Compute orientation of polyline
  Polyline poly(5);
  poly[0] = yield_f_pts[seg_start];
  poly[1] = closest_to_spline;
  poly[2] = yield_f_pts[seg_end];
  poly[3] = trial_pt;
  poly[4] = poly[0];
  
  double orientation = (poly[1].x() - poly[0].x()) * (poly[1].y() + poly[0].y()) + 
                       (poly[2].x() - poly[1].x()) * (poly[2].y() + poly[1].y()) + 
                       (poly[3].x() - poly[2].x()) * (poly[3].y() + poly[2].y()) + 
                       (poly[4].x() - poly[3].x()) * (poly[4].y() + poly[3].y());
  #ifdef DEBUG_CLOSEST_POINT_ELASTIC
    std::cout << "orientation (> 0 = plastic, < 0 = elastic): " << orientation << "\n";
  #endif
  if (orientation > 0.0) {
    return Util::YieldStatus::HAS_YIELDED;
  }

  /*
  Polyline line;
  line.emplace_back(yield_f_pts[seg_start]);
  line.emplace_back(yield_f_pts[seg_end]);

  Point closest_to_line;
  double distSq =
    Vaango::Util::findClosestPoint(trial_pt, line, closest_to_line);
  double distSq_spline = (closest_to_spline - trial_pt).length2();
  double diff_dist     = (distSq > 0) ? std::sqrt(distSq_spline / distSq) - 1
                                  : distSq_spline - distSq;

  if (std::abs(diff_dist) > 1.0e-6) {
    if (diff_dist < 0.0) {
      double dx = closest_to_spline.x() - closest_to_line.x();
      double dy = closest_to_spline.y() - closest_to_line.y();
      double tx = (dx != 0) ? (trial_pt.x() - closest_to_line.x()) / dx : 0.0;
      double ty = (dy != 0) ? (trial_pt.y() - closest_to_line.y()) / dy : 0.0;
      if (!(Util::isInBounds<double>(tx, 0, 1) &&
            Util::isInBounds<double>(ty, 0, 1))) {
        #ifdef DEBUG_CLOSEST_POINT_ELASTIC
        std::cout << "tx = " << tx << " ty = " << ty << "\n";
        std::cout << "closest_index = " << closest_index << "\n";
        std::cout << "indices are (" << seg_start << "," << seg_end << ") from"
                  << "(0," << numPts - 1 << ")\n";
        std::cout << "line = [" << line[0] << ";" << line[1] << "];\n";
        std::cout << "seg = [" << yield_f_pts[seg_start] << ";"
                  << yield_f_pts[seg_start + 1] << ";"
                  << yield_f_pts[seg_start + 2] << "];\n";
        std::cout << "pt = " << trial_pt << ";"
                  << " cspline = " << closest_to_spline << ";"
                  << " cline = " << closest_to_line << ";\n";
        std::cout << std::setprecision(10) << "p-spline = " << distSq_spline
                  << " p-line = " << distSq << " diff = " << diff_dist << "\n";
        std::cout << "orientation distance = plastic\n";
        #endif
        return Util::YieldStatus::HAS_YIELDED;
      }
    }
  }
  #ifdef DEBUG_CLOSEST_POINT_ELASTIC
    std::cout << "orientation distance = elastic\n";
  #endif
  */

  return Util::YieldStatus::IS_ELASTIC;
}

//--------------------------------------------------------------
// Evaluate yield condition max value of sqrtJ2
//--------------------------------------------------------------
double
YieldCond_TabularCap::evalYieldConditionMax(const ModelStateBase* state_input)
{
  const ModelState_TabularCap* state =
    static_cast<const ModelState_TabularCap*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_TabularCap.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  std::array<double, 3> vals = getYieldConditionRange(state->yield_f_pts);
  return vals[2];
}

//--------------------------------------------------------------
// Derivatives needed by return algorithms and Newton iterations

//--------------------------------------------------------------
/*! Compute Derivative with respect to the Cauchy stress (\f$\sigma \f$)
 *  Compute df/dsigma
 *
 *  for the yield function
 *      f := sqrt(J2(s)) - g(p)*Fc(p, X_p) = 0
 *  where
 *      J2 = 1/2 s:s,  s = sigma - p I,  p = 1/3 Tr(sigma)
 *      g(pbar) = table, pbar = -p
 *      X_p = 1/3 X
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
Uintah::Matrix3
YieldCond_TabularCap::df_dsigma(const ModelStateBase* state_input) 
{
  const ModelState_TabularCap* state =
    static_cast<const ModelState_TabularCap*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_TabularCap.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  if (state->yield_f_pts.empty()) {
    std::ostringstream out;
    out << "**ERROR** The yield surface with cap polyline has not been "
           "initialized.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  #ifdef DEBUG_DF_DSIGMA
  std::cout << "p = " << state->I1/3 << " sqrtJ2 = " << state->sqrt_J2
            << " s = " << state->deviatoricStressTensor << "\n";
  #endif

  // Get the bulk and shear moduli and compute sqrt(3/2 K/G)
  double sqrtKG = std::sqrt(1.5 * state->bulkModulus / state->shearModulus);

  // Find the closest point in the transformed 2D stress space
  double p_bar   = -state->I1 / 3;
  double sqrt_J2 = state->sqrt_J2;
  double z = 0.0, rprime = 0.0;
  Vaango::Util::convertToZRprime(sqrtKG, p_bar, sqrt_J2, z, rprime);

  double closest_z = 0.0, closest_rprime = 0.0;
  double tangent_z = 0.0, tangent_rprime = 0.0;
  getClosestPointAndTangent(
    state, z, rprime, closest_z, closest_rprime, tangent_z, tangent_rprime);

  double closest_p_bar = 0.0, closest_sqrt_J2 = 0.0;
  double tangent_p_bar = 0.0, tangent_sqrt_J2 = 0.0;
  Vaango::Util::revertFromZRprime(
    sqrtKG, closest_z, closest_rprime, closest_p_bar, closest_sqrt_J2);
  Vaango::Util::revertFromZRprime(
    sqrtKG, tangent_z, tangent_rprime, tangent_p_bar, tangent_sqrt_J2);

  // Check that the closest point is not at the vertex
  double p_bar_min = d_I1bar_min / 3.0;
  double epsilon = 1.0e-6;
  if (closest_p_bar - epsilon < p_bar_min) {
    return Util::Identity * Util::large_number;
  }

  // Handle vertex (minimum p_bar : assume vertex is in tension)
  if (closest_p_bar < 0 && std::abs(closest_sqrt_J2) < 1.0e-8) {
    //std::cout << std::setprecision(16) << 
    //          "df_dsigma = positive, p_bar = " << p_bar << " closest_p_bar = " << closest_p_bar << "\n";
    return Util::Identity * (Util::large_number);
  }

  // Handle cap (maximum p_bar : assume cap is in compression)
  if (closest_p_bar > 0 && std::abs(closest_sqrt_J2) < 1.0e-8) {
    //std::cout << std::setprecision(16) << 
    //          "df_dsigma = negative, p_bar = " << p_bar << " closest_p_bar = " << closest_p_bar << "\n";
    return Util::Identity * (-Util::large_number);
  }

  #ifdef DEBUG_POINT_ON_TANGENT
    Point closest_zr(closest_z, closest_rprime, 0.0);
    Vector tangent_zr(tangent_z, tangent_rprime, 0.0);
    tangent_zr.normalize();

    Point pt_on_tangent_zr = closest_zr + tangent_zr * 1000;

    double pt_on_tangent_p_bar = 0.0, pt_on_tangent_sqrt_J2 = 0.0;
    Vaango::Util::revertFromZRprime(
      sqrtKG, pt_on_tangent_zr.x(), pt_on_tangent_zr.y(), pt_on_tangent_p_bar, pt_on_tangent_sqrt_J2);

    Point pt_on_tangent(-pt_on_tangent_p_bar, pt_on_tangent_sqrt_J2, 0);
    Point closest(-closest_p_bar, closest_sqrt_J2, 0);
    Vector tangent = pt_on_tangent - closest;
    tangent.normalize();
    Vector normal(-tangent[1], tangent[0], 0.0);
    //std::cout << " pt_tangent = " << tangent.x() << ", " << tangent.y() 
    //          << " tangent = " << tangent_p_bar << "," << tangent_sqrt_J2 
    //          << " normal = " << normal.x() << "," << normal.y() << "\n";

    Point pt_on_normal = closest - normal * 1000;
    //std::cout << "pt_on_tangent = " << pt_on_tangent
    //          << " pt_on_normal = " << pt_on_normal << " closest = " << closest << "\n";

    Matrix3 pt_on_normal_9D = Util::Identity * pt_on_normal.x() +
                              state->deviatoricStressTensor * pt_on_normal.y() / state->sqrt_J2;
    Matrix3 closest_9D = Util::Identity * closest.x() +
                              state->deviatoricStressTensor * closest.y() / state->sqrt_J2;

    //std::cout << "pt_on_normal_9D = " << pt_on_normal_9D << "\nclosest_9D = " << closest_9D << "\n";

    double pt_on_normal_9D_2D_p = pt_on_normal_9D.Trace() / 3.0;
    Matrix3 s_normal = pt_on_normal_9D - Util::Identity * pt_on_normal_9D_2D_p;
    double pt_on_normal_9D_2D_sqrt_J2 = std::sqrt(0.5 * s_normal.Contract(s_normal));
    //std::cout << "pt_on_normal_9D_2D = " << pt_on_normal_9D_2D_p << ", " << pt_on_normal_9D_2D_sqrt_J2 << "\n";

    double closest_9D_2D_p = closest_9D.Trace() / 3.0;
    Matrix3 s = closest_9D - Util::Identity * closest_9D_2D_p;
    double closest_9D_2D_sqrt_J2 = std::sqrt(0.5 * s.Contract(s));
    //std::cout << "closest_9D_2D = " << closest_9D_2D_p << ", " << closest_9D_2D_sqrt_J2 << "\n";

    Matrix3 normal_9D = pt_on_normal_9D - closest_9D;
    normal_9D /= normal_9D.Norm();
    //std::cout << "sqrtJ2 = " << state->sqrt_J2 << " normal_9D = " << normal_9D << "\n";

    Matrix3 pt_on_normal_back = closest_9D + normal_9D * 1000;
    //std::cout << "pt_on_normal_back_9D = " << pt_on_normal_back << "\n";
    double pt_on_normal_back_p = pt_on_normal_back.Trace() / 3.0;
    Matrix3 s_normal_back = pt_on_normal_back - Util::Identity * pt_on_normal_back_p;
    double pt_on_normal_back_sqrt_J2 = std::sqrt(0.5 * s_normal_back.Contract(s_normal_back));
    //std::cout << "pt_on_normal_back = " << pt_on_normal_back_p << ", " << pt_on_normal_back_sqrt_J2 << "\n";
  #endif

  #ifdef DEBUG_DF_DSIGMA
  std::cout << "p_bar = " << p_bar << " sqrtJ2 = " << state->sqrt_J2
            << " closest = " << closest_p_bar << "," << closest_sqrt_J2
            << " tangent = " << tangent_p_bar << "," << tangent_sqrt_J2 <<
            "\n";
  #endif

  // Compute df_dp (cases where tangent_p_bar == 0, at the vertex and cap, have been handled before)
  double dfdp_bar = std::copysign(Util::large_number, tangent_sqrt_J2);
  if (std::abs(tangent_p_bar) > 1.0e-16) {
    dfdp_bar = tangent_sqrt_J2 / tangent_p_bar;
  }
  
  // Compute df_dsqrt(J2)
  double dfdJ2 =
    (closest_sqrt_J2 == 0) ? Util::large_number : 1 / (2 * closest_sqrt_J2);

  #ifdef DEBUG_DF_DSIGMA
  std::cout << "df_dp_bar = " << dfdp_bar << " df_dJ2 = " << dfdJ2 << "\n";
  #endif

  // Compute df_dsigma (dfdp = dfdp_bar for outward normal, -dfdp_bar for inward normal)
  Matrix3 p_term = Util::Identity * (dfdp_bar / 3.0);
  Matrix3 df_dsigma = p_term;

  #ifdef DEBUG_DF_DSIGMA
  std::cout << "p_term = " << p_term << "\n";
  std::cout << "df_dsigma = " << df_dsigma << "\n";
  #endif

  if (state->sqrt_J2 > 0) {
     // Compute stress tensor at closest point
     double s_factor = closest_sqrt_J2 / state->sqrt_J2;
     Matrix3 s_closest = state->deviatoricStressTensor * s_factor;
     Matrix3 s_term = s_closest * (dfdJ2);
     df_dsigma += s_term;

     #ifdef DEBUG_DF_DSIGMA
     std::cout << "s_factor = " << s_factor << "\n";
     std::cout << "s_closest = " << s_closest << "\n";
     std::cout << "s_term = " << s_term << "\n";
     std::cout << "df_dsigma = " << df_dsigma << "\n";
     #endif
  } 

  #ifdef DEBUG_DF_DSIGMA
  std::cout << "s = " << state->deviatoricStressTensor << "\n";
  std::cout << "df_dsigma = " << df_dsigma << "\n";
  #endif

  return df_dsigma;
}

/* Assumes closest point and tangent have already been computed
   elsewhere and the state variables have been updated */
Uintah::Matrix3
YieldCond_TabularCap::df_dsigma(const Matrix3&,
                                const ModelStateBase* state_input)
{
  const ModelState_TabularCap* state =
    static_cast<const ModelState_TabularCap*>(state_input);

  double closest_p_bar = state->closest.x();
  double closest_sqrt_J2 = state->closest.y();
  double tangent_p_bar = state->tangent.x();
  double tangent_sqrt_J2 = state->tangent.y();

  // Check that the closest point is not at the vertex
  double p_bar_min = d_I1bar_min / 3.0;
  double epsilon = 1.0e-6;
  if (closest_p_bar - epsilon < p_bar_min) {
    return Util::Identity * Util::large_number;
  }

  // Handle vertex (minimum p_bar : assume vertex is in tension)
  if (closest_p_bar < 0 && std::abs(closest_sqrt_J2) < 1.0e-8) {
    //std::cout << std::setprecision(16) << 
    //          "df_dsigma = positive, p_bar = " << p_bar << " closest_p_bar = " << closest_p_bar << "\n";
    return Util::Identity * (Util::large_number);
  }

  // Handle cap (maximum p_bar : assume cap is in compression)
  if (closest_p_bar > 0 && std::abs(closest_sqrt_J2) < 1.0e-8) {
    //std::cout << std::setprecision(16) << 
    //          "df_dsigma = negative, p_bar = " << p_bar << " closest_p_bar = " << closest_p_bar << "\n";
    return Util::Identity * (-Util::large_number);
  }

  // Compute df_dp (cases where tangent_p_bar == 0, at the vertex and cap, have been handled before)
  double dfdp_bar = std::copysign(Util::large_number, tangent_sqrt_J2);
  if (std::abs(tangent_p_bar) > 1.0e-16) {
    dfdp_bar = tangent_sqrt_J2 / tangent_p_bar;
  }
  
  // Compute df_dsqrt(J2)
  double dfdJ2 =
    (closest_sqrt_J2 == 0) ? Util::large_number : 1 / (2 * closest_sqrt_J2);

  // Compute df_dsigma (dfdp = dfdp_bar for outward normal, -dfdp_bar for inward normal)
  Matrix3 p_term = Util::Identity * (dfdp_bar / 3.0);
  Matrix3 df_dsigma = p_term;

  if (state->sqrt_J2 > 0) {
     // Compute stress tensor at closest point
     double s_factor = closest_sqrt_J2 / state->sqrt_J2;
     Matrix3 s_closest = state->deviatoricStressTensor * s_factor;
     Matrix3 s_term = s_closest * (dfdJ2);
     df_dsigma += s_term;
  } 

  return df_dsigma;
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
//     X_p = 1/3 X
//
// the derivative is
//     df/dp = -dg/dp * F_c(p, X_p) - g(p) * dF_c/dp
//     dg/dp = dg/dpbar dpbar/dp  = -dg/dpbar
//     dFc/dp = (1/F_c)*(kappa - p)/(kappa - X_p)^2
//--------------------------------------------------------------
double
YieldCond_TabularCap::df_dp(const ModelStateBase* state_input)
{
  const ModelState_TabularCap* state =
    static_cast<const ModelState_TabularCap*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_TabularCap.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  // Get the bulk and shear moduli and compute sqrt(3/2 K/G)
  double sqrtKG = std::sqrt(1.5 * state->bulkModulus / state->shearModulus);

  // Compute dg/dp
  double p_bar   = -state->I1 / 3;
  double sqrt_J2 = state->sqrt_J2;
  double z = 0.0, rprime = 0.0;
  Vaango::Util::convertToZRprime(sqrtKG, p_bar, sqrt_J2, z, rprime);

  double closest_z = 0.0, closest_rprime = 0.0;
  double tangent_z = 0.0, tangent_rprime = 0.0;
  getClosestPointAndTangent(
    state, z, rprime, closest_z, closest_rprime, tangent_z, tangent_rprime);

  double closest_p_bar = 0.0, closest_sqrt_J2 = 0.0;
  double tangent_p_bar = 0.0, tangent_sqrt_J2 = 0.0;
  Vaango::Util::revertFromZRprime(
    sqrtKG, closest_z, closest_rprime, closest_p_bar, closest_sqrt_J2);
  Vaango::Util::revertFromZRprime(
    sqrtKG, tangent_z, tangent_rprime, tangent_p_bar, tangent_sqrt_J2);

  #ifdef DEBUG_CLOSEST_POINT
  std::cout << "p_bar = " << p_bar << " sqrtJ2 = " << state->sqrt_J2
            << " closest = " << closest_p_bar << "," << closest_sqrt_J2
            << " tangent = " << tangent_p_bar << "," << tangent_sqrt_J2 <<
            "\n";
  #endif

  // Check that the closest point is not at the vertex
  double p_bar_min = d_I1bar_min / 3.0;
  double epsilon = 1.0e-6;
  if (closest_p_bar - epsilon < p_bar_min) {
    return Util::large_number;
  }

  // Compute df_dp
  double dg_dpbar =
    (tangent_p_bar == 0) ? Util::large_number : tangent_sqrt_J2 / tangent_p_bar;
  double df_dp = dg_dpbar;

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
YieldCond_TabularCap::df_dq(const ModelStateBase* state_input)
{
  const ModelState_TabularCap* state =
    static_cast<const ModelState_TabularCap*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_TabularCap.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  // Get the bulk and shear moduli and compute sqrt(3/2 K/G)
  double sqrtKG = std::sqrt(1.5 * state->bulkModulus / state->shearModulus);

  // Compute closest point on yield surface
  double p_bar   = -state->I1 / 3;
  double sqrt_J2 = state->sqrt_J2;
  double z = 0.0, rprime = 0.0;
  Vaango::Util::convertToZRprime(sqrtKG, p_bar, sqrt_J2, z, rprime);

  double closest_z = 0.0, closest_rprime = 0.0;
  double tangent_z = 0.0, tangent_rprime = 0.0;
  getClosestPointAndTangent(
    state, z, rprime, closest_z, closest_rprime, tangent_z, tangent_rprime);

  double closest_p_bar = 0.0, closest_sqrt_J2 = 0.0;
  Vaango::Util::revertFromZRprime(
    sqrtKG, closest_z, closest_rprime, closest_p_bar, closest_sqrt_J2);

  #ifdef DEBUG_CLOSEST_POINT
  std::cout << "p_bar = " << p_bar << " sqrtJ2 = " << state->sqrt_J2
            << " closest = " << closest_p_bar << "," << closest_sqrt_J2 << 
            "\n";
  #endif

  double df_dJ2 =
    (closest_sqrt_J2 == 0) ? Util::large_number : 1 / (2 * closest_sqrt_J2);

  return df_dJ2;
}

//--------------------------------------------------------------
// Compute df/deps^p_v
//
//   for the yield function
//     f := sqrt(J2) - g(p)*Fc(p, X_p) = 0
//   where
//     J2 = 1/2 s:s,  s = sigma - p I,  p = 1/3 Tr(sigma)
//     g(pbar) := table,  pbar = -p
//     Fc^2 = 1 - (kappa - p)^2/(kappa - X_p)^2
//     kappa = p_tension - R*(p_tension - X_p)
//     X_p = 1/3 X
//
//   the derivative is
//     df/deps^p_v = df/dX_p*dX_p/deps^p_v
//
//     df/dX_p = -g(p)*dFc/dX_p
//     dX_p/deps^p_v = table interpolation
//
// Requires:  internal variable model
//--------------------------------------------------------------
double
YieldCond_TabularCap::df_depsVol(const ModelStateBase* state_input,
                                 const MPMEquationOfState*,
                                 const ShearModulusModel*)
{
  const ModelState_TabularCap* state =
    static_cast<const ModelState_TabularCap*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_TabularCap.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  // Set up limits
  double p_bar_min = d_I1bar_min / 3.0;
  double X_bar_p   = -state->capX / 3.0;

  // Find location of the center of the cap ellipse
  double kappa_bar =
    p_bar_min + d_yield.capEllipticityRatio * (X_bar_p - p_bar_min);

  // Compute g(p)
  double p_bar = -state->I1 / 3.0;
  DoubleVec1D gg;
  try {
    gg = d_yield.table.interpolate<1>({ { p_bar } });
  } catch (Uintah::InvalidValue& e) {
    std::ostringstream out;
    out << "**ERROR** In compute g(p):"
        << " p_bar = " << p_bar << "\n"
        << e.message();
    throw Uintah::InvalidValue(out.str(), __FILE__, __LINE__);
  }

  // Compute df_dX_p
  double df_dX_p = 0.0;
  if (p_bar > kappa_bar) {
    double denom = X_bar_p - p_bar_min;
    if (denom == 0.0) {
      df_dX_p = Util::large_number;
    } else {
      double numer     = p_bar - kappa_bar;
      double inv_denom = 1.0 / denom;
      double ratio     = numer * inv_denom;
      // double Fc = std::sqrt(1.0 - ratio * ratio);
      double dFc_dX_p =
        ratio / (1 + ratio) * d_yield.capEllipticityRatio * inv_denom;
      df_dX_p = -gg[0] * dFc_dX_p;
    }
  }
  double dX_p_dep_v =
    d_intvar->computeVolStrainDerivOfInternalVariable("capX", state);
  double df_dep_v = df_dX_p * dX_p_dep_v;
  return df_dep_v;
}

/**
 * Function: computeYieldSurfacePolylinePbarSqrtJ2
 *
 * Purpose: Compute a sequence of points representing the yield surface
 *          in pbar-sqrtJ2 space
 *
 * Inputs:
 *  state_old = old state
 *
 * Returns:
 *   std::vector<Point>
 */
Polyline
YieldCond_TabularCap::computeYieldSurfacePolylinePbarSqrtJ2(
  const ModelStateBase* state_input)
{
  const ModelState_TabularCap* state =
    static_cast<const ModelState_TabularCap*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_TabularCap.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  Polyline yield_f_pts;
  computeCapPoints(-state->capX, yield_f_pts);
  return yield_f_pts;
}

/**
 * Function: getUpdatedYieldConditionRange
 *
 * Purpose: Compute range of the yield surface in pbar-sqrtJ2 space
 *          Yield surface changes over time.  So we need to update the range./
 *
 * Inputs:
 *  std::vector<Point>
 *
 * Returns:
 *   std::array<double, 3>  = pbar_min, pbar_max, sqrtJ2_max
 */
std::array<double, 3>
YieldCond_TabularCap::getYieldConditionRange(const Polyline& yield_surface)
{

  if (yield_surface.empty()) {
    std::ostringstream out;
    out << "**ERROR** The yield surface with cap polyline has not been "
           "initialized.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  auto min_max_p = std::minmax_element(
    yield_surface.begin(),
    yield_surface.end(),
    [](const auto& pt1, const auto& p2) { return pt1.x() < p2.x(); });
  double I1bar_min = min_max_p.first->x() * 3.0;
  double I1bar_max = min_max_p.second->x() * 3.0;

  auto max_q = std::max_element(
    yield_surface.begin(),
    yield_surface.end(),
    [](const auto& pt1, const auto& p2) { return pt1.y() < p2.y(); });
  double sqrtJ2_max = max_q->y();
  return std::array<double, 3>(
    { { I1bar_min / 3.0, I1bar_max / 3.0, sqrtJ2_max } });
}

/**
 * Function: getClosestPoint
 * Purpose: Get the point on the yield surface that is closest to a given point
 * (2D)
 * Inputs:
 *  state = current state
 *  z = z-coordinate of point
 *  rprime = r'-coordinate of point
 * Outputs:
 *  cz = z-coordinate of closest point on yield surface
 *  crprime = r'-coordinate of closest point
 */
bool
YieldCond_TabularCap::getClosestPoint(const ModelStateBase* state_input,
                                      const double& z,
                                      const double& rprime,
                                      double& cz,
                                      double& crprime)
{
  const ModelState_TabularCap* state =
    static_cast<const ModelState_TabularCap*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_TabularCap.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  if (state->yield_f_pts.empty()) {
    std::ostringstream out;
    out << "**ERROR** The yield surface with cap polyline has not been "
           "initialized.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  // std::chrono::time_point<std::chrono::system_clock> start, end;
  // start = std::chrono::system_clock::now();
  Point pt(z, rprime, 0.0);
  Point closest(0.0, 0.0, 0.0);
  if (state->yield_f_pts.size() < 5) {
    closest = getClosestPointTable(state, pt);
  } else {
#ifdef USE_NEWTON_CLOSEST_POINT
    Vector tangent(0.0, 0.0, 0.0);
    std::tie(closest, tangent) = getClosestPointSplineNewton(state, pt);
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

bool
YieldCond_TabularCap::getClosestPointAndTangent(
  const ModelStateBase* state_input,
  const double& z,
  const double& rprime,
  double& cz,
  double& crprime,
  double& tz,
  double& trprime)
{
  const ModelState_TabularCap* state =
    static_cast<const ModelState_TabularCap*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_TabularCap.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  if (state->yield_f_pts.empty()) {
    std::ostringstream out;
    out << "**ERROR** The yield surface with cap polyline has not been "
           "initialized.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  Point pt(z, rprime, 0.0);
  Point closest(0.0, 0.0, 0.0);
  Vector tangent(0.0, 0.0, 0.0);

  std::tie(closest, tangent) = getClosestPointSplineNewton(state, pt);

  cz      = closest.x();
  crprime = closest.y();

  tz      = tangent.x();
  trprime = tangent.y();

  return true;
}

bool
YieldCond_TabularCap::getClosestPointAndTangent(
  const ModelStateBase* state_input,
  const Polyline& z_r_table, 
  const Util::PolylineKDTree& z_r_index, 
  const double& z,
  const double& rprime,
  double& cz,
  double& crprime,
  double& tz,
  double& trprime)
{
  const ModelState_TabularCap* state =
    static_cast<const ModelState_TabularCap*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_TabularCap.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  if (state->yield_f_pts.empty()) {
    std::ostringstream out;
    out << "**ERROR** The yield surface with cap polyline has not been "
           "initialized.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  Point pt(z, rprime, 0.0);
  Point closest(0.0, 0.0, 0.0);
  Vector tangent(0.0, 0.0, 0.0);

  std::tie(closest, tangent) = getClosestPointSplineNewtonZR(state, z_r_table, z_r_index, pt);

  cz      = closest.x();
  crprime = closest.y();

  tz      = tangent.x();
  trprime = tangent.y();

  return true;
}

Point
YieldCond_TabularCap::getClosestPointTable(const ModelState_TabularCap* state,
                                           const Point& z_r_pt)
{
  // Get the bulk and shear moduli and compute sqrt(3/2 K/G)
  double sqrtKG = std::sqrt(1.5 * state->bulkModulus / state->shearModulus);

  // Convert tabular data to z-rprime coordinates
  Polyline z_r_table;
  convertToZRprime(sqrtKG, state->yield_f_pts, z_r_table);

  // Find the closest point
  Point z_r_closest;
  Vaango::Util::findClosestPoint(z_r_pt, z_r_table, z_r_closest);

  return z_r_closest;
}

Point
YieldCond_TabularCap::getClosestPointSpline(const ModelState_TabularCap* state,
                                            const Point& z_r_pt)
{
  // Get the bulk and shear moduli and compute sqrt(3/2 K/G)
  double sqrtKG = std::sqrt(1.5 * state->bulkModulus / state->shearModulus);

  // Convert tabular data to z-rprime coordinates
  Polyline z_r_table;
  convertToZRprime(sqrtKG, state->yield_f_pts, z_r_table);

#ifdef DEBUG_CLOSEST_POINT
  std::cout << "z_r_x = (";
  for (auto pt : z_r_table) {
    std::cout << pt.x() << ",";
  }
  std::cout << ")\n";
  std::cout << "z_r_y = (";
  for (auto pt : z_r_table) {
    std::cout << pt.y() << ",";
  }
  std::cout << ")\n";
#endif

  // Find the closest segments
  Polyline z_r_segments;
  std::size_t closest_index = 0;
#ifdef DO_BINARY_CLOSEST_SEGMENT
  closest_index = Vaango::Util::getClosestSegmentsBinarySearch(z_r_pt, z_r_table, 
                                                               z_r_segments);
#else
  if (z_r_table.size() < Vaango::Util::NUM_PTS_KDTREE_SWITCH) {
    closest_index = Vaango::Util::getClosestSegments(z_r_pt, z_r_table, 
                                                     z_r_segments);
  } else {
    closest_index = Vaango::Util::getClosestSegmentsKDTree(z_r_pt, z_r_table, 
                                                           z_r_segments);
  }
#endif

  // Get the yield surface points for the closest segments
  // (Fit quadratic B_spline)
  std::size_t numPts        = z_r_table.size();
  std::size_t ptsPerSegment = 30;
  Polyline z_r_spline;
  auto seg_start = closest_index - 1;
  auto seg_end   = closest_index + 1;
  // if (closest_index < 1) {
  //  seg_start = 0;
  //  seg_end = 1;
  //} else if (closest_index > numPts-3) {
  //  seg_start = numPts-2;
  //  seg_end = numPts-1;
  //}
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

  Vaango::Util::computeOpenUniformQuadraticBSpline(
    z_r_table, seg_start, seg_end, ptsPerSegment, z_r_spline);

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

  if (z_r_spline.empty()) {
    std::ostringstream out;
    out << "**ERROR** Could not fit spline to input yield surface segment: "
        << "indices are (" << seg_start << "," << seg_end << ") from"
        << "(0," << numPts - 1 << ")";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
#endif

  // Find the closest point
  Point z_r_closest;
  Vaango::Util::findClosestPoint(z_r_pt, z_r_spline, z_r_closest);

#ifdef DEBUG_CLOSEST_POINT
  std::cout << "ZRClose = " << std::setprecision(16) << z_r_closest << "\n";
#endif

  return z_r_closest;
}

std::tuple<Point, Vector>
YieldCond_TabularCap::getClosestPointSplineNewton(
  const ModelState_TabularCap* state,
  const Point& z_r_pt)
{
  // Get the bulk and shear moduli and compute sqrt(3/2 K/G)
  double sqrtKG = std::sqrt(1.5 * state->bulkModulus / state->shearModulus);

  #ifdef TIME_POLY_SEARCH
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
  #endif

  // Convert tabular data to z-rprime coordinates
  Polyline z_r_table(state->yield_f_pts.size());
  convertToZRprime(sqrtKG, state->yield_f_pts, z_r_table);

  #ifdef TIME_POLY_SEARCH
    end = std::chrono::system_clock::now();
    std::cout << "ZR poly creation : Time taken = " 
              << std::chrono::duration<double>(end-start).count() << std::endl;
  #endif

  #ifdef DEBUG_CLOSEST_POINT
    std::cout << "z_r_x = (";
    for (auto pt : z_r_table) {
      std::cout << pt.x() << ",";
    }
    std::cout << ")\n";
    std::cout << "z_r_y = (";
    for (auto pt : z_r_table) {
      std::cout << pt.y() << ",";
    }
    std::cout << ")\n";
  #endif

  // Find the closest segments
  Polyline z_r_segments;
  std::size_t closest_index = 0;
#ifdef DO_BINARY_CLOSEST_SEGMENT
  closest_index = Vaango::Util::getClosestSegmentsBinarySearch(z_r_pt, z_r_table, 
                                                               z_r_segments);
#else
  if (z_r_table.size() < Vaango::Util::NUM_PTS_KDTREE_SWITCH) {
    closest_index = Vaango::Util::getClosestSegments(z_r_pt, z_r_table, 
                                                     z_r_segments);
  } else {
    closest_index = Vaango::Util::getClosestSegmentsKDTree(z_r_pt, z_r_table, 
                                                           z_r_segments);
  }
#endif

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
    std::cout << "ZRSeg0 = " << std::setprecision(16) << z_r_table[seg_start]
              << "\n";
    std::cout << "ZRSeg1 = " << std::setprecision(16) << z_r_table[seg_start + 1]
              << "\n";
    std::cout << "ZRSeg2 = " << std::setprecision(16) << z_r_table[seg_start + 2]
              << "\n";
    std::cout << "ZRPoint = " << std::setprecision(16) << z_r_pt << "\n";
    std::cout << "ZRClose = " << std::setprecision(16) << z_r_closest << "\n";
    std::cout << "ZRTangent = " << std::setprecision(16) << z_r_tangent << "\n";
  #endif

  return std::make_tuple(z_r_closest, z_r_tangent);
}

/* Convert yield function data to z_rprime coordinates */
void
YieldCond_TabularCap::convertToZRprime(const double& sqrtKG,
                                       const Polyline& p_q_points,
                                       Polyline& z_r_points) const
{
  // Compute z and r' for the yield surface points
  size_t ii = 0;
  for (const auto& pt : p_q_points) {
    double p_bar   = pt.x();
    double sqrt_J2 = pt.y();
    double z = 0.0, rprime = 0.0;
    Vaango::Util::convertToZRprime(sqrtKG, p_bar, sqrt_J2, z, rprime);
    z_r_points[ii++] = Point(z, rprime, 0);
  }
}

/* The closest point for the situation where the yield function polyline
   has already been converted to z-rprime coordinates in the ModelState */
std::tuple<Point, Vector>
YieldCond_TabularCap::getClosestPointSplineNewtonZR(
  const ModelState_TabularCap* state,
  const Polyline& z_r_table, 
  const Util::PolylineKDTree& z_r_index, 
  const Point& z_r_pt)
{
  // Find the closest segments
  Polyline z_r_segments;
  std::size_t closest_index = 0;
#ifdef DO_BINARY_CLOSEST_SEGMENT
  closest_index = Vaango::Util::getClosestSegmentsBinarySearch(z_r_pt, z_r_table, 
                                                               z_r_segments);
#else
  if (z_r_table.size() < Vaango::Util::NUM_PTS_KDTREE_SWITCH) {
    closest_index = Vaango::Util::getClosestSegments(z_r_pt, z_r_table, 
                                                     z_r_segments);
  } else {
    closest_index = Vaango::Util::getClosestSegmentsKDTree(z_r_pt, z_r_table, 
                                                           z_r_index, z_r_segments);
  }
#endif

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
  std::cout << "ZRSeg0 = " << std::setprecision(16) << z_r_table[seg_start]
            << "\n";
  std::cout << "ZRSeg1 = " << std::setprecision(16) << z_r_table[seg_start + 1]
            << "\n";
  std::cout << "ZRSeg2 = " << std::setprecision(16) << z_r_table[seg_start + 2]
            << "\n";
  std::cout << "ZRPoint = " << std::setprecision(16) << z_r_pt << "\n";
  std::cout << "ZRClose = " << std::setprecision(16) << z_r_closest << "\n";
  std::cout << "ZRTangent = " << std::setprecision(16) << z_r_tangent << "\n";
#endif

  return std::make_tuple(z_r_closest, z_r_tangent);
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
YieldCond_TabularCap::d2f_dp_depsVol(const ModelStateBase* state_input,
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
YieldCond_TabularCap::d2f_dp_depsDev(const ModelStateBase* state_input,
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
YieldCond_TabularCap::d2f_dq_depsVol(const ModelStateBase* state_input,
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
YieldCond_TabularCap::d2f_dq_depsDev(const ModelStateBase* state_input,
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
// Compute df/depse_s
//   df/depse_s =
//
// Requires:  Equation of state, shear modulus model
//--------------------------------------------------------------
double
YieldCond_TabularCap::df_depsDev(const ModelStateBase* state_input,
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

/*! Derivative with respect to the \f$xi\f$ where \f$\xi = s \f$
    where \f$s\f$ is deviatoric part of Cauchy stress */
Uintah::Matrix3
YieldCond_TabularCap::df_dxi(const Matrix3& stress,
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
YieldCond_TabularCap::df_dsigmaDev_dbeta(const Matrix3& stress,
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
YieldCond_TabularCap::computeElasPlasTangentModulus(
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
      << "called by models that use the Tabular yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return;
}

void
YieldCond_TabularCap::computeTangentModulus(const TangentModulusTensor& Ce,
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
operator<<(std::ostream& out, const YieldCond_TabularCap& yc)
{
  DoubleVec1D pvals, qvals;
  try {
    pvals =
      yc.d_yield.table.getIndependentVarData("Pressure", IndexKey(0, 0, 0, 0));
  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << std::endl;
  }

  try {
    qvals =
      yc.d_yield.table.getDependentVarData("SqrtJ2", IndexKey(0, 0, 0, 0));
  } catch (Uintah::InvalidValue e) {
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
