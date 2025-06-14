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

#include <CCA/Components/MPM/ConstitutiveModel/Utilities/YieldCondUtils.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCond_ArenaMixture.h>
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

const double YieldCond_ArenaMixture::sqrt_two       = std::sqrt(2.0);
const double YieldCond_ArenaMixture::sqrt_three     = std::sqrt(3.0);
const double YieldCond_ArenaMixture::one_sqrt_three = 1.0 / sqrt_three;

YieldCond_ArenaMixture::YieldCond_ArenaMixture(Uintah::ProblemSpecP& ps)
{
  ps->require("vol_frac.phase1", d_volfrac[0]); // Volume fractions
  d_volfrac[1] = 1.0 - d_volfrac[0];

  // Nonlinear Drucker-Prager parameters
  ps->require("PEAKI1.phase1",
              d_yieldParam[0].PEAKI1); // Shear Limit Surface Parameter
  ps->require("FSLOPE.phase1",
              d_yieldParam[0].FSLOPE); // Shear Limit Surface Parameter
  ps->require("STREN.phase1",
              d_yieldParam[0].STREN); // Shear Limit Surface Parameter
  ps->require("YSLOPE.phase1",
              d_yieldParam[0].YSLOPE); // Shear Limit Surface Parameter
  ps->getWithDefault(
    "PEAKI1_failed.phase1", d_yieldParam[0].PEAKI1_failed, 1.0e-5);
  ps->getWithDefault("FSLOPE_failed.phase1",
                     d_yieldParam[0].FSLOPE_failed,
                     0.5 * d_yieldParam[0].FSLOPE);
  ps->getWithDefault("STREN_failed.phase1",
                     d_yieldParam[0].STREN_failed,
                     0.1 * d_yieldParam[0].STREN);
  ps->getWithDefault(
    "YSLOPE_failed.phase1", d_yieldParam[0].YSLOPE_failed, 1.0e-5);

  ps->require("PEAKI1.phase2",
              d_yieldParam[1].PEAKI1); // Shear Limit Surface Parameter
  ps->require("FSLOPE.phase2",
              d_yieldParam[1].FSLOPE); // Shear Limit Surface Parameter
  ps->require("STREN.phase2",
              d_yieldParam[1].STREN); // Shear Limit Surface Parameter
  ps->require("YSLOPE.phase2",
              d_yieldParam[1].YSLOPE); // Shear Limit Surface Parameter
  ps->getWithDefault(
    "PEAKI1_failed.phase2", d_yieldParam[1].PEAKI1_failed, 1.0e-5);
  ps->getWithDefault("FSLOPE_failed.phase2",
                     d_yieldParam[1].FSLOPE_failed,
                     0.5 * d_yieldParam[1].FSLOPE);
  ps->getWithDefault("STREN_failed.phase2",
                     d_yieldParam[1].STREN_failed,
                     0.1 * d_yieldParam[1].STREN);
  ps->getWithDefault(
    "YSLOPE_failed.phase2", d_yieldParam[1].YSLOPE_failed, 1.0e-5);

  // Non-associativity parameters
  ps->require("BETA.phase1",
              d_nonAssocParam[0].BETA); // Nonassociativity Parameter
  ps->require("BETA.phase2",
              d_nonAssocParam[1].BETA); // Nonassociativity Parameter

  // Cap parameters: CR = (peakI1-kappa)/(peakI1-X)
  ps->require("CR.phase1", d_capParam[0].CR); // Cap Shape Parameter
  ps->require("CR.phase2", d_capParam[1].CR); // Cap Shape Parameter

  // Duvaut-Lion rate parameters
  ps->require("T1.phase1", d_rateParam[0].T1); // Rate dependence parameter
  ps->require("T2.phase1", d_rateParam[0].T2); // Rate dependence parameter
  ps->require("T1.phase2", d_rateParam[1].T1); // Rate dependence parameter
  ps->require("T2.phase2", d_rateParam[1].T2); // Rate dependence parameter

  // Check the input parameters
  checkInputParameters();

  // Compute the model parameters from the input parameters
  computeModelParameters(1.0);

  // Now optionally get the variablity information for each parameter
  // Mean values
  double PEAKI1 = d_volfrac[0] * d_yieldParam[0].PEAKI1 +
                  d_volfrac[1] * d_yieldParam[1].PEAKI1;
  double FSLOPE = d_volfrac[0] * d_yieldParam[0].FSLOPE +
                  d_volfrac[1] * d_yieldParam[1].FSLOPE;
  double STREN =
    d_volfrac[0] * d_yieldParam[0].STREN + d_volfrac[1] * d_yieldParam[1].STREN;
  double YSLOPE = d_volfrac[0] * d_yieldParam[0].YSLOPE +
                  d_volfrac[1] * d_yieldParam[1].YSLOPE;
  double BETA = d_volfrac[0] * d_nonAssocParam[0].BETA +
                d_volfrac[1] * d_nonAssocParam[1].BETA;
  double CR = d_volfrac[0] * d_capParam[0].CR + d_volfrac[1] * d_capParam[1].CR;
  double T1 =
    d_volfrac[0] * d_rateParam[0].T1 + d_volfrac[1] * d_rateParam[1].T1;
  double T2 =
    d_volfrac[0] * d_rateParam[0].T2 + d_volfrac[1] * d_rateParam[1].T2;
  std::string weibullDist;
  ps->getWithDefault("weibullDist_PEAKI1", weibullDist, std::to_string(PEAKI1));
  d_weibull_PEAKI1.WeibullParser(weibullDist);
  proc0cout << d_weibull_PEAKI1 << std::endl;

  ps->getWithDefault("weibullDist_FSLOPE", weibullDist, std::to_string(FSLOPE));
  d_weibull_FSLOPE.WeibullParser(weibullDist);
  proc0cout << d_weibull_FSLOPE << std::endl;

  ps->getWithDefault("weibullDist_STREN", weibullDist, std::to_string(STREN));
  d_weibull_STREN.WeibullParser(weibullDist);
  proc0cout << d_weibull_STREN << std::endl;

  ps->getWithDefault("weibullDist_YSLOPE", weibullDist, std::to_string(YSLOPE));
  d_weibull_YSLOPE.WeibullParser(weibullDist);
  proc0cout << d_weibull_YSLOPE << std::endl;

  ps->getWithDefault("weibullDist_BETA", weibullDist, std::to_string(BETA));
  d_weibull_BETA.WeibullParser(weibullDist);
  proc0cout << d_weibull_BETA << std::endl;

  ps->getWithDefault("weibullDist_CR", weibullDist, std::to_string(CR));
  d_weibull_CR.WeibullParser(weibullDist);
  proc0cout << d_weibull_CR << std::endl;

  ps->getWithDefault("weibullDist_T1", weibullDist, std::to_string(T1));
  d_weibull_T1.WeibullParser(weibullDist);
  proc0cout << d_weibull_T1 << std::endl;

  ps->getWithDefault("weibullDist_T2", weibullDist, std::to_string(T2));
  d_weibull_T2.WeibullParser(weibullDist);
  proc0cout << d_weibull_T2 << std::endl;

  // Initialize local labels for parameter variability
  initializeLocalMPMLabels();
}

YieldCond_ArenaMixture::YieldCond_ArenaMixture(const YieldCond_ArenaMixture* yc)
{
  for (int ii = 0; ii < 2; ii++) {
    d_volfrac[ii] = yc->d_volfrac[ii];

    d_modelParam[ii]    = yc->d_modelParam[ii];
    d_yieldParam[ii]    = yc->d_yieldParam[ii];
    d_nonAssocParam[ii] = yc->d_nonAssocParam[ii];
    d_capParam[ii]      = yc->d_capParam[ii];
    d_rateParam[ii]     = yc->d_rateParam[ii];
  }

  // Copy parameter variability information
  d_weibull_PEAKI1 = yc->d_weibull_PEAKI1;
  d_weibull_FSLOPE = yc->d_weibull_FSLOPE;
  d_weibull_STREN  = yc->d_weibull_STREN;
  d_weibull_YSLOPE = yc->d_weibull_YSLOPE;
  d_weibull_BETA   = yc->d_weibull_BETA;
  d_weibull_CR     = yc->d_weibull_CR;
  d_weibull_T1     = yc->d_weibull_T1;
  d_weibull_T2     = yc->d_weibull_T2;

  // Initialize local labels for parameter variability
  initializeLocalMPMLabels();
}

YieldCond_ArenaMixture::~YieldCond_ArenaMixture()
{
  VarLabel::destroy(pPEAKI1Label);
  VarLabel::destroy(pPEAKI1Label_preReloc);
  VarLabel::destroy(pFSLOPELabel);
  VarLabel::destroy(pFSLOPELabel_preReloc);
  VarLabel::destroy(pSTRENLabel);
  VarLabel::destroy(pSTRENLabel_preReloc);
  VarLabel::destroy(pYSLOPELabel);
  VarLabel::destroy(pYSLOPELabel_preReloc);

  VarLabel::destroy(pBETALabel);
  VarLabel::destroy(pBETALabel_preReloc);

  VarLabel::destroy(pCRLabel);
  VarLabel::destroy(pCRLabel_preReloc);

  VarLabel::destroy(pT1Label);
  VarLabel::destroy(pT1Label_preReloc);
  VarLabel::destroy(pT2Label);
  VarLabel::destroy(pT2Label_preReloc);
}

void
YieldCond_ArenaMixture::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP yield_ps = ps->appendChild("yield_condition");
  yield_ps->setAttribute("type", "arena_mixture");

  yield_ps->appendElement("vol_frac.phase1", d_volfrac[0]);

  yield_ps->appendElement("FSLOPE.phase1", d_yieldParam[0].FSLOPE);
  yield_ps->appendElement("PEAKI1.phase1", d_yieldParam[0].PEAKI1);
  yield_ps->appendElement("STREN.phase1", d_yieldParam[0].STREN);
  yield_ps->appendElement("YSLOPE.phase1", d_yieldParam[0].YSLOPE);
  yield_ps->appendElement("FSLOPE_failed.phase1",
                          d_yieldParam[0].FSLOPE_failed);
  yield_ps->appendElement("PEAKI1_failed.phase1",
                          d_yieldParam[0].PEAKI1_failed);
  yield_ps->appendElement("STREN_failed.phase1", d_yieldParam[0].STREN_failed);
  yield_ps->appendElement("YSLOPE_failed.phase1",
                          d_yieldParam[0].YSLOPE_failed);

  yield_ps->appendElement("BETA.phase1", d_nonAssocParam[0].BETA);

  yield_ps->appendElement("CR.phase1", d_capParam[0].CR);

  yield_ps->appendElement("T1.phase1", d_rateParam[0].T1);
  yield_ps->appendElement("T2.phase1", d_rateParam[0].T2);

  yield_ps->appendElement("FSLOPE.phase2", d_yieldParam[1].FSLOPE);
  yield_ps->appendElement("PEAKI1.phase2", d_yieldParam[1].PEAKI1);
  yield_ps->appendElement("STREN.phase2", d_yieldParam[1].STREN);
  yield_ps->appendElement("YSLOPE.phase2", d_yieldParam[1].YSLOPE);
  yield_ps->appendElement("FSLOPE_failed.phase2",
                          d_yieldParam[1].FSLOPE_failed);
  yield_ps->appendElement("PEAKI1_failed.phase2",
                          d_yieldParam[1].PEAKI1_failed);
  yield_ps->appendElement("STREN_failed.phase2", d_yieldParam[1].STREN_failed);
  yield_ps->appendElement("YSLOPE_failed.phase2",
                          d_yieldParam[1].YSLOPE_failed);

  yield_ps->appendElement("BETA.phase2", d_nonAssocParam[1].BETA);

  yield_ps->appendElement("CR.phase2", d_capParam[1].CR);

  yield_ps->appendElement("T1.phase2", d_rateParam[1].T1);
  yield_ps->appendElement("T2.phase2", d_rateParam[1].T2);

  yield_ps->appendElement("weibullDist_PEAKI1", d_weibull_PEAKI1.getWeibDist());
  yield_ps->appendElement("weibullDist_FSLOPE", d_weibull_FSLOPE.getWeibDist());
  yield_ps->appendElement("weibullDist_STREN", d_weibull_STREN.getWeibDist());
  yield_ps->appendElement("weibullDist_YSLOPE", d_weibull_YSLOPE.getWeibDist());

  yield_ps->appendElement("weibullDist_BETA", d_weibull_BETA.getWeibDist());

  yield_ps->appendElement("weibullDist_CR", d_weibull_CR.getWeibDist());

  yield_ps->appendElement("weibullDist_T1", d_weibull_T1.getWeibDist());
  yield_ps->appendElement("weibullDist_T2", d_weibull_T2.getWeibDist());
}

//--------------------------------------------------------------
// Check that the input parameters are reasonable
//--------------------------------------------------------------
void
YieldCond_ArenaMixture::checkInputParameters()
{
  std::ostringstream warn;
  if (d_yieldParam[0].PEAKI1 < 0.0 || d_yieldParam[0].PEAKI1_failed < 0.0) {
    warn << "Phase 1: PEAKI1 must be nonnegative. PEAKI1 = "
         << d_yieldParam[0].PEAKI1 << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_yieldParam[0].FSLOPE < 0.0 || d_yieldParam[0].FSLOPE_failed < 0.0) {
    warn << "Phase 1: FSLOPE must be nonnegative. FSLOPE = "
         << d_yieldParam[0].FSLOPE << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_yieldParam[0].FSLOPE < d_yieldParam[0].YSLOPE ||
      d_yieldParam[0].FSLOPE_failed < d_yieldParam[0].YSLOPE_failed) {
    warn << "Phase 1: FSLOPE must be greater than YSLOPE. FSLOPE = "
         << d_yieldParam[0].FSLOPE << ", YSLOPE = " << d_yieldParam[0].YSLOPE
         << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_nonAssocParam[0].BETA <= 0.0) {
    warn << "Phase 1: BETA (nonassociativity factor) must be positive. BETA = "
         << d_nonAssocParam[0].BETA << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_capParam[0].CR >= 1 || d_capParam[0].CR <= 0.0) {
    warn << "Phase 1: CR must be 0<CR<1. CR = " << d_capParam[0].CR
         << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_rateParam[0].T1 < 0.0) {
    warn << "Phase 1: T1 must be nonnegative. T1 = " << d_rateParam[0].T1
         << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_rateParam[0].T2 < 0.0) {
    warn << "Phase 1: T2 must be nonnegative. T2 = " << d_rateParam[0].T2
         << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if ((d_rateParam[0].T1 > 0.0 || d_rateParam[0].T2 > 0.0) !=
      (d_rateParam[0].T1 > 0.0 && d_rateParam[0].T2 > 0.0)) {
    warn
      << "Phase 1: For rate dependence both T1 and T2 must be positive. T1 = "
      << d_rateParam[0].T1 << ", T2 = " << d_rateParam[0].T2 << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  if (d_yieldParam[1].PEAKI1 < 0.0 || d_yieldParam[1].PEAKI1_failed < 0.0) {
    warn << "Phase 2: PEAKI1 must be nonnegative. PEAKI1 = "
         << d_yieldParam[1].PEAKI1 << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_yieldParam[1].FSLOPE < 0.0 || d_yieldParam[1].FSLOPE_failed < 0.0) {
    warn << "Phase 2: FSLOPE must be nonnegative. FSLOPE = "
         << d_yieldParam[1].FSLOPE << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_yieldParam[1].FSLOPE < d_yieldParam[1].YSLOPE ||
      d_yieldParam[1].FSLOPE_failed < d_yieldParam[1].YSLOPE_failed) {
    warn << "Phase 2: FSLOPE must be greater than YSLOPE. FSLOPE = "
         << d_yieldParam[1].FSLOPE << ", YSLOPE = " << d_yieldParam[1].YSLOPE
         << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_nonAssocParam[1].BETA <= 0.0) {
    warn << "Phase 2: BETA (nonassociativity factor) must be positive. BETA = "
         << d_nonAssocParam[1].BETA << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_capParam[1].CR >= 1 || d_capParam[1].CR <= 0.0) {
    warn << "Phase 2: CR must be 0<CR<1. CR = " << d_capParam[1].CR
         << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_rateParam[1].T1 < 0.0) {
    warn << "Phase 2: T1 must be nonnegative. T1 = " << d_rateParam[1].T1
         << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_rateParam[1].T2 < 0.0) {
    warn << "Phase 2: T2 must be nonnegative. T2 = " << d_rateParam[1].T2
         << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if ((d_rateParam[1].T1 > 0.0 || d_rateParam[1].T2 > 0.0) !=
      (d_rateParam[1].T1 > 0.0 && d_rateParam[1].T2 > 0.0)) {
    warn
      << "Phase 2: For rate dependence both T1 and T2 must be positive. T1 = "
      << d_rateParam[1].T1 << ", T2 = " << d_rateParam[1].T2 << std::endl;
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
YieldCond_ArenaMixture::computeModelParameters(double)
{
  computeModelParameters(0);
  computeModelParameters(1);
}

void
YieldCond_ArenaMixture::computeModelParameters(int phase)
{
  double FSLOPE = d_yieldParam[phase].FSLOPE, // Slope at I1=PEAKI1
    STREN       = d_yieldParam[phase].STREN,  // Value of rootJ2 at I1=0
    YSLOPE      = d_yieldParam[phase].YSLOPE, // High pressure slope
    PEAKI1      = d_yieldParam[phase].PEAKI1; // Value of I1 at strength=0
  double FSLOPE_failed =
           d_yieldParam[phase].FSLOPE_failed,         // Slope at I1=PEAKI1
    STREN_failed  = d_yieldParam[phase].STREN_failed, // Value of rootJ2 at I1=0
    YSLOPE_failed = d_yieldParam[phase].YSLOPE_failed, // High pressure slope
    PEAKI1_failed =
      d_yieldParam[phase].PEAKI1_failed; // Value of I1 at strength=0

  std::vector<double> limitParameters =
    computeModelParameters(PEAKI1, FSLOPE, STREN, YSLOPE);
  std::vector<double> limitParameters_failed = computeModelParameters(
    PEAKI1_failed, FSLOPE_failed, STREN_failed, YSLOPE_failed);

  d_modelParam[phase].a1        = limitParameters[0];
  d_modelParam[phase].a2        = limitParameters[1];
  d_modelParam[phase].a3        = limitParameters[2];
  d_modelParam[phase].a4        = limitParameters[3];
  d_modelParam[phase].a1_failed = limitParameters_failed[0];
  d_modelParam[phase].a2_failed = limitParameters_failed[1];
  d_modelParam[phase].a3_failed = limitParameters_failed[2];
  d_modelParam[phase].a4_failed = limitParameters_failed[3];
}

std::vector<double>
YieldCond_ArenaMixture::computeModelParameters(const double& PEAKI1,
                                               const double& FSLOPE,
                                               const double& STREN,
                                               const double& YSLOPE) const
{
  double a1, a2, a3, a4;
  if (FSLOPE > 0.0 && PEAKI1 >= 0.0 && STREN == 0.0 &&
      YSLOPE == 0.0) { // ----------------------------------------------Linear
                       // Drucker Prager
    a1 = PEAKI1 * FSLOPE;
    a2 = 0.0;
    a3 = 0.0;
    a4 = FSLOPE;
  } else if (FSLOPE == 0.0 && PEAKI1 == 0.0 && STREN > 0.0 &&
             YSLOPE ==
               0.0) { // -------------------------------------------------------
                      // Von Mises
    a1 = STREN;
    a2 = 0.0;
    a3 = 0.0;
    a4 = 0.0;
  } else if (FSLOPE > 0.0 && YSLOPE == 0.0 && STREN > 0.0 &&
             PEAKI1 ==
               0.0) { // -------------------------------------------------------
                      // 0 PEAKI1 to vonMises
    a1 = STREN;
    a2 = FSLOPE / STREN;
    a3 = STREN;
    a4 = 0.0;
  } else if (FSLOPE > YSLOPE && YSLOPE > 0.0 && STREN > YSLOPE * PEAKI1 &&
             PEAKI1 >=
               0.0) { // -------------------------------------------------------
                      // Nonlinear Drucker-Prager
    a1 = STREN;
    a2 = (FSLOPE - YSLOPE) / (STREN - YSLOPE * PEAKI1);
    a3 = (STREN - YSLOPE * PEAKI1) * exp(-a2 * PEAKI1);
    a4 = YSLOPE;
  } else {
    // Bad inputs, throw exception:
    std::ostringstream warn;
    warn << "Bad input parameters for shear limit surface. "
         << "FSLOPE = " << FSLOPE << ", YSLOPE = " << YSLOPE
         << ", PEAKI1 = " << PEAKI1 << ", STREN = " << STREN << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  std::vector<double> limitParameters = { a1, a2, a3, a4 };
  return limitParameters;
}

//--------------------------------------------------------------
// Evaluate yield condition
//
// f := J2 - Ff^2*Fc^2 = 0
// where
//     J2 = 1/2 s:s,  s = sigma - p I,  p = 1/3 Tr(sigma)
//     I1_eff = 3*(p + pbar_w)
//     X_eff  = X + 3*pbar_w
//     kappa = I1_peak - CR*(I1_peak - X_eff)
//     Ff := a1 - a3*exp(a2*I1_eff) - a4*I1_eff
//     Fc^2 := 1 - (kappa - I1_eff)^2/(kappa - X_eff)^2
//
// Returns:
//   hasYielded = -1.0 (if elastic)
//              =  1.0 (otherwise)
//--------------------------------------------------------------
std::pair<double, Util::YieldStatus>
YieldCond_ArenaMixture::evalYieldCondition(const ModelStateBase* state_input)
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

  // Get the particle specific internal variables from the model state
  double PEAKI1 = state->yieldParams.at("PEAKI1");
  double FSLOPE = state->yieldParams.at("FSLOPE");
  double STREN  = state->yieldParams.at("STREN");
  double YSLOPE = state->yieldParams.at("YSLOPE");
  double CR     = state->yieldParams.at("CR");

  std::vector<double> limitParameters =
    computeModelParameters(PEAKI1, FSLOPE, STREN, YSLOPE);
  double a1 = limitParameters[0];
  double a2 = limitParameters[1];
  double a3 = limitParameters[2];
  double a4 = limitParameters[3];

  // Get the local vars from the model state
  double X_eff = state->capX + 3.0 * state->pbar_w;

  // Cauchy stress invariants: I1_eff = 3*(p + pbar_w), J2 = q^2/3
  double I1_eff  = state->I1_eff;
  double sqrt_J2 = state->sqrt_J2;

  // --------------------------------------------------------------------
  // *** SHEAR LIMIT FUNCTION (Ff) ***
  // --------------------------------------------------------------------
  double Ff = a1 - a3 * exp(a2 * I1_eff) - a4 * I1_eff;

  // --------------------------------------------------------------------
  // *** Branch Point (Kappa) ***
  // --------------------------------------------------------------------
  double kappa = PEAKI1 - CR * (PEAKI1 - X_eff); // Branch Point

  // --------------------------------------------------------------------
  // *** COMPOSITE YIELD FUNCTION ***
  // --------------------------------------------------------------------
  // Evaluate Composite Yield Function F(I1) = Ff(I1)*fc(I1) in each region.
  // The elseif statements have nested if statements, which is not equivalent
  // to them having a single elseif(A&&B&&C)
  //---------------------------------------------------(I1<X)
  if (I1_eff < X_eff) {
    // std::cout << " I1_eff < X_eff " << I1_eff << "," << X_eff << std::endl;
    return std::make_pair(1.0, Util::YieldStatus::HAS_YIELDED);
  }

  // **Elliptical Cap Function: (fc)**
  // fc = sqrt(1.0 - Pow((Kappa-I1mZ)/(Kappa-X)),2.0);
  // faster version: fc2 = fc^2
  // **WARNING** p3 is the maximum achievable volumetric plastic strain in
  // compresson
  // so if a value of 0 has been specified this indicates the user
  // wishes to run without porosity, and no cap function is used, i.e. fc=1
  if ((X_eff < I1_eff) && (I1_eff < kappa)) { // ---------------(X<I1<kappa)

    double kappaRatio = (kappa - I1_eff) / (kappa - X_eff);
    double fc2        = 1.0 - kappaRatio * kappaRatio;
    if (sqrt_J2 * sqrt_J2 > Ff * Ff * fc2) {
      // std::cout << " X_eff < I1_eff " << I1_eff << "," << X_eff << std::endl;
      // std::cout << " I1_eff < kappa " << I1_eff << "," << kappa << std::endl;
      // std::cout << " J2 < Ff^2*Fc^2 " << sqrt_J2 << "," << Ff << ", " << fc2
      // << std::endl;
      return std::make_pair(1.0, Util::YieldStatus::HAS_YIELDED);
    }
  } else { // --------- X >= I1 or kappa <= I1

    if (I1_eff <= PEAKI1) { // ----- (kappa <= I1 <= PEAKI1)
      if (sqrt_J2 > Ff) {
        // std::cout << " I1_eff < PEAKI1 " << I1_eff << "," << PEAKI1 <<
        // std::endl;
        // std::cout << " sqrt(J2) > Ff " << sqrt_J2 << "," << Ff << std::endl;
        return std::make_pair(1.0, Util::YieldStatus::HAS_YIELDED);
      }
    } else { // I1 > PEAKI1
      // std::cout << " I1_eff > PEAKI1 " << I1_eff << "," << PEAKI1 <<
      // std::endl;
      return std::make_pair(1.0, Util::YieldStatus::HAS_YIELDED);
    }
  }

  return std::make_pair(-1.0, Util::YieldStatus::IS_ELASTIC);
}

double
YieldCond_ArenaMixture::computeYieldFunction(const ModelStateBase* state_input) const
{
  const ModelState_Arena* state =
    static_cast<const ModelState_Arena*>(state_input);

  // Get the particle specific internal variables from the model state
  double PEAKI1 = state->yieldParams.at("PEAKI1");
  double FSLOPE = state->yieldParams.at("FSLOPE");
  double STREN  = state->yieldParams.at("STREN");
  double YSLOPE = state->yieldParams.at("YSLOPE");
  double CR     = state->yieldParams.at("CR");

  std::vector<double> limitParameters =
    computeModelParameters(PEAKI1, FSLOPE, STREN, YSLOPE);
  double a1 = limitParameters[0];
  double a2 = limitParameters[1];
  double a3 = limitParameters[2];
  double a4 = limitParameters[3];

  // Get the local vars from the model state
  double X_eff = state->capX + 3.0 * state->pbar_w;

  // Cauchy stress invariants: I1_eff = 3*(p + pbar_w), J2 = q^2/3
  double I1_eff  = state->I1_eff;
  double sqrt_J2 = state->sqrt_J2;

  // --------------------------------------------------------------------
  // *** SHEAR LIMIT FUNCTION (Ff) ***
  // --------------------------------------------------------------------
  double Ff = a1 - a3 * exp(a2 * I1_eff) - a4 * I1_eff;

  // --------------------------------------------------------------------
  // *** Branch Point (Kappa) ***
  // --------------------------------------------------------------------
  double epsilon    = std::numeric_limits<double>::epsilon();
  double kappa      = PEAKI1 - CR * (PEAKI1 - X_eff); // Branch Point
  double kappaRatio = (kappa - I1_eff) / (kappa - X_eff + epsilon);
  double Fc_sq      = 1.0 - kappaRatio * kappaRatio;

  // --------------------------------------------------------------------
  // *** COMPOSITE YIELD FUNCTION ***
  // --------------------------------------------------------------------
  double f = sqrt_J2 * sqrt_J2 - Ff * Ff * Fc_sq;

  return f;
}

//--------------------------------------------------------------
// Evaluate yield condition max  value of sqrtJ2
//--------------------------------------------------------------
double
YieldCond_ArenaMixture::evalYieldConditionMax(const ModelStateBase* state_input)
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

  // Get the particle specific internal variables from the model state
  // Store in a local struct
  d_local.PEAKI1 = state->yieldParams.at("PEAKI1");
  d_local.FSLOPE = state->yieldParams.at("FSLOPE");
  d_local.STREN  = state->yieldParams.at("STREN");
  d_local.YSLOPE = state->yieldParams.at("YSLOPE");
  d_local.BETA   = state->yieldParams.at("BETA");
  d_local.CR     = state->yieldParams.at("CR");

  std::vector<double> limitParameters = computeModelParameters(
    d_local.PEAKI1, d_local.FSLOPE, d_local.STREN, d_local.YSLOPE);
  d_local.a1 = limitParameters[0];
  d_local.a2 = limitParameters[1];
  d_local.a3 = limitParameters[2];
  d_local.a4 = limitParameters[3];

  // Get the plastic internal variables from the model state
  double pbar_w = state->pbar_w;
  double X_eff  = state->capX + 3.0 * pbar_w;

  // Compute kappa
  double kappa = d_local.PEAKI1 - d_local.CR * (d_local.PEAKI1 - X_eff);

  // Number of points
  int num_points = 10;

  // Set up I1 values
  // double I1eff_min = 0.99999*X_eff;
  // double I1eff_max = 0.99999*d_local.PEAKI1;
  // std::vector<double> I1_eff_vec;
  // Vaango::Util::linspace(I1eff_min, I1eff_max, num_points, I1_eff_vec);
  double rad       = 0.5 * (d_local.PEAKI1 - X_eff);
  double cen       = 0.5 * (d_local.PEAKI1 + X_eff);
  double theta_min = 0.0;
  double theta_max = M_PI;
  std::vector<double> theta_vec;
  Vaango::Util::linspace(theta_min, theta_max, num_points, theta_vec);
  double J2_max = std::numeric_limits<double>::min();
  // for (auto I1_eff : I1_eff_vec) {
  for (auto theta : theta_vec) {

    double I1_eff = cen + rad * std::cos(theta);

    // Compute F_f
    double Ff = d_local.a1 - d_local.a3 * std::exp(d_local.a2 * I1_eff) -
                d_local.a4 * (I1_eff);
    double Ff_sq = Ff * Ff;

    // Compute Fc
    double Fc_sq = 1.0;
    if ((I1_eff < kappa) && (X_eff < I1_eff)) {
      double ratio = (kappa - I1_eff) / (kappa - X_eff);
      Fc_sq        = 1.0 - ratio * ratio;
    }

    // Compute J2
    J2_max = std::max(J2_max, Ff_sq * Fc_sq);
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
Matrix3
YieldCond_ArenaMixture::df_dsigma(const ModelStateBase* state) 
{
  return df_dsigma(Vaango::Util::Identity, state);
}

Matrix3
YieldCond_ArenaMixture::df_dsigma(const Matrix3&,
                                  const ModelStateBase* state_input)
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

  double dfdp  = df_dp(state_input);
  double dfdJ2 = df_dq(state_input);

  Matrix3 One;
  One.Identity();
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
YieldCond_ArenaMixture::df_dp(const ModelStateBase* state_input)
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

  // Get the particle specific internal variables from the model state
  double PEAKI1 = state->yieldParams.at("PEAKI1");
  double FSLOPE = state->yieldParams.at("FSLOPE");
  double STREN  = state->yieldParams.at("STREN");
  double YSLOPE = state->yieldParams.at("YSLOPE");
  double CR     = state->yieldParams.at("CR");

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
  double kappa_X_eff  = kappa - X_eff;
  double kappaRatio   = kappa_I1_eff / kappa_X_eff;
  double Fc_sq        = 1.0 - kappaRatio * kappaRatio;

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
YieldCond_ArenaMixture::df_dq(const ModelStateBase*)
{
  /*
  const ModelState_Arena* state =
    static_cast<const ModelState_Arena*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Arena.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

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
YieldCond_ArenaMixture::getInternalPoint(
  const ModelStateBase* state_old_input,
  const ModelStateBase* state_trial_input)
{
  const ModelState_Arena* state_old =
    static_cast<const ModelState_Arena*>(state_old_input);
  const ModelState_Arena* state_trial =
    static_cast<const ModelState_Arena*>(state_trial_input);
  /*
  if ((!state_old) || (!state_trial)) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Arena.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

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
  double upperI1         = PEAKI1;
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
YieldCond_ArenaMixture::getClosestPoint(const ModelStateBase* state_input,
                                        const double& px,
                                        const double& py,
                                        double& cpx,
                                        double& cpy)
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
YieldCond_ArenaMixture::getClosestPointGeometricBisect(
  const ModelState_Arena* state,
  const Uintah::Point& z_r_pt,
  Uintah::Point& z_r_closest)
{
  // Get the particle specific internal variables from the model state
  // Store in a local struct
  d_local.PEAKI1 = state->yieldParams.at("PEAKI1");
  d_local.FSLOPE = state->yieldParams.at("FSLOPE");
  d_local.STREN  = state->yieldParams.at("STREN");
  d_local.YSLOPE = state->yieldParams.at("YSLOPE");
  d_local.BETA   = state->yieldParams.at("BETA");
  d_local.CR     = state->yieldParams.at("CR");

  std::vector<double> limitParameters = computeModelParameters(
    d_local.PEAKI1, d_local.FSLOPE, d_local.STREN, d_local.YSLOPE);
  d_local.a1 = limitParameters[0];
  d_local.a2 = limitParameters[1];
  d_local.a3 = limitParameters[2];
  d_local.a4 = limitParameters[3];

  // Get the plastic internal variables from the model state
  double pbar_w = state->pbar_w;
  double X_eff  = state->capX + 3.0 * pbar_w;

  // Compute kappa
  double I1_diff = d_local.PEAKI1 - X_eff;
  double kappa   = d_local.PEAKI1 - d_local.CR * I1_diff;

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
  double eta_mid   = 0.5 * (eta_lo + eta_hi);

  // Do bisection
  int iters        = 1;
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
    getYieldSurfacePointsAll_RprimeZ(
      X_eff, kappa, sqrtKG, I1eff_min, I1eff_max, num_points, z_r_points);

    // Find the closest point
    Vaango::Util::findClosestPoint(z_r_pt, z_r_points, z_r_closest);

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
      eta_hi    = eta_mid;
    } else {
      I1eff_min = I1eff_mid;
      eta_lo    = eta_mid;
    }

    I1eff_mid = 0.5 * (I1eff_min + I1eff_max);
    eta_mid   = 0.5 * (eta_lo + eta_hi);

    // Distance to old closest point
    if (iters > 10 && (z_r_closest - z_r_closest_old).length2() < 1.0e-16) {
      break;
    }
    z_r_closest_old = z_r_closest;

    ++iters;
  }

  return;
}

void
YieldCond_ArenaMixture::getClosestPointAlgebraicBisect(
  const ModelState_Arena* state,
  const Uintah::Point& z_r_pt,
  Uintah::Point& z_r_closest)
{
  // Get the particle specific internal variables from the model state
  // Store in a local struct
  d_local.PEAKI1 = state->yieldParams.at("PEAKI1");
  d_local.FSLOPE = state->yieldParams.at("FSLOPE");
  d_local.STREN  = state->yieldParams.at("STREN");
  d_local.YSLOPE = state->yieldParams.at("YSLOPE");
  d_local.BETA   = state->yieldParams.at("BETA");
  d_local.CR     = state->yieldParams.at("CR");

  std::vector<double> limitParameters = computeModelParameters(
    d_local.PEAKI1, d_local.FSLOPE, d_local.STREN, d_local.YSLOPE);
  d_local.a1 = limitParameters[0];
  d_local.a2 = limitParameters[1];
  d_local.a3 = limitParameters[2];
  d_local.a4 = limitParameters[3];

  // Get the plastic internal variables from the model state
  double pbar_w = state->pbar_w;
  double X_eff  = state->capX + 3.0 * pbar_w;

  // Compute kappa
  double kappa = d_local.PEAKI1 - d_local.CR * (d_local.PEAKI1 - X_eff);

  // Compute factor
  double beta_KG_fac =
    d_local.BETA * std::sqrt(3.0 * state->bulkModulus / state->shearModulus);

  // Set up I1 and z_eff limits
  double I1eff_min = X_eff;
  double I1eff_max = d_local.PEAKI1;

  // Get the trial point
  double zeff_trial   = z_r_pt.x();
  double rprime_trial = z_r_pt.y();

  // Set up lambda to calculate g(z_eff)
  auto gfun = [this, &kappa, &X_eff, &beta_KG_fac, &rprime_trial, &zeff_trial](double I1eff) {

    // Compute F_f
    double a3_exp_a2_I1 = d_local.a3 * std::exp(d_local.a2 * I1eff);
    double Ff           = d_local.a1 - a3_exp_a2_I1 - d_local.a4 * (I1eff);

    // Compute dFf_dzeff
    double dFf_dzeff = -sqrt_three * (d_local.a2 * a3_exp_a2_I1 + d_local.a4);

    // Compute Fc and dFc_dzeff
    double Fc        = 1.0;
    double dFc_dzeff = 0.0;
    if ((I1eff < kappa) && (X_eff < I1eff)) {
      double ratio = (kappa - I1eff) / (kappa - X_eff);
      // TODO: Add check for negative values of 1 - ratio^2
      Fc        = std::sqrt(1.0 - ratio * ratio);
      dFc_dzeff = sqrt_three * ratio / (Fc * (kappa - X_eff));
    }

    // Compute g(zeff)
    double zeff = I1eff * one_sqrt_three;
    double gval = (zeff_trial - zeff) +
                  beta_KG_fac * (rprime_trial - beta_KG_fac * Ff * Fc) *
                    (Fc * dFf_dzeff + Ff * dFc_dzeff);

    // Compute r'
    double rprime = beta_KG_fac * Ff * Fc;

    return std::vector<double>{ gval, zeff, rprime };

  };

  // First check the end points
  std::vector<double> gfun_min = gfun(I1eff_min);
  std::vector<double> gfun_max = gfun(I1eff_max);
  double gmin                  = gfun_min[0];
  double gmax                  = gfun_max[0];

  if (std::signbit(gmax) == std::signbit(gmin)) {
    std::cout << "gmin = " << gmin << " gmax = " << gmax << std::endl;
    std::cout << " g(z_min) and g(z_max) have the same sign."
              << " Doing geometric bisection." << std::endl;
    getClosestPointGeometricBisect(state, z_r_pt, z_r_closest);
    return;
  }

  // The first point is the closest point
  if (std::abs(gmin) < std::numeric_limits<double>::min()) {
    z_r_closest.x(gfun_min[1]);
    z_r_closest.y(gfun_min[2]);
    return;
  }

  // The last point is the closest point
  if (std::abs(gmax) < std::numeric_limits<double>::min()) {
    z_r_closest.x(gfun_max[1]);
    z_r_closest.y(gfun_max[2]);
    return;
  }

  // Set up bisection
  double TOLERANCE =
    std::min(1.0e-10, 1.0e-16 * std::abs(I1eff_max - I1eff_min));
  int MAX_ITER = (int)std::ceil(std::log2((I1eff_max - I1eff_min) / TOLERANCE));
  int iter     = 0;
  bool isSuccess = false;

  double I1eff_mid = 0.0, zeff_mid = 0.0, rprime_mid = 0.0;
  while (iter < MAX_ITER) {

    iter++;

    I1eff_mid = 0.5 * (I1eff_min + I1eff_max);

    std::vector<double> gfun_mid = gfun(I1eff_mid);
    double gmid                  = gfun_mid[0];
    zeff_mid                     = gfun_mid[1];
    rprime_mid                   = gfun_mid[2];

    // Check g(zeff = 0) or (zeff_max - zeff_min)/2 < TOLERANCE
    if ((std::abs(gmid) < std::numeric_limits<double>::min()) ||
        (0.5 * (I1eff_max - I1eff_min) < TOLERANCE)) {
      isSuccess = true;
      break;
    }

    std::vector<double> gfun_min = gfun(I1eff_min);
    double gmin                  = gfun_min[0];

    if (std::signbit(gmid) == std::signbit(gmin)) {
      I1eff_min = I1eff_mid;
    } else {
      I1eff_max = I1eff_mid;
    }
  }

  if (isSuccess) {
    z_r_closest.x(zeff_mid);
    z_r_closest.y(rprime_mid);
  } else {
    getClosestPointGeometricBisect(state, z_r_pt, z_r_closest);
  }

#ifdef DEBUG_YIELD_BISECTION
  // Compute g for several values of I
  int num_points = 20;
  double rad     = 0.5 * (d_local.PEAKI1 - X_eff);
  double cen     = 0.5 * (d_local.PEAKI1 + X_eff);
  std::vector<double> theta_vec;
  Vaango::Util::linspace(0.0, M_PI, num_points, theta_vec);

  std::vector<double> gvec;
  std::vector<double> I1vec;
  for (auto theta : theta_vec) {
    double I1_eff = std::max(cen + rad * std::cos(theta), X_eff);
    I1vec.push_back(I1_eff);
    gvec.push_back(gfun(I1_eff)[0]);
  }
  std::cout << "I1vec = [";
  for (auto& I1 : I1vec) {
    std::cout << I1 << " ";
  }
  std::cout << "];" << std::endl;
  std::cout << "gvec = [";
  for (auto& g : gvec) {
    std::cout << g << " ";
  }
  std::cout << "];" << std::endl;
  std::cout << "plot(I1vec, gvec);" << std::endl;

  // Get the yield surface points
  std::vector<Uintah::Point> z_r_points;
  double sqrtKG = std::sqrt(1.5 * state->bulkModulus / state->shearModulus);
  I1eff_min     = 0.999999 * X_eff;
  I1eff_max     = 0.999999 * d_local.PEAKI1;
  getYieldSurfacePointsAll_RprimeZ(
    X_eff, kappa, sqrtKG, I1eff_min, I1eff_max, num_points, z_r_points);
  // Compute distances
  std::vector<double> distSq;
  for (auto& pt : z_r_points) {
    distSq.push_back((z_r_pt - pt).length2());
  }

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
  std::cout << "plot(z_r_pt(1), z_r_pt(2), 'ko');" << std::endl;
  std::cout << "plot(z_r_closest(1), z_r_closest(2), 'gx');" << std::endl;
  std::cout
    << "plot([z_r_pt(1) z_r_closest(1)],[z_r_pt(2) z_r_closest(2)], 'r--');"
    << std::endl;
  std::cout << "z_r_distSq = [";
  for (auto& dist : distSq) {
    std::cout << dist << " ";
  }
  std::cout << "];" << std::endl;
  std::cout << "plot(z_r_yield_z, z_r_distSq); hold on;" << std::endl;

#endif

  return;
}

/* Get the points on the yield surface */
void
YieldCond_ArenaMixture::getYieldSurfacePointsAll_RprimeZ(
  const double& X_eff,
  const double& kappa,
  const double& sqrtKG,
  const double& I1eff_min,
  const double& I1eff_max,
  const int& num_points,
  std::vector<Uintah::Point>& z_r_vec)
{
  // Compute z_eff and r'
  computeZeff_and_RPrime(
    X_eff, kappa, sqrtKG, I1eff_min, I1eff_max, num_points, z_r_vec);

  return;
}

/* Get the points on two segments the yield surface */
void
YieldCond_ArenaMixture::getYieldSurfacePointsSegment_RprimeZ(
  const double& X_eff,
  const double& kappa,
  const double& sqrtKG,
  const Uintah::Point& start_point,
  const Uintah::Point& end_point,
  const int& num_points,
  std::vector<Uintah::Point>& z_r_poly)
{

  // Find the start I1 and end I1 values of the segments
  // **TODO** make sure that the start and end points are differenet
  double z_effStart  = start_point.x();
  double z_effEnd    = end_point.x();
  double I1_effStart = sqrt_three * z_effStart;
  double I1_effEnd   = sqrt_three * z_effEnd;

  // Compute z_eff and r'
  computeZeff_and_RPrime(
    X_eff, kappa, sqrtKG, I1_effStart, I1_effEnd, num_points, z_r_poly);

  return;
}

/*! Compute a vector of z_eff, r' values given a range of I1_eff values */
void
YieldCond_ArenaMixture::computeZeff_and_RPrime(
  const double& X_eff,
  const double& kappa,
  const double& sqrtKG,
  const double& I1eff_min,
  const double& I1eff_max,
  const int& num_points,
  std::vector<Uintah::Point>& z_r_vec)
{
  // Set up points
  double rad       = 0.5 * (d_local.PEAKI1 - X_eff);
  double cen       = 0.5 * (d_local.PEAKI1 + X_eff);
  double theta_max = std::acos(std::max((I1eff_min - cen) / rad, -1.0));
  double theta_min = std::acos(std::min((I1eff_max - cen) / rad, 1.0));
  std::vector<double> theta_vec;
  Vaango::Util::linspace(theta_min, theta_max, num_points, theta_vec);

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
      Fc_sq        = std::max(1.0 - ratio * ratio, 0.0);
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
YieldCond_ArenaMixture::d2f_dp_depsVol([[maybe_unused]] const ModelStateBase* state_input,
                                       [[maybe_unused]] const MPMEquationOfState* eos,
                                       const ShearModulusModel*)
{
  std::ostringstream out;
  out << "**ERROR** d2f_dp_depsVol should not be called by "
      << " models that use the SoilMix yield criterion.";
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
YieldCond_ArenaMixture::d2f_dp_depsDev([[maybe_unused]] const ModelStateBase* state_input,
                                       [[maybe_unused]] const MPMEquationOfState* eos,
                                       const ShearModulusModel*)
{
  std::ostringstream out;
  out << "**ERROR** d2f_dp_depsDev should not be called by "
      << " models that use the SoilMix yield criterion.";
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
YieldCond_ArenaMixture::d2f_dq_depsVol([[maybe_unused]] const ModelStateBase* state_input,
                                       const MPMEquationOfState*,
                                       [[maybe_unused]] const ShearModulusModel* shear)
{
  std::ostringstream out;
  out << "**ERROR** d2f_dq_depsVol should not be called by "
      << " models that use the SoilMix yield criterion.";
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
YieldCond_ArenaMixture::d2f_dq_depsDev([[maybe_unused]] const ModelStateBase* state_input,
                                       const MPMEquationOfState*,
                                       [[maybe_unused]] const ShearModulusModel* shear)
{
  std::ostringstream out;
  out << "**ERROR** d2f_dq_depsDev should not be called by "
      << " models that use the SoilMix yield criterion.";
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
YieldCond_ArenaMixture::df_depsVol([[maybe_unused]] const ModelStateBase* state_input,
                                   [[maybe_unused]] const MPMEquationOfState* eos,
                                   [[maybe_unused]] const ShearModulusModel* shear)
{
  std::ostringstream out;
  out << "**ERROR** df_depsVol should not be called by "
      << " models that use the SoilMix yield criterion.";
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
YieldCond_ArenaMixture::df_depsDev([[maybe_unused]] const ModelStateBase* state_input,
                                   [[maybe_unused]] const MPMEquationOfState* eos,
                                   [[maybe_unused]] const ShearModulusModel* shear)
{
  std::ostringstream out;
  out << "**ERROR** df_depsVol should not be called by "
      << " models that use the SoilMix yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return 0.0;
}

// Evaluate yield condition (s = deviatoric stress
//                           p = state->p)
double
YieldCond_ArenaMixture::evalYieldCondition(const Uintah::Matrix3&,
                                           [[maybe_unused]] const ModelStateBase* state_input)
{
  std::ostringstream out;
  out << "**ERROR** evalYieldCondition with a Matrix3 argument should not be "
         "called by "
      << " models that use the SoilMix yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return 0.0;
}

//--------------------------------------------------------------
// Other derivatives

/*! Derivative with respect to the \f$xi\f$ where \f$\xi = s \f$
    where \f$s\f$ is deviatoric part of Cauchy stress */
Matrix3
YieldCond_ArenaMixture::df_dxi([[maybe_unused]] const Matrix3& stress,
                               const ModelStateBase*)

{
  std::ostringstream out;
  out << "**ERROR** df_dxi with a Matrix3 argument should not be "
      << "called by models that use the SoilMix yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return Matrix3(0.0);
}

/* Derivative with respect to \f$ s \f$ and \f$ \beta \f$ */
std::pair<Matrix3, Matrix3>
YieldCond_ArenaMixture::df_dsigmaDev_dbeta([[maybe_unused]] const Matrix3& stress,
                                           const ModelStateBase*)
{
  std::ostringstream out;
  out << "**ERROR** df_dsigmaDev_dbeta with a Matrix3 argument should not be "
      << "called by models that use the SoilMix yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);

  return std::make_pair(Matrix3(0.0), Matrix3(0.0));
}

//--------------------------------------------------------------
// Tangent moduli
void
YieldCond_ArenaMixture::computeElasPlasTangentModulus(
  [[maybe_unused]] const TangentModulusTensor& Ce,
  [[maybe_unused]] const Matrix3& sigma,
  [[maybe_unused]] double sigY,
  [[maybe_unused]] double dsigYdep,
  [[maybe_unused]] double porosity,
  [[maybe_unused]] double,
  [[maybe_unused]] TangentModulusTensor& Cep)
{
  std::ostringstream out;
  out << "**ERROR** computeElasPlasTangentModulus with a Matrix3 argument "
         "should not be "
      << "called by models that use the SoilMix yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return;
}

void
YieldCond_ArenaMixture::computeTangentModulus([[maybe_unused]] const TangentModulusTensor& Ce,
                                              [[maybe_unused]] const Matrix3& f_sigma,
                                              [[maybe_unused]] double f_q1,
                                              [[maybe_unused]] double h_q1,
                                              [[maybe_unused]] TangentModulusTensor& Cep)
{
  std::ostringstream out;
  out << "**ERROR** coputeTangentModulus with a Matrix3 argument should not be "
      << "called by models that use the SoilMix yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return;
}

/**
 *  This is used to scale and update the yield parameters
 */
void
YieldCond_ArenaMixture::updateLocalVariables(
  ParticleSubset* pset,
  DataWarehouse* old_dw,
  DataWarehouse* new_dw,
  constParticleVariable<double>& pCoherence_old,
  const ParticleVariable<double>& pCoherence_new)
{
  constParticleVariable<double> pPEAKI1_old, pFSLOPE_old, pSTREN_old,
    pYSLOPE_old;
  constParticleVariable<double> pBETA_old, pCR_old, pT1_old, pT2_old;
  old_dw->get(pPEAKI1_old, pPEAKI1Label, pset);
  old_dw->get(pFSLOPE_old, pFSLOPELabel, pset);
  old_dw->get(pSTREN_old, pSTRENLabel, pset);
  old_dw->get(pYSLOPE_old, pYSLOPELabel, pset);
  old_dw->get(pBETA_old, pBETALabel, pset);
  old_dw->get(pCR_old, pCRLabel, pset);
  old_dw->get(pT1_old, pT1Label, pset);
  old_dw->get(pT2_old, pT2Label, pset);

  ParticleVariable<double> pPEAKI1_new, pFSLOPE_new, pSTREN_new, pYSLOPE_new;
  ParticleVariable<double> pBETA_new, pCR_new, pT1_new, pT2_new;
  new_dw->allocateAndPut(pPEAKI1_new, pPEAKI1Label_preReloc, pset);
  new_dw->allocateAndPut(pFSLOPE_new, pFSLOPELabel_preReloc, pset);
  new_dw->allocateAndPut(pSTREN_new, pSTRENLabel_preReloc, pset);
  new_dw->allocateAndPut(pYSLOPE_new, pYSLOPELabel_preReloc, pset);
  new_dw->allocateAndPut(pBETA_new, pBETALabel_preReloc, pset);
  new_dw->allocateAndPut(pCR_new, pCRLabel_preReloc, pset);
  new_dw->allocateAndPut(pT1_new, pT1Label_preReloc, pset);
  new_dw->allocateAndPut(pT2_new, pT2Label_preReloc, pset);

  double PEAKI1_failed = d_volfrac[0] * d_yieldParam[0].PEAKI1_failed +
                         d_volfrac[1] * d_yieldParam[1].PEAKI1_failed;
  double FSLOPE_failed = d_volfrac[0] * d_yieldParam[0].FSLOPE_failed +
                         d_volfrac[1] * d_yieldParam[1].FSLOPE_failed;
  double STREN_failed = d_volfrac[0] * d_yieldParam[0].STREN_failed +
                        d_volfrac[1] * d_yieldParam[1].STREN_failed;
  double YSLOPE_failed = d_volfrac[0] * d_yieldParam[0].YSLOPE_failed +
                         d_volfrac[1] * d_yieldParam[1].YSLOPE_failed;
  for (int idx : *pset) {
    // Get the coherence values
    double coher_old = pCoherence_old[idx];
    double coher_new = pCoherence_new[idx];

    // Compute intact values of the parameters
    double PEAKI1_intact =
      (pPEAKI1_old[idx] - (1.0 - coher_old) * PEAKI1_failed) / coher_old;
    double FSLOPE_intact =
      (pFSLOPE_old[idx] - (1.0 - coher_old) * FSLOPE_failed) / coher_old;
    double YSLOPE_intact =
      (pYSLOPE_old[idx] - (1.0 - coher_old) * YSLOPE_failed) / coher_old;
    double STREN_intact =
      (pSTREN_old[idx] - (1.0 - coher_old) * STREN_failed) / coher_old;

    // Compute the damaged values of the parameters
    pPEAKI1_new[idx] =
      coher_new * PEAKI1_intact + (1.0 - coher_new) * PEAKI1_failed;
    pFSLOPE_new[idx] =
      coher_new * FSLOPE_intact + (1.0 - coher_new) * FSLOPE_failed;
    pSTREN_new[idx] =
      coher_new * STREN_intact + (1.0 - coher_new) * STREN_failed;
    pYSLOPE_new[idx] =
      coher_new * YSLOPE_intact + (1.0 - coher_new) * YSLOPE_failed;

    // Copy the other parameters
    pBETA_new[idx] = pBETA_old[idx];
    pCR_new[idx]   = pCR_old[idx];
    pT1_new[idx]   = pT1_old[idx];
    pT2_new[idx]   = pT2_old[idx];
  }
}
