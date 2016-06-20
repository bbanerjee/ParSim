/*
 * The MIT License
 *
 * Copyright (c) 2015-2016 Parresia Research Limited, New Zealand
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


#include <CCA/Components/MPM/ConstitutiveModel/Models/YieldCond_MasonSand.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <cmath>

using namespace Vaango;

const double YieldCond_MasonSand::sqrt_three = std::sqrt(3.0);
const double YieldCond_MasonSand::one_sqrt_three = 1.0/sqrt_three;

YieldCond_MasonSand::YieldCond_MasonSand(Uintah::ProblemSpecP& ps)
{
  // Nonlinear Drucker-Prager parameters
  ps->require("PEAKI1", d_yieldParam.PEAKI1);  // Shear Limit Surface Parameter
  ps->require("FSLOPE", d_yieldParam.FSLOPE);  // Shear Limit Surface Parameter
  ps->require("STREN",  d_yieldParam.STREN);   // Shear Limit Surface Parameter
  ps->require("YSLOPE", d_yieldParam.YSLOPE);  // Shear Limit Surface Parameter
  ps->getWithDefault("PEAKI1_failed", d_yieldParam.PEAKI1_failed, 1.0e-5);
  ps->getWithDefault("FSLOPE_failed", d_yieldParam.FSLOPE_failed, 0.5*d_yieldParam.FSLOPE);
  ps->getWithDefault("STREN_failed",  d_yieldParam.STREN_failed,  0.1*d_yieldParam.STREN);
  ps->getWithDefault("YSLOPE_failed", d_yieldParam.YSLOPE_failed, 1.0e-5);

  // Non-associativity parameters
  ps->require("BETA",   d_nonAssocParam.BETA); // Nonassociativity Parameter

  // Cap parameters: CR = (peakI1-kappa)/(peakI1-X)
  ps->require("CR", d_capParam.CR);            // Cap Shape Parameter 

  // Duvaut-Lion rate parameters
  ps->require("T1", d_rateParam.T1);           // Rate dependence parameter
  ps->require("T2", d_rateParam.T2);           // Rate dependence parameter
                                          
  // Check the input parameters
  checkInputParameters();

  // Compute the model parameters from the input parameters
  computeModelParameters(1.0);

  // Now optionally get the variablity information for each parameter
  std::string weibullDist;
  ps->getWithDefault("weibullDist_PEAKI1", weibullDist, std::to_string(d_yieldParam.PEAKI1));
  d_weibull_PEAKI1.WeibullParser(weibullDist);
  proc0cout << d_weibull_PEAKI1 << std::endl;

  ps->getWithDefault("weibullDist_FSLOPE", weibullDist, std::to_string(d_yieldParam.FSLOPE));
  d_weibull_FSLOPE.WeibullParser(weibullDist);
  proc0cout << d_weibull_FSLOPE << std::endl;

  ps->getWithDefault("weibullDist_STREN", weibullDist, std::to_string(d_yieldParam.STREN));
  d_weibull_STREN.WeibullParser(weibullDist);
  proc0cout << d_weibull_STREN << std::endl;

  ps->getWithDefault("weibullDist_YSLOPE", weibullDist, std::to_string(d_yieldParam.YSLOPE));
  d_weibull_YSLOPE.WeibullParser(weibullDist);
  proc0cout << d_weibull_YSLOPE << std::endl;

  ps->getWithDefault("weibullDist_BETA", weibullDist, std::to_string(d_nonAssocParam.BETA));
  d_weibull_BETA.WeibullParser(weibullDist);
  proc0cout << d_weibull_BETA << std::endl;

  ps->getWithDefault("weibullDist_CR", weibullDist, std::to_string(d_capParam.CR));
  d_weibull_CR.WeibullParser(weibullDist);
  proc0cout << d_weibull_CR << std::endl;

  ps->getWithDefault("weibullDist_T1", weibullDist, std::to_string(d_rateParam.T1));
  d_weibull_T1.WeibullParser(weibullDist);
  proc0cout << d_weibull_T1 << std::endl;

  ps->getWithDefault("weibullDist_T2", weibullDist, std::to_string(d_rateParam.T2));
  d_weibull_T2.WeibullParser(weibullDist);
  proc0cout << d_weibull_T2 << std::endl;

  // Initialize local labels for parameter variability
  initializeLocalMPMLabels();
}
         
YieldCond_MasonSand::YieldCond_MasonSand(const YieldCond_MasonSand* yc)
{
  d_modelParam = yc->d_modelParam; 
  d_yieldParam = yc->d_yieldParam; 
  d_nonAssocParam = yc->d_nonAssocParam; 
  d_capParam = yc->d_capParam; 
  d_rateParam = yc->d_rateParam; 

  // Copy parameter variability information
  d_weibull_PEAKI1 = yc->d_weibull_PEAKI1;
  d_weibull_FSLOPE = yc->d_weibull_FSLOPE;
  d_weibull_STREN = yc->d_weibull_STREN;
  d_weibull_YSLOPE = yc->d_weibull_YSLOPE;
  d_weibull_BETA = yc->d_weibull_BETA;
  d_weibull_CR = yc->d_weibull_CR;
  d_weibull_T1 = yc->d_weibull_T1;
  d_weibull_T2 = yc->d_weibull_T2;

  // Initialize local labels for parameter variability
  initializeLocalMPMLabels();
}
         
YieldCond_MasonSand::~YieldCond_MasonSand()
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
YieldCond_MasonSand::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP yield_ps = ps->appendChild("plastic_yield_condition");
  yield_ps->setAttribute("type", "mason_sand");

  yield_ps->appendElement("FSLOPE", d_yieldParam.FSLOPE);
  yield_ps->appendElement("PEAKI1", d_yieldParam.PEAKI1);
  yield_ps->appendElement("STREN",  d_yieldParam.STREN);
  yield_ps->appendElement("YSLOPE", d_yieldParam.YSLOPE);
  yield_ps->appendElement("FSLOPE_failed", d_yieldParam.FSLOPE_failed);
  yield_ps->appendElement("PEAKI1_failed", d_yieldParam.PEAKI1_failed);
  yield_ps->appendElement("STREN_failed",  d_yieldParam.STREN_failed);
  yield_ps->appendElement("YSLOPE_failed", d_yieldParam.YSLOPE_failed);

  yield_ps->appendElement("BETA",   d_nonAssocParam.BETA);

  yield_ps->appendElement("CR", d_capParam.CR);

  yield_ps->appendElement("T1", d_rateParam.T1);
  yield_ps->appendElement("T2", d_rateParam.T2);

  yield_ps->appendElement("weibullDist_PEAKI1", d_weibull_PEAKI1.getWeibDist());
  yield_ps->appendElement("weibullDist_FSLOPE", d_weibull_FSLOPE.getWeibDist());
  yield_ps->appendElement("weibullDist_STREN",  d_weibull_STREN.getWeibDist());
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
YieldCond_MasonSand::checkInputParameters()
{
  std::ostringstream warn;
  if (d_yieldParam.PEAKI1 <0.0 || d_yieldParam.PEAKI1_failed < 0.0) {
    warn << "PEAKI1 must be nonnegative. PEAKI1 = " << d_yieldParam.PEAKI1 << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_yieldParam.FSLOPE<0.0 || d_yieldParam.FSLOPE_failed < 0.0) {
    warn << "FSLOPE must be nonnegative. FSLOPE = " << d_yieldParam.FSLOPE << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_yieldParam.FSLOPE < d_yieldParam.YSLOPE ||
       d_yieldParam.FSLOPE_failed < d_yieldParam.YSLOPE_failed) {
    warn << "FSLOPE must be greater than YSLOPE. FSLOPE = " << d_yieldParam.FSLOPE
         << ", YSLOPE = " << d_yieldParam.YSLOPE << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_nonAssocParam.BETA <= 0.0) {
    warn << "BETA (nonassociativity factor) must be positive. BETA = "
         << d_nonAssocParam.BETA << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_capParam.CR >= 1 || d_capParam.CR <= 0.0) {
    warn << "CR must be 0<CR<1. CR = " << d_capParam.CR << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_rateParam.T1 < 0.0) {
    warn << "T1 must be nonnegative. T1 = "<< d_rateParam.T1 << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_rateParam.T2 < 0.0) {
    warn << "T2 must be nonnegative. T2 = "<< d_rateParam.T2 << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if ( (d_rateParam.T1 > 0.0 || d_rateParam.T2 > 0.0)
       != (d_rateParam.T1 > 0.0 && d_rateParam.T2 > 0.0) ) {
    warn << "For rate dependence both T1 and T2 must be positive. T1 = "
         << d_rateParam.T1 << ", T2 = " << d_rateParam.T2 << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
}

//--------------------------------------------------------------
// Compute the model parameters a1, a2, a3, a4, beta from the 
// input parameters FSLOPE, PEAKI1, STREN, SLOPE, BETA_nonassoc
// The shear limit surface is defined in terms of the a1,a2,a3,a4 parameters, but
// the user inputs are the more intuitive set of FSLOPE. YSLOPE, STREN, and PEAKI1.
//
// Note: This routine computes the a_i parameters from the user inputs.  The code was
// originally written by R.M. Brannon, with modifications by M.S. Swan.
//--------------------------------------------------------------
void 
YieldCond_MasonSand::computeModelParameters(double)
{
  double  FSLOPE = d_yieldParam.FSLOPE,  // Slope at I1=PEAKI1
          STREN  = d_yieldParam.STREN,   // Value of rootJ2 at I1=0
          YSLOPE = d_yieldParam.YSLOPE,  // High pressure slope
          PEAKI1 = d_yieldParam.PEAKI1;  // Value of I1 at strength=0
  double  FSLOPE_failed = d_yieldParam.FSLOPE_failed,  // Slope at I1=PEAKI1
          STREN_failed  = d_yieldParam.STREN_failed,   // Value of rootJ2 at I1=0
          YSLOPE_failed = d_yieldParam.YSLOPE_failed,  // High pressure slope
          PEAKI1_failed = d_yieldParam.PEAKI1_failed;  // Value of I1 at strength=0

  std::vector<double> limitParameters = 
    computeModelParameters(PEAKI1, FSLOPE, STREN, YSLOPE);
  std::vector<double> limitParameters_failed = 
    computeModelParameters(PEAKI1_failed, FSLOPE_failed, STREN_failed, YSLOPE_failed);

  d_modelParam.a1 = limitParameters[0];
  d_modelParam.a2 = limitParameters[1];
  d_modelParam.a3 = limitParameters[2];
  d_modelParam.a4 = limitParameters[3];
  d_modelParam.a1_failed = limitParameters_failed[0];
  d_modelParam.a2_failed = limitParameters_failed[1];
  d_modelParam.a3_failed = limitParameters_failed[2];
  d_modelParam.a4_failed = limitParameters_failed[3];
}
  
std::vector<double> 
YieldCond_MasonSand::computeModelParameters(const double& PEAKI1,
                                            const double& FSLOPE,
                                            const double& STREN,
                                            const double& YSLOPE)
{
  double a1, a2, a3, a4;
  if (FSLOPE > 0.0 && PEAKI1 >= 0.0 && STREN == 0.0 && YSLOPE == 0.0)
  {// ----------------------------------------------Linear Drucker Prager
    a1 = PEAKI1*FSLOPE;
    a2 = 0.0;
    a3 = 0.0;
    a4 = FSLOPE;
  } 
  else if (FSLOPE == 0.0 && PEAKI1 == 0.0 && STREN > 0.0 && YSLOPE == 0.0)
  { // ------------------------------------------------------- Von Mises
    a1 = STREN;
    a2 = 0.0;
    a3 = 0.0;
    a4 = 0.0;
  }
  else if (FSLOPE > 0.0 && YSLOPE == 0.0 && STREN > 0.0 && PEAKI1 == 0.0)
  { // ------------------------------------------------------- 0 PEAKI1 to vonMises
    a1 = STREN;
    a2 = FSLOPE/STREN;
    a3 = STREN;
    a4 = 0.0;
  }
  else if (FSLOPE > YSLOPE && YSLOPE > 0.0 && STREN > YSLOPE*PEAKI1 && PEAKI1 >= 0.0)
  { // ------------------------------------------------------- Nonlinear Drucker-Prager
    a1 = STREN;
    a2 = (FSLOPE-YSLOPE)/(STREN-YSLOPE*PEAKI1);
    a3 = (STREN-YSLOPE*PEAKI1)*exp(-a2*PEAKI1);
    a4 = YSLOPE;
  }
  else
  {
    // Bad inputs, throw exception:
    std::ostringstream warn;
    warn << "Bad input parameters for shear limit surface. "
         << "FSLOPE = " << FSLOPE
         << ", YSLOPE = " << YSLOPE
         << ", PEAKI1 = " << PEAKI1
         << ", STREN = " << STREN << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  std::vector<double> limitParameters = {a1, a2, a3, a4};
  return limitParameters;
}

//--------------------------------------------------------------
// Evaluate yield condition (q = state->q
//                           p = state->p)
// f := J2 - Ff^2*Fc^2 = 0
// where
//     J2 = q^2/3
//     I1 = 3*(p - pbar_w)
//     Ff := a1 - a3*exp(a2*I1) - a4*I1 
//     Fc^2 := 1 - (kappa - 3*p)^2/(kappa - X)^2
//     kappa = I1_0 - CR*(I1_0 - X)
//
// Returns:
//   hasYielded = -1.0 (if elastic)
//              =  1.0 (otherwise)
//--------------------------------------------------------------
double 
YieldCond_MasonSand::evalYieldCondition(const ModelStateBase* state_input)
{
  const ModelState_MasonSand* state = dynamic_cast<const ModelState_MasonSand*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_MasonSand.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

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
  double X_eff = state->capX + 3.0*state->pbar_w;
  double kappa = state->kappa;

  // Initialize hasYielded to -1
  double hasYielded = -1.0;

  // Cauchy stress invariants: I1_eff = 3*(p + pbar_w), J2 = q^2/3
  double I1_eff = state->I1_eff;
  double sqrt_J2 = state->sqrt_J2;

  // --------------------------------------------------------------------
  // *** SHEAR LIMIT FUNCTION (Ff) ***
  // --------------------------------------------------------------------
  double Ff = a1 - a3*exp(a2*I1_eff) - a4*I1_eff;

  // --------------------------------------------------------------------
  // *** Branch Point (Kappa) ***
  // --------------------------------------------------------------------
  kappa = PEAKI1 - CR*(PEAKI1 - X_eff); // Branch Point

  // --------------------------------------------------------------------
  // *** COMPOSITE YIELD FUNCTION ***
  // --------------------------------------------------------------------
  // Evaluate Composite Yield Function F(I1) = Ff(I1)*fc(I1) in each region.
  // The elseif statements have nested if statements, which is not equivalent
  // to them having a single elseif(A&&B&&C)
  if (I1_eff < X_eff) {//---------------------------------------------------(I1<X)
    hasYielded = 1.0;
    //std::cout << " I1_eff < X_eff " << I1_eff << "," << X_eff << std::endl;
    return hasYielded;
  }

  // **Elliptical Cap Function: (fc)**
  // fc = sqrt(1.0 - Pow((Kappa-I1mZ)/(Kappa-X)),2.0);
  // faster version: fc2 = fc^2
  // **WARNING** p3 is the maximum achievable volumetric plastic strain in compresson
  // so if a value of 0 has been specified this indicates the user
  // wishes to run without porosity, and no cap function is used, i.e. fc=1
  if ((X_eff < I1_eff) && (I1_eff < kappa)) {// ---------------(X<I1<kappa)

    double kappaRatio = (kappa - I1_eff)/(kappa - X_eff);
    double fc2 = 1.0 - kappaRatio*kappaRatio;
    if (sqrt_J2*sqrt_J2 > Ff*Ff*fc2 ) {
      //std::cout << " X_eff < I1_eff " << I1_eff << "," << X_eff << std::endl;
      //std::cout << " I1_eff < kappa " << I1_eff << "," << kappa << std::endl;
      //std::cout << " J2 < Ff^2*Fc^2 " << sqrt_J2 << "," << Ff << ", " << fc2 << std::endl;
      hasYielded = 1.0;
    }
  } else { // --------- X >= I1 or kappa <= I1

    if (I1_eff <= PEAKI1) { // ----- (kappa <= I1 <= PEAKI1)
      if (sqrt_J2 > Ff) {
        //std::cout << " I1_eff < PEAKI1 " << I1_eff << "," << PEAKI1 << std::endl;
        //std::cout << " sqrt(J2) > Ff " << sqrt_J2 << "," << Ff << std::endl;
        hasYielded = 1.0;
      }
    } else { // I1 > PEAKI1 
      //std::cout << " I1_eff > PEAKI1 " << I1_eff << "," << PEAKI1 << std::endl;
      hasYielded = 1.0;
    }
  }

  return hasYielded;
}

//--------------------------------------------------------------
// Evaluate yield condition max (q = state->q
//                               p = state->p)
//--------------------------------------------------------------
double 
YieldCond_MasonSand::evalYieldConditionMax(const ModelStateBase* )
{
  std::ostringstream out;
  out << "**ERROR** evalYieldConditionMax should not be called by "
      << " models that use the MasonSand yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
         
  return 0.0;
}

//--------------------------------------------------------------
// Derivatives needed by return algorithms and Newton iterations

//--------------------------------------------------------------
// Compute df/dp  where p = volumetric stress = 1/3 Tr(sigma)
//   df/dp = derivative of the yield function wrt p
//
// f := J2 - Ff^2*Fc^2 = 0
// where
//     J2 = q^2/3
//     I1 = 3*(p - pbar_w)
//     Ff := a1 - a3*exp(a2*I1) - a4*I1 
//     Fc^2 := 1 - (kappa - 3*p)^2/(kappa - X)^2
//     kappa = I1_0 - CR*(I1_0 - X)
//
// df/dp = 6*Ff*(a2*a3*exp(a2*I1) + a4)*Fc^2 - 
//             6*Ff^2*(kappa - I1)/(kappa - X)^2
//--------------------------------------------------------------
double 
YieldCond_MasonSand::computeVolStressDerivOfYieldFunction(const ModelStateBase* state_input)
{
  const ModelState_MasonSand* state = dynamic_cast<const ModelState_MasonSand*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_MasonSand.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

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
  double X_eff = state->capX + 3.0*state->pbar_w;
  double kappa = state->kappa;

  // Cauchy stress invariants: I1 = 3*p, J2 = q^2/3
  double I1_eff = state->I1_eff;

  // --------------------------------------------------------------------
  // *** SHEAR LIMIT FUNCTION (Ff) ***
  // --------------------------------------------------------------------
  double Ff = a1 - a3*exp(a2*I1_eff) - a4*I1_eff;

  // --------------------------------------------------------------------
  // *** Branch Point (Kappa) ***
  // --------------------------------------------------------------------
  kappa = PEAKI1 - CR*(PEAKI1 - X_eff); // Branch Point

  // --------------------------------------------------------------------
  // **Elliptical Cap Function: (fc)**
  // --------------------------------------------------------------------
  double kappaRatio = (kappa - I1_eff)/(kappa - X_eff);
  double Fc_sq = 1.0 - kappaRatio*kappaRatio;

  // --------------------------------------------------------------------
  // Derivatives
  // --------------------------------------------------------------------
  // term1 = 6*Ff*(a2*a3*exp(a2*I1) + a4)*Fc^2  
  double term1 = 6.0*Ff*Fc_sq*(a2*a3*exp(a2*I1_eff) + a4);
  // term2 = 6*Ff^2*(kappa - I1)/(kappa - X)^2
  double term2 = 6.0*Ff*Ff*kappaRatio/(kappa - X_eff);

  return (term1 - term2);
}

//--------------------------------------------------------------
// Compute df/dq  
// f := J2 - Ff^2*Fc^2 = 0
// where
//     J2 = q^2/3
// df/dq = 2q/3
//--------------------------------------------------------------
double 
YieldCond_MasonSand::computeDevStressDerivOfYieldFunction(const ModelStateBase* state_input)
{
  const ModelState_MasonSand* state = dynamic_cast<const ModelState_MasonSand*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_MasonSand.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  return (2.0/std::sqrt(3.0))*state->sqrt_J2;
}

//--------------------------------------------------------------
// Compute d/depse_v(df/dp)
//   df/dp = 6*Ff*(a2*a3*exp(3*a2*p) + a4)*Fc^2 - 
//             6*Ff^2*(kappa - I1)/(kappa - X)^2
//   d/depse_v(df/dp) = 
//
// Requires:  Equation of state and internal variable
//--------------------------------------------------------------
double
YieldCond_MasonSand::computeVolStrainDerivOfDfDp(const ModelStateBase* state_input,
                                                 const PressureModel* eos,
                                                 const ShearModulusModel* ,
                                                 const InternalVariableModel* )
{
  std::ostringstream out;
  out << "**ERROR** computeVolStrainDerivOfDfDp should not be called by "
      << " models that use the MasonSand yield criterion.";
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
YieldCond_MasonSand::computeDevStrainDerivOfDfDp(const ModelStateBase* state_input,
                                                 const PressureModel* eos,
                                                 const ShearModulusModel* ,
                                                 const InternalVariableModel* )
{
  std::ostringstream out;
  out << "**ERROR** computeDevStrainDerivOfDfDp should not be called by "
      << " models that use the MasonSand yield criterion.";
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
YieldCond_MasonSand::computeVolStrainDerivOfDfDq(const ModelStateBase* state_input,
                                                 const PressureModel* ,
                                                 const ShearModulusModel* shear,
                                                 const InternalVariableModel* )
{
  std::ostringstream out;
  out << "**ERROR** computeVolStrainDerivOfDfDq should not be called by "
      << " models that use the MasonSand yield criterion.";
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
YieldCond_MasonSand::computeDevStrainDerivOfDfDq(const ModelStateBase* state_input,
                                                 const PressureModel* ,
                                                 const ShearModulusModel* shear,
                                                 const InternalVariableModel* )
{
  std::ostringstream out;
  out << "**ERROR** computeDevStrainDerivOfDfDq should not be called by "
      << " models that use the MasonSand yield criterion.";
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
YieldCond_MasonSand::computeVolStrainDerivOfYieldFunction(const ModelStateBase* state_input,
                                                          const PressureModel* eos,
                                                          const ShearModulusModel* shear,
                                                          const InternalVariableModel* )
{
  std::ostringstream out;
  out << "**ERROR** computeVolStrainDerivOfYieldFunction should not be called by "
      << " models that use the MasonSand yield criterion.";
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
YieldCond_MasonSand::computeDevStrainDerivOfYieldFunction(const ModelStateBase* state_input,
                                                          const PressureModel* eos,
                                                          const ShearModulusModel* shear,
                                                          const InternalVariableModel* )
{
  std::ostringstream out;
  out << "**ERROR** computeVolStrainDerivOfYieldFunction should not be called by "
      << " models that use the MasonSand yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
         
  return 0.0;
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
YieldCond_MasonSand::getInternalPoint(const ModelStateBase* state_old_input,
                                      const ModelStateBase* state_trial_input)
{
  const ModelState_MasonSand* state_old = 
    dynamic_cast<const ModelState_MasonSand*>(state_old_input);
  const ModelState_MasonSand* state_trial = 
    dynamic_cast<const ModelState_MasonSand*>(state_trial_input);
  if ((!state_old) || (!state_trial)) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_MasonSand.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  // Compute effective trial stress
  double  I1_eff_trial = state_trial->I1_eff - state_trial->pbar_w + state_old->pbar_w;

  // Get the particle specific internal variables from the model state
  double PEAKI1 = state_old->yieldParams.at("PEAKI1");

  // It may be better to use an interior point at the center of the yield surface, rather than at 
  // pbar_w, in particular when PEAKI1=0.  Picking the midpoint between PEAKI1 and X would be 
  // problematic when the user has specified some no porosity condition (e.g. p0=-1e99)
  double I1_eff_interior = 0.0;
  double upperI1 = PEAKI1;
  if (I1_eff_trial < upperI1) {
    if (I1_eff_trial > state_old->capX + 3.0*state_old->pbar_w) { // Trial is above yield surface
      I1_eff_interior = state_trial->I1_eff;
    } else { // Trial is past X, use yield midpoint as interior point
      I1_eff_interior = -3.0*state_old->pbar_w + 0.5*(PEAKI1 + state_old->capX + 3.0*state_old->pbar_w);
    }
  } else { // I1_trial + pbar_w >= I1_peak => Trial is past vertex
    double lTrial = sqrt(I1_eff_trial*I1_eff_trial + state_trial->sqrt_J2*state_trial->sqrt_J2);
    double lYield = 0.5*(PEAKI1 - state_old->capX - 3.0*state_old->pbar_w);
    I1_eff_interior = -3.0*state_old->pbar_w + upperI1 - std::min(lTrial, lYield);
  }
  
  return I1_eff_interior;
}

/**
 * Function: getClosestPoint
 *
 * Purpose: Get the point on the yield surface that is closest to a given point (2D)
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
YieldCond_MasonSand::getClosestPoint(const ModelStateBase* state_input,
                                     const double& px, const double& py,
                                     double& cpx, double& cpy)
{
  const ModelState_MasonSand* state = 
    dynamic_cast<const ModelState_MasonSand*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_MasonSand.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  // Save the input point
  point_type pt(px, py);

  // Get the yield surface points
  int num_points = 1000;
  std::vector<point_type> poly_points = getYieldSurfacePointsAll_RprimeZ(state, num_points);

  // Find two yield surface segments that are closest to input point
  std::vector<point_type> segments = getClosestSegments(pt, poly_points);
  cpx = boost::geometry::get<0>(segments[1]);
  cpy = boost::geometry::get<1>(segments[1]);

  // Discretize the closest segments
  std::vector<point_type> segment_points = 
     getYieldSurfacePointsSegment_RprimeZ(state, segments[0], segments[2], num_points);

  // Find the closest point
  point_type cpt = findClosestPoint(pt, segment_points);
  cpx = boost::geometry::get<0>(cpt);
  cpy = boost::geometry::get<1>(cpt);

  /*
  if (cpx == boost::numeric::bounds<double>::highest()) {
    return false;
  }
  */

  return true;
}

/* Get the points on the yield surface */
std::vector<point_type> 
YieldCond_MasonSand::getYieldSurfacePointsAll_RprimeZ(const ModelState_MasonSand* state,
                                                      const int& num_points)
{

  // Get the particle specific internal variables from the model state
  double PEAKI1 = state->yieldParams.at("PEAKI1");

  // Get the plastic internal variables from the model state
  double pbar_w = state->pbar_w;
  double X_eff = state->capX + 3.0*pbar_w;

   // Set up I1 values
  double I1eff_min = 0.99999*X_eff;
  double I1eff_max = 0.99999*PEAKI1;

  // Compute z_eff and r'
  std::vector<double> z_eff_vec, rprime_vec;
  computeZeff_and_RPrime(state, I1eff_min, I1eff_max, num_points, z_eff_vec, 
                         rprime_vec);

  // Create a point_type vector
  std::vector<point_type> polyline = polylineFromRelectedPoints(z_eff_vec, rprime_vec);

  // Add the first point to close the polygon
  polyline.push_back(point_type(*(z_eff_vec.begin()), *(rprime_vec.begin())));
  //std::cout << std::endl;

  return polyline;
}

/* Get the points on the yield surface */
std::vector<point_type> 
YieldCond_MasonSand::getYieldSurfacePointsSegment_RprimeZ(const ModelState_MasonSand* state,
                                                          const point_type& start_point,
                                                          const point_type& end_point,
                                                          const int& num_points)
{

  // Find the start I1 and end I1 values of the segments
  // **TODO** make sure that the start and end points are differenet
  double z_effStart = boost::geometry::get<0>(start_point);
  double z_effEnd = boost::geometry::get<0>(end_point);
  double I1_effStart = std::sqrt(3.0)*z_effStart;
  double I1_effEnd = std::sqrt(3.0)*z_effEnd;

  // Compute z_eff and r'
  std::vector<double> z_eff_vec, rprime_vec;
  computeZeff_and_RPrime(state, I1_effStart, I1_effEnd, num_points, z_eff_vec, 
                         rprime_vec);

  // Create a point_type vector
  std::vector<point_type> polyline = polylineFromRelectedPoints(z_eff_vec, rprime_vec);

  return polyline;
}

/*! Compute a vector of z_eff, r' values given a range of I1_eff values */
void
YieldCond_MasonSand::computeZeff_and_RPrime(const ModelState_MasonSand* state,
                                            const double& I1eff_min,
                                            const double& I1eff_max,
                                            const int& num_points,
                                            std::vector<double>& z_eff_vec,
                                            std::vector<double>& rprime_vec) 
{
  // Get the particle specific internal variables from the model state
  double PEAKI1 = state->yieldParams.at("PEAKI1");
  double FSLOPE = state->yieldParams.at("FSLOPE");
  double STREN  = state->yieldParams.at("STREN");
  double YSLOPE = state->yieldParams.at("YSLOPE");
  double BETA   = state->yieldParams.at("BETA");
  double CR     = state->yieldParams.at("CR");

  std::vector<double> limitParameters = 
    computeModelParameters(PEAKI1, FSLOPE, STREN, YSLOPE);
  double a1 = limitParameters[0];
  double a2 = limitParameters[1];
  double a3 = limitParameters[2];
  double a4 = limitParameters[3];

  // Get the plastic internal variables from the model state
  double X_eff = state->capX + 3.0*state->pbar_w;
  double kappa =  PEAKI1 - CR*(PEAKI1 - X_eff);

  // Get the bulk and shear moduli and compute sqrt(3/2 K/G)
  double sqrtKG = std::sqrt(1.5*state->bulkModulus/state->shearModulus);

  // Set up I1 values
  std::vector<double> I1_eff_vec = linspace(I1eff_min, I1eff_max, num_points);
  for (auto I1_eff : I1_eff_vec) {

    // Compute F_f
    double Ff = a1 - a3*std::exp(a2*I1_eff) - a4*(I1_eff);
    double Ff_sq = Ff*Ff;

    // Compute Fc
    double Fc_sq = 1.0;
    if ((I1_eff < kappa) && (X_eff < I1_eff)) {
      double ratio = (kappa - I1_eff)/(kappa - X_eff);
      Fc_sq = 1.0 - ratio*ratio;
    }

    // Compute J2
    double J2 = Ff_sq*Fc_sq;
    z_eff_vec.push_back(I1_eff/std::sqrt(3.0));
    rprime_vec.push_back(BETA*std::sqrt(2.0*J2)*sqrtKG);
  }

  return;
}

/* linspace function */
std::vector<double> 
YieldCond_MasonSand::linspace(double start, double end, int num)
{
  double delta = (end - start) / (double)num;

  std::vector<double> linspaced;
  for (int i=0; i < num+1; ++i) {
    linspaced.push_back(start + delta * (double) i);
  }
  return linspaced;
}

/*! Create a polyline containing the original and reflected points
    after reflecting the r' values about the z_eff axis */
std::vector<point_type>
YieldCond_MasonSand::polylineFromRelectedPoints(const std::vector<double> z_eff_vec,
                                                const std::vector<double> rprime_vec)
{
  // Create a point_type vector
  // Iterate forward and add points to the positive r' side of the yield polygon
  std::vector<point_type> polyline;
  auto z_iter = z_eff_vec.begin();
  auto r_iter = rprime_vec.begin();
  while (z_iter != z_eff_vec.end() || r_iter != rprime_vec.end()) {
    double z_eff = *z_iter;
    double r = *r_iter;
    //std::cout << "(" << z_eff << "," << r << ")";
    polyline.push_back(point_type(z_eff, r));
    ++z_iter; ++r_iter;
  }
  //std::cout << std::endl;

  // Iterate in reverse and add points to the negative r' side of the yield polygon
  auto rev_z_iter = z_eff_vec.rbegin();
  auto rev_r_iter = rprime_vec.rbegin();
  while (rev_z_iter != z_eff_vec.rend() || rev_r_iter != rprime_vec.rend()) {
    double z_eff = *rev_z_iter;
    double r = *rev_r_iter;
    //std::cout << "(" << z_eff << "," << -r << ")";
    polyline.push_back(point_type(z_eff, -r));
    ++rev_z_iter; ++rev_r_iter;
  }

  return polyline;
}

/* Find two yield surface segments that are closest to input point */
std::vector<point_type> 
YieldCond_MasonSand::getClosestSegments(const point_type& pt, 
                                        const std::vector<point_type> poly)
{
  // Set up the first segment to start from the end of the polygon
  // **TODO** Make sure that the second to last point is being chosen because the
  //          polygon has been closed
  point_type p_prev = *(poly.rbegin()+1);

  // Set up the second segment to start from the beginning of the polygon
  auto iterNext = poly.begin();
  ++iterNext;
  point_type p_next = *iterNext;
  point_type min_p_prev, min_p, min_p_next;

  double min_dSq = boost::numeric::bounds<double>::highest();

  // Loop through the polygon
  for (auto poly_pt : poly) {

    //std::cout << "Poly pt = " << boost::geometry::dsv(poly_pt)
    //          << " Next = " << boost::geometry::dsv(*iterNext) << std::endl;

    // Compute distance sq
    double dSq = boost::geometry::distance(pt, poly_pt);
    if (dSq < min_dSq) {
      min_dSq = dSq;
      min_p = poly_pt;
      min_p_prev = p_prev;
      min_p_next = p_next;
    }

    // Since the polygon is closed, ignore the last point
    ++iterNext;
    if (iterNext == poly.end()) {
      break;
    }
   
    // Update prev and next
    p_prev = poly_pt;
    p_next = *iterNext; 
  }
 
  // Return the three points
  std::vector<point_type> segments;
  segments.push_back(min_p_prev);
  segments.push_back(min_p);
  segments.push_back(min_p_next);
  //std::cout << "Closest segments = " 
  //          << boost::geometry::dsv(min_p_prev)
  //          << boost::geometry::dsv(min_p)
  //          << boost::geometry::dsv(min_p_next) << std::endl;

  return segments;

}

/* Get the closest point on the yield surface */
point_type 
YieldCond_MasonSand::findClosestPoint(const point_type& p, const std::vector<point_type>& poly)
{
  // Get point coordinates
  double xx = boost::geometry::get<0>(p);
  double yy = boost::geometry::get<1>(p);

  std::vector<point_type> XP;

  // Loop through the segments of the polyline
  auto iterStart = poly.begin();
  auto iterEnd   = poly.end();
  auto iterNext = iterStart;
  ++iterNext;
  for ( ; iterNext != iterEnd; ++iterStart, ++iterNext) {
    double xstart = boost::geometry::get<0>(*iterStart);
    double ystart = boost::geometry::get<1>(*iterStart);
    double xnext = boost::geometry::get<0>(*iterNext);
    double ynext = boost::geometry::get<1>(*iterNext);

    // segments that connect the vertices
    double xab = xnext - xstart;
    double yab = ynext - ystart;

    // segment length (squared)
    double abSq = xab*xab + yab*yab;

    // find the projection of point p = (x,y) on each segment
    double xpa = xx - xstart;
    double ypa = yy - ystart;

    // find t = (p - a)/(b - a);
    double pa_dot_ab = xpa*xab + ypa*yab;
    double tt = pa_dot_ab/abSq;

    // Find projction point
    if (tt < 0.0) {
      XP.push_back(point_type(xstart, ystart));
    } else if (tt > 1.0) {
      XP.push_back(point_type(xnext, ynext));
    } else {
      //std::cout << " tt = " << tt << " xp = " <<  xstart + tt*xab << " yp = " << ystart + tt*yab
      //          << std::endl;
      double xp = xstart + tt*xab;
      double yp = ystart + tt*yab;
      XP.push_back(point_type(xp, yp));
    }
  }

  /*
  if (!(XP.size() > 0)) {
    std::cout << "No closest point" << std::endl;
    return point_type(boost::numeric::bounds<double>::highest(), 0.0);
  }
  */

  point_type min_p;
  double min_d = boost::numeric::bounds<double>::highest();
  BOOST_FOREACH(point_type const& xp, XP) {
    double d = boost::geometry::comparable_distance(p, xp);
    if (d < min_d) {
      min_d = d;
      min_p = xp;
    }
  }

  //std::cout << "Closest: " << boost::geometry::dsv(min_p) << std::endl
  //          << "At: " << boost::geometry::distance(p, min_p) << std::endl;
  return min_p;
}

//--------------------------------------------------------------
// Other yield condition functions

// Evaluate the yield function.
double 
YieldCond_MasonSand::evalYieldCondition(const double p,
                                        const double q,
                                        const double dummy0,
                                        const double dummy1,
                                        double& dummy2)
{
  std::ostringstream out;
  out << "**ERROR** Deprecated evalYieldCondition with double arguments. "
      << " Should not be called by models that use the MasonSand yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
         
  return 0.0;
}

// Evaluate yield condition (s = deviatoric stress
//                           p = state->p)
double 
YieldCond_MasonSand::evalYieldCondition(const Uintah::Matrix3& ,
                                        const ModelStateBase* state_input)
{
  std::ostringstream out;
  out << "**ERROR** evalYieldCondition with a Matrix3 argument should not be called by "
      << " models that use the MasonSand yield criterion.";
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
YieldCond_MasonSand::evalDerivOfYieldFunction(const Uintah::Matrix3& sig,
                                              const double p_c,
                                              const double ,
                                              Uintah::Matrix3& derivative)
{
  std::ostringstream out;
  out << "**ERROR** evalDerivOfYieldCondition with a Matrix3 argument should not be "
      << "called by models that use the MasonSand yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
         
  return;
}

// Compute df/ds  where s = deviatoric stress
//    df/ds = 
void 
YieldCond_MasonSand::evalDevDerivOfYieldFunction(const Uintah::Matrix3& sigDev,
                                                 const double ,
                                                 const double ,
                                                 Uintah::Matrix3& derivative)
{
  std::ostringstream out;
  out << "**ERROR** evalDerivOfYieldCondition with a Matrix3 argument should not be "
      << "called by models that use the MasonSand yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
         
  return;
}

/*! Derivative with respect to the Cauchy stress (\f$\sigma \f$) */
void 
YieldCond_MasonSand::eval_df_dsigma(const Matrix3& sig,
                                    const ModelStateBase* state_input,
                                    Matrix3& df_dsigma)
{
  std::ostringstream out;
  out << "**ERROR** eval_df_dsigma with a Matrix3 argument should not be "
      << "called by models that use the MasonSand yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
         
  return;
}

/*! Derivative with respect to the \f$xi\f$ where \f$\xi = s \f$  
    where \f$s\f$ is deviatoric part of Cauchy stress */
void 
YieldCond_MasonSand::eval_df_dxi(const Matrix3& sigDev,
                                 const ModelStateBase* ,
                                 Matrix3& df_ds)
         
{
  std::ostringstream out;
  out << "**ERROR** eval_df_dxi with a Matrix3 argument should not be "
      << "called by models that use the MasonSand yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return;
}

/* Derivative with respect to \f$ s \f$ and \f$ \beta \f$ */
void 
YieldCond_MasonSand::eval_df_ds_df_dbeta(const Matrix3& sigDev,
                                           const ModelStateBase*,
                                           Matrix3& df_ds,
                                           Matrix3& df_dbeta)
{
  std::ostringstream out;
  out << "**ERROR** eval_df_ds_df_dbeta with a Matrix3 argument should not be "
      << "called by models that use the MasonSand yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return;
}

/*! Derivative with respect to the plastic strain (\f$\epsilon^p \f$) */
double 
YieldCond_MasonSand::eval_df_dep(const Matrix3& ,
                                 const double& dsigy_dep,
                                 const ModelStateBase* )
{
  std::ostringstream out;
  out << "**ERROR** eval_df_dep with a Matrix3 argument should not be "
      << "called by models that use the MasonSand yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return 0.0;
}

/*! Derivative with respect to the porosity (\f$\epsilon^p \f$) */
double 
YieldCond_MasonSand::eval_df_dphi(const Matrix3& ,
                                  const ModelStateBase* )
{
  std::ostringstream out;
  out << "**ERROR** eval_df_dphi with a Matrix3 argument should not be "
      << "called by models that use the MasonSand yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return 0.0;
}

/*! Compute h_alpha  where \f$d/dt(ep) = d/dt(gamma)~h_{\alpha}\f$ */
double 
YieldCond_MasonSand::eval_h_alpha(const Matrix3& ,
                                    const ModelStateBase* )
{
  std::ostringstream out;
  out << "**ERROR** eval_h_alpha with a Matrix3 argument should not be "
      << "called by models that use the MasonSand yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return 1.0;
}

/*! Compute h_phi  where \f$d/dt(phi) = d/dt(gamma)~h_{\phi}\f$ */
double 
YieldCond_MasonSand::eval_h_phi(const Matrix3& ,
                                  const double& ,
                                  const ModelStateBase* )
{
  std::ostringstream out;
  out << "**ERROR** eval_h_phi with a Matrix3 argument should not be "
      << "called by models that use the MasonSand yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return 0.0;
}

//--------------------------------------------------------------
// Tangent moduli
void 
YieldCond_MasonSand::computeElasPlasTangentModulus(const TangentModulusTensor& Ce,
                                                   const Matrix3& sigma, 
                                                   double sigY,
                                                   double dsigYdep,
                                                   double porosity,
                                                   double ,
                                                   TangentModulusTensor& Cep)
{
  std::ostringstream out;
  out << "**ERROR** computeElasPlasTangentModulus with a Matrix3 argument should not be "
      << "called by models that use the MasonSand yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return;
}

void 
YieldCond_MasonSand::computeTangentModulus(const TangentModulusTensor& Ce,
                                           const Matrix3& f_sigma, 
                                           double f_q1,
                                           double h_q1,
                                           TangentModulusTensor& Cep)
{
  std::ostringstream out;
  out << "**ERROR** coputeTangentModulus with a Matrix3 argument should not be "
      << "called by models that use the MasonSand yield criterion.";
  throw InternalError(out.str(), __FILE__, __LINE__);
  return;
}


/**
 *  This is used to scale and update the yield parameters 
 */
void 
YieldCond_MasonSand::updateLocalVariables(ParticleSubset* pset,
                                          DataWarehouse* old_dw,
                                          DataWarehouse* new_dw,
                                          constParticleVariable<double>& pCoherence_old,
                                          const ParticleVariable<double>& pCoherence_new)
{
  constParticleVariable<double> pPEAKI1_old, pFSLOPE_old, pSTREN_old, pYSLOPE_old; 
  constParticleVariable<double> pBETA_old, pCR_old, pT1_old, pT2_old;
  old_dw->get(pPEAKI1_old, pPEAKI1Label,    pset);
  old_dw->get(pFSLOPE_old, pFSLOPELabel,    pset);
  old_dw->get(pSTREN_old,  pSTRENLabel,     pset);
  old_dw->get(pYSLOPE_old, pYSLOPELabel,    pset);
  old_dw->get(pBETA_old,   pBETALabel,      pset);
  old_dw->get(pCR_old,     pCRLabel,        pset);
  old_dw->get(pT1_old,     pT1Label,        pset);
  old_dw->get(pT2_old,     pT2Label,        pset);

  ParticleVariable<double> pPEAKI1_new, pFSLOPE_new, pSTREN_new, pYSLOPE_new; 
  ParticleVariable<double> pBETA_new, pCR_new, pT1_new, pT2_new;
  new_dw->allocateAndPut(pPEAKI1_new, pPEAKI1Label_preReloc,    pset);
  new_dw->allocateAndPut(pFSLOPE_new, pFSLOPELabel_preReloc,    pset);
  new_dw->allocateAndPut(pSTREN_new,  pSTRENLabel_preReloc,     pset);
  new_dw->allocateAndPut(pYSLOPE_new, pYSLOPELabel_preReloc,    pset);
  new_dw->allocateAndPut(pBETA_new,   pBETALabel_preReloc,      pset);
  new_dw->allocateAndPut(pCR_new,     pCRLabel_preReloc,        pset);
  new_dw->allocateAndPut(pT1_new,     pT1Label_preReloc,        pset);
  new_dw->allocateAndPut(pT2_new,     pT2Label_preReloc,        pset);

  double PEAKI1_failed = d_yieldParam.PEAKI1_failed;
  double FSLOPE_failed = d_yieldParam.FSLOPE_failed;
  double STREN_failed = d_yieldParam.STREN_failed;
  double YSLOPE_failed = d_yieldParam.YSLOPE_failed;
  for (auto iter = pset->begin(); iter != pset->end(); iter++) {
    particleIndex idx = *iter;

    // Get the coherence values
    double coher_old = pCoherence_old[idx];
    double coher_new = pCoherence_new[idx];

    // Compute intact values of the parameters
    double PEAKI1_intact = (pPEAKI1_old[idx] - (1.0 - coher_old)*PEAKI1_failed)/coher_old;
    double FSLOPE_intact = (pFSLOPE_old[idx] - (1.0 - coher_old)*FSLOPE_failed)/coher_old;
    double YSLOPE_intact = (pYSLOPE_old[idx] - (1.0 - coher_old)*YSLOPE_failed)/coher_old;
    double STREN_intact = (pSTREN_old[idx] - (1.0 - coher_old)*STREN_failed)/coher_old;

    // Compute the damaged values of the parameters
    pPEAKI1_new[idx]    = coher_new*PEAKI1_intact + (1.0 - coher_new)*PEAKI1_failed;
    pFSLOPE_new[idx]    = coher_new*FSLOPE_intact + (1.0 - coher_new)*FSLOPE_failed;
    pSTREN_new[idx]     = coher_new*STREN_intact + (1.0 - coher_new)*STREN_failed;
    pYSLOPE_new[idx]    = coher_new*YSLOPE_intact + (1.0 - coher_new)*YSLOPE_failed;

    // Copy the other parameters
    pBETA_new[idx]      = pBETA_old[idx];
    pCR_new[idx]        = pCR_old[idx];
    pT1_new[idx]        = pT1_old[idx];
    pT2_new[idx]        = pT2_old[idx];
  }
}

