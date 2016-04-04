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

YieldCond_MasonSand::YieldCond_MasonSand(Uintah::ProblemSpecP& ps,
                                         InternalVariableModel* intvar)
{
  d_intvar = intvar;

  // Nonlinear Drucker-Prager parameters
  ps->require("PEAKI1", d_yieldParam.PEAKI1);  // Shear Limit Surface Parameter
  ps->require("FSLOPE", d_yieldParam.FSLOPE);  // Shear Limit Surface Parameter
  ps->require("STREN",  d_yieldParam.STREN);   // Shear Limit Surface Parameter
  ps->require("YSLOPE", d_yieldParam.YSLOPE);  // Shear Limit Surface Parameter

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
  double coherence = 1.0;
  computeModelParameters(coherence);

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
  d_intvar = yc->d_intvar;

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

  VarLabel::destroy(pCoherenceLabel);
  VarLabel::destroy(pCoherenceLabel_preReloc);
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
  if (d_yieldParam.PEAKI1 <0.0 ) {
    warn << "PEAKI1 must be nonnegative. PEAKI1 = " << d_yieldParam.PEAKI1 << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_yieldParam.FSLOPE<0.0) {
    warn << "FSLOPE must be nonnegative. FSLOPE = " << d_yieldParam.FSLOPE << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_yieldParam.FSLOPE < d_yieldParam.YSLOPE) {
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
YieldCond_MasonSand::computeModelParameters(double coherence)
{
  double  FSLOPE = d_yieldParam.FSLOPE,            // Slope at I1=PEAKI1
          STREN  = d_yieldParam.STREN,             // Value of rootJ2 at I1=0
          YSLOPE = d_yieldParam.YSLOPE,            // High pressure slope
          PEAKI1 = coherence*d_yieldParam.PEAKI1;  // Value of I1 at strength=0

  std::vector<double> limitParameters = 
    computeModelParameters(PEAKI1, FSLOPE, STREN, YSLOPE);

  d_modelParam.a1 = limitParameters[0];
  d_modelParam.a2 = limitParameters[1];
  d_modelParam.a3 = limitParameters[2];
  d_modelParam.a4 = limitParameters[3];
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
    a3 = (STREN-YSLOPE*PEAKI1)*exp(-d_modelParam.a2*PEAKI1);
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
//     I1 = 3*(p - zeta)
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
    throw SCIRun::InternalError(out.str(), __FILE__, __LINE__);
  }

  // Get the particle specific internal variables from the model state
  // ** WARNING ** the sequence is hardcoded
  double PEAKI1 = state->yieldParams[0];
  double FSLOPE = state->yieldParams[1];
  double STREN  = state->yieldParams[2];
  double YSLOPE = state->yieldParams[3];
  //double BETA   = state->yieldParams[4];
  double CR     = state->yieldParams[5];
  //double T1     = state->yieldParams[6];
  //double T2     = state->yieldParams[7];
  //double Coher  = state->yieldParams[8];

  std::vector<double> limitParameters = 
    computeModelParameters(PEAKI1, FSLOPE, STREN, YSLOPE);
  double a1 = limitParameters[0];
  double a2 = limitParameters[1];
  double a3 = limitParameters[2];
  double a4 = limitParameters[3];

  // Get the local vars from the model state
  double zeta = state->zeta;
  double kappa = state->kappa;
  double capX = state->capX;

  // Initialize hasYielded to -1
  double hasYielded = -1.0;

  // Cauchy stress invariants: I1 = 3*p, J2 = q^2/3
  double I1 = state->I1;
  double sqrt_J2 = state->sqrt_J2;

  // Shift stress to evaluate yield condition where zeta is the isotropic back stress
  double I1_eff = (I1 - zeta);

  // --------------------------------------------------------------------
  // *** SHEAR LIMIT FUNCTION (Ff) ***
  // --------------------------------------------------------------------
  double Ff = a1 - a3*exp(a2*I1_eff) - a4*I1_eff;

  // --------------------------------------------------------------------
  // *** Branch Point (Kappa) ***
  // --------------------------------------------------------------------
  kappa = PEAKI1 - CR*(PEAKI1 - capX); // Branch Point

  // --------------------------------------------------------------------
  // *** COMPOSITE YIELD FUNCTION ***
  // --------------------------------------------------------------------
  // Evaluate Composite Yield Function F(I1) = Ff(I1)*fc(I1) in each region.
  // The elseif statements have nested if statements, which is not equivalent
  // to them having a single elseif(A&&B&&C)
  if (I1_eff < capX) {//---------------------------------------------------(I1<X)
    hasYielded = 1.0;
    return hasYielded;
  }

  // **Elliptical Cap Function: (fc)**
  // fc = sqrt(1.0 - Pow((Kappa-I1mZ)/(Kappa-X)),2.0);
  // faster version: fc2 = fc^2
  // **WARNING** p3 is the maximum achievable volumetric plastic strain in compresson
  // so if a value of 0 has been specified this indicates the user
  // wishes to run without porosity, and no cap function is used, i.e. fc=1
  if ((capX < I1_eff) && (I1_eff < kappa)) {// ---------------(X<I1<kappa)

    double kappaRatio = (kappa - I1_eff)/(kappa - capX);
    double fc2 = 1.0 - kappaRatio*kappaRatio;
    if (sqrt_J2*sqrt_J2 > Ff*Ff*fc2 ) {
      hasYielded = 1.0;
    }
  } else { // --------- X >= I1 or kappa <= I1

    if (I1_eff <= PEAKI1) { // ----- (kappa <= I1 <= PEAKI1)
      if (sqrt_J2 > Ff) {
        hasYielded = 1.0;
      }
    } else { // I1 > PEAKI1 
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
//     I1 = 3*(p - zeta)
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
    throw SCIRun::InternalError(out.str(), __FILE__, __LINE__);
  }

  // Get the particle specific internal variables from the model state
  // ** WARNING ** the sequence is hardcoded
  double PEAKI1 = state->yieldParams[0];
  double FSLOPE = state->yieldParams[1];
  double STREN  = state->yieldParams[2];
  double YSLOPE = state->yieldParams[3];
  //double BETA   = state->yieldParams[4];
  double CR     = state->yieldParams[5];
  //double T1     = state->yieldParams[6];
  //double T2     = state->yieldParams[7];
  //double Coher  = state->yieldParams[8];

  std::vector<double> limitParameters = 
    computeModelParameters(PEAKI1, FSLOPE, STREN, YSLOPE);
  double a1 = limitParameters[0];
  double a2 = limitParameters[1];
  double a3 = limitParameters[2];
  double a4 = limitParameters[3];

  // Get the plastic internal variables from the model state
  double zeta = state->zeta;
  double kappa = state->kappa;
  double capX = state->capX;

  // Cauchy stress invariants: I1 = 3*p, J2 = q^2/3
  double I1 = state->I1;

  // Shift stress to evaluate yield condition where zeta is the isotropic back stress
  double I1_eff = (I1 - zeta);

  // --------------------------------------------------------------------
  // *** SHEAR LIMIT FUNCTION (Ff) ***
  // --------------------------------------------------------------------
  double Ff = a1 - a3*exp(a2*I1_eff) - a4*I1_eff;

  // --------------------------------------------------------------------
  // *** Branch Point (Kappa) ***
  // --------------------------------------------------------------------
  kappa = PEAKI1 - CR*(PEAKI1 - capX); // Branch Point

  // --------------------------------------------------------------------
  // **Elliptical Cap Function: (fc)**
  // --------------------------------------------------------------------
  double kappaRatio = (kappa - I1_eff)/(kappa - capX);
  double Fc_sq = 1.0 - kappaRatio*kappaRatio;

  // --------------------------------------------------------------------
  // Derivatives
  // --------------------------------------------------------------------
  // term1 = 6*Ff*(a2*a3*exp(a2*I1) + a4)*Fc^2  
  double term1 = 6.0*Ff*Fc_sq*(a2*a3*exp(a2*I1_eff) + a4);
  // term2 = 6*Ff^2*(kappa - I1)/(kappa - X)^2
  double term2 = 6.0*Ff*Ff*kappaRatio/(kappa - capX);

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
    throw SCIRun::InternalError(out.str(), __FILE__, __LINE__);
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
                                                 const InternalVariableModel* intvar)
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
                                                          const InternalVariableModel* intvar)
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
    throw SCIRun::InternalError(out.str(), __FILE__, __LINE__);
  }

  // Compute effective trial stress
  double  I1eff_trial = state_trial->I1 - state_old->zeta;

  // Get the particle specific internal variables from the model state
  // ** WARNING ** the sequence is hardcoded
  double PEAKI1 = state_old->yieldParams[0];
  double coher  = state_old->yieldParams[8];

  // It may be better to use an interior point at the center of the yield surface, rather than at 
  // zeta, in particular when PEAKI1=0.  Picking the midpoint between PEAKI1 and X would be 
  // problematic when the user has specified some no porosity condition (e.g. p0=-1e99)
  double I1_interior = 0.0;
  double upperI1 = coher*PEAKI1;
  if (I1eff_trial < upperI1) {
    if (I1eff_trial > state_old->capX) { // Trial is above yield surface
      I1_interior = state_trial->I1;
    } else { // Trial is past X, use yield midpoint as interior point
      I1_interior = state_old->zeta + 0.5*(coher*PEAKI1 + state_old->capX);
    }
  } else { // I1_trial - zeta >= coher*I1_peak => Trial is past vertex
    double lTrial = sqrt(I1eff_trial*I1eff_trial + state_trial->sqrt_J2*state_trial->sqrt_J2);
    double lYield = 0.5*(coher*PEAKI1 - state_old->capX);
    I1_interior = state_old->zeta + upperI1 - std::min(lTrial, lYield);
  }
  
  return I1_interior;
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
    throw SCIRun::InternalError(out.str(), __FILE__, __LINE__);
  }

  // Save the input point
  point_type pt(px, py);

  // Get the yield surface points
  int num_points = 1000;
  std::vector<point_type> poly_points = getYieldSurfacePoints(state, num_points);

  // Find two yield surface segments that are closest to input point
  std::vector<point_type> segment_points = getClosestSegments(pt, poly_points);
  cpx = boost::geometry::get<0>(segment_points[1]);
  cpy = boost::geometry::get<1>(segment_points[1]);

  // Find the closest point
  point_type cpt = findClosestPoint(pt, segment_points);
  cpx = boost::geometry::get<0>(cpt);
  cpy = boost::geometry::get<1>(cpt);

  if (cpx == boost::numeric::bounds<double>::highest()) {
    return false;
  }
  return true;
}

/* Get the points on the yield surface */
std::vector<point_type> 
YieldCond_MasonSand::getYieldSurfacePoints(const ModelState_MasonSand* state,
                                           const int& num_points)
{

  // Get the particle specific internal variables from the model state
  // ** WARNING ** the sequence is hardcoded
  double PEAKI1 = state->yieldParams[0];
  double FSLOPE = state->yieldParams[1];
  double STREN  = state->yieldParams[2];
  double YSLOPE = state->yieldParams[3];
  double BETA   = state->yieldParams[4];
  double CR     = state->yieldParams[5];

  std::vector<double> limitParameters = 
    computeModelParameters(PEAKI1, FSLOPE, STREN, YSLOPE);
  double a1 = limitParameters[0];
  double a2 = limitParameters[1];
  double a3 = limitParameters[2];
  double a4 = limitParameters[3];

  // Get the plastic internal variables from the model state
  double zeta = state->zeta;
  double capX = state->capX;
  double kappa =  PEAKI1 - CR*(PEAKI1 - capX);

  // Get the bulk and shear moduli and compute sqrt(3/2 K/G)
  double sqrtKG = std::sqrt(1.5*state->bulkModulus/state->shearModulus);

   // Set up I1 values
  std::vector<double> I1_vec = linspace(0.99999*capX+zeta, 0.99999*PEAKI1+zeta, num_points);
  std::vector<double> J2_vec;
  for (auto I1 : I1_vec) {

    // Compute F_f
    double I1_minus_zeta = I1 - zeta;
    double Ff = a1 - a3*std::exp(a2*I1_minus_zeta) - a4*(I1_minus_zeta);
    double Ff_sq = Ff*Ff;

    // Compute Fc
    double Fc_sq = 1.0;
    if ((I1_minus_zeta < kappa) && (capX < I1_minus_zeta)) {
      double ratio = (kappa - I1_minus_zeta)/(kappa - capX);
      Fc_sq = 1.0 - ratio*ratio;
    }

    // Compute J2
    double J2 = Ff_sq*Fc_sq;
    J2_vec.push_back(J2);
  }

  // Convert I1 vs J2 to r' vs. z
  std::vector<double> z_vec;
  for (auto I1 : I1_vec) {
    z_vec.push_back(I1/std::sqrt(3.0));
  }

  std::vector<double> rprime_vec;
  for (auto J2 : J2_vec) {
    rprime_vec.push_back(BETA*std::sqrt(2.0*J2)*sqrtKG);
  }

  // Create a point_type vector
  std::vector<point_type> polyline;
  auto z_iter = z_vec.begin();
  auto r_iter = rprime_vec.begin();
  while (z_iter != z_vec.end() || r_iter != rprime_vec.end()) {
    double z = *z_iter;
    double r = *r_iter;
    //std::cout << "(" << z << "," << r << ")";
    polyline.push_back(point_type(z, r));
    ++z_iter; ++r_iter;
  }
  //std::cout << std::endl;

  auto rev_z_iter = z_vec.rbegin();
  auto rev_r_iter = rprime_vec.rbegin();
  while (rev_z_iter != z_vec.rend() || rev_r_iter != rprime_vec.rend()) {
    double z = *rev_z_iter;
    double r = *rev_r_iter;
    //std::cout << "(" << z << "," << -r << ")";
    polyline.push_back(point_type(z, -r));
    ++rev_z_iter; ++rev_r_iter;
  }
  polyline.push_back(point_type(*(z_vec.begin()), *(rprime_vec.begin())));
  //std::cout << std::endl;

  return polyline;
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

/* Find two yield surface segments that are closest to input point */
std::vector<point_type> 
YieldCond_MasonSand::getClosestSegments(const point_type& pt, 
                                        const std::vector<point_type> poly)
{
  point_type p_prev = *(poly.rbegin());
  auto iterNext = poly.begin();
  ++iterNext;
  point_type p_next = *iterNext;
  point_type min_p_prev, min_p, min_p_next;

  double min_dSq = boost::numeric::bounds<double>::highest();
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

  if (!(XP.size() > 0)) {
    std::cout << "No closest point" << std::endl;
    return point_type(boost::numeric::bounds<double>::highest(), 0.0);
  }
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


