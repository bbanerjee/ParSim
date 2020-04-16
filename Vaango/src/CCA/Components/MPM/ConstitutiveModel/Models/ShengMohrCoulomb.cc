/*
 * The MIT License
 *
 * Copyright (c) 1997-2019 Center for the Simulation of Accidental Fires and
 * Explosions (CSAFE), and  Scientific Computing and Imaging Institute (SCI),
 * University of Utah.
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
 *
 * License for the specific language governing rights and limitations under
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#include "ShengMohrCoulomb.h"

#include <Core/Exceptions/InvalidValue.h>
#include <Core/Util/DebugStream.h>

#include <cmath>

using namespace Uintah;

static DebugStream dbg("Sheng", false);
static DebugStream dbg_unloading("Sheng_unloading", false);

extern double USE_ERROR;
extern double STEP_MAX, STEP_MIN, ERROR_DEF, USE_ERROR_STEP, MIN_DIVISION_SIZE;
extern double ADDTOLYIELD, CHGEPSINTOL;
extern int    ALGORITHM_TYPE, USE_NICE_SCHEME;

/**
 * Constructors
 */
ShengMohrCoulomb::ShengMohrCoulomb()
{
  double G = 10000.0;                                   
  double K = 20000.0;                                   
  double cohesion = 0;
  double phi = 30;
  setModelParameters(G, K, cohesion, phi, phi);
  d_integration.setDefaults(d_yield);
}

ShengMohrCoulomb::ShengMohrCoulomb(double G, double K, double cohesion, double phi,
                                   double psi)
{
  setModelParameters(G, K, cohesion, phi, psi);
  d_integration.setDefaults(d_yield);
}

/**
 * Initialization methods
 */
void
ShengMohrCoulomb::setModelParameters(double G, double K, double cohesion, double phi,
                                     double psi)
{
  d_elastic.set(G, K);
  d_yield.set(cohesion, phi);
  d_potential.set(psi);
  if (phi != psi) {
    d_nonAssociated = true;
  } else {
    d_nonAssociated = false;
  }
}

void
ShengMohrCoulomb::setIntegrationParameters(int maxIterPegasus,
                                           double integrationTolerance, 
                                           double betaFactor, 
                                           double yieldLocTolerance,
                                           SolutionAlgorithm solutionAlgorithm,
                                           ToleranceMethod toleranceMethod,
                                           DriftCorrection driftCorrection) 
{
  d_int.d_maxIter        = maxIterPegasus;
  d_int.d_integrationTol = integrationTolerance; 
  d_int.d_yieldTol       = yieldLocTolerance;
  d_int.d_betaFactor     = betaFactor; 

  d_int.d_solutionAlgorithm = solutionAlgorithm;
  d_int.d_tolMethod         = toleranceMethod;
  d_int.d_driftCorrection   = driftCorrection;
}

/**
 * Integration
 */
void
ShengMohrCoulomb::integrate(const Vector6& strainIncrement, 
                            const StateMohrCoulomb& initialState)
{
  Vector7 purelyElasticStrain, purelyPlasticStrain;

  if (dbg.active()) {
    dbg << "Strain increment:" << strainIncrement << "\n";
  }

  // checking if on the YL (Returns true if on the YL or outside, in
  // the plastic part)
  bool onTheYieldLocus = checkYieldNormalized(initialState); 

  // Initial point copied to the final point, later
  // Final point updated with the elastic stress
  StateMohrCoulomb finalState;
  calcElastic(strainIncrement, initialState, finalState); 

  // checkYieldNormalized returns true  if in the plastic part of the YL. If on the YL
  // (within Tolerance INT_TOL returns false)
  bool elasticPlastic = checkYieldNormalized(finalState); 

  if (dbg.active()) {
    dbg << "p = " << finalState.meanStress() << " q = "<< finalState.shearStress()
        << "  Yield? " << std::boolalpha << purelyElastic << "\n";
  }

  bool unloading = false;
  if (elasticPlastic) {
    if (onTheYieldLocus) {

      if (dbg_unloading.active()) {
        // checking if there is any unloading part; in finalState purely elastic values; 
        // we care only for correct change of suction
        unloading = checkGradient(initialState, finalState); 
      }

      if (unloading) {

        if (dbg_unloading.active()) {
          dbg_unloading << "\n\n Elasto-Plastic unloading=" << unloading << "\n";
        }

        // find elastic part
        findIntersectionUnloading(strainIncrement, initialState,
                                  purelyElasticStrain,
                                  purelyPlasticStrain); 

        // calculating elastic part, updated in the same point.
        calcElastic(purelyElasticStrain, initialState, finalState); 

        // Update specific volume
        double spVol = computeNu(finalState.stress, finalState.state, finalState.suction());
        finalState.specificVolume(spVol);

        calculatePlastic(purelyPlasticStrain, finalState); 

      } else {

        calculatePlastic(strainIncrement, finalState);

      }

    } else {
      // not on the yield locus, finding intersection and subsequently
      // calculate elastic and plastic stress increment
      findIntersection(strainIncrement, initialState, 
                       purelyElasticStrain,
                       purelyPlasticStrain); 

      calcElastic(purelyElasticStrain, initialState, finalState);

      double spVol = computeNu(finalState.stress, finalState.state, finalState.suction());
      finalState.specificVolume(spVol);

      calculatePlastic(purelyPlasticStrain, &finalState);
    }
  }
  // for (int i=0; i<3; i++) cout<<strainIncrement[i]<<"  "<<"\n";
  // for (int i=0; i<3; i++) cout<<finalState.stress[i]<<"  "<<"\n";

  double spVol = computeNu(finalState.stress, finalState.state, finalState.suction());
  finalState.setSpecificVolume(spVol);
  finalState.Copy(initialState);
  // getchar();
}

/**
 * Check yield
 *  check of the standard yield surface
 *  Note: value of function is normalised by the mean stress+2*cohesion!
 */
bool
ShengMohrCoulomb::checkYieldNormalized(const StateMohrCoulomb& state) const
{
  double yieldFnValue = computeYieldNormalized(state.stress);
  if (yieldFnValue > (-d_int.d_yieldTol)) {
    return true;
  }
  return false;
}

double
ShengMohrCoulomb::computeYieldNormalized(const Vector6& stress) const
{
  double meanStress = firstInvariant(stress) / 3.0;
  double shearStress = std::sqrt(3.0 * secondDevInvariant(stress));
  double J3 = thirdDevInvariant(stress);
  double M = (6.0 * d_yield.d_sin_phi) / (3.0 - d_yield.d_sin_phi);
  if (shearStress > TINY) {
    double factor = -27.0 * J3 / (2.0 * shearStress * shearStress * shearStress);
    if (!std::isfinite(factor)) {
      factor = 1.0;
      std::cout << "**WARNING** Factor in checkYield Normalised not finite. set to 1\n"
                   "  alpha4 = " << d_yield.d_alpha4 << " J3 = " << J3
                   " shearStress = " << shearStress << "\n";
    } else if (factor > 1) {
      factor = 1;
    } else if (factor < -1) {
      factor = -1;
    } else {
      factor = 1 + d_yield.d_alpha4 - (1 - d_yield.d_alpha4) * factor;
    }
    double scaledAlpha = d_yield.d_alpha / std::pow(0.5 * factor, 0.25);
    if (scaledAlpha > -1.0 && scaledAlpha < 1.0) {
      M *= M * scaledAlpha;
    } 
  } else {
    M *= d_yield.d_alpha / std::pow(0.5 * (1 + d_yield.d_alpha4), 0.25);
  }

  double cohesion = d_yield.d_cohesion;
  double yieldFnValue = shearStress / M - 2.0 * d_yield.d_cohesion / M - meanStress;

  // normalisation
  yieldFnValue /= (std::abs(meanStress) + 2.0 * d_yield.d_cohesion);

  if (dbg.active()) {
    std::cout << "Check Yield: Mean Stress = " << meanStress 
              << " Shear Stress=" << shearStress 
              << " cohesion = "<< d_yield.d_cohesion
              << " M = "<< M
              << " Yield Function = " << yieldFnValue << "\n";
  }
  return yieldFnValue;
}

/**
 * Calc elastic state 
 */
void
ShengMohrCoulomb::calcElastic(const Vector7& strainInc, 
                              const StateMohrCoulomb& initialState,
                              StateMohrCoulomb& finalState) const
{
  finalState = initialState;

  Vector6 stressInc = calcStressIncElast(initialState.specificVolume(), 
                                         initialState.stress, 
                                         initialState.strain, 
                                         strainInc); 

  std::cout << "Stress increment is: " << stressInc << "\n";

  Vector6 plasticStrainInc = Vector6::Zero();
  double p0StarInc = 0.0;
  finalState.update(plasticStrainInc, strainInc, stressInc, p0StarInc);
}

Vector6
ShengMohrCoulomb::calcStressIncElast(double nu0, 
                                     const Vector6& s0, 
                                     const Vector7& eps0,
                                     const Vector7& deps) 
{
  double K = d_elastic.d_K;
  double G = d_elastic.d_G;

  double K43G = K + 4.0 * G / 3.0;
  double K23G = K - 2.0 * G / 3.0;

  Vector6 ds;
  ds(0) =  K43G * deps(0) + K23G * (deps(1) + deps(2));
  ds(1) =  K43G * deps(1) + K23G * (deps(0) + deps(2));
  ds(2) =  K43G * deps(2) + K23G * (deps(0) + deps(1));
  ds(3) = 2 * G * deps(3);
  ds(4) = 2 * G * deps(4);
  ds(5) = 2 * G * deps(5);

  return ds;
}

/*
 *  Check the gradient between initial and final state
 */
bool
ShengMohrCoulomb::checkGradient(const StateMohrCoulomb& initialState, 
                                const StateMohrCoulomb& finalState) const
{
  Vector7 strainInc = finalState->strain - initialState->strain;

  double max = 0.0;
  for (int i = 0; i < 6; i++) {
    if (std::abs(strainInc(i)) > max) {
      max = std::abs(strainInc(i));
    }
  }
  for (int i = 0; i < 6; i++) {
    strainInc(i) /= (max * 10E-10);
  }
  strainInc(6) = finalState->suction() - initialState->suction();

  // The normalisation is important to catch the unloading due to shear stress
  // 12 13 23. As the q is always positive, large unloading - loading values 
  // of 12 13 or 23 component lead to larger q, the shear stress change is 
  // taken as positive and there is no unloading. This fix, though not the 
  // most elegant, should work.
  Matrix67 tangentElastic = calculateElasticTangentMatrix(initialState);

  Vector6 stressInc = tangentElastic * strainInc;

  if (dbg.active()) {
    dbg << "Strain inc = " << strainInc << "\n";
    dbg << "TangentElastic: "<< tangentElastic << "\n";
    dbg << "Stress inc = " << stressInc << "\n";
    dbg << "Suction Increment = " << strainInc(7) << "\n";
  }

  Vector6 df; // not initialised; will contain values of derivatives of the
              // yield locus function
  double cosinus = findGradient(initialState->state, initialState->stress, 
                                stressInc, df, 0, 0);

  if (cosinus > -d_int.d_yieldTol) {
    return false;
  } 

  if (dbg.active()) {
    dbg << "Unloading occuring... cosinus = " << cosinus << "\n";
  }

  return true; //(negative cosinus means unloading occurs)
}

/*
 *  Compute elastic tangent stiffness
 */
Matrix67
ShengMohrCoulomb::calculateElasticTangentMatrix(const StateMohrCoulomb& state) const
{
  double K = d_elastic.d_K;
  double G = d_elastic.d_G;
  double K43G = K + 4.0 * G / 3.0;
  double K23G = K - 2.0 * G / 3.0;

  Matrix67 DEP = Matrix67::Zeros();
  DEP(0, 0) =  K43G;
  DEP(0, 1) =  K23G;
  DEP(0, 2) =  K23G; 
  DEP(1, 0) =  K23G;
  DEP(1, 1) =  K43G;
  DEP(1, 2) =  K23G; 
  DEP(2, 0) =  K23G;
  DEP(2, 1) =  K23G;
  DEP(2, 2) =  K43G;
  DEP(3, 3) =  2.0 * G;
  DEP(4, 4) =  2.0 * G;
  DEP(5, 5) =  2.0 * G;
  return DEP;
}

/**
  This procedure finds a gradient to the yield locus at given point.
  It is used later to determine whether we have partly elastic unloading or only
  plastic step
  The gradient may be found only when we are on "main" plastic yield locus;
  Otherwise results may be erroneous.

  Parameters: state[], s[] (stress), ds[] (stress increment), suction,

  Gradient = gradient[1], gradient[2] etc - there the result will be stored

  First the C(s) part must be calculated, then the gradient.

  returns cosinus of the angle...
 */
double
ShengMohrCoulomb::findGradient(const Vector3& state, const Vector6& s, const Vector6& ds, 
                               Vector6s& df_dSigma, double suction, double dsuction) const
{
  double I1 = firstInvariant(s);
  double I2 = secondInvariant(s);
  double I3 = thirdInvariant(s);
  double J2 = secondDevInvariant(s);
  double J3 = thirdDevInvariant(s);

  double shearStress = sqrt(3.0 * J2);
  double factor = -27.0 * J3 / (2.0 * shearStress * shearStress * shearStress);
  if (!std::isfinite(factor)) {
    factor = 1.0;
    std::cout << "**WARNING** Factor in findGradient not finite. set to 1\n"
                 "  alpha4 = " << alpha4 << " J3 = " << J3
                 " shearStress = " << shearStress << "\n";
  } else if (factor > 1) {
    factor = 1;
  } else if (factor < -1) {
    factor = -1;
  } 
  factor = 1 + alpha4 - (1 - alpha4) * factor;

  double factor025 = std::pow(factor, 0.25);
  double factor075 = std::pow(factor, -0.75);
  double M = (3 - d_yield.d_sin_phi) / (6 * d_yield.d_alpha * d_yield.d_sin_phi);

  Vector6 dI1_dSig = Vector6::Zeros();
  dI1_dSig(0) =  1.0;
  dI1_dSig(1) =  1.0;
  dI1_dSig(2) =  1.0; //{1,1,1,0,0,0}

  Vector6 dI2_dSig = Vector6::Zeros();
  dI2_dSig(0) =  s(1) + s(2);
  dI2_dSig(1) =  s(0) + s(2);
  dI2_dSig(2) =  s(0) + s(1);
  dI2_dSig(3) =  -2 * s(3);
  dI2_dSig(4) =  -2 * s(4);
  dI2_dSig(5) =  -2 * s(5);

  if (dbg.active()) {
    dbg << setprecision(15) << dI2_dSig << "\n";
  }

  Vector6 dI3_dSig = Vector6::Zeros(); 
  dI3_dSig(0) =  s(1) * s(2) - s(5) * s(5);
  dI3_dSig(1) =  s(0) * s(2) - s(4) * s(4);
  dI3_dSig(2) =  s(0) * s(1) - s(3) * s(3);
  dI3_dSig(3) =  2 * s(5) * s(4) - 2 * s(2) * s(3);
  dI3_dSig(4) =  2 * s(3) * s(5) - 2 * s(1) * s(4);
  dI3_dSig(5) =  2 * s(3) * s(4) - 2 * s(0) * s(5);

  if (dbg.active()) {
    dbg << setprecision(15) << dI3_dSig << "\n";
  }

  Vector6 dJ2_dSig = Vector6::Zeros();
  dJ2_dSig(0) =  (2 * s(0) - s(1) - s(2)) / 3.0;
  dJ2_dSig(1) =  (2 * s(1) - s(0) - s(2)) / 3.0;
  dJ2_dSig(2) =  (2 * s(2) - s(0) - s(1)) / 3.0;
  dJ2_dSig(3) =  2 * s(3);
  dJ2_dSig(4) =  2 * s(4);
  dJ2_dSig(5) =  2 * s(5);

  Vector6 dq_dSig = dJ2_dSig;
  dq_dSig *= std::sqrt(3.0) * 0.5 / std::sqrt(J2);

  Vector6 dJ3_dSig = dI1_dSig;
  dJ3_dSig *= (2.0 / 9.0 * I1 * I1 - I2 / 3.0);

  if (dbg.active()) {
    dbg << setprecision(15) << dJ3_dSig << "\n";
  }

  dJ3_dSig += (dI2_dSig * (-I1 / 3.0) + dI3_dSig);

  if (dbg.active()) {
    dbg << setprecision(15) << dJ3_dSig << "\n";
  }

  df_dSigma = dJ2_dSig;
  df_dSigma *= (0.21022410391343 * M * factor025 / shearStress);

  for (int i = 0; i < 6; i++) {
    if (!std::isfinite(df_dSigma(i))) {
      df_dSigma(i) =  0.0;
    }
  }

  if (dbg.active()) {
    dbg << "dF/dJ2 part\n";
    dbg << setprecision(15) << df_dSigma << "\n";
  }

  Vector6 temp = dJ3_dSig;
  temp *= (-1.419012700740643 * (d_yield.alpha4 - 1) * M * sqrt(J2) * factor075 /
             (std::sqrt(3.0) * shearStress * shearStress * shearStress) +
           0.94600846716043 * (d_yield.alpha4 - 1) * d_yield.cohesion * M * factor075 /
             (shearStress * shearStress * shearStress));

  if (dbg.active()) {
    dbg << "dF/dJ3 part\n";
    dbg << setprecision(15) << temp << "\n";
  }

  df_dSigma += temp;;

  for (int i = 0; i < 6; i++) {
    if (!std::isfinite(df_dSigma(i))) {
      df_dSigma(i) =  0.0;
    }
  }

  temp = dq_dSig;
  temp *= (-2.838025401481287 * (d_yield.alpha4 - 1) * d_yield.cohesion * M * J3 * factor075 /
             (shearStress * shearStress * shearStress * shearStress) +
           4.257038102221929 * (d_yield.alpha4 - 1) * M * std::sqrt(J2) * J3 * factor075 /
             (std::sqrt(3.0) * shearStress * shearStress * shearStress * shearStress));

  if (dbg.active()) {
    dbg << "dF/dq part\n";
    dbg << setprecision(15) << temp << "\n";
  }

  df_dSigma += temp;

  df_dsigma += dI1_dSig * (-1.0 / 3.0);

  for (int i = 0; i < 6; i++) {
    if (!std::isfinite(df_dSigma(i))) {
      df_dSigma(i) =  0.0;
    }
  }

  if (dbg.active()) {
    dbg << "df_dSigma = " << df_dSigma << "\n";
  }

  // calculate length - gradient - total length of vector ==1#
  double dFlength = 0, dslength = 0;
  for (int i = 0; i < 6; i++) {
    dFlength += df_dSigma(i) * df_dSigma(i); // calculated length
    dslength += ds(i) * ds(i);
  }
  dslength = std::sqrt(dslength);
  dFlength = std::sqrt(dFlength);

  // calculate cosinus of the theta angle...
  double cosin = 0;
  for (int i = 0; i < 6; i++) {
    cosin += df_dsigma(i) * ds(i) / (dslength * dFlength); 
  }

  if (dbg.active()) {
    if (cosin < -d_int.d_yieldTol) {
      dbg << "Check if no problem in the findGradient\n";
      dbg << "cosin = " << cosin << "\n";
    }
  }

  return cosin; 
}

/**
 Find intersection with yield surface during unloading
 */
void
ShengMohrCoulomb::findIntersectionUnloading(const Vector7& strainIncrement,
                                            const StateMohrCoulomb& initialState,
                                            Vector7& elasticStrainInc,
                                            Vector7& plasticStrainInc) const
{
  double alpha = findYieldAlpha(initialState.state, initialState.stress,
                                initialState.strain, strainIncrement);
  
  elasticStrainInc = strainIncrement * alpha;
  plasticStrainInc = strainIncrement - elasticStrainInc;
}

void
ShengMohrCoulomb::findIntersection(const Vector7& strainIncrement,
                                   const StateMohrCoulomb& initialState,
                                   Vector7& elasticStrainInc,
                                   Vector7& plasticStrainInc)
{
  double alpha = findYieldModified(initialState.state, initialState.stress,
                                   initialState.strain, strainIncrement);
  elasticStrainInc = strainIncrement[i] * alpha;
  plasticStrainInc = strainIncrement - elasticStrainInc;
}

/**
  Main purpose of this procedure is to find the intersection between yield
  surface and the stress vector. This is a special version of the procedure 
  used in case of elastic-plastic unloading - when we have plastic state
  at the beginning of the step and then elastic unloading and plastic again. We
  need to know at what value of epsilon (strain) we are on the surface - where 
  plastic yielding will begin.

  Parameters: state vector, initial stress s0, initial strain eps0, strain
  increment deps, a- to which value of alfa will be written

  Changed - a

  Algorithm - this is a Pegasus algorithm; we are looking for alfa such, that
  value of yield function F==0.
  This alfa is found and written to a.
*/
double
ShengMohrCoulomb::findYieldAlpha(const Vector3& state, 
                                 const Vector6& s0, 
                                 const Vector7& eps0,
                                 const Vector7& deps) const 
{
  double f0 = computeYieldNormalized(s0);
  double yieldTol_old = d_int.d_yieldTol;
  double delta = 0.0;
  if (f0 < 0) {
    delta = 0;
    d_int.d_yieldTol = -0.9 * f0;
  } else {
    if (dbg.active()) {
      dbg << "findYieldAlpha : case not coded for yet\n";
    }
    delta = (f0 + d_yield.d_yieldTol) / 2;
    d_yield.d_yieldTol = 0.9 * (d_yield.d_yieldTol - f0) / 2;
  }

  if (dbg.active()) {
    dbg << "Procedure Find Yield Yield ... s = " << s0 << "\n";
    dbg << "Procedure Find Yield Yield ... dEps = " << dep << "\n";
    dbg << "Initial normalised Yield Locus... f0 =" << f0 << "\n";
    dbg << "Procedure Find Yield Yield ... delta =" << delta << "\n";
    dbg << "Call of this procedure is generally rare and most likely improper\n";
    dbg << "Parameters are set to:\n";
    dbg << "Maximum number of iteration d_maxIter: "<< d_int.d_maxIter << "\n";;
    dbg << "Value of d_alfaCheck - when the additional iter. is to "
           "enter: " << d_int.d_alfaCheck << "\n";
    dbg << "Value of change of alfa in the step: "<< d_int.d_alfaChange << "\n";;
    dbg << "alfa old/alfa ratio: " << d_int.d_alfaRatio << "\n";
  }

  double alfa0 = 0;
  double alfa1 = 1;
  Vector7 epsIni = deps * alfa0;
  Vector7 epsFin = deps * alfa1;

  Vector6 sIni = s0;
  Vector6 sFin = calcStressIncElast(state(2), s0, eps0, epsFin);
  sFin += s0;

  f0 = computeYieldNormalized(sIni);
  f0 -= delta;

  double f1 = computeYieldNormalized(sFin);
  f1 -= delta;

  if (dbg.active()) {
    dbg << "Procedure findYieldAlpha. Initial values of f0 = " << f0 
        << " f1 = " << f1 << "\n";
    dbg << "f0 should lie on the yield locus within tolerance: "
        << yieldTol_old << "\n";
  }

  double alphaYield = 0.0;
  double alfa = 0.5;
  double alfa_old = 0;
  for (int iter = 0; iter < d_maxIter; iter++) {
    bool problems = false;
    alfa_old = alfa;

    alfa = f0 / (f0 - f1);
    Vector7 epsAlfa = deps * (alfa0 + alfa * (alfa1 - alfa0));

    // calculate stress increment for current alfa
    Vector6 sAlfa = calcStressIncElast(state(2), s0, eps0, epsAlfa); 

    // update stress
    sAlfa += s0;

    // calculate yield function for current alfa
    double fAlfa = computeYieldNormalized(sAlfa);
    fAlfa -= delta;

    if (dbg.active()) {
      dbg << "In iteration " << iter << " alfa = " << alfa << " and f = "<< fAlfa << "\n";
    }

    // if fAlfa within tolerance, we have the solution
    if (std::abs(fAlfa) < d_int.d_yieldTol) {

      alphaYield = alfa0 + alfa * (alfa1 - alfa0); 

      if (dbg.active()) {
        dbg << "Solution in findYieldAlpha procedure was found after "<< iter 
            << " iterations." << "\n";
        dbg << "Value of the yield function is equal to: "<< fAlfa << "\n";
        dbg << "Modified value of tolerance is: " << d_int.d_yieldTol << "\n";
        dbg << "Value of delta is: " << delta << "\n";
        dbg << "Alfa is equal to: " << alphaYield << "\n";
      }

      if (iter > 50) {
        std::cout << "**WARNING** Large number of iterations!!! Solution is however correct... "
             << "\n";
      }

      d_yieldTol = yieldTol_old;

      if (dbg.active()) {
        dbg << "Yield tolerance is set back to: " << d_int.d_yieldTol << "\n";
      }

      return alphaYield;
    }

    if (fAlfa > 0) {

      // if fAlfa > 0 - we are yielding - max alfa set to current alfa

      if (((1.0 - alfa) < d_int.d_alfaCheck) && 
          ((1.0 - alfa_old) / (1.0 - alfa)) < d_int.d_alfaRatio) {
        problems = true;
      }

      alfa1 = alfa0 + alfa * (alfa1 - alfa0);
      f1 = fAlfa; 

      if (problems) {
        if (dbg.active()) {
          dbg << "Problematic iteration entered !!!"<<"\n";
        }

        alfa = alfa1 - d_int.d_alfaChange * (alfa1 - alfa0);
        epsAlfa = deps * (alfa0 + alpfa * (alfa1 - alfa0));

        sAlpha = calcStressIncElast(state(2), s0, eps0, epsAlfa);
        sAlpha += s0;

        fAlfa = computeYieldNormalized(sAlfa);
        fAlfa -= delta;

        if (fAlfa > 0) {
          // if fAlfa > 0 - we are yielding - max alfa set to current alfa
          alfa1 = alfa0 + alfa * (alfa1 - alfa0);
          f1 = fAlfa; 
        } else {
          // if fAlfa < 0 - we are elastic - minimum alfa is set to current alfa
          alfa0 += alfa * (alfa1 - alfa0);
          f0 = fAlfa; 
        }
      }

    } else {

      // if fAlfa < 0 - we are elastic - minimum alfa is set to current alfa

      if ((alfa < d_alfaCheck) && (alfa_old / alfa) < d_alfaRatio) {
        problems = true;
      }
      alfa0 = alfa0 + alfa * (alfa1 - alfa0);
      f0 = fAlfa; 

      if (problems) {
        if (dbg.active()) {
          dbg << "Problematic iteration entered !!!"<<"\n";
        }

        alfa = alfa0 + d_int.d_alfaChange * (alfa1 - alfa0);
        epsAlfa = deps * (alfa0 + alpfa * (alfa1 - alfa0));

        sAlpha = calcStressIncElast(state(2), s0, eps0, epsAlfa);
        sAlpha += s0;

        fAlfa = computeYieldNormalized(sAlfa);
        fAlfa -= delta;

        if (fAlfa > 0) {
          // if fAlfa > 0 - we are yielding - max alfa set to current alfa
          alfa1 = alfa0 + alfa * (alfa1 - alfa0);
          f1 = fAlfa; 
        } else {
          // if fAlfa < 0 - we are elastic - minimum alfa is set to current alfa
          alfa0 += alfa * (alfa1 - alfa0);
          f0 = fAlfa; 
        }
      }
    }
  }

  d_yieldTol = yieldTol_old;

  // if we are here, we must have perforemed to many iterations...
  std::ostringstream err;
  err << "**ERROR** in procedure findYieldAlpha" << "\n";
  err << "  After " << d_maxIter << " iterations crossing point not found" << "\n";
  err << "  This is likely to cause incorrect results... Results obtained "
         "  should not be taken too seriously..." << "\n";
  throw InvalidValue(err.str(), __FILE__, __LINE__);

  return -999.0;
}

/**
  Main purpose of this procedure is to find the intersection between yield
  surface and the stress vector. We need to know at what value of epsilon 
  (strain) we are on the surface - where plastic yielding will begin.

  Parameters: state vector, initial stress s0, initial strain eps0, strain
  increment deps, a- to which value of alfa will be written

  Changed - a

  Algorithm - this is a Pegasus algorithm; we are looking for alfa such, that
  value of yield function F==0.
  This alfa is found and written to a.
*/

double
ShengMohrCoulomb::findYieldModified(const Vector3& state, 
                                    const Vector6& s0, 
                                    const Vector7& eps0,
                                    const Vector7& deps) const
{

  if (dbg.active()) {
    dbg << "Parameters are set to:\n";
    dbg << "Maximum number of iteration d_maxIter: "<< d_int.d_maxIter << "\n";;
    dbg << "Value of d_alfaCheck - when the additional iter. is to "
           "enter: " << d_int.d_alfaCheck << "\n";
    dbg << "Value of change of alfa in the step: "<< d_int.d_alfaChange << "\n";;
    dbg << "alfa old/alfa ratio: " << d_int.d_alfaRatio << "\n";
  }

  double alfa0 = 0;
  double alfa1 = 1;
  Vector7 epsIni = deps * alfa0;
  Vector7 epsFin = deps * alfa1;

  Vector6 sIni = s0;
  Vector6 sFin = calcStressIncElast(state(2), s0, eps0, epsFin);
  sFin += s0;

  double f0 = computeYieldNormalized(s0);
  double f1 = computeYieldNormalized(sFin);

  if (dbg.active()) {
    dbg << "Procedure findYieldModified. Initial values of f0 = " << f0 
        << " f1 = " << f1 << "\n";
    dbg << "Value of f0 should be negative, and value of f1 should be positive.\n";
    dbg << "Values should be larger than tolerance for yield: "
        << d_int.d_yieldTol << "\n";
  }

  double alfa = 0.5;
  double alfa_old = 0;

  for (int iter = 0; iter < d_maxIter; iter++) {

    bool problems = false;
    alfa_old = alfa;
    alfa = f0 / (f0 - f1);

    Vector7 epsAlfa = deps * (alfa0 + alfa * (alfa1 - alfa0));
    Vector6 sAlfa = calcStressIncElast(state(2), s0, eps0, epsAlfa);
    sAlfa += s0;

    double fAlfa = computeYieldNormalized(sAlfa);

    if (dbg.active()) {
      dbg << "In iteration " << iter << " alfa = " << alfa << " and f = " << fAlfa;
      dbg << " alfa0 = " << alfa0 << "\n";
    }

    // if the difference is below numerical the accuracy, we have the solution
    if ((alfa1 - alfa0) < TINY) {
      fAlfa = 0.0; 
    }

    if (std::abs(fAlfa) < d_yieldTol) {

      // if fAlfa within tolerance, we have the solution
      double alphaYield = alfa0 + alfa * (alfa1 - alfa0);

      if (iter > 50) {
        std::cout << "**WARNING** Large number of iterations in findYieldModified procedure!!! "
                     "Solution is however correct..." << "\n";
      }
      if (dbg.active()) {
        dbg << "Solution in findYieldModified procedure was found after " 
            << iter << " iterations."<<"\n";
        dbg << "solution is: " << alphaYield << "\n";
      }
      return alphaYield;
    }

    if (fAlfa > 0) {

      // if fAlfa >0 - we are yielding - max alfa set to current alfa
      fAlfa = computeYieldNormalized(sAlfa);
      if (((1 - alfa) < d_int.d_alfaCheck) && 
          ((1 - alfa_old) / (1 - alfa)) < d_int.d_alfaRatio) {
        problems = true;
      }

      alfa1 = alfa0 + alfa * (alfa1 - alfa0);
      f1 = fAlfa; 

      if (problems) {

        if (dbg.active()) {
          dbg << "Problematic iteration entered !!!"<<"\n";
        }
        alfa = alfa1 - d_int.d_alfaChange * (alfa1 - alfa0);
        epsAlfa = deps * (alfa0 + alfa * (alfa1 - alfa0));

        sAlfa = calcStressIncElast(state(2), s0, eps0, epsAlfa);
        sAlfa += s0;

        falfa = computeYieldNormalized(sAlfa);

        if (fAlfa > 0) {
          // if fAlfa > 0 - we are yielding - max alfa set to current alfa
          alfa1 = alfa0 + alfa * (alfa1 - alfa0);
          f1 = fAlfa; 
        } else {
          // if fAlfa < 0 - we are elastic - minimum alfa is set to current alfa
          alfa0 = alfa0 + alfa * (alfa1 - alfa0);
          f0 = fAlfa; 
        }
      }

    } else {

      // if fAlfa <0 - we are elastic - minimum alfa is set to current alfa
      fAlfa = computeYieldNormalized(sAlfa);
      if ((alfa < d_int.d_alfaCheck) && (alfa_old / alfa) < d_int.d_alfaRatio) {
        problems = true;
      }
      alfa0 = alfa0 + alfa * (alfa1 - alfa0);
      f0 = fAlfa; 

      if (problems) {

        if (dbg.active()) {
          dbg << "Problematic iteration entered !!!"<<"\n";
        }

        alfa = alfa0 + d_alfaChange * (alfa1 - alfa0);
        epsAlfa = deps * (alfa0 + alfa * (alfa1 - alfa0));

        sAlfa = calcStressIncElast(state(2), s0, eps0, epsAlfa);
        sAlfa += s0;

        fAlfa = computeYieldNormalized(sAlfa);
        if (fAlfa > 0) {
          // if fAlfa > 0 - we are yielding - max alfa set to current alfa
          alfa1 = alfa0 + alfa * (alfa1 - alfa0);
          f1 = fAlfa; 
        } else {
          // if fAlfa < 0 - we are elastic - minimum alfa is set to current alfa
          alfa0 = alfa0 + alfa * (alfa1 - alfa0);
          f0 = fAlfa; 
        }
      }
    }
  }

  // if we are here, we must have perforemed to many iterations...
  std::ostringstream err;
  err << "Error in procedure findYieldModified" << "\n";
  err << "After " << d_maxIter << " iterations crossing point not found" << "\n";
  err << "alphamin = " << alfa0 << " alphamax = " << alfa1
       << " dalpha = " << alfa1 - alfa0;
  err << "Yield Function value Min=" << f0 << " Max=" << f1 << "\n";
  err << "Stress:" << s0 << "\n";
  err << "Strain:" << deps << "\n";
  err << "G: " << G << " K: " << K << " cohesion: " << cohesion << " phi: " << phi
       << "\n";
  err << "This is likely to cause incorrect results... Results obtained "
          "should not be taken too seriously..."
       << "\n";
  throw InvalidValue(err.str(), __FILE__, __LINE__);

  return -999.0;
}

double
ShengMohrCoulomb::computeNu(const Matrix6& s, const Matrix3& state, double suction)
{
  // does nothing for SMC
  return 1.0;
}

void
ShengMohrCoulomb::findElStrGrad(double nu0, double* s0, double* eps0,
                                double* deps, double* ds)
{
  // to be easily used
  double K43G = K + 4.0 * G / 3.0;
  double K23G = K - 2.0 * G / 3.0;

  ds[0] = K43G * deps[0] + K23G * (deps[1] + deps[2]);
  ds[1] = K43G * deps[1] + K23G * (deps[0] + deps[2]);
  ds[2] = K43G * deps[2] + K23G * (deps[0] + deps[1]);
  ds[3] = 2 * G * deps[3];
  ds[4] = 2 * G * deps[4];
  ds[5] = 2 * G * deps[5];

  // for (int i=0; i<6; i++) cout<<"Delta stress i="<<ds[i]<<"\n";
}

bool
ShengMohrCoulomb::computeYieldNormalized(StateMohrCoulomb* point)
{
  // check of the standard yield surface
  // Note: value of function is not normalised by the mean stress

  double meanStress = Point->getmeanStress();
  double shearStress = Point->shearStress();
  double J3 = Point->thirdDevInvariant();
  double M = (6 * sin_phi) / (3 - sin_phi);
  if (shearStress > TINY) {
    double Factor = -27 * J3 / (2 * shearStress * shearStress * shearStress);
    if (Factor > 1)
      Factor = 1;
    if (Factor < -1)
      Factor = -1;
    if (!finite(Factor)) {
      Factor = 1.0;
      cout << "Factor in checkYield not finite. set to 1" << "\n";
      cout << "alpha4=" << alpha4 << " J3=" << J3
           << " shearStress=" << shearStress << "\n";
    }
    Factor = 1 + alpha4 - (1 - alpha4) * Factor;

    if ((alpha / pow(0.5 * Factor, 0.25) < 1.0) &&
        (alpha / pow(0.5 * Factor, 0.25) > -1.0))
      M = M * alpha / pow(0.5 * Factor, 0.25);
  } else
    M = M * alpha / pow(0.5 * (1 + alpha4), 0.25);

  double Value = shearStress / M - 2 * cohesion / M - meanStress;
  /*
      cout<<"check Yield: Mean Stress="<<meanStress;
          cout<<" Shear Stress="<<shearStress;
          cout<<" cohesion="<<cohesion;
          cout<<" M="<<M;
          cout<<" Yield Function="<<Value<<"\n";
  */
  if (Value > (-d_yieldTol))
    return true;
  else
    return false;
}

double
ShengMohrCoulomb::computeYieldFunction(StateMohrCoulomb* point)
{
  // Note: value of function is not normalised by the mean stress+2*cohesion!

  double meanStress = Point->getmeanStress();
  double shearStress = Point->shearStress();
  double J3 = Point->thirdDevInvariant();
  double M = (6 * sin_phi) / (3 - sin_phi);
  if (shearStress > TINY) {
    double Factor = -27 * J3 / (2 * shearStress * shearStress * shearStress);
    if (Factor > 1)
      Factor = 1;
    if (Factor < -1)
      Factor = -1;
    if (!finite(Factor)) {
      Factor = 1.0;
      cout << "Factor in computeYieldFunction [parameters *Point] is not "
              "finite. set to 1"
           << "\n";
      cout << "alpha4=" << alpha4 << " J3=" << J3
           << " shearStress=" << shearStress << "\n";
    }
    Factor = 1 + alpha4 - (1 - alpha4) * Factor;
    if ((alpha / pow(0.5 * Factor, 0.25) < 1.0) &&
        (alpha / pow(0.5 * Factor, 0.25) > -1.0))
      M = M * alpha / pow(0.5 * Factor, 0.25);
  } else
    M = M * alpha / pow(0.5 * (1 + alpha4), 0.25);

  double Value = shearStress / M - 2 * cohesion / M - meanStress;

  // normalisation
  // Value=Value/(meanStress+2*cohesion);

  return Value;
}

double
ShengMohrCoulomb::computeYieldFunctionNN(StateMohrCoulomb* point)
{
  // Note: value of function is normalised by the mean stress+2*cohesion!

  double meanStress = Point->getmeanStress();
  double shearStress = Point->shearStress();
  double J3 = Point->thirdDevInvariant();
  double M = (6 * sin_phi) / (3 - sin_phi);
  if (shearStress > TINY) {
    double Factor = -27 * J3 / (2 * shearStress * shearStress * shearStress);
    if (Factor > 1)
      Factor = 1;
    if (Factor < -1)
      Factor = -1;
    if (!finite(Factor)) {
      Factor = 1.0;
      cout << "Factor in computeYieldFunctionNN [parameters *Point] is not "
              "finite. set to 1"
           << "\n";
      cout << "alpha4=" << alpha4 << " J3=" << J3
           << " shearStress=" << shearStress << "\n";
    }
    Factor = 1 + alpha4 - (1 - alpha4) * Factor;
    if ((alpha / pow(0.5 * Factor, 0.25) < 1.0) &&
        (alpha / pow(0.5 * Factor, 0.25) > -1.0))
      M = M * alpha / pow(0.5 * Factor, 0.25);
  } else
    M = M * alpha / pow(0.5 * (1 + alpha4), 0.25);

  double Value = shearStress / M - 2 * cohesion / M - meanStress;

  return Value;
}


int
ShengMohrCoulomb::calcPlasticFaster(StateMohrCoulomb point, double* epStrain,
                                    BBMMatrix* dSigma, double* plasticStrain,
                                    double* dpZeroStar, double fValue,
                                    double* dS, double* dLambda)
{
  calcPlastic(Point, epStrain, dSigma, plasticStrain, dpZeroStar, fValue, dS,
              dLambda);
  return 0;
  // so we had the increase of strenghtening parameter...
  // that's all?
}

int
ShengMohrCoulomb::calcPlasticPQ(StateMohrCoulomb point, double* epStrain,
                                BBMMatrix* dSigma, double* plasticStrain,
                                double* dpZeroStar, double fValue, double* dS,
                                double* dLambda)
{
  calcPlastic(Point, epStrain, dSigma, plasticStrain, dpZeroStar, fValue, dS,
              dLambda);
  return 0;
  // so we had the increase of strenghtening parameter...
  // that's all?
}

int
ShengMohrCoulomb::calcPlastic(StateMohrCoulomb point, double* epStrain,
                              BBMMatrix* dSigma, double* plasticStrain,
                              double* dpZeroStar, double fValue, double* dS,
                              double* dLambda)
{
  // at the beginning we need to calculate the derivatives a, b, c, d, g, p...)
  BBMMatrix A(6, 1); // dF/dsigma
  BBMMatrix DENOMINATOR(1, 1);

  BBMMatrix GG(6, 1);  // dG/dsigma, in case of associated flow
                       // (nonAssociated==false) same as A
  BBMMatrix MM(6, 1);  // will be vector(1,1,1,0,0,0)T
  BBMMatrix DEL(6, 6); // D elastic matrix...
  BBMMatrix DEPS(6, 1);

  for (int i = 1; i < 6; i++)
    DEPS.PutElement(i, 1, epStrain[i - 1]); // increase of epsilon copied
  // we do not have g, as a == g in associated flow rule.
  if (!Point.checkIfFinite()) {
    cout << "Error in the calcPlastic Procedure. point internal values are not "
            "finite. Press any key."
         << "\n";
    getchar();
  }

  double I1 = Point.firstInvariant();
  double I2 = Point.secondInvariant();
  double J2 = Point.secondDevInvariant();
  if (fabs(J2) < TINY)
    J2 = TINY;
  double J3 = Point.thirdDevInvariant();
  double shearStress = Point.shearStress();
  if (fabs(shearStress) < TINY)
    shearStress = TINY;
  double Factor = -27 * J3 / (2 * shearStress * shearStress * shearStress);
  if (Factor > 1)
    Factor = 1;
  if (Factor < -1)
    Factor = -1;
  if (!finite(Factor)) {
    Factor = 1.0;
    cout << "Factor in calcPlastic is not finite. set to 1" << "\n";
    cout << "alpha4=" << alpha4 << " J3=" << J3
         << " shearStress=" << shearStress << "\n";
  }
  Factor = 1 + alpha4 - (1 - alpha4) * Factor;
  double Factor025 = pow(Factor, 0.25);
  double Factor075 = pow(Factor, -0.75);
  double M = (3 - sin_phi) /
             (alpha * sin_phi); // not M but useful in derivatives, see file
  double Mpsi = (3 - sin_psi) / (alpha * sin_psi);

  BBMMatrix dJ2_dSig(6, 1), dJ3_dSig(6, 1), dq_dSig(6, 1), dI1_dSig(6, 1),
    dI2_dSig(6, 1), dI3_dSig(6, 1), TEMP(6, 1), NUMERATOR(6, 1);

  // derivatives of the invariants
  double s[6] = { Point.stress[0], Point.stress[1], Point.stress[2],
                  Point.stress[3], Point.stress[4], Point.stress[5] };

  dI1_dSig(0, 0) =  1.0);
  dI1_dSig(1, 0) =  1.0);
  dI1_dSig(2, 0) =  1.0); //{1,1,1,0,0,0}

  dI2_dSig(0, 0) =  s[1] + s[2]);
  dI2_dSig(1, 0) =  s[0] + s[2]);
  dI2_dSig(2, 0) =  s[0] + s[1]);
  dI2_dSig(3, 0) =  -2 * s[3]);
  dI2_dSig(4, 0) =  -2 * s[4]);
  dI2_dSig(5, 0) =  -2 * s[5]);

  // dI2_dSig.PrintPrecise();

  dI3_dSig(0, 0) =  s[1] * s[2] - s[5] * s[5]);
  dI3_dSig(1, 0) =  s[0] * s[2] - s[4] * s[4]);
  dI3_dSig(2, 0) =  s[0] * s[1] - s[3] * s[3]);
  dI3_dSig(3, 0) =  2 * s[5] * s[4] - 2 * s[2] * s[3]);
  dI3_dSig(4, 0) =  2 * s[3] * s[5] - 2 * s[1] * s[4]);
  dI3_dSig(5, 0) =  2 * s[3] * s[4] - 2 * s[0] * s[5]);

  // dI3_dSig.PrintPrecise();

  dJ2_dSig(0, 0) =  (2 * s[0] - s[1] - s[2]) / 3.0);
  dJ2_dSig(1, 0) =  (2 * s[1] - s[0] - s[2]) / 3.0);
  dJ2_dSig(2, 0) =  (2 * s[2] - s[0] - s[1]) / 3.0);
  dJ2_dSig(3, 0) =  2 * s[3]);
  dJ2_dSig(4, 0) =  2 * s[4]);
  dJ2_dSig(5, 0) =  2 * s[5]);

  dJ2_dSig.Copy(&dq_dSig);
  dq_dSig.Multiply(sqrt(3.0) * 0.5 / sqrt(J2), &dq_dSig);

  dI1_dSig.Copy(&dJ3_dSig);
  dJ3_dSig.Multiply(2.0 / 9.0 * I1 * I1 - I2 / 3.0, &dJ3_dSig);
  // dJ3_dSig.PrintPrecise();
  dI2_dSig.Copy(&TEMP);
  TEMP.Multiply(-I1 / 3.0, &TEMP);
  dJ3_dSig.Add(&TEMP, &dJ3_dSig);
  dJ3_dSig.Add(&dI3_dSig, &dJ3_dSig);
  // dJ3_dSig.PrintPrecise();
  // finished dJ3_dSig

  // above needs to be *CORRECTED* as some errors in ij i!=j components is
  // likely (2x too big)

  dJ2_dSig.Copy(&A);
  A.Multiply(0.21022410391343 * M * Factor025 / shearStress, &A);

  for (int i = 1; i < 7; i++) {
    if (!finite(A.getElement(i, 1))) {
      cout << "calcPlastic(1) DF/dSigma not finite at A(" << i
           << "). set to 1.0. Results may be incorrect. Press any key." << "\n";
      getchar();
      A.PutElement(i, 1, 0.0);
    }
  }

  // cout<<"dF/dJ2 part"<<"\n";
  // A.PrintPrecise();
  dJ3_dSig.Copy(&TEMP);
  TEMP.Multiply(-1.419012700740643 * (alpha4 - 1) * M * sqrt(J2) * Factor075 /
                    (sqrt(3.0) * shearStress * shearStress * shearStress) +
                  0.94600846716043 * (alpha4 - 1) * cohesion * M * Factor075 /
                    (shearStress * shearStress * shearStress),
                &TEMP);
  // cout<<"dF/dJ3 part"<<"\n";
  // TEMP.PrintPrecise();
  A.Add(&TEMP, &A);

  for (int i = 1; i < 7; i++) {
    if (!finite(A.getElement(i, 1))) {
      cout << "calcPlastic(2) DF/dSigma not finite at A(" << i
           << "). set to 1.0. Results may be incorrect. Press any key." << "\n";
      getchar();
      A.PutElement(i, 1, 0.0);
    }
  }

  dq_dSig.Copy(&TEMP);
  TEMP.Multiply(
    -2.838025401481287 * (alpha4 - 1) * cohesion * M * J3 * Factor075 /
        (shearStress * shearStress * shearStress * shearStress) +
      4.257038102221929 * (alpha4 - 1) * M * sqrt(J2) * J3 * Factor075 /
        (sqrt(3.0) * shearStress * shearStress * shearStress * shearStress),
    &TEMP);
  // cout<<"dF/dq part"<<"\n";
  // TEMP.PrintPrecise();
  A.Add(&TEMP, &A);
  dI1_dSig.Copy(&TEMP);
  TEMP.Multiply(-1 / 3.0, &TEMP);
  A.Add(&TEMP, &A);

  for (int i = 1; i < 7; i++) {
    if (!finite(A.getElement(i, 1))) {
      cout << "calcPlastic(3) DF/dSigma not finite at A(" << i
           << "). set to 1.0. Results may be incorrect. Press any key." << "\n";
      getchar();
      A.PutElement(i, 1, 0.0);
    }
  }

  /*A.PrintPrecise ();
  //FINISHED dF/dSigma

  StateMohrCoulomb CopyPoint;
  Point.Copy(&CopyPoint);
  double Yield1,Yield2,dSs=0.0001;
  Yield1=computeYieldFunction(&CopyPoint);
  CopyPoint.stress[0]+=dSs;
  Yield2=computeYieldFunctionNN(&CopyPoint);
  cout<<"DSigma="<<dSs<<"\n";
  cout<<"dF/dSigma11="<<(Yield2-Yield1)/dSs<<"\n";
  CopyPoint.stress[0]-=dSs;

  CopyPoint.stress[1]+=dSs;
  Yield2=computeYieldFunctionNN(&CopyPoint);
  cout<<"dF/dSigma22="<<(Yield2-Yield1)/dSs<<"\n";
  CopyPoint.stress[1]-=dSs;

  CopyPoint.stress[2]+=dSs;
  Yield2=computeYieldFunctionNN(&CopyPoint);
  cout<<"dF/dSigma33="<<(Yield2-Yield1)/dSs<<"\n";
  CopyPoint.stress[2]-=dSs;

  dSs=10*dSs;
  CopyPoint.stress[0]+=dSs;
  Yield2=computeYieldFunctionNN(&CopyPoint);
  cout<<"DSigma="<<dSs<<"\n";
  cout<<"dF/dSigma11="<<(Yield2-Yield1)/dSs<<"\n";
  CopyPoint.stress[0]-=dSs;

  CopyPoint.stress[1]+=dSs;
  Yield2=computeYieldFunctionNN(&CopyPoint);
  cout<<"dF/dSigma22="<<(Yield2-Yield1)/dSs<<"\n";
  CopyPoint.stress[1]-=dSs;

  CopyPoint.stress[2]+=dSs;
  Yield2=computeYieldFunctionNN(&CopyPoint);
  cout<<"dF/dSigma33="<<(Yield2-Yield1)/dSs<<"\n";
  CopyPoint.stress[2]-=dSs;

  dSs=10*dSs;
  CopyPoint.stress[0]+=dSs;
  Yield2=computeYieldFunctionNN(&CopyPoint);
  cout<<"DSigma="<<dSs<<"\n";
  cout<<"dF/dSigma11="<<(Yield2-Yield1)/dSs<<"\n";
  CopyPoint.stress[0]-=dSs;

  CopyPoint.stress[1]+=dSs;
  Yield2=computeYieldFunctionNN(&CopyPoint);
  cout<<"dF/dSigma22="<<(Yield2-Yield1)/dSs<<"\n";
  CopyPoint.stress[1]-=dSs;

  CopyPoint.stress[2]+=dSs;
  Yield2=computeYieldFunctionNN(&CopyPoint);
  cout<<"dF/dSigma33="<<(Yield2-Yield1)/dSs<<"\n";
  CopyPoint.stress[2]-=dSs;

  dSs=10*dSs;
  CopyPoint.stress[0]+=dSs;
  Yield2=computeYieldFunctionNN(&CopyPoint);
  cout<<"DSigma="<<dSs<<"\n";
  cout<<"dF/dSigma11="<<(Yield2-Yield1)/dSs<<"\n";
  CopyPoint.stress[0]-=dSs;

  CopyPoint.stress[1]+=dSs;
  Yield2=computeYieldFunctionNN(&CopyPoint);
  cout<<"dF/dSigma22="<<(Yield2-Yield1)/dSs<<"\n";
  CopyPoint.stress[1]-=dSs;

  CopyPoint.stress[2]+=dSs;
  Yield2=computeYieldFunctionNN(&CopyPoint);
  cout<<"dF/dSigma33="<<(Yield2-Yield1)/dSs<<"\n";
  CopyPoint.stress[2]-=dSs;

  dSs=10*dSs;
  CopyPoint.stress[0]+=dSs;
  Yield2=computeYieldFunctionNN(&CopyPoint);
  cout<<"DSigma="<<dSs<<"\n";
  cout<<"dF/dSigma11="<<(Yield2-Yield1)/dSs<<"\n";
  CopyPoint.stress[0]-=dSs;

  CopyPoint.stress[1]+=dSs;
  Yield2=computeYieldFunctionNN(&CopyPoint);
  cout<<"dF/dSigma22="<<(Yield2-Yield1)/dSs<<"\n";
  CopyPoint.stress[1]-=dSs;

  CopyPoint.stress[2]+=dSs;
  Yield2=computeYieldFunctionNN(&CopyPoint);
  cout<<"dF/dSigma33="<<(Yield2-Yield1)/dSs<<"\n";
  CopyPoint.stress[2]-=dSs;

  dSs=10*dSs;
  CopyPoint.stress[0]+=dSs;
  Yield2=computeYieldFunctionNN(&CopyPoint);
  cout<<"DSigma="<<dSs<<"\n";
  cout<<"dF/dSigma11="<<(Yield2-Yield1)/dSs<<"\n";
  CopyPoint.stress[0]-=dSs;

  CopyPoint.stress[1]+=dSs;
  Yield2=computeYieldFunctionNN(&CopyPoint);
  cout<<"dF/dSigma22="<<(Yield2-Yield1)/dSs<<"\n";
  CopyPoint.stress[1]-=dSs;

  CopyPoint.stress[2]+=dSs;
  Yield2=computeYieldFunctionNN(&CopyPoint);
  cout<<"dF/dSigma33="<<(Yield2-Yield1)/dSs<<"\n";
  CopyPoint.stress[2]-=dSs;
  */

  dJ2_dSig.Copy(&GG);
  GG.Multiply(0.21022410391343 * Mpsi * Factor025 / shearStress, &GG);
  dJ3_dSig.Copy(&TEMP);
  TEMP.Multiply(-1.419012700740643 * (alpha4 - 1) * Mpsi * sqrt(J2) *
                    Factor075 /
                    (sqrt(3.0) * shearStress * shearStress * shearStress) +
                  0.94600846716043 * (alpha4 - 1) * cohesion * Mpsi *
                    Factor075 / (shearStress * shearStress * shearStress),
                &TEMP);
  GG.Add(&TEMP, &GG);
  dq_dSig.Copy(&TEMP);
  TEMP.Multiply(
    -2.838025401481287 * (alpha4 - 1) * cohesion * Mpsi * J3 * Factor075 /
        (shearStress * shearStress * shearStress * shearStress) +
      4.257038102221929 * (alpha4 - 1) * Mpsi * sqrt(J2) * J3 * Factor075 /
        (sqrt(3.0) * shearStress * shearStress * shearStress * shearStress),
    &TEMP);
  GG.Add(&TEMP, &GG);
  dI1_dSig.Copy(&TEMP);
  TEMP.Multiply(
    -1 / 3.0,
    &TEMP); // First correction: sign here, but most likely not enough!!!
  GG.Add(&TEMP, &GG);

  // FINISHED dQ/dSigma

  double K43G = K + 4.0 * G / 3.0;
  double K23G = K - 2.0 * G / 3.0;

  DEL(0, 0) =  K43G);
  DEL(0, 1) =  K23G);
  DEL(0, 2) =  K23G); // rest of the line are zeros and rightly so
  DEL(1, 0) =  K23G);
  DEL(1, 1) =  K43G);
  DEL(1, 2) =  K23G); // yes, the matrix is symmetrical, but it is
                              // faster to put this 3 additional elements
  DEL.PutElement(
    3, 1,
    K23G); // than just mirror all, including zeros, which are there already
  DEL.PutElement(3, 2, K23G);
  DEL.PutElement(3, 3, K43G);
  DEL.PutElement(4, 4, 2.0 * G);
  DEL.PutElement(5, 5, 2.0 * G);
  DEL.PutElement(6, 6, 2.0 * G); // rest of the matrix is filled with zeros...

  // getting lambda and Dep

  A.Transpose(&NUMERATOR);
  NUMERATOR.Multiply(&DEL, &NUMERATOR); // NUMERATOR=aT*Del -->Numerator of
                                        // Lambda without multiplication by
                                        // dEpsilon
  NUMERATOR.Multiply(&GG, &DENOMINATOR);
  NUMERATOR.Multiply(&DEPS, &NUMERATOR); // NUMERATOR=aT*Del*dEps -->Numerator
                                         // of Lambda without multiplication by
                                         // dEpsilon

  if (USE_NICE_SCHEME > 0) {
    // NUMERATOR.PutElement(1,1,NUMERATOR.getElement(1,1)+computeYieldFunction(&Point));
    NUMERATOR(0, 0) =  NUMERATOR.getElement(1, 1) + fValue);
    NUMERATOR(0, 0) =  NUMERATOR.getElement(1, 1) +
                                 (*dLambda) * DENOMINATOR.getElement(1, 1));
    *dLambda = fValue / DENOMINATOR.getElement(1, 1);
    DEL.Multiply(&GG, &TEMP);
    TEMP.Multiply(*dLambda, &TEMP);
    dS[0] = -TEMP.getElement(1, 1);
    dS[1] = -TEMP.getElement(2, 1);
    dS[2] = -TEMP.getElement(3, 1);
    dS[3] = -TEMP.getElement(4, 1);
    dS[4] = -TEMP.getElement(5, 1);
    dS[5] = -TEMP.getElement(6, 1);
  }

  DEL.Multiply(&GG, &TEMP);
  TEMP.Multiply(&NUMERATOR, &TEMP);
  if (DENOMINATOR.getElement(1, 1) < TINY)
    cout << "Denominator of plastic multiplier is very small. Some error may "
            "arise and results may be incorrect"
         << "\n";

  TEMP.Multiply(1 / DENOMINATOR.getElement(1, 1), &TEMP);
  DEL.Multiply(&DEPS, dSigma);
  dSigma->Substract(&TEMP, dSigma);
  for (int i = 1; i < 7; i++) {
    if (!finite(dSigma->getElement(i, 1))) {
      cout << "calcPlastic: dSigma not finite: set to zero. Results may be "
              "erronous."
           << "\n";
      dSigma->PutElement(i, 1, 0.0);
    }
  }
  *dpZeroStar = 1.0;
  // Stress increment computed and in dSigma matrix

  return 0;

  // that's all?
}

void
ShengMohrCoulomb::getTangentMatrix(StateMohrCoulomb* point, BBMMatrix* DEP)
{
  /*
  This procedure returns tangent matrix. If given point lies on the yield locus,
  the matrix is elastoplastic.
  When the point lies inside the yield locus, the matrix is elastic.
  Point is supposed to have all the required data about describing its state,
  like stress state, specific volume nu and the
  hardening parameter p0*.

  Input: point, Matrix (7x6 - six 6 rows, seven 7 columns) that will be filled
  with coefficients corresponding to the tangent matrix
  Initialise BBMMatrix (6,7);
  */

  DEP->Resize(6, 7);

  if (computeYieldNormalized(Point)) {
    // calculate elasto-plastic matrix
    calculateElastoPlasticTangentMatrix(Point, DEP);
    // DEP->Print();

    // calculateElasticTangentMatrix (Point, DEP);
    // cout<<"Elasto-Plastic Matrix used"<<"\n";
  } else {
    // calculate elastic matrix
    calculateElasticTangentMatrix(Point, DEP);
  }
}

void
ShengMohrCoulomb::getTangentMatrixPQ(StateMohrCoulomb* point, BBMMatrix* DEP)
{
  cout << "getTangentMatrixPQ is unsupported for Mohr-Coulomb Model" << "\n";
  getchar();
  // calculateElastoPlasticTangentMatrix (Point,DEP);
}

void
ShengMohrCoulomb::calculateElasticTangentMatrixPQ(StateMohrCoulomb* point,
                                                  BBMMatrix* DEP)
{
  cout
    << "calculateElasticTangentMatrixPQ is unsupported for Mohr-Coulomb Model"
    << "\n";
  getchar();
  // calculateElasticTangentMatrix (Point, DEP);
}

void
ShengMohrCoulomb::calculateElastoPlasticTangentMatrix(StateMohrCoulomb* point,
                                                      BBMMatrix* DEP)
{
  /*************************************
  Here the ElastoPlastic tangent matrix is calculated.
  For more theoretical explanation what is being done below please see
  'Evaluation of the tangent elastic matrix Del & elastoplastic matrix Dep'
  (evaluation of the tangent matrix D.doc)

  This procedure is quite similar to procedure calcPlastic where the stress
  increment using
  the tangent elasto-plastic matrix is calculated. The matrix itself, however,
  is not explicitly
  calculated there.


          // at the beginning we need to calculate the derivatives a, b, c, d,
  g, p...)
  BBMMatrix A (6,1); //dF/dsigma
  BBMMatrix DENOMINATOR (1,1);

  BBMMatrix GG (6,1); //dG/dsigma, in case of associated flow
  (nonAssociated==false) same as A
  BBMMatrix MM (6,1); // will be vector(1,1,1,0,0,0)T
  BBMMatrix DEL (6,6); //D elastic matrix...

  // we do not have g, as a == g in associated flow rule.

  double I1=Point->firstInvariant();
  double I2=Point->secondInvariant();
  double J2=Point->secondDevInvariant ();
  double J3=Point->thirdDevInvariant ();
  double shearStress=Point->shearStress ();
  double Factor=-27*J3/(2*shearStress*shearStress*shearStress);
  if (Factor>1) Factor=1;
  if (Factor<-1) Factor=-1;
  if (!finite(Factor))
                  {
                       Factor=1.0;
                       cout<<"Factor in calculate Elasto-Plastic Tangent Matrix
  is not finite. set to 1"<<"\n";
                       cout<<"alpha4="<<alpha4<<" J3="<<J3<<"
  shearStress="<<shearStress<<"\n";
                  }
  Factor=1+alpha4-(1-alpha4)*Factor;
  double Factor025=pow(Factor,0.25);
  double Factor075=pow(Factor,-0.75);
  double M=(3-sin_phi)/(6*alpha*sin_phi);
  double Mpsi=(3-sin_psi)/(6*alpha*sin_psi);



  BBMMatrix dJ2_dSig (6,1), dJ3_dSig (6,1), dq_dSig (6,1), dI1_dSig (6,1), dI2_dSig
  (6,1), dI3_dSig (6,1), TEMP(6,1), NUMERATOR(6,1);

  //derivatives of the invariants
  double
  s[6]={Point->stress[0],Point->stress[1],Point->stress[2],Point->stress[3],Point->stress[4],Point->stress[5]};

  dI1_dSig.PutElement (1,1,1.0); dI1_dSig.PutElement (2,1,1.0);dI1_dSig.PutElement
  (3,1,1.0); //{1,1,1,0,0,0}

  dI2_dSig.PutElement (1,1,s[1]+s[2]);dI2_dSig.PutElement
  (2,1,s[0]+s[2]);dI2_dSig.PutElement (3,1,s[0]+s[1]);
  dI2_dSig.PutElement (4,1,-2*s[3]);dI2_dSig.PutElement
  (5,1,-2*s[4]);dI2_dSig.PutElement (6,1,-2*s[5]);

  dI3_dSig.PutElement (1,1,s[1]*s[2]-s[5]*s[5]);dI3_dSig.PutElement
  (2,1,s[0]*s[2]-s[4]*s[4]);dI3_dSig.PutElement (3,1,s[0]*s[1]-s[3]*s[3]);
  dI3_dSig.PutElement (4,1,2*s[5]*s[4]-2*s[2]*s[3]);dI3_dSig.PutElement
  (5,1,2*s[3]*s[5]-2*s[1]*s[4]);dI3_dSig.PutElement
  (6,1,2*s[3]*s[4]-2*s[0]*s[5]);

  dJ2_dSig.PutElement (1,1,(2*s[0]-s[1]-s[2])/3.0);dJ2_dSig.PutElement
  (2,1,(2*s[1]-s[0]-s[2])/3.0);dJ2_dSig.PutElement (3,1,(2*s[2]-s[0]-s[1])/3.0);
  dJ2_dSig.PutElement (4,1,2*s[3]);dJ2_dSig.PutElement
  (5,1,2*s[4]);dJ2_dSig.PutElement (6,1,2*s[5]);

  dJ2_dSig.Copy(&dq_dSig);
  dq_dSig.Multiply(sqrt(3.0)*0.5/sqrt(J2),&dq_dSig);

  dI1_dSig.Copy(&dJ3_dSig);
  dJ3_dSig.Multiply(6.0/27.0*I1*I1-I2/3.0,&dJ3_dSig);
  dI2_dSig.Copy(&TEMP);
  TEMP.Multiply(-I1/3.0,&TEMP);
  dJ3_dSig.Add(&TEMP,&dJ3_dSig);
  dJ3_dSig.Add(&dI3_dSig,&dJ3_dSig);
  //finished dJ3_dSig

  //above needs to be *CORRECTED* as some errors in ij i!=j components is likely
  (2x too big)

  dJ2_dSig.Copy(&A);
  A.Multiply (0.21022410391343*M*Factor025/shearStress,&A);
  dJ3_dSig.Copy(&TEMP);
  TEMP.Multiply
  (1.419012700740643*(alpha4-1)*M*sqrt(J2)*Factor075/(sqrt(3.0)*shearStress*shearStress*shearStress)-
                             0.94600846716043*(alpha4-1)*cohesion*M*Factor075/(shearStress*shearStress*shearStress),&TEMP);
  A.Add (&TEMP,&A);
  dq_dSig.Copy(&TEMP);
  TEMP.Multiply
  (2.838025401481287*(alpha4-1)*cohesion*M*J3*Factor075/(shearStress*shearStress*shearStress*shearStress)-
                             4.257038102221929*(alpha4-1)*M*sqrt(J2)*J3*Factor075/(sqrt(3.0)*shearStress*shearStress*shearStress*shearStress),&TEMP);
  A.Add (&TEMP,&A);
  dI1_dSig.Copy(&TEMP);
  TEMP.Multiply(-1/3.0, &TEMP);
  A.Add (&TEMP,&A);

  //FINISHED dF/dSigma

  //FINISHED dF/dSigma

  dJ2_dSig.Copy(&GG);
  GG.Multiply (0.21022410391343*Mpsi*Factor025/shearStress,&GG);
  dJ3_dSig.Copy(&TEMP);
  TEMP.Multiply
  (1.419012700740643*(alpha4-1)*Mpsi*sqrt(J2)*Factor075/(sqrt(3.0)*shearStress*shearStress*shearStress)-
                             0.94600846716043*(alpha4-1)*cohesion*Mpsi*Factor075/(shearStress*shearStress*shearStress),&TEMP);
  GG.Add (&TEMP,&GG);
  dq_dSig.Copy(&TEMP);
  TEMP.Multiply
  (2.838025401481287*(alpha4-1)*cohesion*Mpsi*J3*Factor075/(shearStress*shearStress*shearStress*shearStress)-
                             4.257038102221929*(alpha4-1)*Mpsi*sqrt(J2)*J3*Factor075/(sqrt(3.0)*shearStress*shearStress*shearStress*shearStress),&TEMP);
  GG.Add (&TEMP,&GG);
  dI1_dSig.Copy(&TEMP);
  TEMP.Multiply(-1/3.0, &TEMP);  //First correction: sign here, but most likely
  not enough!!!
  GG.Add (&TEMP,&GG);

  //FINISHED dQ/dSigma

  double K43G=K+4.0*G/3.0;
  double K23G=K-2.0*G/3.0;

  DEL.PutElement (1,1,K43G);
  DEL.PutElement (1,2,K23G);
  DEL.PutElement (1,3,K23G); //rest of the line are zeros and rightly so
  DEL.PutElement (2,1,K23G);
  DEL.PutElement (2,2,K43G);
  DEL.PutElement (2,3,K23G); //yes, the matrix is symmetrical, but it is faster
  to put this 3 additional elements
  DEL.PutElement (3,1,K23G); //than just mirror all, including zeros, which are
  there already
  DEL.PutElement (3,2,K23G);
  DEL.PutElement (3,3,K43G);
  DEL.PutElement (4,4,2.0*G);
  DEL.PutElement (5,5,2.0*G);
  DEL.PutElement (6,6,2.0*G); //rest of the matrix is filled with zeros...

  //getting lambda and Dep

  A.Transpose(&NUMERATOR);
  NUMERATOR.Multiply(&DEL,&NUMERATOR); //NUMERATOR=aT*Del -->Numerator of Lambda
  without multiplication by dEpsilon
  NUMERATOR.Multiply(&GG,&DENOMINATOR);

  DEL.Multiply(&GG,&TEMP);
  TEMP.Multiply(&NUMERATOR,&TEMP);
  TEMP.Multiply(1/DENOMINATOR.getElement(1,1),&TEMP);
  DEL.Substract(&TEMP,&DEL);
  for (int i=1;i<7; i++) for (int j=1; j<7; j++)
  DEP->PutElement(i,j,DEL.getElement(i,j));
  //so if 6x7 is needed, it is preserved for compatibility purposes

  //Stress increment computed and in dSigma matrix

  //DEP->Print();
  //getchar(); */
  // FINAL MATRIX PUT TOGETHER, 6 rows, 7 columns
}



void
ShengMohrCoulomb::getDerivative(double meanStress, double shearStress,
                                double suction, double PZero, double* state,
                                double* deriv)
{
  cout << "getDerivative is not implemented for Mohr-Coulomb model." << "\n";
  getchar();
  /*	double SuctionPressure=k*suction;
          deriv[0]=M*M*(2*meanStress+SuctionPressure-PZero);  // dF/dp
          deriv[1]=2*shearStress;	//dF/dq
          //now we have to calculate df/ds... this is the most difficult (see
     word file)
          //p0* is == state [0]
          //pc is a parameter
          double P0Star=state[0];
          double denominator=LambdaZero*((1-r)*exp(-1*Beta*suction)+r);
          double numerator=LambdaZero-KappaP;
          deriv[2]=-1*Beta*(1-r)*exp(-1*Beta*suction);
          //cout<<"1:"<<deriv[2]<<"\n";
          deriv[2]=deriv[2]*-1*numerator/(denominator*denominator);
          //cout<<"2:"<<deriv[2]<<"\n";
          deriv[2]=deriv[2]*pc*log(P0Star/pc);
          //cout<<"3:"<<deriv[2]<<"\n";
          deriv[2]=deriv[2]*pow((P0Star/pc),numerator/denominator); //it's
     finished with dp0/ds
          //cout<<"4:"<<deriv[2]<<"\n";
          deriv[2]=M*M*(meanStress+SuctionPressure)*deriv[2];
          //cout<<"5:"<<deriv[2]<<"\n";
          deriv[2]=-1*M*M*k*(PZero-meanStress)-deriv[2];	//final result

  */
  // double ShengMohrCoulomb::findGradientPQ (double * state, double * s, double
  // *ds, double * dF, double suction, double dsuction)
}

double
ShengMohrCoulomb::findGradientPQ(StateMohrCoulomb* point, double* ds, double* dF,
                                 double dsuction)
{
  cout << "findGradientPQ is unsupported for Mohr-Coulomb Model" << "\n";
  /*
  double PZero, SuctionPressure, dmeanStress, dshearStress, fValue, dFlength,
  dslength, cosin;
  double P0Star=Point->p0Star();
  double meanStress=Point->getmeanStress();
  double shearStress=Point->shearStress();
  double Suction=Point->suction();

          //check of the second yield surface:
  if (Point->yieldSuction()-Suction<-d_suctionTol)
          {
                  cout<<"yield of the suction surface occurred"<<"\n";
                  if (dsuction>0)
                  {
                          cout <<"Unable to find gradient; Whole step
  plastic"<<"\n";
                          return -2;
                  //this check may be unused...
                  }
          }
  //check of the standard yield surface

          //cout<<"meanStress="<<meanStress<<"\n";
          //cout<<"shearStress="<<shearStress<<"\n";

          PZero=LambdaZero*((1-r)*exp(-1*Beta*Suction)+r);
      PZero=(LambdaZero-KappaP)/(PZero-KappaP);
          PZero=pc*pow((P0Star/pc),PZero);
          SuctionPressure=k*Suction;
  //	cout<<"suction="<<suction<<"\n";
  //	cout<<"PZero="<<PZero<<"\n";
          fValue=shearStress*shearStress-M*M*(meanStress+SuctionPressure)*(PZero-meanStress);
          fValue=fValue/((PZero+SuctionPressure)*(PZero+SuctionPressure));
  //value of Yield function calculated and normalized
          if (fValue<-d_yieldTol)
          {
                  cout<<"!!!there is no yield at the beginning !!!"<<"\n";
                  cout<<"fValue is="<<fValue<<"\n";
                  cout<<"!!!find gradient procedure
  terminated!!!"<<"\n"<<"\n"<<"\n"<<"\n";
                  return -3;
                  // this check may be disabled later
          }

          dF[0]=M*M*(2*meanStress+SuctionPressure-PZero);  // dF/dp
          dF[1]=2*shearStress;	//dF/dq
          //now we have to calculate df/ds... this is the most difficult (see
  word file)
          //p0* is == state [0]
          //pc is a parameter

          double denominator=LambdaZero*((1-r)*exp(-1*Beta*Suction)+r);
          double numerator=LambdaZero-KappaP;
          dF[2]=-1*Beta*(1-r)*exp(-1*Beta*Suction);
          //cout<<"1:"<<dF[2]<<"\n";
          dF[2]=dF[2]*-1*numerator/(denominator*denominator);
          //cout<<"2:"<<dF[2]<<"\n";
          dF[2]=dF[2]*pc*log(P0Star/pc);
          //cout<<"3:"<<dF[2]<<"\n";
          dF[2]=dF[2]*pow((P0Star/pc),numerator/denominator); //it's finished
  with dp0/ds
          //cout<<"4:"<<dF[2]<<"\n";
          dF[2]=M*M*(meanStress+SuctionPressure)*dF[2];
          //cout<<"5:"<<dF[2]<<"\n";
          dF[2]=-1*M*M*k*(PZero-meanStress)-dF[2];	//final result
          //cout<<dF[2]<<"\n";
          //calculate changes in p and q...;

          dmeanStress=(ds[0]+ds[1]+ds[2])/3;
          double dds[6];
          for (int i=0; i<6; i++) dds[i]=ds[i]+Point->stress[i];


          dshearStress=(dds[0]-dds[1])*(dds[0]-dds[1])+(dds[0]-dds[2])*(dds[0]-dds[2])+(dds[1]-dds[2])*(dds[1]-dds[2]);
          dshearStress=dshearStress+6*(dds[3]*dds[3]+dds[4]*dds[4]+dds[5]*dds[5]);
          dshearStress=dshearStress/2;
          dshearStress=sqrt(dshearStress);
          //we must not calculate change of shear stress only basing on change
  of stress; it must be done this way... it's not linear

          dshearStress=dshearStress-shearStress;
          //this line should be removed later
          //dF[2]=0;
          //dF[3]=0;dF[4]=0;dF[5]=0;dF[6]=0; //don't take suction into account
  at the moment

          //dmeanStress=ds[0];
          //dshearStress=ds[1];

          dFlength=dF[0]*dF[0]+dF[1]*dF[1]+dF[2]*dF[2]; //calculated length
          dFlength=sqrt(dFlength);
          dslength=dshearStress*dshearStress+dmeanStress*dmeanStress+dsuction*dsuction;
          dslength=sqrt(dslength);
          for (int i=0; i<3; i++) dF[i]=dF[i]/dFlength;
          cosin=(dF[0]*dmeanStress+dF[1]*dshearStress+dF[2]*dsuction)/dslength;
  //it should be d vector multiplied by gradient
          cout<<"Mean Stress="<<meanStress<<"  Shear Stress="<<shearStress<< "
  Suction="<<Suction<<"\n";
          cout<<"dMean Stress="<<dmeanStress<<"  dShear
  Stress="<<dshearStress<<" dSuction="<<dsuction<<"\n";
          for (int i=0; i<3; i++) cout<<"value of F["<<i<<"] is:"<<dF[i]<<"\n";
          cout<<"df length is:"<<dFlength<<"\n";
          cout<<"cosinus is:"<<cosin<<"\n";
          getchar (); */
  //	return cosin;
  return 0;
}

void
ShengMohrCoulomb::correctDrift(StateMohrCoulomb* Point)
{

  /*
  This procedure should correct the drift, as described in the word file. The
  method of drift correction is based on
  the Potts&Zdravkovic book; It is however slightly changes so to use it in
  unsaturated soil model. Detailed
  description of the algorithm is presented in the word file. It is noted
  however, that in current version the algorithm
  does calculate all the derivatives in the forbidden space outside the yield
  locus which may be cause of concerns.

  Input: point to correct
  Output: Updated point, corrected using values of stress at the end
  */

  // Need dF/Ds, D, dF/dP0*, dP0*/DEpsPl, m, dF/dS
  // at the beginning we need to calculate the derivatives a, b, c, d, g, p...)

  // cout<<"Correct Drift Procedure entered!"<<"\n";

  BBMMatrix A(6, 1); // dF/dsigma
  BBMMatrix DENOMINATOR(1, 1);

  BBMMatrix GG(6, 1);  // dG/dsigma, in case of associated flow
                       // (nonAssociated==false) same as A
  BBMMatrix MM(6, 1);  // will be vector(1,1,1,0,0,0)T
  BBMMatrix DEL(6, 6); // D elastic matrix...
  BBMMatrix DEPS(6, 1);
  BBMMatrix dSigma(6, 1);

  double DSigma[6], epStrain[6], zeros[7];
  for (int i = 0; i < 7; i++)
    zeros[i] = 0;

  for (int i = 1; i < 6; i++)
    DEPS.PutElement(i, 1, epStrain[i - 1]); // increase of epsilon copied
  // we do not have g, as a == g in associated flow rule.

  bool correct;

  double I1 = Point->firstInvariant();
  double I2 = Point->secondInvariant();
  double J2 = Point->secondDevInvariant();
  if (fabs(J2) < TINY)
    J2 = TINY;
  double J3 = Point->thirdDevInvariant();
  double shearStress = Point->shearStress();
  if (fabs(shearStress) < TINY)
    shearStress = TINY;
  double Factor = -27 * J3 / (2 * shearStress * shearStress * shearStress);
  if (Factor > 1)
    Factor = 1;
  if (Factor < -1)
    Factor = -1;
  if (!finite(Factor)) {
    Factor = 1.0;
    cout << "Factor in Correct Drift is not finite. set to 1" << "\n";
    cout << "alpha4=" << alpha4 << " J3=" << J3
         << " shearStress=" << shearStress << "\n";
  }
  Factor = 1 + alpha4 - (1 - alpha4) * Factor;
  double Factor025 = pow(Factor, 0.25);
  double Factor075 = pow(Factor, -0.75);
  double M = (3 - sin_phi) / (alpha * sin_phi);
  double Mpsi = (3 - sin_psi) / (alpha * sin_psi);
  double fValue;

  int numberIter = 0;

  do {
    computeYieldNormalized(Point->state, point->stress, point->suction(),
               &fValue); // 20 Feb 2006, preparations for the drift correction

    if ((fValue / (Point->getmeanStress() + 2 * cohesion) < -d_yieldTol) ||
        (fValue / (Point->getmeanStress() + 2 * cohesion) > d_yieldTol))
      correct = TRUE;
    else
      correct = FALSE;
    if (correct == TRUE) {
      numberIter++;
      // cout<<"Drift Correction, Iteration="<<numberIter<<" Function
      // Value="<<fValue<<"\n";
      // CORRECT FOR DRIFT
      // HERE THE DRIFT WILL BE CORRECTED BY USING THE D MATRIX FROM THE
      // FORBIDDEN SPACE
      // ALTHOUGH BECAUSE THE DRIFT WILL BE CHECKED AGAIN, IT SHOULDN'T POSE
      // MUCH PROBLEM.

      BBMMatrix dJ2_dSig(6, 1), dJ3_dSig(6, 1), dq_dSig(6, 1), dI1_dSig(6, 1),
        dI2_dSig(6, 1), dI3_dSig(6, 1), TEMP(6, 1), NUMERATOR(6, 1);

      // derivatives of the invariants
      double s[6] = { Point->stress[0], Point->stress[1], Point->stress[2],
                      Point->stress[3], Point->stress[4], Point->stress[5] };

      dI1_dSig(0, 0) =  1.0);
      dI1_dSig(1, 0) =  1.0);
      dI1_dSig(2, 0) =  1.0); //{1,1,1,0,0,0}

      dI2_dSig(0, 0) =  s[1] + s[2]);
      dI2_dSig(1, 0) =  s[0] + s[2]);
      dI2_dSig(2, 0) =  s[0] + s[1]);
      dI2_dSig(3, 0) =  -2 * s[3]);
      dI2_dSig(4, 0) =  -2 * s[4]);
      dI2_dSig(5, 0) =  -2 * s[5]);

      // dI2_dSig.PrintPrecise();

      dI3_dSig(0, 0) =  s[1] * s[2] - s[5] * s[5]);
      dI3_dSig(1, 0) =  s[0] * s[2] - s[4] * s[4]);
      dI3_dSig(2, 0) =  s[0] * s[1] - s[3] * s[3]);
      dI3_dSig(3, 0) =  2 * s[5] * s[4] - 2 * s[2] * s[3]);
      dI3_dSig(4, 0) =  2 * s[3] * s[5] - 2 * s[1] * s[4]);
      dI3_dSig(5, 0) =  2 * s[3] * s[4] - 2 * s[0] * s[5]);

      // dI3_dSig.PrintPrecise();

      dJ2_dSig(0, 0) =  (2 * s[0] - s[1] - s[2]) / 3.0);
      dJ2_dSig(1, 0) =  (2 * s[1] - s[0] - s[2]) / 3.0);
      dJ2_dSig(2, 0) =  (2 * s[2] - s[0] - s[1]) / 3.0);
      dJ2_dSig(3, 0) =  2 * s[3]);
      dJ2_dSig(4, 0) =  2 * s[4]);
      dJ2_dSig(5, 0) =  2 * s[5]);

      dJ2_dSig.Copy(&dq_dSig);
      dq_dSig.Multiply(sqrt(3.0) * 0.5 / sqrt(J2), &dq_dSig);

      dI1_dSig.Copy(&dJ3_dSig);
      dJ3_dSig.Multiply(2.0 / 9.0 * I1 * I1 - I2 / 3.0, &dJ3_dSig);
      // dJ3_dSig.PrintPrecise();
      dI2_dSig.Copy(&TEMP);
      TEMP.Multiply(-I1 / 3.0, &TEMP);
      dJ3_dSig.Add(&TEMP, &dJ3_dSig);
      dJ3_dSig.Add(&dI3_dSig, &dJ3_dSig);
      // dJ3_dSig.PrintPrecise();
      // finished dJ3_dSig

      // above needs to be *CORRECTED* as some errors in ij i!=j components is
      // likely (2x too big)

      dJ2_dSig.Copy(&A);
      A.Multiply(0.21022410391343 * M * Factor025 / shearStress, &A);

      for (int i = 1; i < 7; i++) {
        if (!finite(A.getElement(i, 1))) {
          A.PutElement(i, 1, 0.0);
        }
      }

      // cout<<"dF/dJ2 part"<<"\n";
      // A.PrintPrecise();
      dJ3_dSig.Copy(&TEMP);
      TEMP.Multiply(-1.419012700740643 * (alpha4 - 1) * M * sqrt(J2) *
                        Factor075 /
                        (sqrt(3.0) * shearStress * shearStress * shearStress) +
                      0.94600846716043 * (alpha4 - 1) * cohesion * M *
                        Factor075 / (shearStress * shearStress * shearStress),
                    &TEMP);
      // cout<<"dF/dJ3 part"<<"\n";
      // TEMP.PrintPrecise();
      A.Add(&TEMP, &A);

      for (int i = 1; i < 7; i++) {
        if (!finite(A.getElement(i, 1))) {
          A.PutElement(i, 1, 0.0);
        }
      }

      dq_dSig.Copy(&TEMP);
      TEMP.Multiply(
        -2.838025401481287 * (alpha4 - 1) * cohesion * M * J3 * Factor075 /
            (shearStress * shearStress * shearStress * shearStress) +
          4.257038102221929 * (alpha4 - 1) * M * sqrt(J2) * J3 * Factor075 /
            (sqrt(3.0) * shearStress * shearStress * shearStress * shearStress),
        &TEMP);
      // cout<<"dF/dq part"<<"\n";
      // TEMP.PrintPrecise();
      A.Add(&TEMP, &A);
      dI1_dSig.Copy(&TEMP);
      TEMP.Multiply(-1 / 3.0, &TEMP);
      A.Add(&TEMP, &A);

      for (int i = 1; i < 7; i++) {
        if (!finite(A.getElement(i, 1))) {
          A.PutElement(i, 1, 1.0);
        }
      }

      double meanStress = Point->getmeanStress();
      if (meanStress < -cohesion * M) {
        // tension cut-off plane
        A(0, 0) =  1.0 / 3.0);
        A(1, 0) =  1.0 / 3.0);
        A(2, 0) =  1.0 / 3.0);
        A(3, 0) =  0.0);
        A(4, 0) =  0.0);
        A(5, 0) =  0.0);
      }

      /*A.PrintPrecise ();
      //FINISHED dF/dSigma

      StateMohrCoulomb CopyPoint;
      Point->Copy(&CopyPoint);
      double Yield1,Yield2,dSs=0.0001;
      Yield1=computeYieldFunction(&CopyPoint);
      CopyPoint.stress[0]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      cout<<"DSigma="<<dSs<<"\n";
      cout<<"dF/dSigma11="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[0]-=dSs;

      CopyPoint.stress[1]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      cout<<"dF/dSigma22="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[1]-=dSs;

      CopyPoint.stress[2]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      cout<<"dF/dSigma33="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[2]-=dSs;

      dSs=10*dSs;
      CopyPoint.stress[0]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      cout<<"DSigma="<<dSs<<"\n";
      cout<<"dF/dSigma11="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[0]-=dSs;

      CopyPoint.stress[1]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      cout<<"dF/dSigma22="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[1]-=dSs;

      CopyPoint.stress[2]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      cout<<"dF/dSigma33="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[2]-=dSs;

      dSs=10*dSs;
      CopyPoint.stress[0]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      cout<<"DSigma="<<dSs<<"\n";
      cout<<"dF/dSigma11="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[0]-=dSs;

      CopyPoint.stress[1]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      cout<<"dF/dSigma22="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[1]-=dSs;

      CopyPoint.stress[2]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      cout<<"dF/dSigma33="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[2]-=dSs;

      dSs=10*dSs;
      CopyPoint.stress[0]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      cout<<"DSigma="<<dSs<<"\n";
      cout<<"dF/dSigma11="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[0]-=dSs;

      CopyPoint.stress[1]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      cout<<"dF/dSigma22="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[1]-=dSs;

      CopyPoint.stress[2]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      cout<<"dF/dSigma33="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[2]-=dSs;

      dSs=10*dSs;
      CopyPoint.stress[0]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      cout<<"DSigma="<<dSs<<"\n";
      cout<<"dF/dSigma11="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[0]-=dSs;

      CopyPoint.stress[1]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      cout<<"dF/dSigma22="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[1]-=dSs;

      CopyPoint.stress[2]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      cout<<"dF/dSigma33="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[2]-=dSs;

      dSs=10*dSs;
      CopyPoint.stress[0]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      cout<<"DSigma="<<dSs<<"\n";
      cout<<"dF/dSigma11="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[0]-=dSs;

      CopyPoint.stress[1]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      cout<<"dF/dSigma22="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[1]-=dSs;

      CopyPoint.stress[2]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      cout<<"dF/dSigma33="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[2]-=dSs;
      */
      dJ2_dSig.Copy(&GG);
      GG.Multiply(0.21022410391343 * Mpsi * Factor025 / shearStress, &GG);
      dJ3_dSig.Copy(&TEMP);
      TEMP.Multiply(-1.419012700740643 * (alpha4 - 1) * Mpsi * sqrt(J2) *
                        Factor075 /
                        (sqrt(3.0) * shearStress * shearStress * shearStress) +
                      0.94600846716043 * (alpha4 - 1) * cohesion * Mpsi *
                        Factor075 / (shearStress * shearStress * shearStress),
                    &TEMP);
      GG.Add(&TEMP, &GG);
      dq_dSig.Copy(&TEMP);
      TEMP.Multiply(
        -2.838025401481287 * (alpha4 - 1) * cohesion * Mpsi * J3 * Factor075 /
            (shearStress * shearStress * shearStress * shearStress) +
          4.257038102221929 * (alpha4 - 1) * Mpsi * sqrt(J2) * J3 * Factor075 /
            (sqrt(3.0) * shearStress * shearStress * shearStress * shearStress),
        &TEMP);
      GG.Add(&TEMP, &GG);
      dI1_dSig.Copy(&TEMP);
      TEMP.Multiply(
        -1 / 3.0,
        &TEMP); // First correction: sign here, but most likely not enough!!!
      GG.Add(&TEMP, &GG);

      // FINISHED dQ/dSigma

      if (meanStress < -cohesion * M) {
        // tension cut-off plane
        GG(0, 0) =  1.0 / 3.0);
        GG(1, 0) =  1.0 / 3.0);
        GG(2, 0) =  1.0 / 3.0);
        GG(3, 0) =  0.0);
        GG(4, 0) =  0.0);
        GG(5, 0) =  0.0);
      }

      double K43G = K + 4.0 * G / 3.0;
      double K23G = K - 2.0 * G / 3.0;

      DEL(0, 0) =  K43G);
      DEL(0, 1) =  K23G);
      DEL(0, 2) =  K23G); // rest of the line are zeros and rightly so
      DEL(1, 0) =  K23G);
      DEL(1, 1) =  K43G);
      DEL(1, 2) =  K23G); // yes, the matrix is symmetrical, but it is
                                  // faster to put this 3 additional elements
      DEL(2, 0) =  K23G); // than just mirror all, including zeros,
                                  // which are there already
      DEL.PutElement(3, 2, K23G);
      DEL.PutElement(3, 3, K43G);
      DEL.PutElement(4, 4, 2.0 * G);
      DEL.PutElement(5, 5, 2.0 * G);
      DEL.PutElement(6, 6,
                     2.0 * G); // rest of the matrix is filled with zeros...

      // getting lambda and Dep

      A.Transpose(&NUMERATOR);
      NUMERATOR.Multiply(&DEL, &NUMERATOR); // NUMERATOR=aT*Del -->Numerator of
                                            // Lambda without multiplication by
                                            // dEpsilon
      NUMERATOR.Multiply(&GG, &DENOMINATOR);

      double Lambda = fValue / DENOMINATOR.getElement(1, 1);

      A.Multiply(Lambda, &TEMP); // delta epsilon plastic= -delta epsilon
                                 // elastic

      // cout<<"Delta Epsilon Plastic:"<<"\n";
      // TEMP.Print();

      for (int i = 1; i < 7; i++)
        epStrain[i - 1] = TEMP.getElement(i, 1); // epsilon pl change

      DEL.Multiply(&GG, &TEMP);
      TEMP.Multiply(-Lambda, &dSigma);
      // final result for stress change, with the negative sign, so the stress
      // should be ADDED
      // be ADDED to get the right result of corrected stresses.
      for (int i = 0; i < 6; i++)
        DSigma[i] = dSigma.getElement(i + 1, 1);
      Point->Update(epStrain, zeros, DSigma, 0);
    }
    if (numberIter > 10) {
      // cout<<"Drift Correction Procedure failed"<<"\n";
      correct = FALSE;
    }
  } while (correct == TRUE);

  /*


bool correct;
double fValue, dpZeroStar, Lambda, DSigma[6], epStrain[6], zeros[7];
for (int i=0; i<7; i++) zeros[i]=0;

double temp, PZero, meanStress, shearStress, Suction, PZeroStar, SpecificVolume, LambdaS, Fraction;
// we do not have g, as a == g in associated flow rule.
	do
	{
	checkYield (PointCopy.state, pointCopy.stress, pointCopy.suction(), &fValue); //20 Feb 2006, preparations for the drift correction

		if ((fValue<-d_yieldTol)||(fValue>d_yieldTol)) correct=TRUE; else correct=FALSE;
		if (correct==TRUE)
		{
				// CORRECT FOR DRIFT
				//HERE THE DRIFT WILL BE CORRECTED BY USING THE D MATRIX FROM THE FORBIDDEN SPACE
				//ALTHOUGH BECAUSE THE DRIFT WILL BE CHECKED AGAIN, IT SHOULDN'T POSE MUCH PROBLEM.
				BBMMatrix A (6,1); //dF/dsigma
				BBMMatrix P (1,1); //dF/dP0*
				BBMMatrix GG(6,1);	//dG/dsigma
				BBMMatrix PEP (1,1); //dp0* /depsvpl
				BBMMatrix MM (6,1); // will be vector(1,1,1,0,0,0)T
				BBMMatrix DEL (6,6); //D elastic matrix...
				BBMMatrix DPZEROSTAR (1,1);
				BBMMatrix TEMP (1,1);
				BBMMatrix dSigma (6,1);
				MM.PutElement (1,1,1);
				MM.PutElement (2,1,1);
				MM.PutElement (3,1,1); //rest is zero as initialized


				SpecificVolume=PointCopy.getSpecVol();  //specific volume need to be used the right one
				//cout<<"Specific Volume:"<<SpecificVolume<<"\n";
				PZeroStar=PointCopy.p0Star ();

				//convention - matrices names are made from CAPITALIZED letters


				meanStress=PointCopy.getmeanStress();
				shearStress=PointCopy.shearStress();
				Suction=PointCopy.suction();

				LambdaS=(1-r)*exp(-Beta*Suction)+r;
				LambdaS=LambdaS*LambdaZero;
				Fraction=(LambdaZero-KappaP)/(LambdaS- KappaP);
				PZero=pc*pow(PZeroStar/pc,Fraction);  //get.calculated.pzero;

				//cout<<"PZero = "<<PZero<<"\n";
				//cout<<"PZeroStar = "<<PZeroStar<<"\n";
				//cout<<"p = "<<meanStress<<"\n";
				//cout<<"q = "<<shearStress<<"\n";
				//cout<<"s = "<<Suction<<"\n";

				temp=2*PointCopy.stress[0]-PointCopy.stress[1]-PointCopy.stress[2]+M*M/3*(2*meanStress+k*Suction-PZero);
				A.PutElement(1,1,temp);
				temp=2*PointCopy.stress[1]-PointCopy.stress[0]-PointCopy.stress[2]+M*M/3*(2*meanStress+k*Suction-PZero);
				A.PutElement(2,1,temp);
				temp=2*PointCopy.stress[2]-PointCopy.stress[0]-PointCopy.stress[1]+M*M/3*(2*meanStress+k*Suction-PZero);
				A.PutElement(3,1,temp);
				temp=6*PointCopy.stress[3];
				A.PutElement(4,1,temp);
				temp=6*PointCopy.stress[4];
				A.PutElement(5,1,temp);
				temp=6*PointCopy.stress[5];
				A.PutElement(6,1,temp);
				//cout<<"A:"<<"\n"; A.Print();
				//dF/dsigma - inserted into A

				if (nonAssociated)
				{
					temp=alfa*(2*Point->stress[0]-Point->stress[1]-Point->stress[2])+M*M/3*(2*meanStress+k*Suction-PZero);
					GG.PutElement(1,1,temp);
					temp=alfa*(2*Point->stress[1]-Point->stress[0]-Point->stress[2])+M*M/3*(2*meanStress+k*Suction-PZero);
					GG.PutElement(2,1,temp);
					temp=alfa*(2*Point->stress[2]-Point->stress[0]-Point->stress[1])+M*M/3*(2*meanStress+k*Suction-PZero);
					GG.PutElement(3,1,temp);
					temp=6*alfa*Point->stress[3];
					GG.PutElement(4,1,temp);
					temp=6*alfa*Point->stress[4];
					GG.PutElement(5,1,temp);
					temp=6*alfa*Point->stress[5];
					GG.PutElement(6,1,temp);
				}
				else */ /*A.Copy (&GG);


                                //d

                                temp=0;
                                temp=-M*M*(meanStress+k*Suction)*Fraction*pow((PZeroStar/pc),
                                Fraction-1);

                                P.PutElement (1,1,temp);

                                //cout<<"P:"<<"\n"; P.Print();

                                temp=PZeroStar*SpecificVolume/(LambdaZero-KappaP);
                                PEP.PutElement (1,1,temp); //dP0* /depsplv
                                //cout<<"dpZeroStar/Depsvpl:"<<temp<<"\n";
                                //DEL... elastic matrix... values of K. Here we
                                need values of K at the point...
                                We need to have the K - bulk modulus of the soil
                                calculated, and then it is
                                possible to fill into the DEL Matrix...
                                So, the way of doing it will be repeated,
                                algorithm as before... , in procedure find
                                stress elast, but, this time
                                it will be made inside and the results will be
                                put into the matrix.
                                */

  /*
                                  // checks whether Mean stress is large enough
  to hold K in right range.
                                  if ((meanStress<d_minMeanStress)&&(meanStress>(-d_minMeanStress)))
                                  {
                                          cout<<"WARNING !!! Mean stress too
  low. Mean stress is adjusted to d_minMeanStress value!!!"<<"\n";
                                          meanStress=d_minMeanStress;
                                  }

                                  K=meanStress*SpecificVolume/KappaP;  //tangent
  bulk modulus K=p*specificVol/KappaP, from eq dSpecVol=-KappaP*ln(p/pzero);
                                  //cout<<"K="<<K<<"\n";

                                  // ****************************** K need
  correcting, but with const. suction seems to be ok
  ******************************

                                  // Tiny set in stdafx to 1e-12
                                  // calculate helpful variables:
                                  double K43G=K+4.0*G/3.0;
                                  double K23G=K-2.0*G/3.0;


                                  // Fill in the matrix:

                                  DEL.PutElement (1,1,K43G);
                                  DEL.PutElement (1,2,K23G);
                                  DEL.PutElement (1,3,K23G); //rest of the line
  are zeros and rightly so
                                  DEL.PutElement (2,1,K23G);
                                  DEL.PutElement (2,2,K43G);
                                  DEL.PutElement (2,3,K23G); //yes, the matrix
  is symmetrical, but it is faster to put this 3 additional elements
                                  DEL.PutElement (3,1,K23G); //than just mirror
  all, including zeros, which are there already
                                  DEL.PutElement (3,2,K23G);
                                  DEL.PutElement (3,3,K43G);
                                  DEL.PutElement (4,4,2*G);
                                  DEL.PutElement (5,5,2*G);
                                  DEL.PutElement (6,6,2*G); //rest of the matrix
  is filled with zeros...

                                  A.Transpose (&TEMP);
                                  TEMP.Multiply (&DEL,&TEMP);
                                  TEMP.Multiply (&GG, &TEMP);

                                  temp=TEMP.getElement (1,1);  //first part of
  the denominator
                                  //cout<<"First Part of Denominator done...
  temp="<<temp<<"\n";


                                  P.Multiply (&PEP,&TEMP);
                                  MM.Transpose (&MM);
                                  TEMP.Multiply (&MM,&TEMP);	//MM is
  transposed
                                  TEMP.Multiply (&GG,&TEMP);

                                  temp=temp+TEMP.getElement (1,1); //'end of the
  denominator

                                  //cout<<"Denominator="<<temp<<"\n";

                                  Lambda=fValue*(PZero+Suction*k)*(PZero+Suction*k)/temp;
  //because we need the value, not the reduced value...

                                  //cout<<"Lambda="<<Lambda<<"\n";

                                  A.Multiply (Lambda, &TEMP); //delta epsilon
  plastic= -delta epsilon elastic

                                  //cout<<"Delta Epsilon Plastic:"<<"\n";
                                  //TEMP.Print();

                                  for (int i=1; i<7; i++) epStrain
  [i-1]=TEMP.getElement (i,1); //epsilon pl change
                                  temp=epStrain[0]+epStrain[1]+epStrain[2];
  //DepsilonV
                                  dpZeroStar=PEP.getElement(1,1)*temp;
                                  //cout<<"dpZeroStar="<<*dpZeroStar<<"\n";

                                  DEL.Multiply (&GG, &TEMP);
                                  TEMP.Multiply (-Lambda, &dSigma);
                          //final result for stress change, with the negative
  sign, so the stress should be ADDED
                          //be ADDED to get the right result of corrected
  stresses.

                          //cout<<"Delta Sigma="<<"\n";
                          //dSigma->Print();


                          //cout<<"Press any key (end of Correct Drift
  Procedure)"<<"\n";
                          //getchar();
                          for (int i=0; i<6; i++) DSigma[i]=dSigma.getElement
  (i+1,1);
                          PointCopy.Update (epStrain, zeros, DSigma,
  dpZeroStar);

                  }
          }
          while (correct==TRUE);

  PointCopy.Copy(Point);
                          //this finishes the algorithm */
}

void
ShengMohrCoulomb::correctDriftBeg(StateMohrCoulomb* Point, StateMohrCoulomb* PointOld)
{

  /*
  This procedure should correct the drift, as described in the word file. The
  method of drift correction is based on
  the Potts&Zdravkovic book; It is however slightly changes so to use it in
  unsaturated soil model. Detailed
  description of the algorithm is presented in the word file. It is noted
  however, that in current version the algorithm
  does calculate all the derivatives in the forbidden space outside the yield
  locus which may be cause of concerns.

  Input: point to correct
  Output: Updated point, corrected using values of stress at the beginning
  */

  // Need dF/Ds, D, dF/dP0*, dP0*/DEpsPl, m, dF/dS
  // at the beginning we need to calculate the derivatives a, b, c, d, g, p...)

  // cout<<"Correct Drift Procedure entered!"<<"\n";

  BBMMatrix A(6, 1); // dF/dsigma
  BBMMatrix DENOMINATOR(1, 1);

  BBMMatrix GG(6, 1);  // dG/dsigma, in case of associated flow
                       // (nonAssociated==false) same as A
  BBMMatrix MM(6, 1);  // will be vector(1,1,1,0,0,0)T
  BBMMatrix DEL(6, 6); // D elastic matrix...
  BBMMatrix DEPS(6, 1);
  BBMMatrix dSigma(6, 1);

  double DSigma[6], epStrain[6], zeros[7];
  for (int i = 0; i < 7; i++)
    zeros[i] = 0;

  for (int i = 1; i < 6; i++)
    DEPS.PutElement(i, 1, epStrain[i - 1]); // increase of epsilon copied
  // we do not have g, as a == g in associated flow rule.

  bool correct;

  double I1 = PointOld->firstInvariant();
  double I2 = PointOld->secondInvariant();
  double J2 = PointOld->secondDevInvariant();
  if (fabs(J2) < TINY)
    J2 = TINY;
  double J3 = PointOld->thirdDevInvariant();
  double shearStress = PointOld->shearStress();
  if (fabs(shearStress) < TINY)
    shearStress = TINY;
  double Factor = -27 * J3 / (2 * shearStress * shearStress * shearStress);
  if (Factor > 1)
    Factor = 1;
  if (Factor < -1)
    Factor = -1;
  if (!finite(Factor)) {
    Factor = 1.0;
    cout << "Factor in Correct Drift Beg [parameters *Point] is not finite. "
            "set to 1"
         << "\n";
    cout << "alpha4=" << alpha4 << " J3=" << J3
         << " shearStress=" << shearStress << "\n";
  }
  Factor = 1 + alpha4 - (1 - alpha4) * Factor;
  double Factor025 = pow(Factor, 0.25);
  double Factor075 = pow(Factor, -0.75);
  double M = (3 - sin_phi) / (6 * alpha * sin_phi);
  double Mpsi = (3 - sin_psi) / (6 * alpha * sin_psi);
  double fValue;

  BBMMatrix dJ2_dSig(6, 1), dJ3_dSig(6, 1), dq_dSig(6, 1), dI1_dSig(6, 1),
    dI2_dSig(6, 1), dI3_dSig(6, 1), TEMP(6, 1), NUMERATOR(6, 1);

  // derivatives of the invariants
  double s[6] = {
    PointOld->stress[0], PointOld->stress[1], PointOld->stress[2],
    PointOld->stress[3], PointOld->stress[4], PointOld->stress[5]
  };

  dI1_dSig(0, 0) =  1.0);
  dI1_dSig(1, 0) =  1.0);
  dI1_dSig(2, 0) =  1.0); //{1,1,1,0,0,0}

  dI2_dSig(0, 0) =  s[1] + s[2]);
  dI2_dSig(1, 0) =  s[0] + s[2]);
  dI2_dSig(2, 0) =  s[0] + s[1]);
  dI2_dSig(3, 0) =  -2 * s[3]);
  dI2_dSig(4, 0) =  -2 * s[4]);
  dI2_dSig(5, 0) =  -2 * s[5]);

  // dI2_dSig.PrintPrecise();

  dI3_dSig(0, 0) =  s[1] * s[2] - s[5] * s[5]);
  dI3_dSig(1, 0) =  s[0] * s[2] - s[4] * s[4]);
  dI3_dSig(2, 0) =  s[0] * s[1] - s[3] * s[3]);
  dI3_dSig(3, 0) =  2 * s[5] * s[4] - 2 * s[2] * s[3]);
  dI3_dSig(4, 0) =  2 * s[3] * s[5] - 2 * s[1] * s[4]);
  dI3_dSig(5, 0) =  2 * s[3] * s[4] - 2 * s[0] * s[5]);

  // dI3_dSig.PrintPrecise();

  dJ2_dSig(0, 0) =  (2 * s[0] - s[1] - s[2]) / 3.0);
  dJ2_dSig(1, 0) =  (2 * s[1] - s[0] - s[2]) / 3.0);
  dJ2_dSig(2, 0) =  (2 * s[2] - s[0] - s[1]) / 3.0);
  dJ2_dSig(3, 0) =  2 * s[3]);
  dJ2_dSig(4, 0) =  2 * s[4]);
  dJ2_dSig(5, 0) =  2 * s[5]);

  dJ2_dSig.Copy(&dq_dSig);
  dq_dSig.Multiply(sqrt(3.0) * 0.5 / sqrt(J2), &dq_dSig);

  dI1_dSig.Copy(&dJ3_dSig);
  dJ3_dSig.Multiply(6.0 / 27.0 * I1 * I1 - I2 / 3.0, &dJ3_dSig);
  dI2_dSig.Copy(&TEMP);
  TEMP.Multiply(-I1 / 3.0, &TEMP);
  dJ3_dSig.Add(&TEMP, &dJ3_dSig);
  dJ3_dSig.Add(&dI3_dSig, &dJ3_dSig);
  // finished dJ3_dSig

  // above needs to be *CORRECTED* as some errors in ij i!=j components is
  // likely (2x too big)

  dJ2_dSig.Copy(&A);
  A.Multiply(0.21022410391343 * M * Factor025 / shearStress, &A);
  dJ3_dSig.Copy(&TEMP);
  TEMP.Multiply(-1.419012700740643 * (alpha4 - 1) * M * sqrt(J2) * Factor075 /
                    (sqrt(3.0) * shearStress * shearStress * shearStress) +
                  0.94600846716043 * (alpha4 - 1) * cohesion * M * Factor075 /
                    (shearStress * shearStress * shearStress),
                &TEMP);
  A.Add(&TEMP, &A);
  dq_dSig.Copy(&TEMP);
  TEMP.Multiply(
    -2.838025401481287 * (alpha4 - 1) * cohesion * M * J3 * Factor075 /
        (shearStress * shearStress * shearStress * shearStress) +
      4.257038102221929 * (alpha4 - 1) * M * sqrt(J2) * J3 * Factor075 /
        (sqrt(3.0) * shearStress * shearStress * shearStress * shearStress),
    &TEMP);
  A.Add(&TEMP, &A);
  dI1_dSig.Copy(&TEMP);
  TEMP.Multiply(-1 / 3.0, &TEMP);
  A.Add(&TEMP, &A);

  // FINISHED dF/dSigma

  dJ2_dSig.Copy(&GG);
  GG.Multiply(0.21022410391343 * Mpsi * Factor025 / shearStress, &GG);
  dJ3_dSig.Copy(&TEMP);
  TEMP.Multiply(-1.419012700740643 * (alpha4 - 1) * Mpsi * sqrt(J2) *
                    Factor075 /
                    (sqrt(3.0) * shearStress * shearStress * shearStress) +
                  0.94600846716043 * (alpha4 - 1) * cohesion * Mpsi *
                    Factor075 / (shearStress * shearStress * shearStress),
                &TEMP);
  GG.Add(&TEMP, &GG);
  dq_dSig.Copy(&TEMP);
  TEMP.Multiply(
    -2.838025401481287 * (alpha4 - 1) * cohesion * Mpsi * J3 * Factor075 /
        (shearStress * shearStress * shearStress * shearStress) +
      4.257038102221929 * (alpha4 - 1) * Mpsi * sqrt(J2) * J3 * Factor075 /
        (sqrt(3.0) * shearStress * shearStress * shearStress * shearStress),
    &TEMP);
  GG.Add(&TEMP, &GG);
  dI1_dSig.Copy(&TEMP);
  TEMP.Multiply(
    -1 / 3.0,
    &TEMP); // First correction: sign here, but most likely not enough!!!
  GG.Add(&TEMP, &GG);

  // FINISHED dQ/dSigma

  double K43G = K + 4.0 * G / 3.0;
  double K23G = K - 2.0 * G / 3.0;

  DEL(0, 0) =  K43G);
  DEL(0, 1) =  K23G);
  DEL(0, 2) =  K23G); // rest of the line are zeros and rightly so
  DEL(1, 0) =  K23G);
  DEL(1, 1) =  K43G);
  DEL(1, 2) =  K23G); // yes, the matrix is symmetrical, but it is
                              // faster to put this 3 additional elements
  DEL.PutElement(
    3, 1,
    K23G); // than just mirror all, including zeros, which are there already
  DEL.PutElement(3, 2, K23G);
  DEL.PutElement(3, 3, K43G);
  DEL.PutElement(4, 4, 2.0 * G);
  DEL.PutElement(5, 5, 2.0 * G);
  DEL.PutElement(6, 6, 2.0 * G); // rest of the matrix is filled with zeros...

  // getting lambda and Dep

  A.Transpose(&NUMERATOR);
  NUMERATOR.Multiply(&DEL, &NUMERATOR); // NUMERATOR=aT*Del -->Numerator of
                                        // Lambda without multiplication by
                                        // dEpsilon
  NUMERATOR.Multiply(&GG, &DENOMINATOR);

  int numberIter = 0;

  do {
    numberIter++;
    computeYieldNormalized(Point->state, point->stress, point->suction(),
               &fValue); // 20 Feb 2006, preparations for the drift correction

    if ((fValue / (Point->getmeanStress() + 2 * cohesion) < -d_yieldTol) ||
        (fValue / (Point->getmeanStress() + 2 * cohesion) > d_yieldTol))
      correct = TRUE;
    else
      correct = FALSE;
    if (correct == TRUE) {
      // CORRECT FOR DRIFT

      double Lambda = fValue / DENOMINATOR.getElement(1, 1);
      A.Multiply(Lambda, &TEMP); // delta epsilon plastic= -delta epsilon
                                 // elastic

      // cout<<"Delta Epsilon Plastic:"<<"\n";
      // TEMP.Print();

      for (int i = 1; i < 7; i++)
        epStrain[i - 1] = TEMP.getElement(i, 1); // epsilon pl change
      DEL.Multiply(&GG, &TEMP);
      TEMP.Multiply(-Lambda, &dSigma);
      // final result for stress change, with the negative sign, so the stress
      // should be ADDED
      // be ADDED to get the right result of corrected stresses.
      for (int i = 0; i < 6; i++)
        DSigma[i] = dSigma.getElement(i + 1, 1);
      Point->Update(epStrain, zeros, DSigma, 0);
    }
    if (numberIter > 10) {
      // cout<<"Drift Correction Procedure failed"<<"\n";
      correct = FALSE;
    }
  } while (correct == TRUE);

  /*
StateMohrCoulomb pointCopy, pointEnd;
PointOld->Copy(&PointCopy);
Point->Copy(&PointEnd);

bool correct;
double fValue, dpZeroStar, Lambda, DSigma[6], epStrain[6], zeros[7];
for (int i=0; i<7; i++) zeros[i]=0;

double temp, PZero, meanStress, shearStress, Suction, PZeroStar, SpecificVolume, LambdaS, Fraction;
// we do not have g, as a == g in associated flow rule.
	do
	{
	checkYield (PointEnd.state, pointEnd.stress, pointEnd.suction(), &fValue); //20 Feb 2006, preparations for the drift correction

	if (fabs(fValue)>d_yieldTol) correct=TRUE; else correct=FALSE;
		if (correct==TRUE)
		{
				// CORRECT FOR DRIFT
				//HERE THE DRIFT WILL BE CORRECTED BY USING THE D MATRIX FROM THE FORBIDDEN SPACE
				//ALTHOUGH BECAUSE THE DRIFT WILL BE CHECKED AGAIN, IT SHOULDN'T POSE MUCH PROBLEM.
				BBMMatrix A (6,1); //dF/dsigma
				BBMMatrix P (1,1); //dF/dP0*
				BBMMatrix GG (6,1);
				BBMMatrix PEP (1,1); //dp0* /depsvpl
				BBMMatrix MM (6,1); // will be vector(1,1,1,0,0,0)T
				BBMMatrix DEL (6,6); //D elastic matrix...
				BBMMatrix DPZEROSTAR (1,1);
				BBMMatrix TEMP (1,1);
				BBMMatrix dSigma (6,1);
				MM.PutElement (1,1,1);
				MM.PutElement (2,1,1);
				MM.PutElement (3,1,1); //rest is zero as initialized



				SpecificVolume=PointCopy.getSpecVol();  //specific volume need to be used the right one
				//cout<<"Specific Volume:"<<SpecificVolume<<"\n";
				PZeroStar=PointCopy.p0Star ();

				//convention - matrices names are made from CAPITALIZED letters


				meanStress=PointCopy.getmeanStress();
				shearStress=PointCopy.shearStress();
				Suction=PointCopy.suction();

				LambdaS=(1-r)*exp(-Beta*Suction)+r;
				LambdaS=LambdaS*LambdaZero;
				Fraction=(LambdaZero-KappaP)/(LambdaS- KappaP);
				PZero=pc*pow(PZeroStar/pc,Fraction);  //get.calculated.pzero;

				//cout<<"PZero = "<<PZero<<"\n";
				//cout<<"PZeroStar = "<<PZeroStar<<"\n";
				//cout<<"p = "<<meanStress<<"\n";
				//cout<<"q = "<<shearStress<<"\n";
				//cout<<"s = "<<Suction<<"\n";

				temp=2*PointCopy.stress[0]-PointCopy.stress[1]-PointCopy.stress[2]+M*M/3*(2*meanStress+k*Suction-PZero);
				A.PutElement(1,1,temp);
				temp=2*PointCopy.stress[1]-PointCopy.stress[0]-PointCopy.stress[2]+M*M/3*(2*meanStress+k*Suction-PZero);
				A.PutElement(2,1,temp);
				temp=2*PointCopy.stress[2]-PointCopy.stress[0]-PointCopy.stress[1]+M*M/3*(2*meanStress+k*Suction-PZero);
				A.PutElement(3,1,temp);
				temp=6*PointCopy.stress[3];
				A.PutElement(4,1,temp);
				temp=6*PointCopy.stress[4];
				A.PutElement(5,1,temp);
				temp=6*PointCopy.stress[5];
				A.PutElement(6,1,temp);
				//cout<<"A:"<<"\n"; A.Print();
				//dF/dsigma - inserted into A

				if (nonAssociated)
				{
					temp=alfa*(2*Point->stress[0]-Point->stress[1]-Point->stress[2])+M*M/3*(2*meanStress+k*Suction-PZero);
					GG.PutElement(1,1,temp);
					temp=alfa*(2*Point->stress[1]-Point->stress[0]-Point->stress[2])+M*M/3*(2*meanStress+k*Suction-PZero);
					GG.PutElement(2,1,temp);
					temp=alfa*(2*Point->stress[2]-Point->stress[0]-Point->stress[1])+M*M/3*(2*meanStress+k*Suction-PZero);
					GG.PutElement(3,1,temp);
					temp=6*alfa*Point->stress[3];
					GG.PutElement(4,1,temp);
					temp=6*alfa*Point->stress[4];
					GG.PutElement(5,1,temp);
					temp=6*alfa*Point->stress[5];
					GG.PutElement(6,1,temp);
				}
				else */ /* A.Copy (&GG);

                                //d

                                temp=0;
                                temp=-M*M*(meanStress+k*Suction)*Fraction*pow((PZeroStar/pc),
                                Fraction-1);

                                P.PutElement (1,1,temp);

                                //cout<<"P:"<<"\n"; P.Print();

                                temp=PZeroStar*SpecificVolume/(LambdaZero-KappaP);
                                PEP.PutElement (1,1,temp); //dP0* /depsplv
                                //cout<<"dpZeroStar/Depsvpl:"<<temp<<"\n";
                                //DEL... elastic matrix... values of K. Here we
                                need values of K at the point...
                                We need to have the K - bulk modulus of the soil
                                calculated, and then it is
                                possible to fill into the DEL Matrix...
                                So, the way of doing it will be repeated,
                                algorithm as before... , in procedure find
                                stress elast, but, this time
                                it will be made inside and the results will be
                                put into the matrix.
                                */

  /*
                                  // checks whether Mean stress is large enough
  to hold K in right range.
                                  if ((meanStress<d_minMeanStress)&&(meanStress>(-d_minMeanStress)))
                                  {
                                          cout<<"WARNING !!! Mean stress too
  low. Mean stress is adjusted to d_minMeanStress value!!!"<<"\n";
                                          meanStress=d_minMeanStress;
                                  }

                                  K=meanStress*SpecificVolume/KappaP;  //tangent
  bulk modulus K=p*specificVol/KappaP, from eq dSpecVol=-KappaP*ln(p/pzero);
                                  //cout<<"K="<<K<<"\n";

                                  // ****************************** K need
  correcting, but with const. suction seems to be ok
  ******************************

                                  // Tiny set in stdafx to 1e-12
                                  // calculate helpful variables:
                                  double K43G=K+4.0*G/3.0;
                                  double K23G=K-2.0*G/3.0;


                                  // Fill in the matrix:

                                  DEL.PutElement (1,1,K43G);
                                  DEL.PutElement (1,2,K23G);
                                  DEL.PutElement (1,3,K23G); //rest of the line
  are zeros and rightly so
                                  DEL.PutElement (2,1,K23G);
                                  DEL.PutElement (2,2,K43G);
                                  DEL.PutElement (2,3,K23G); //yes, the matrix
  is symmetrical, but it is faster to put this 3 additional elements
                                  DEL.PutElement (3,1,K23G); //than just mirror
  all, including zeros, which are there already
                                  DEL.PutElement (3,2,K23G);
                                  DEL.PutElement (3,3,K43G);
                                  DEL.PutElement (4,4,2*G);
                                  DEL.PutElement (5,5,2*G);
                                  DEL.PutElement (6,6,2*G); //rest of the matrix
  is filled with zeros...

                                  A.Transpose (&TEMP);
                                  TEMP.Multiply (&DEL,&TEMP);
                                  TEMP.Multiply (&GG, &TEMP);

                                  temp=TEMP.getElement (1,1);  //first part of
  the denominator
                                  //cout<<"First Part of Denominator done...
  temp="<<temp<<"\n";


                                  P.Multiply (&PEP,&TEMP);
                                  MM.Transpose (&MM);
                                  TEMP.Multiply (&MM,&TEMP);	//MM is
  transposed
                                  TEMP.Multiply (&GG,&TEMP);

                                  temp=temp+TEMP.getElement (1,1); //'end of the
  denominator

                                  //cout<<"Denominator="<<temp<<"\n";

                                  Lambda=fValue*(PZero+Suction*k)*(PZero+Suction*k)/temp;
  //because we need the value, not the reduced value...

                                  //cout<<"Lambda="<<Lambda<<"\n";

                                  A.Multiply (Lambda, &TEMP); //delta epsilon
  plastic= -delta epsilon elastic

                                  //cout<<"Delta Epsilon Plastic:"<<"\n";
                                  //TEMP.Print();

                                  for (int i=1; i<7; i++) epStrain
  [i-1]=TEMP.getElement (i,1); //epsilon pl change
                                  temp=epStrain[0]+epStrain[1]+epStrain[2];
  //DepsilonV
                                  dpZeroStar=PEP.getElement(1,1)*temp;
                                  //cout<<"dpZeroStar="<<*dpZeroStar<<"\n";

                                  DEL.Multiply (&GG, &TEMP);
                                  TEMP.Multiply (-Lambda, &dSigma);
                          //final result for stress change, with the negative
  sign, so the stress should be ADDED
                          //be ADDED to get the right result of corrected
  stresses.

                          //cout<<"Delta Sigma="<<"\n";
                          //dSigma->Print();


                          //cout<<"Press any key (end of Correct Drift
  Procedure)"<<"\n";
                          //getchar();
                          for (int i=0; i<6; i++) DSigma[i]=dSigma.getElement
  (i+1,1);
                          PointEnd.Update (epStrain, zeros, DSigma, dpZeroStar);

                  }
          }
          while (correct==TRUE);

  PointEnd.Copy(Point);
  //this finishes the algorithm */
}



void
ShengMohrCoulomb::integrate(double* strainIncrement, double suctionIncrement,
                            StateMohrCoulomb* initialState, double* stressIncrement,
                            double P0StarIncrement,
                            double* plasticStrainIncrement)
{
}

void
ShengMohrCoulomb::integrateConst(double* strainIncrement,
                                 StateMohrCoulomb* initialState, int stepNo, int method)

// method does not work at the moment
{
  StateMohrCoulomb finalState;
  bool purelyElastic;
  bool onTheYieldLocus;
  bool unloading;
  double purelyElasticStrain[7], purelyPlasticStrain[7];

  onTheYieldLocus = computeYieldNormalized(initialState); // checking if on the YL (Returns
                                              // true if on the YL or outside,
                                              // in the plastic part)
  calcElastic(strainIncrement, initialState,
              &finalState); // Initial point copied to the final point, later
                            // Final point updated with the elastic stress
  purelyElastic = !(checkYieldNormalized(&finalState)); // checkYieldNormalized returns true
                                                  // if in the plastic part of
                                                  // the YL. If on the YL
                                                  // (within Tolerance INT_TOL
                                                  // returns false)
  // cout<<purelyElastic<<"\n";

  if (purelyElastic) {
    // finalState.Copy(initialState);
    // update with what we have - all updated after the procedure
    // cout<<"Purely Elastic Step"<<"\n";
  } else // we have elasto-plastic step
  {
    if (onTheYieldLocus) {
      unloading =
        checkGradient(initialState, &finalState); // checking if there is any
                                                  // unloading part; in
                                                  // finalState purely elastic
                                                  // values; we care only for
                                                  // correct change of suction
      // cout<<"\n"<<"\n"<<"unloading="<<unloading<<"\n";
      // getchar();
      if (unloading) {
        cout << "\n"
             << "\n"
             << "Elasto-Plastic unloading=" << unloading << "\n";
        findIntersectionUnloading(strainIncrement, initialState,
                                  purelyElasticStrain,
                                  purelyPlasticStrain); // finding elastic part
        calcElastic(purelyElasticStrain, initialState,
                    &finalState); // calculating elastic part, updated in the
                                  // same point.
        calculatePlasticConst(purelyPlasticStrain, &finalState,
                              stepNo); // and subsequently plastic part
        // finalState.Copy(initialState);
      } else {
        calculatePlasticConst(strainIncrement, &finalState, stepNo);
      }
    } else {
      findIntersection(strainIncrement, initialState, purelyElasticStrain,
                       purelyPlasticStrain); // not on the yield locus, finding
                                             // intersection and subsequently
                                             // calculate elastic and plastic
                                             // stress increment
      calcElastic(purelyElasticStrain, initialState, &finalState);
      calculatePlasticConst(purelyPlasticStrain, &finalState, stepNo);
      // finalState.Copy(initialState);
    }
  }
  // for (int i=0; i<3; i++) cout<<strainIncrement[i]<<"  "<<"\n";
  // for (int i=0; i<3; i++) cout<<finalState.stress[i]<<"  "<<"\n";
  finalState.setSpecificVolume(
    computeNu(finalState.stress, finalState.state, finalState.suction()));
  finalState.Copy(initialState);
}

void
ShengMohrCoulomb::findElStrGradPQ(double nu0, double* s0, double* eps0,
                                  double* deps, double* ds)
{
  /*
  Assumed that into s0 is the vector p,q inserted s0[0]=p. s0[1]=q
  Output - into ds the p, q is written
  */
  /*
          double K, dEpsSV, dEpsSig, dEpsV, dEpsq;

          dEpsq=(deps[0]+deps[1])*(deps[0]+deps[1])+(deps[0]+deps[2])*(deps[0]+deps[2])+(deps[1]+deps[2])*(deps[1]+deps[2]);
          dEpsq=dEpsq+6*(deps[3]*deps[3]+deps[4]*deps[4]+deps[5]*deps[5]);
          dEpsq=sqrt(2*dEpsq);
          dEpsq=dEpsq/3;	//distortional strain computed

          K=s0[0]*nu0/KappaP;
          dEpsV=(deps[0]+deps[1]+deps[2])/3;

          dEpsSV=(KappaS/nu0)*log((eps0[6]+deps[6]+PAtmos)/(eps0[6]+PAtmos));
          dEpsSig=dEpsV-dEpsSV;	//volumetric strain due to change of stress
     computed
          ds[0]=K*dEpsSig;
          ds[1]=G*dEpsq;	//p & q computed  */
  // !!! all stresses computed using tangent K modulus- only for gradient
  // procedures !!!
}

double
ShengMohrCoulomb::calculatePZero(StateMohrCoulomb* point)
{
  /*	double PZero=LambdaZero*((1-r)*exp(-1*Beta*Point->suction())+r);
          PZero=(LambdaZero-KappaP)/(PZero-KappaP);
          PZero=pc*pow((Point->p0Star()/pc),PZero);
          return PZero;*/
  return 0;
}

/*
bool ShengMohrCoulomb::checkYield (double *state, double *s, double suction)
{


        Purpose of this routine is to calculate value of yield function to
determine, whether we have yielding
        or not.

        Arguments: *state - table of state parameters; in fact we need only p*.
                                s- stress state
                                suction - suction
        Returns: FALSE when NOT yielding, TRUE when yielding.
        */
/*	double PZero, SuctionPressure, meanStress, shearStress, fValue;
        //check of the second yield surface:
        if (state[1]-suction<-d_suctionTol) return true;
        //check of the standard yield surface
        meanStress=(s[0]+s[1]+s[2])/3;
        shearStress=(s[0]-s[1])*(s[0]-s[1])+(s[0]-s[2])*(s[0]-s[2])+(s[1]-s[2])*(s[1]-s[2]);
        shearStress=shearStress+6*(s[3]*s[3]+s[4]*s[4]+s[5]*s[5]);
        shearStress=shearStress/2;
        shearStress=sqrt (shearStress);  //Naylor Pande

        PZero=LambdaZero*((1-r)*exp(-1*Beta*suction)+r);
        PZero=(LambdaZero-KappaP)/(PZero-KappaP);
        PZero=pc*pow((state[0]/pc),PZero);
        SuctionPressure=k*suction;
        fValue=shearStress*shearStress-M*M*(meanStress+SuctionPressure)*(PZero-meanStress);
        fValue=fValue/((PZero+SuctionPressure)*(PZero+SuctionPressure)); //value
of Yield function calculated and normalized
        if (fValue<d_yieldTol) return false;
         else return true;

}
*/

void
ShengMohrCoulomb::findYieldOriginal(double* state, double* s0, double* eps0,
                                    double* deps, double* a)
{

  /*
  Main purpose of this procedure is to find the intersection between yield
  surface and the stress vector. We need to know
  at what value of epsilon (strain) we are on the surface - where plastic
  yielding will begin.

  Parameters: state vector, initial stress s0, initial strain eps0, strain
  increment deps, a- to which value of alfa will be written

  Changed - a

  Algorithm - this is a Pegasus algorithm; we are looking for alfa such, that
  value of yield function F==0.
  This alfa is found and written to a.

  */

  double F0, F1, fAlfa, epsini[7], epsfini[7], epsAlfa[7], sini[6], sfini[6],
    sAlfa[6], alfa0, alfa1, alfa;
  alfa0 = 0;
  alfa1 = 1;

  for (int i = 0; i < 6; i++)
    sini[i] = 0;
  for (int i = 0; i < 7; i++) {
    epsini[i] = deps[i] * alfa0;
    epsfini[i] = deps[i] * alfa1;
  }

  // check of the second yield surface:
  if ((eps0[6] + deps[6]) >
      state[1]) // so suction is larger then maximum experienced suction
  {
    alfa0 = (state[1] - eps0[6]) / deps[6];
    for (int i = 0; i < 7; i++) {
      epsini[i] = deps[i] * alfa0;
    }
    calcStressIncElast(state[2], s0, eps0, epsini, sini);
    for (int i = 0; i < 6; i++)
      sini[i] = sini[i] + s0[i];
    computeYieldNormalized(state, sini, eps0[6] + epsini[6], &F0);
    // main F checked
    if (F0 < d_yieldTol) {
      *a = alfa0; // so we have elastic state in other case,
      return;
    } else {
      alfa1 = alfa0; // we have the upper limit of yiels
      alfa0 = 0;
    }
    // the second surface is the most important - and it is violated...
    // we have to find alfa depending on first surface - which is violated first
  }

  // now because the Main yield surface is valid, the whole check is made...
  // however we start with alfa0 different if the suction was greater than
  // maksimum suction...

  for (int i = 0; i < 7; i++) {
    epsini[i] = deps[i] * alfa0;
    epsfini[i] = deps[i] * alfa1;
  }

  calcStressIncElast(state[2], s0, eps0, epsfini, sfini);

  for (int i = 0; i < 6; i++) {
    sfini[i] = sfini[i] + s0[i]; // otherwise we have sfini already calculated
    sini[i] = s0[i];
  }
  computeYieldNormalized(state, sini, eps0[6] + epsini[6], &F0);
  computeYieldNormalized(state, sfini, eps0[6] + epsfini[6], &F1);

  // cout <<"F0="<<F0<<"\n"<<"F1="<<F1<<"\n";

  for (int iter = 0; iter < d_maxIter; iter++) {
    alfa = F0 / (F0 - F1);
    for (int i = 0; i < 7; i++) {
      epsAlfa[i] = alfa0 * deps[i] + alfa * (alfa1 - alfa0) * deps[i];
    }

    calcStressIncElast(state[2], s0, eps0, epsAlfa,
                    sAlfa); // calculated stress increment for current alfa
    // state[2]=computeNu (sAlfa, state, eps0[6]+epsAlfa[6]);
    for (int i = 0; i < 6; i++)
      sAlfa[i] = sAlfa[i] + s0[i]; // update stress
    computeYieldNormalized(state, sAlfa, eps0[6] + epsAlfa[6],
               &fAlfa); // calculated yield function for current alfa
    // cout<<"In iteration "<<iter<<" alfa="<<alfa<<" and F="<<fAlfa<<"\n";

    if ((fAlfa > -d_yieldTol) && (fAlfa < d_yieldTol)) {
      *a = alfa0 +
           alfa *
             (alfa1 - alfa0); // if fAlfa within tolerance, we have the solution
      cout << "Solution in findYieldOriginal procedure was found after " << iter
           << " iterations." << "\n";
      if (iter > 50) {
        getchar();
        cout << "Large number of iterations!!! Solution is however correct... "
                "Press any key..."
             << "\n";
      }
      return;
    }
    if (fAlfa > 0) {
      alfa1 = alfa0 + alfa * (alfa1 - alfa0);
      F1 = fAlfa; // if fAlfa >0 - we are yielding - max alfa set to current
                  // alfa
    } else {
      alfa0 = alfa0 + alfa * (alfa1 - alfa0);
      F0 = fAlfa; // if fAlfa <0 - we are elastic - minimum alfa is set to
                  // current alfa
    }
  }
  *a = -1;
  // if we are here, we must have perforemed to many iterations...
  // cout<<"Error in procedure findYieldOriginal"<<"\n";
  // cout<<"After "<<d_maxIter<<" iterations crossing point not found"<<"\n";
  // cout<<"This is likely to cause incorrect results... Results obtained should
  // not be taken too seriously..."<<"\n";
  *a = -1; // set value to error...
}


void
ShengMohrCoulomb::moveYieldaBit(double* state, double* s0, double* ds,
                                double* eps0, double* deps, double* gradient,
                                double F0)
/*

This function is to move a point a little bit from just outside the yield locus,
just inside the yield locus.
This movement should be made within given tolerance. This - of course - raise
some additional problems...
Especially about tolerance levels...

Algorithm is modified to what it was originally thought.

*/
{
  double ddeps[7], dds[6];
  double ddF, x;
  F0 = F0 + ADDTOLYIELD * d_yieldTol;
  for (int i = 0; i < 7; i++)
    ddeps[i] = -CHGEPSINTOL * eps0[i];
  calcStressIncElast(state[2], s0, eps0, ddeps, dds);
  for (int i = 0; i < 6; i++)
    dds[i] = dds[i] + s0[i];
  computeYieldNormalized(state, dds, eps0[6] + ddeps[6], &ddF);
  ddF = ddF + ADDTOLYIELD * d_yieldTol; // so we added the same to both values...
  // now x from doc file is calculated
  x = (F0 * CHGEPSINTOL) / (F0 - ddF);
  for (int i = 0; i < 7; i++)
    ddeps[i] = -x * eps0[i];
  calcStressIncElast(state[2], s0, eps0, ddeps, dds);
  for (int i = 0; i < 6; i++)
    dds[i] = dds[i] + s0[i];
  computeYieldNormalized(state, dds, eps0[6] + ddeps[6], &ddF);
  if (ddF > 0)
    cout << "Error in procedure moveYieldaBit; Yield Function ddF is greater "
            "then zero...";
  else {
    for (int i = 0; i < 6; i++) {
      deps[i] = deps[i] - ddeps[i];
      eps0[i] = eps0[i] + ddeps[i];
      s0[i] = dds[i];
    }
    deps[6] = deps[6] - ddeps[6]; // generally speaking ddeps <0 and adding
                                  // means substracting...
    eps0[6] = eps0[6] + ddeps[6];
  }
  return;
}

void
ShengMohrCoulomb::read()
{
  // reads data from file "ShengMohrCoulomb.dta"
  ifstream infile("ShengMohrCoulomb.dta", ios_base::in);

  // file opened
  string s;
  int slength = 0, index = 0, line = 0;
  // double temp=0;

  do {

    getline(infile, s, '\n');
    line++;
    // cout << s <<" Prep: In line no:"<<line<<"\n";
    // getchar ();
    if (!infile.good()) {
      cout << "Wrong Data File";
      break;
    }
  } while (s != "***Start of data***");
  // I ignore file until "Start of data";

  // Reading data - 5 material parameters
  double storage[12];
  for (int j = 0; j < 4; j++) {

    getline(infile, s, '\n'); // read whole line
    line++;
    // cout << s <<"In line no:"<<line<<"\n";
    // getchar ();
    bool notcomment = true;
    if (s == "")
      notcomment = false;
    if ((s[0] == '/') && (s[1] == '/'))
      notcomment = false; // check whether not a comment line or empty line

    if (notcomment) {
      slength = s.length(); // get length of line
      index = s.find(";");  // find where is a ; char
      if (index != 0) {
        s.erase(index, slength - index); // delete all after ;
        storage[j] = atof(s.c_str());    // converse to double
      } else
        cout << "No ; in line:" << line << " May cause errors."
             << "\n"; // warn about lack of ; in line
    } else
      j--;
    if (!infile.good())
      break;
  }

  // Moving data from storage to object variables

  setModelParameters(storage[0], storage[1], storage[2],
                     storage[3]); // G, K, cohesion, Friction Angle

  /*	KappaP=storage[1];
          KappaS=storage[2];
          PAtmos=storage[3];
          pc=storage[4];
          k=storage[5];
          r=storage[6];
          Beta=storage[7];
          LambdaZero=storage[8];
          NZero=storage[9];
          M=storage[10];
          if (storage [11]==0) nonAssociated=false; else nonAssociated=true;
  if (nonAssociated)
          {
                  alfa=(M*(M-9)*(M-3))/(9*(6-M)*(1-KappaP/LambdaZero));
                  cout<<"Non associated flow rule used. Value of
  alfa:"<<alfa<<"\n";
          }
  else {
          alfa=1;
          }

          // finished  */
  infile.close(); // close file
  // all done
}

void
ShengMohrCoulomb::write()
{
  /*
          cout<<"Model Parameters:"<<"\n";
          cout<<"Shear Modulus G="<<G<<"\n";
          cout<<"Kappa p="<<KappaP<<"\n";
          cout<<"Kappa s="<<KappaS<<"\n";
          cout<<"Pc="<<pc<<"\n";
          cout<<"k="<<k<<"\n";
          cout<<"r="<<r<<"\n";
          cout<<"Beta="<<Beta<<"\n";
          cout<<"Lambda (0)="<<LambdaZero<<"\n";
          cout<<"N (0)="<<NZero<<"\n";
          cout<<"M="<<M<<" deg"<<"\n";

          cout<<"Computed parameters:"<<"\n";
          cout<<"K="<<K<<"\n";
          cout<<"Ks="<<Ks<<"\n";
          cout<<"Kp="<<Kp<<"\n";*/
}

void
ShengMohrCoulomb::paintLocus(double* state, double suction, int Max)
{
  // FILE *stream;
  // stream = fopen( "yieldshape.dta", "w" );

  /* Reassign "stderr" to "freopen.out": */
  // locus starts from -s0 and finishes at P0.
  // double minP, maxP, PZero, difference;
  /*
  //cout<<Max<<"\n"; // in first line we put how many points we have...
  fprintf( stream, "%d\n", Max );
  minP=-k*suction;
  PZero=LambdaZero*((1-r)*exp(-1*Beta*suction)+r);
  PZero=(LambdaZero-KappaP)/(PZero-KappaP);
  PZero=pc*pow((state[0]/pc),PZero);
  maxP=PZero;
  //cout<<minP<<"\n";	//minimum
  //cout<<maxP<<"\n";	//maximum, used to set the view...
  fprintf( stream, "%f\n",minP );
  fprintf( stream, "%f\n",maxP );
  difference=maxP-minP;	//this is the max difference...
  //cout<<difference<<"\n";
  double p,q;

  for (int i=0; i<Max; i++)
  {
          p=i*difference/(Max-1)+minP;
          q=M*sqrt((p-minP)*(maxP-p));

          fprintf( stream, "%f\n", p);
          fprintf( stream, "%f\n", q);
   //point written to the file
  }
  */
  // fclose( stream ); //closing stream...
}

double
ShengMohrCoulomb::calculatePlasticConst(double* purelyPlasticStrain,
                                        StateMohrCoulomb* point, int stepNo)
{
  double time, stressIncrAbs[7], RelativeError;
  int numberIter = stepNo;
  int methodSteps = 6;

  double methodOrder = 5;
  double A[8][8]; // matrix must be this size to be used in the doRungeutta
                  // method
  bool errorEstimate =
    false; // we give a 4th order solution, not an error estimate

  /*A - matrix with coefficients, B - error estimate, BRes - result
   * coefficients, C - x coefficients. */
  double C[6] = { 0.0, 0.2, 0.3, 0.6, 1, 0.875 };
  double BRes[6] = { 37.0 / 378.0,  0, 250.0 / 621.0,
                     125.0 / 594.0, 0, 512.0 / 1771.0 };
  double B[6] = { 2825.0 / 27648.0, 0,   18575.0 / 48384.0, 13525.0 / 55296.0,
                  277.0 / 14336.0,  0.25 };

  A[0][0] = 0;

  A[1][0] = 0.2;

  A[2][0] = 0.075;
  A[2][1] = 0.225;

  A[3][0] = 0.3;
  A[3][1] = -0.9;
  A[3][2] = 1.2;

  A[4][0] = -11.0 / 54.0;
  A[4][1] = 2.5;
  A[4][2] = -70.0 / 27.0;
  A[4][3] = 35.0 / 27.0;

  A[5][0] = 1631.0 / 55296.0;
  A[5][1] = 175.0 / 512.0;
  A[5][2] = 575.0 / 13824.0;
  A[5][3] = 44275.0 / 110592.0;
  A[5][4] = 253.0 / 4096.0;

  time = doRungeuttaEqualStep(A, B, BRes, C, point, purelyPlasticStrain,
                             stressIncrAbs, &RelativeError, numberIter,
                             methodOrder, methodSteps, errorEstimate);
  return time;
}

double
ShengMohrCoulomb::calculatePlastic(double* purelyPlasticStrain, StateMohrCoulomb* point)
{
  double time = 0, stressIncrAbs[7];
  int numberIter;

  switch (d_solutionAlgorithm) {
    case 1: {
      // calculate using Modified Euler scheme
      time =
        plasticRKME221(Point, purelyPlasticStrain, stressIncrAbs, &numberIter);
    } break;
    case 2: {
      // calculate using 3rd order RK scheme (Nystrom)
      time =
        plasticRK332(Point, purelyPlasticStrain, stressIncrAbs, &numberIter);
    } break;
    case 3: {
      // calculate using 3rd order RK scheme (Bogacki - Shampine)
      time =
        plasticRKBog432(Point, purelyPlasticStrain, stressIncrAbs, &numberIter);
    } break;
    case 4: {
      // calculate using 4th order RK scheme
      time =
        plasticRK543(Point, purelyPlasticStrain, stressIncrAbs, &numberIter);
    } break;
    case 5: {
      // calculate using 5th order RK scheme (England)
      time =
        plasticRKEng654(Point, purelyPlasticStrain, stressIncrAbs, &numberIter);
    } break;
    case 6: {
      // calculate using 5th order RK scheme (Cash - Karp)
      time =
        plasticRKCK654(Point, purelyPlasticStrain, stressIncrAbs, &numberIter);
    } break;
    case 7: {
      // calculate using 5th order RK scheme (Dormand - Prince)
      time =
        plasticRKDP754(Point, purelyPlasticStrain, stressIncrAbs, &numberIter);
    } break;
    case 8: {
      // calculate using 5th order RK scheme (Bogacki - Shampine)
      time = plasticRKErr8544(Point, purelyPlasticStrain, stressIncrAbs,
                              &numberIter);
    } break;
    case 9: {
      // calculate using extrapolation method (Bulirsch - Stoer)
      time =
        plasticExtrapol(Point, purelyPlasticStrain, stressIncrAbs, &numberIter);
    } break;
    default: {
      cout
        << "Unknown Solution Algorithm. Value of d_solutionAlgorithm variable "
           "is set to:"
        << d_solutionAlgorithm << "\n";
      cout
        << "Acceptable values are ints from 1 to 9. Procedure calculate "
           "Plastic exits without any calculations done. Please press any key "
           "to continue"
        << "\n";
      getchar();
    } break;
  }

  return time;
}

double
ShengMohrCoulomb::plasticEuler(StateMohrCoulomb* point, double* epStrain,
                               double* absStress, int numberIterations)
{

  BBMMatrix dSigma(6, 1);
  double dpZeroStar;
  double DSigma[6];
  double CurrentStrain[7], plasticStrain[6];
  double PZeroStarTot = 0;
  double dSigmaDrift[6], dLambdaDrift;
  StateMohrCoulomb OldPoint;
  Point->Copy(&OldPoint);

  vector<double> stressInc;
  vector<double> strainInc;
  vector<double> PZeroLambda;
  vector<double>::iterator Iter;

  // cout<<"Euler Procedure Entered"<<"\n";
  // cout<<"Number of iterations in Euler Algorithm:"<<numberIterations<<"\n";

  clock_t StartTime, EndTime;
  StartTime = clock();

  for (int i = 0; i < 7; i++)
    absStress[i] = 0;
  for (int i = 0; i < 7; i++)
    CurrentStrain[i] = epStrain[i] / numberIterations;
  for (int loop = 0; loop < numberIterations; loop++) {
    double fValue = 0;
    if (USE_NICE_SCHEME > 0)
      fValue = computeYieldFunction(Point);
    calcPlastic(*Point, CurrentStrain, &dSigma, plasticStrain, &dpZeroStar,
                fValue, dSigmaDrift, &dLambdaDrift);
    for (int i = 0; i < 6; i++)
      DSigma[i] = dSigma.getElement(i + 1, 1);
    Point->Update(plasticStrain, CurrentStrain, DSigma, dpZeroStar);

    PZeroLambda.push_back(dpZeroStar);
    PZeroLambda.push_back(plasticStrain[0]);
    PZeroStarTot = PZeroStarTot + dpZeroStar;

    for (int i = 0; i < 6; i++)
      absStress[i] = absStress[i] + fabs(DSigma[i]);
    absStress[6] = absStress[6] + fabs(dpZeroStar);

    /*for (int i=0; i<6; i++)
            {
            stressInc.push_back (DSigma[i]);  //is ok, as updated total stress?
            strainInc.push_back (CurrentStrain [i]);
            }
    strainInc.push_back (CurrentStrain [6]);*/
  }

  EndTime = clock();

  /*
  cout<<"calculation took:"<<double(EndTime-StartTime)/CLOCKS_PER_SEC<<"
  s."<<"\n";
  cout<<"The d_integrationTol parameter is equal to:"<<d_integrationTol<<"\n";
  cout<<"Total number of steps done:"<<numberIterations<<"\n";

  cout<<"Total PZeroStar change: "<<PZeroStarTot<<"\n";
  checkYield (Point->state, point->stress, point->strain[6],&fValue);
  cout<<"Yield Function value is:"<<fValue<<"\n";

  cout<<"Over the whole step change of stress is:"<<"\n";
  for (int i=0; i<6; i++)
  {
          cout<<"s["<<i<<"]="<<Point->stress[i]-OldPoint.stress[i]<<"\n";
  }

  cout<<"Initial specific volume="<<OldPoint.getSpecVol()<<"\n";
  cout<<"Final specific volume="<<Point->getSpecVol()<<"\n";
  cout<<"Change of specific volume is equal
  to="<<(OldPoint.getSpecVol()-Point->getSpecVol())<<"\n";
  cout<<"Change of mean stress p is equal
  to:"<<Point->getmeanStress()-OldPoint.getmeanStress()<<"\n";
  cout<<"Change of shear stress q is equal
  to:"<<Point->shearStress()-OldPoint.shearStress()<<"\n";



  FILE * StressFile;
  FILE * StrainFile;
  FILE * PZeroFile;
  StressFile = fopen( "stressIncEuler.dta", "w" );

  for (Iter=stressInc.begin() ; Iter!=stressInc.end(); )
  {
  for (int i=0; i<6; i++)
          {
          fprintf( StressFile, "%.20f , ",*Iter);
          Iter++;
          }
  fprintf( StressFile, "\n");
  }
  fclose (StressFile);

  StrainFile = fopen( "strainIncEuler.dta", "w" );

  for (Iter=strainInc.begin() ; Iter!=strainInc.end(); )
  {
  for (int i=0; i<7; i++)
          {
          fprintf( StrainFile, "%.20f , ",*Iter);
          Iter++;
          }
  fprintf( StrainFile, "\n");
  }
  fclose (StrainFile);


  PZeroFile = fopen( "PZeroIncEuler.dta", "w" );
  for (Iter=PZeroLambda.begin() ; Iter!=PZeroLambda.end(); )
  {
  for (int i=0; i<2; i++)
          {
          fprintf( PZeroFile, "%.25f , ",*Iter);
          Iter++;
          }
  fprintf( PZeroFile, "\n");
  }
  fclose (PZeroFile);



  cout<<"Press any key..."<<"\n";
  getchar();
  */

  return (EndTime - StartTime);
}

double
ShengMohrCoulomb::doRungeutta(double A[][8], double* B, double* BRes, double* C,
                             StateMohrCoulomb* point, double* epStrain,
                             double* absStress, int* numberIter,
                             double methodOrder, int methodSteps,
                             bool errorEstimate)
/*
This procedure calculate any Runge - Kutta pair, given the coefficient of
stress in A for each x used to calculate values in C, where B gives the
coefficient to calculate
error estimate (or lower order solution) and BRes to calculate result. The
procedure will integrate whole step of strain
re-using the initial derivative for rejected steps.
methodOrder contain order of the method (required for substep prediction)
methodSteps contain the number of stages to get the result
errorEstimate if true - the B table contains error estimate. If false, it
contains 4th order solution
*/
{
  BBMMatrix dSigma(6, 1);
  StateMohrCoulomb MidPoint[8], OldPoint, TrialPoint;

  double CriticalStepSize = 1E-4; // HACK HACK HACK - prevents algorithm to have
                                  // more than 1e4 steps; cost: accuracy
                                  // reduction in some unusual cases, but the
                                  // whole thing will keep on going

  double DSigma[8][6], DSigmaTemp[7], Result[7];
  double dpZeroStar[8], dpZeroStarTemp, plasticStrainTemp[6],
    plasticStrain[8][6];
  double Error[7], ReUseRes[7], RError, NewStepSize, TotalSize, StepLength,
    Temp, methodPower;
  double dSigmaDrift[6], dSigmadriftConst[6], dSigmaDriftOldConst[6],
    dLambdaDrift;

  double Frequency = 100000 / methodOrder; // how often display info about steps
  bool ReUseStep = false;

  StepLength = 1;
  TotalSize = 0;
  NewStepSize = 1;
  methodPower = pow(2.0, methodOrder) * d_integrationTol;

  for (int i = 0; i < methodOrder; i++) {
    for (int j = 0; j < 6; j++)
      DSigma[i][j] = 0;
    dpZeroStar[i] = 0;
    dSigmaDrift[i] = 0;
  }

  bool Finished = false, StepAccepted = false;
  double SubstepStrain[7], CurrentStrain[7];
  double MicroStep = 0;
  double StepAccuracycheck;

  for (int i = 0; i < 7; i++) {
    SubstepStrain[i] = epStrain[i];
    absStress[i] = 0;
  }
  int stepNo = 0;
  vector<double> stressInc;
  vector<double> strainInc;
  vector<double> PZeroLambda;
  vector<double>::iterator Iter;

  clock_t StartTime, EndTime;
  StartTime = clock();

  NewStepSize = 0;
  for (int i = 0; i < 6; i++)
    NewStepSize = NewStepSize + fabs(SubstepStrain[i]);
  NewStepSize = 0.01 / (NewStepSize * methodOrder);
  // cout<<"NewStepSize="<<NewStepSize<<"\n";
  if (NewStepSize > 1)
    NewStepSize = 1;

  do {
    StepAccepted = false;
    stepNo++;

    /*if (stepNo>1)
{
	if (NewStepSize>10) NewStepSize=10;
	if (NewStepSize<0.1) NewStepSize=0.1;
}*/ // limiting step increase/decrease does not improve results/ enables faster
    // convergence...

    StepLength = StepLength * NewStepSize; // size of a step
    if (StepLength < CriticalStepSize)
      StepLength = CriticalStepSize; // HACK here, it is a fairly large value
    if ((StepLength + TotalSize) > 1)
      StepLength =
        1 - TotalSize; // check whether the step not exceed the whole increment

    for (int i = 0; i < 7; i++) {
      SubstepStrain[i] =
        StepLength * epStrain[i]; // strain increment in current step
      CurrentStrain[i] = 0;
    }
    RError = 0;

    // cout<<"Step Length="<<StepLength<<"\n";
    // cout<<"Current strain [0]="<<CurrentStrain[0]<<"\n";

    for (int i = 0; i < methodSteps; i++)
      point->Copy(&MidPoint[i]); // point is unchanged in  procedure

    // Below the main R-K loop to calculate the value of intermediate stresses;
    // values stored in DSigmaTemp[][]
    // ReUseStep=false;
    if (ReUseStep) {
      for (int i = 0; i < 6; i++) {
        DSigma[0][i] = ReUseRes[i];
        plasticStrain[0][i] = plasticStrainTemp[i];
      }
      dpZeroStar[0] = ReUseRes[6];
      // add line about plastic strain...
    } else {
      switch (ALGORITHM_TYPE) {
        case 0: {
          calcPlastic(MidPoint[0], SubstepStrain, &dSigma, plasticStrainTemp,
                      &dpZeroStar[0], computeYieldFunctionNN(Point),
                      dSigmaDrift, &dLambdaDrift);

          /*	BBMMatrix DEPTEST (6,7), SIGTEST (6,1), EPSTEST (7,1);
                  for (int i=1; i<8; i++) EPSTEST.PutElement (i,1,SubstepStrain
             [i-1]);
                  calculateElastoPlasticTangentMatrix (&MidPoint[0],&DEPTEST);
                  DEPTEST.Multiply (&EPSTEST,&SIGTEST);
                  SIGTEST.Print ();
                  dSigma.Print();
                  getchar (); */
        } break;

        case 1:
          calcPlasticFaster(MidPoint[0], SubstepStrain, &dSigma,
                            plasticStrainTemp, &dpZeroStar[0],
                            computeYieldFunctionNN(Point), dSigmaDrift,
                            &dLambdaDrift);
          break;

        case 2:
          calcPlasticPQ(MidPoint[0], SubstepStrain, &dSigma, plasticStrainTemp,
                        &dpZeroStar[0], computeYieldFunctionNN(Point),
                        dSigmaDrift, &dLambdaDrift);
          /*dSigma.Print();
          cout<<"dpZeroStar="<<dpZeroStar[0]<<"\n";
          calcPlastic (MidPoint[0], SubstepStrain, &dSigma, plasticStrainTemp,
          &dpZeroStar[0]);
          dSigma.Print();
          cout<<"dpZeroStar="<<dpZeroStar[0]<<"\n"; */

          break;

        default:
          calcPlasticFaster(MidPoint[0], SubstepStrain, &dSigma,
                            plasticStrainTemp, &dpZeroStar[0],
                            computeYieldFunctionNN(Point), dSigmaDrift,
                            &dLambdaDrift);
          break;
      }

      for (int i = 0; i < 6; i++) {
        DSigma[0][i] = dSigma.getElement(i + 1, 1);
        plasticStrain[0][i] = plasticStrainTemp[i];
      }
    }

    if (USE_NICE_SCHEME == 3)
      for (int i = 0; i < 6; i++)
        dSigmadriftConst[i] = dSigmaDrift[i];
    else
      for (int i = 0; i < 6; i++)
        dSigmadriftConst[i] = 0.0;
    if (USE_NICE_SCHEME == 4)
      for (int i = 0; i < 6; i++) {
        dSigmadriftConst[i] = dSigmaDrift[i] + dSigmaDriftOldConst[i];
        dSigmaDriftOldConst[i] = dSigmaDrift[i]; //+0.5*dSigmaDriftOldConst[i];
      }

    for (int rkloop = 1; rkloop < methodSteps; rkloop++) {
      for (int i = 0; i < 6; i++) {
        DSigmaTemp[i] = 0;
        plasticStrainTemp[i] = 0;
      }
      dpZeroStarTemp = 0;
      for (int i = 0; i < 7; i++)
        CurrentStrain[i] =
          C[rkloop] *
          SubstepStrain[i]; // set the beginning point of the procedure
      for (int i = 0; i < rkloop; i++) {
        for (int j = 0; j < 6; j++) {
          DSigmaTemp[j] = DSigmaTemp[j] + A[rkloop][i] * DSigma[i][j];
          plasticStrainTemp[j] =
            plasticStrainTemp[j] + A[rkloop][i] * plasticStrain[i][j];
        }
        dpZeroStarTemp = dpZeroStarTemp + A[rkloop][i] * dpZeroStar[i];
      }
      MidPoint[rkloop].Update(plasticStrainTemp, CurrentStrain, DSigmaTemp,
                              dpZeroStarTemp);
      // double dummy;
      // StateMohrCoulomb TempPoint;
      // BBMMatrix TEMPMATRIX (6,7), TEMPEPSILON(7,1);
      // SubstepStrain[6]=0;
      // for (int i=1; i<8; i++) TEMPEPSILON.PutElement(i,1,SubstepStrain[i-1]);
      // MidPoint[rkloop].Copy(&TempPoint);
      double YieldFunctionValue;

      switch (USE_NICE_SCHEME) {
        case 0:
          break;

        case 1: {
          // nice scheme with yield function correction
          YieldFunctionValue = computeYieldFunctionNN(Point);
          dLambdaDrift = 0;
        } break;

        case 2: {
          // nice scheme with carrying lambda
          YieldFunctionValue = 0;
        } break;

        case 3: {
          // nice scheme with carrying stresses only
          YieldFunctionValue = 0;
          dLambdaDrift = 0;
        } break;

        case 4: {
          // predictive nice scheme
          cout
            << "WARNING: Predictive NICE scheme is not accurate and may lead "
               "to large errors. Please use for debug purposes only"
            << "\n";
          YieldFunctionValue = 0;
          dLambdaDrift = 0;
        } break;

        default: {
          YieldFunctionValue = computeYieldFunctionNN(Point);
          dLambdaDrift = 0;
        }
      }

      switch (ALGORITHM_TYPE) {
        case 0:
          calcPlastic(MidPoint[rkloop], SubstepStrain, &dSigma,
                      plasticStrainTemp, &dpZeroStar[rkloop],
                      YieldFunctionValue, dSigmaDrift, &dLambdaDrift);
          break;

        case 1:
          calcPlasticFaster(MidPoint[rkloop], SubstepStrain, &dSigma,
                            plasticStrainTemp, &dpZeroStar[rkloop],
                            YieldFunctionValue, dSigmaDrift, &dLambdaDrift);
          break;

        case 2:
          calcPlasticPQ(MidPoint[rkloop], SubstepStrain, &dSigma,
                        plasticStrainTemp, &dpZeroStar[rkloop],
                        YieldFunctionValue, dSigmaDrift, &dLambdaDrift);
          break;

        default:
          calcPlastic(MidPoint[rkloop], SubstepStrain, &dSigma,
                      plasticStrainTemp, &dpZeroStar[rkloop],
                      YieldFunctionValue, dSigmaDrift, &dLambdaDrift);
          break;
      }

      // calcPlastic (MidPoint[rkloop], SubstepStrain, &dSigma,
      // plasticStrainTemp, &dpZeroStar[rkloop]);
      // calculateElastoPlasticTangentMatrix(&TempPoint,&TEMPMATRIX);
      // TEMPMATRIX.Multiply(&TEMPEPSILON,&TEMPMATRIX);
      /*	for (int i=1; i<7; i++)
              {
                      dummy=fabs(dSigma.getElement(i,1)-TEMPMATRIX.getElement(i,1));
                      if (fabs(dSigma.getElement(i,1))>TINY)
         dummy=dummy/fabs(dSigma.getElement(i,1));
                      if (dummy>0.01)
                      {
                              cout<<"Problems with the alternative
         matrix."<<i<<" Change of stress is:"<<dSigma.getElement(i,1)<<"\n";
                              cout<<"calculated change of stress with the DEP
         matrix is:"<<TEMPMATRIX.getElement(i,1)<<"\n";
                              getchar();
                      }
                      else cout<<"*";
              }


      */
      for (int i = 0; i < 6; i++) {
        DSigma[rkloop][i] = dSigma.getElement(i + 1, 1) + dSigmadriftConst[i];
        plasticStrain[rkloop][i] = plasticStrainTemp[i];
      }
    }

    // needed: result, error

    for (int i = 0; i < 6; i++) {
      Result[i] = 0;
      plasticStrainTemp[i] = 0;
      for (int j = 0; j < methodSteps; j++) {
        Result[i] = Result[i] + BRes[j] * DSigma[j][i];
        plasticStrainTemp[i] =
          plasticStrainTemp[i] + BRes[j] * plasticStrain[j][i];
      }
    }
    Result[6] = 0;
    for (int j = 0; j < methodSteps; j++)
      Result[6] = Result[6] + BRes[j] * dpZeroStar[j];

    for (int i = 0; i < 7; i++)
      Error[i] = 0;

    for (int i = 0; i < methodSteps; i++) {
      for (int j = 0; j < 6; j++)
        Error[j] = Error[j] + B[i] * DSigma[i][j];
      Error[6] = Error[6] + B[i] * dpZeroStar[i];
    }
    if (!errorEstimate)
      for (int i = 0; i < 7; i++)
        Error[i] = Error[i] - Result[i]; // error estimate calculated in case we
                                         // have lower order solution instead of
                                         // error estimate

    // check the error norm

    switch (int(d_tolMethod)) {
      case 0: {
        // RError=checkNorm (DSigmaTemp, dpZeroStarTemp, point, Error);
        // //returns RError
        RError = checkNorm(Result, Result[6], point, Error); // returns RError
      } break;
      case 1: {
        // SLOAN NORM
        // RError=checkNormSloan (DSigmaTemp, dpZeroStarTemp, point, Error);
        // //returns RError
        RError = checkNormSloan(Result, Result[6], point, Error); // returns
                                                                  // RError
      } break;
      default: {
        cout << "ERROR !!!! Improper d_tolMethod in increment.dta" << "\n";
        getchar();
      }
    }

    for (int i = 0; i < 7; i++)
      if (!finite(Result[i])) {
        cout << "Results not a number. Correcting issue, but results may be "
                "incorrect..."
             << "\n";
        Result[i] = 0;
        RError = TINY;
        // if (RError<methodPower) RError=methodPower;
      }

    // if ((Point->getmeanStress()+(Result[0]+Result[1]+Result[2]))/3<0) if
    // (RError<methodPower) RError=methodPower;
    // if ((Point->p0Star()+Result[6])<0) if (RError<methodPower)
    // RError=methodPower;

    if (RError < d_integrationTol) {
      StepAccepted = true;
    } else {
      StepAccepted = false;
      if (StepLength <= CriticalStepSize) {
        StepAccepted = true;
        // StepLength=1E-10;
        // HACK here - you need to adjust all the check in the procedure to the
        // same value
      }
    }

    if (RError < TINY)
      RError = TINY;
    if (d_tolMethod == 0)
      NewStepSize =
        d_betaFactor * pow(d_integrationTol / RError, (1 / (methodOrder - 1.0)));
    else
      NewStepSize =
        d_betaFactor * pow(d_integrationTol / RError, (1 / methodOrder));

    // cout<<d_betaFactor;
    // cout<<"What is going on????"<<"\n";

    if (!StepAccepted) {
      // here we should take care about correct re - usage of the first
      // evaluation of derivative
      ReUseStep = true;
      for (int i = 0; i < 6; i++) {
        ReUseRes[i] = DSigma[0][i] * NewStepSize;
        plasticStrainTemp[i] = plasticStrain[0][i] * NewStepSize;
      }
      ReUseRes[6] = dpZeroStar[0] * NewStepSize;
    } else {
      // here we need to update all the point data.
      Point->Copy(&OldPoint);
      Point->Update(plasticStrainTemp, SubstepStrain, Result, Result[6]);
      Point->Copy(&TrialPoint);
      // this value is not used in any calculations, it is just to show at the
      // end of the step the p0* increase

      if (d_driftCorrection == 3)
        correctDrift(Point);
      if (d_driftCorrection == 2)
        correctDriftBeg(
          Point,
          &OldPoint); // value of OldPoint copied before updating the point

      // re - evaluate the error in the point:

      for (int i = 0; i < 6; i++) {
        Error[i] = Error[i] + Point->stress[i] - TrialPoint.stress[i];
      }
      Error[6] = Error[6] + Point->p0Star() - TrialPoint.p0Star();
      // error vector updated, norm should be re-evaluated:

      switch (int(d_tolMethod)) {
        case 0: {
          // Temp=checkNorm (DSigmaTemp, dpZeroStarTemp, point, Error);
          // //returns RError
          Temp = checkNorm(Result, Result[6], point, Error); // returns RError
          // if (Temp>RError) RError=Temp;
        } break;
        case 1: {
          // SLOAN NORM
          // Temp=checkNormSloan (DSigmaTemp, dpZeroStarTemp, point, Error);
          // //returns RError
          Temp = checkNormSloan(Result, Result[6], point, Error); // returns
                                                                  // RError
          // if (Temp>RError) RError=Temp;
        } break;
        default: {
          cout << "ERROR !!!! Improper d_tolMethod in increment.dta" << "\n";
          getchar();
        }
      }
      if (!finite(RError))
        if (RError < methodPower)
          RError = methodPower;

      if (d_tolMethod == 0)
        NewStepSize =
          d_betaFactor * pow(d_integrationTol / RError, (1 / (methodOrder - 1.0)));
      else
        NewStepSize =
          d_betaFactor * pow(d_integrationTol / RError, (1 / methodOrder));

      for (int i = 0; i < 6; i++)
        dSigmaDriftOldConst[i] = dSigmaDriftOldConst[i] * NewStepSize;
      // cout<<dSigmaDriftOldConst[0]<<"\n"; getchar();

      if (RError < d_integrationTol)
        StepAccepted = true;
      else {
        StepAccepted = false;
        if (StepLength <= CriticalStepSize) {
          StepAccepted = true;
          // StepLength=1E-6;
        }
        // HACK HERE: this is relatively large value
      }

      if (!StepAccepted) {
        ReUseStep = true;
        for (int i = 0; i < 6; i++) {
          ReUseRes[i] = DSigma[0][i] * NewStepSize;
          plasticStrainTemp[i] = plasticStrain[0][i] * NewStepSize;
        }
        ReUseRes[6] = dpZeroStar[0] * NewStepSize;
        OldPoint.Copy(Point);
        Temp = double(stepNo) / Frequency;
        if (modf(Temp, &Temp) == 0) {
          cout << "Step number:" << stepNo << "\n";
          cout << "Total size done is: " << TotalSize
               << " of whole step. Current StepLength=" << StepLength << "\n";
          cout << "Stress state" << Point->stress[0] << ' ' << Point->stress[1]
               << ' ' << Point->stress[2] << ' ' << Point->stress[3] << ' '
               << Point->stress[4] << ' ' << Point->stress[5] << ' ' << "\n";
        }
      } else {
        // this may be done only after successful re-evaluation of the substep
        for (int i = 0; i < 6; i++)
          absStress[i] =
            absStress[i] + fabs(Point->stress[i] - OldPoint.stress[i]);
        absStress[6] =
          absStress[6] + fabs(Point->p0Star() - OldPoint.p0Star());
        ReUseStep = false;
        StepAccuracycheck = TotalSize;
        MicroStep = MicroStep + StepLength;
        TotalSize = TotalSize + MicroStep; // total part of step done updated
        MicroStep = MicroStep - (TotalSize - StepAccuracycheck);
        if (TotalSize >= 1)
          Finished = true;
        Temp = double(stepNo) / Frequency;
        if (modf(Temp, &Temp) == 0) {
          cout << "Step number:" << stepNo << "\n";
          cout << "Total size done is: " << TotalSize
               << " of whole step. Current StepLength=" << StepLength << "\n";
          cout << "Stress state" << Point->stress[0] << ' ' << Point->stress[1]
               << ' ' << Point->stress[2] << ' ' << Point->stress[3] << ' '
               << Point->stress[4] << ' ' << Point->stress[5] << ' ' << "\n";
        }
        /*
                //Debug
                strainInc.push_back (StepLength);
                for (int i=0; i<6; i++) stressInc.push_back (Point->stress[i]);
                stressInc.push_back (Point->p0Star());
                //End Debug */
      }
    }

  } while (!Finished);
  EndTime = clock();
  *numberIter = stepNo;
  /*
  //debug
  FILE * ResultsFile;
  ResultsFile = fopen( "ResultsRK.dta", "w" );

  for (Iter=stressInc.begin();Iter!=stressInc.end();)
  {
          for (int i=0; i<7; i++)
                          {
                          fprintf( ResultsFile, "%.15g , ",*Iter);
                          Iter++;
                          }
          fprintf( ResultsFile, "\n");
  }
  fclose (ResultsFile);

  ResultsFile = fopen( "StrainRK.dta", "w" );
  for (Iter=strainInc.begin();Iter!=strainInc.end();)
  {
          fprintf( ResultsFile, "%.15g \n",*Iter);
          Iter++;
  }
  fclose (ResultsFile);
  //end debug */

  return (EndTime - StartTime);
}

double
ShengMohrCoulomb::doRungeuttaEqualStep(double A[][8], double* B, double* BRes,
                                      double* C, StateMohrCoulomb* point,
                                      double* epStrain, double* absStress,
                                      double* RelError, int numberIter,
                                      double methodOrder, int methodSteps,
                                      bool errorEstimate)
/*
This procedure calculate any Runge - Kutta pair, given the coefficient of
stress in A for each x used to calculate values in C, where B gives the
coefficient to calculate
error estimate (or lower order solution) and BRes to calculate result. The
procedure will integrate whole step of strain
re-using the initial derivative for rejected steps.
methodOrder contain order of the method (required for substep prediction)
methodSteps contain the number of stages to get the result
errorEstimate if true - the B table contains error estimate. If false, it
contains 4th order solution
*/

{
  BBMMatrix dSigma(6, 1);
  StateMohrCoulomb MidPoint[8], OldPoint, TrialPoint;

  double DSigma[8][6], DSigmaTemp[7], Result[7];
  double dpZeroStar[8], dpZeroStarTemp, plasticStrain[8][6],
    plasticStrainTemp[6];
  double Error[7], RError, TotRError, TotalSize, StepLength, Temp, methodPower;
  double Frequency = 15000 / methodOrder; // how often display info about steps

  for (int i = 0; i < methodOrder; i++) {
    for (int j = 0; j < 6; j++)
      DSigma[i][j] = 0;
    dpZeroStar[i] = 0;
  }

  // bool Finished=false, StepAccepted=false;
  double SubstepStrain[7], CurrentStrain[7];
  double MicroStep = 0;
  double StepAccuracycheck;
  double dStressDrift[6], dLambdaDrift;

  int stepNo = 0;
  vector<double> stressInc;
  vector<double> strainInc;
  vector<double> PZeroLambda;
  vector<double>::iterator Iter;

  StepLength = 1.0 / numberIter;

  for (int i = 0; i < 7; i++) {
    CurrentStrain[i] = 0;
    SubstepStrain[i] =
      StepLength * epStrain[i]; // strain increment in all steps (equally sized)
    absStress[i] = 0;
  }
  TotalSize = 0; // NewStepSize=1;
  methodPower = pow(2.0, methodOrder) * d_integrationTol;
  // StepAccepted=true;
  RError = 0;
  TotRError = 0;

  for (int loop = 0; loop < numberIter; loop++) {

    stepNo++;

    /*if (stepNo>1)
{
	if (NewStepSize>10) NewStepSize=10;
	if (NewStepSize<0.1) NewStepSize=0.1;
}*/ // limiting step increase/decrease does not improve results/ enables faster
    // convergence...

    // cout<<"Step Length="<<StepLength<<"\n";
    // cout<<"Current strain [0]="<<CurrentStrain[0]<<"\n";

    for (int i = 0; i < methodSteps; i++)
      point->Copy(&MidPoint[i]); // point is unchanged in calcPlastic procedure

    // Below the main R-K loop to calculate the value of intermediate stresses;
    // values stored in DSigmaTemp[][]
    // ReUseStep=false;

    for (int rkloop = 0; rkloop < methodSteps; rkloop++) {
      for (int i = 0; i < 6; i++) {
        DSigmaTemp[i] = 0;
        plasticStrainTemp[i] = 0;
      }
      dpZeroStarTemp = 0;
      for (int i = 0; i < 7; i++)
        CurrentStrain[i] =
          C[rkloop] *
          SubstepStrain[i]; // set the beginning point of the procedure
      for (int i = 0; i < rkloop; i++) {
        for (int j = 0; j < 6; j++) {
          DSigmaTemp[j] = DSigmaTemp[j] + A[rkloop][i] * DSigma[i][j];
          plasticStrainTemp[j] =
            plasticStrainTemp[j] + A[rkloop][i] * plasticStrain[i][j];
        }
        dpZeroStarTemp = dpZeroStarTemp + A[rkloop][i] * dpZeroStar[i];
      }
      MidPoint[rkloop].Update(plasticStrainTemp, CurrentStrain, DSigmaTemp,
                              dpZeroStarTemp);
      calcPlastic(MidPoint[rkloop], SubstepStrain, &dSigma, plasticStrainTemp,
                  &dpZeroStar[rkloop], computeYieldFunctionNN(Point),
                  dStressDrift, &dLambdaDrift);
      for (int i = 0; i < 6; i++) {
        DSigma[rkloop][i] = dSigma.getElement(i + 1, 1);
        plasticStrain[rkloop][i] = plasticStrainTemp[i];
      }
    }

    // needed: result, error

    for (int i = 0; i < 6; i++) {
      Result[i] = 0;
      plasticStrainTemp[i] = 0;
      for (int j = 0; j < methodSteps; j++) {
        Result[i] = Result[i] + BRes[j] * DSigma[j][i];
        plasticStrainTemp[i] =
          plasticStrainTemp[i] + BRes[j] * plasticStrain[j][i];
      }
    }
    Result[6] = 0;
    for (int j = 0; j < methodSteps; j++)
      Result[6] = Result[6] + BRes[j] * dpZeroStar[j];

    for (int i = 0; i < 7; i++)
      Error[i] = 0;

    for (int i = 0; i < methodSteps; i++) {
      for (int j = 0; j < 6; j++)
        Error[j] = Error[j] + B[i] * DSigma[i][j];
      Error[6] = Error[6] + B[i] * dpZeroStar[i];
    }
    if (!errorEstimate)
      for (int i = 0; i < 7; i++)
        Error[i] = Error[i] - Result[i]; // error estimate calculated in case we
                                         // have lower order solution instead of
                                         // error estimate

    // check the error norm

    switch (int(d_tolMethod)) {
      case 0: {
        RError = checkNorm(Result, Result[6], point, Error); // returns RError
      } break;
      case 1: {
        // SLOAN NORM
        RError = checkNormSloan(Result, Result[6], point, Error); // returns
                                                                  // RError
        // cout<<"Runge - Kutta constant size procedure. RError="<<RError<<"\n";
        // cout<<"Error[0]="<<Error[0]<<"  Result[0]="<<Result[0]<<"\n";
      } break;
      default: {
        cout << "ERROR !!!! Improper d_tolMethod in increment.dta" << "\n";
        getchar();
      }
    }

    // cout<<"Procedure R-K constant. RError="<<RError<<"\n";

    for (int i = 0; i < 7; i++)
      if (!finite(Result[i])) {
        Result[i] = 0;
        if (RError < methodPower)
          RError = methodPower;
      }

    if ((Point->getmeanStress() + (Result[0] + Result[1] + Result[2])) / 3 <
        0) {
      // cout<<"Mean Stress less then 0!!!
      // Result:"<<((Result[0]+Result[1]+Result[2])/3)<<"  Mean
      // stress:"<<Point->getmeanStress()<<"\n";

      if (RError < methodPower)
        RError = methodPower;
    }
    if ((Point->p0Star() + Result[6]) < 0) {
      // cout<<"P Zero Star less then 0!!!"<<"\n";
      if (RError < methodPower)
        RError = methodPower;
    }
    // here we need to update all the point data.
    Point->Copy(&OldPoint);
    Point->Update(plasticStrainTemp, SubstepStrain, Result, Result[6]);
    Point->Copy(&TrialPoint);
    // this value is not used in any calculations, it is just to show at the end
    // of the step the p0* increase

    // cout<<"Procedure R-K constant. RError="<<RError<<"\n";

    if (d_driftCorrection == 3)
      correctDrift(Point);
    if (d_driftCorrection == 2)
      correctDriftBeg(
        Point,
        &OldPoint); // value of OldPoint copied before updating the point

    // re - evaluate the error in the point:

    for (int i = 0; i < 6; i++) {
      Error[i] = Error[i] + Point->stress[i] - TrialPoint.stress[i];
    }
    Error[6] = Error[6] + Point->p0Star() - TrialPoint.p0Star();
    // error vector updated, norm should be re-evaluated:

    switch (int(d_tolMethod)) {
      case 0: {
        Temp = checkNorm(Result, Result[6], point, Error); // returns RError
        if (Temp > RError)
          RError = Temp;
      } break;
      case 1: {
        // SLOAN NORM
        Temp =
          checkNormSloan(Result, Result[6], point, Error); // returns RError
        if (Temp > RError)
          RError = Temp;
      } break;
      default: {
        cout << "ERROR !!!! Improper d_tolMethod in increment.dta" << "\n";
        getchar();
      }
    }
    if (!finite(RError))
      if (RError < methodPower)
        RError = methodPower;

    // cout<<"Procedure R-K constant. RError="<<RError<<"\n";
    TotRError = TotRError + RError;
    // cout<<"Procedure R-K constant. Total RError="<<TotRError<<"\n";

    // this is done anyway
    for (int i = 0; i < 6; i++)
      absStress[i] = absStress[i] + fabs(Point->stress[i] - OldPoint.stress[i]);
    absStress[6] = absStress[6] + fabs(Point->p0Star() - OldPoint.p0Star());
    StepAccuracycheck = TotalSize;
    MicroStep = MicroStep + StepLength;
    TotalSize = TotalSize + MicroStep; // total part of step done updated
    MicroStep = MicroStep - (TotalSize - StepAccuracycheck);
    // if (TotalSize>=1) Finished=true;
    Temp = double(stepNo) / Frequency;
    /*if (modf(Temp,&Temp)==0)
            {
                    cout<<"Step number:"<<stepNo<<"\n";
                    cout<<"Total size done is: "<<TotalSize<<" of whole step.
       Current StepLength="<<StepLength<<"\n";
             }*/
    /*
            //Debug
            strainInc.push_back (StepLength);
            for (int i=0; i<6; i++) stressInc.push_back (Point->stress[i]);
            stressInc.push_back (Point->p0Star());
            //End Debug */
  }

  *RelError = TotRError / numberIter;

  /*
  //debug
  FILE * ResultsFile;
  ResultsFile = fopen( "ResultsRK.dta", "w" );

  for (Iter=stressInc.begin();Iter!=stressInc.end();)
  {
          for (int i=0; i<7; i++)
                          {
                          fprintf( ResultsFile, "%.15g , ",*Iter);
                          Iter++;
                          }
          fprintf( ResultsFile, "\n");
  }
  fclose (ResultsFile);

  ResultsFile = fopen( "StrainRK.dta", "w" );
  for (Iter=strainInc.begin();Iter!=strainInc.end();)
  {
          fprintf( ResultsFile, "%.15g \n",*Iter);
          Iter++;
  }
  fclose (ResultsFile);
  //end debug */

  return (0);
}

double
ShengMohrCoulomb::doRungeuttaExtrapol(double A[][8], double* B, double* BRes,
                                     double* C, StateMohrCoulomb* point,
                                     double* epStrain, double* absStress,
                                     int* numberIter, double methodOrder,
                                     int methodSteps, bool errorEstimate)
{
  /*
  this procedure tend to use Runge Kutta scheme. In case the step is not
  accepted, and the substep size is just a bit smaller
  than current one , instead of cancelling the step, it goes into extrapolation
  procedure, just to save the substep done.
  */

  BBMMatrix dSigma(6, 1);
  StateMohrCoulomb TrialPoint, OldPoint;
  double WORTH_EXTRAPOL = 3;

  // double DSigma[8][6];
  // double dpZeroStar[8];
  double RError, NewStepSize, TotalSize, StepLength, Temp; // methodPower;
  double Frequency = 15000 / methodOrder; // how often display info about steps
  // bool ReUseStep=false;

  StepLength = 1;
  TotalSize = 0;
  NewStepSize = 1;
  // methodPower=pow(2.0,methodOrder)*d_integrationTol;

  // for (int i=0; i<methodOrder; i++) {
  // for (int j=0; j<6; j++) DSigma[i][j]=0;
  // dpZeroStar[i]=0;	}

  bool Finished = false, StepAccepted = false;
  double SubstepStrain[7]; // CurrentStrain[7];
  double MicroStep = 0;
  double StepAccuracycheck;

  for (int i = 0; i < 7; i++) {
    SubstepStrain[i] = epStrain[i];
    absStress[i] = 0;
  }
  int stepNo = 0;
  vector<double> stressInc;
  vector<double> strainInc;
  vector<double> PZeroLambda;
  vector<double>::iterator Iter;

  clock_t StartTime, EndTime;
  StartTime = clock();

  // assuming that maximum tolerable step is 0.5% for order 2 Initial step size
  // is being calculated:
  NewStepSize = 0;
  for (int i = 0; i < 6; i++)
    NewStepSize = NewStepSize + fabs(SubstepStrain[i]);
  NewStepSize = 0.01 / (NewStepSize * methodOrder);
  // cout<<"NewStepSize="<<NewStepSize<<"\n";
  if (NewStepSize > 1)
    NewStepSize = 1;
  // getchar(); */

  do {
    StepAccepted = false;

    /*if (stepNo>1)
{
	if (NewStepSize>10) NewStepSize=10;
	if (NewStepSize<0.1) NewStepSize=0.1;
}*/ // limiting step increase/decrease does not improve results/ enables faster
    // convergence...

    StepLength = StepLength * NewStepSize; // size of a step
    if ((StepLength + TotalSize) > 1)
      StepLength =
        1 - TotalSize; // check whether the step not exceed the whole increment

    for (int i = 0; i < 7; i++) {
      SubstepStrain[i] =
        StepLength * epStrain[i]; // strain increment in current step
                                  // CurrentStrain[i]=0;
    }
    RError = 0;

    Point->Copy(&TrialPoint);
    Point->Copy(&OldPoint);
    doRungeuttaEqualStep(A, B, BRes, C, &TrialPoint, SubstepStrain, absStress,
                        &RError, 2, methodOrder, methodSteps, errorEstimate);

    stepNo = stepNo + 2;
    // cout<<"Step Length="<<StepLength<<"\n";
    // cout<<"Current strain [0]="<<SubstepStrain[0]<<"\n";
    // cout<<"RError="<<RError<<"\n";
    // getchar();

    if (RError < d_integrationTol) {
      StepAccepted = true;
    } else {
      StepAccepted = false;
      if (StepLength < 1e-20)
        StepAccepted = true;
    }

    if (d_tolMethod == 0)
      NewStepSize =
        d_betaFactor * pow(d_integrationTol / RError, (1 / (methodOrder - 1.0)));
    else
      NewStepSize =
        d_betaFactor * pow(d_integrationTol / RError, (1 / methodOrder));

    if (!StepAccepted) {
      // check the extrapolation and go into it, if it looks sensible. otherwise
      // reject the step...
      if (NewStepSize < WORTH_EXTRAPOL) {
        // double Result[7];
        int TempNumber = 0;
        // for (int i=0; i<6; i++)
        // Result[i]=TrialPoint.stress[i]-OldPoint.stress[i];
        // Result[6]=TrialPoint.p0Star()-OldPoint.p0Star();
        TrialPoint.Copy(&OldPoint);
        Point->Copy(&TrialPoint);
        doRKExtrapolation(
          A, B, BRes, C, &TrialPoint, SubstepStrain, absStress, &OldPoint,
          &RError, &TempNumber, methodOrder, methodSteps,
          errorEstimate); // Extrapolate and finally accept the step
        stepNo = stepNo + TempNumber;
        StepAccepted = true;
        Point->Copy(&OldPoint);
      } else
        ;
      /*	//here we should take care about correct re - usage of the first evaluation of derivative
	ReUseStep=true;
	for (int i=0; i<6; i++) ReUseRes[i]=DSigma[0][i]*NewStepSize;
	ReUseRes[6]=dpZeroStar[0]*NewStepSize; */ // No step re-using so far
      // reject the step.
    }

    if (StepAccepted) {
      // here we need to update all the point data. and nothing else...
      for (int i = 0; i < 6; i++)
        absStress[i] =
          absStress[i] + fabs(OldPoint.stress[i] - TrialPoint.stress[i]);
      absStress[6] =
        absStress[6] + fabs(Point->p0Star() - TrialPoint.p0Star());
      TrialPoint.Copy(Point);
      // Point->Update(0,SubstepStrain,DSigmaTemp,OldPoint.p0Star()-TrialPoint.p0Star());
      // drift is already corrected
      // ReUseStep=false;
      StepAccuracycheck = TotalSize;
      MicroStep = MicroStep + StepLength;
      TotalSize = TotalSize + MicroStep; // total part of step done updated
      MicroStep = MicroStep - (TotalSize - StepAccuracycheck);
      if (TotalSize >= 1)
        Finished = true;
      Temp = double(stepNo) / Frequency;
      if (modf(Temp, &Temp) == 0) {
        cout << "Step number:" << stepNo << "\n";
        cout << "Total size done is: " << TotalSize
             << " of whole step. Current StepLength=" << StepLength << "\n";
      }
      /*
              //Debug
              strainInc.push_back (StepLength);
              for (int i=0; i<6; i++) stressInc.push_back (Point->stress[i]);
              stressInc.push_back (Point->p0Star());
              //End Debug */
    }

  } while (!Finished);
  EndTime = clock();
  *numberIter = stepNo;
  /*
  //debug
  FILE * ResultsFile;
  ResultsFile = fopen( "ResultsRK.dta", "w" );

  for (Iter=stressInc.begin();Iter!=stressInc.end();)
  {
          for (int i=0; i<7; i++)
                          {
                          fprintf( ResultsFile, "%.15g , ",*Iter);
                          Iter++;
                          }
          fprintf( ResultsFile, "\n");
  }
  fclose (ResultsFile);

  ResultsFile = fopen( "StrainRK.dta", "w" );
  for (Iter=strainInc.begin();Iter!=strainInc.end();)
  {
          fprintf( ResultsFile, "%.15g \n",*Iter);
          Iter++;
  }
  fclose (ResultsFile);
  //end debug */
  return (EndTime - StartTime);
}

double
ShengMohrCoulomb::plasticRKErr8544(StateMohrCoulomb* point, double* epStrain,
                                   double* absStress, int* numberIter)
{
  /* the procedure uses the embedded Runge - Kutta integration scheme with
  Adaptive Stepsize Control
  the constants are as proposed by Bogacki and Shampine (1996), An efficient R-K
  (4,5) pair, Computers Math Applic, Vol 32 No 6 pp 15-28
  with FSAL feauture the method allows for getting the error estimate and
  calculating value in one go
  It is arguably better than any other 5(4) RK pair; It has double 4th order
  error estimate
  */

  double A[8][8]; // matrix for aij components for RK method, as in the Fortran
                  // source www.netlib.org/ode/rksuite

  A[0][0] = 0;
  A[1][0] = 1.0 / 6.0;
  A[2][0] = 2.0 / 27.0;
  A[2][1] = 4.0 / 27.0;
  A[3][0] = 183.0 / 1372.0;
  A[3][1] = -162.0 / 343.0;
  A[3][2] = 1053.0 / 1372.0;
  A[4][0] = 68.0 / 297.0;
  A[4][1] = -4.0 / 11.0;
  A[4][2] = 42.0 / 143.0;
  A[4][3] = 1960.0 / 3861.0;
  A[5][0] = 597.0 / 22528.0;
  A[5][1] = 81.0 / 352.0;
  A[5][2] = 63099.0 / 585728.0;
  A[5][3] = 58653.0 / 366080.0;
  A[5][4] = 4617.0 / 20480.0;
  A[6][0] = 174197.0 / 959244.0;
  A[6][1] = -30942.0 / 79937.0;
  A[6][2] = 8152137.0 / 19744439.0;
  A[6][3] = 666106.0 / 1039181.0;
  A[6][4] = -29421.0 / 29068.0;
  A[6][5] = 482048.0 / 414219.0;
  A[7][0] = 587.0 / 8064.0;
  A[7][1] = 0.0;
  A[7][2] = 4440339.0 / 15491840.0;
  A[7][3] = 24353.0 / 124800.0;
  A[7][4] = 387.0 / 44800.0;
  A[7][5] = 2152.0 / 5985.0;
  A[7][6] = 7267.0 / 94080.0;

  double B[8]; //  The coefficients B[*] refer to the formula of order 4.

  B[0] = 2479.0 / 34992.0;
  B[1] = 0.0;
  B[2] = 123.0 / 416.0;
  B[3] = 612941.0 / 3411720.0;
  B[4] = 43.0 / 1440.0;
  B[5] = 2272.0 / 6561.0;
  B[6] = 79937.0 / 1113912.0;
  B[7] = 3293.0 / 556956.0;

  /*  The coefficients E(*) refer to an estimate of the local error based on
  C  the first formula of order 4.  It is the difference of the fifth order
  C  result, here located in A(8,*), and the fourth order result.  By
  C  construction both ErrorCoef[1] and ErrorCoef[6] are zero. */

  double
    ErrorCoef[7]; // first error estimate, does not require knowing the result

  ErrorCoef[0] = -3.0 / 1280.0;
  ErrorCoef[1] = 0.0;
  ErrorCoef[2] = 6561.0 / 632320.0;
  ErrorCoef[3] = -343.0 / 20800.0;
  ErrorCoef[4] = 243.0 / 12800.0;
  ErrorCoef[5] = -1.0 / 95.0;
  ErrorCoef[6] = 0.0;

  double C[8]; // ci matrix, parameters for x

  C[0] = 0.0;
  C[1] = 1.0 / 6.0;
  C[2] = 2.0 / 9.0;
  C[3] = 3.0 / 7.0;
  C[4] = 2.0 / 3.0;
  C[5] = 3.0 / 4.0;
  C[6] = 1.0;
  C[7] = 1.0;

  double BRes[8];
  for (int i = 0; i < 7; i++)
    BRes[i] = A[7][i];
  BRes[7] = 0; // as this scheme is the FSAL scheme

  // All the RK matrices are put from the Fortran code, the indices are reduced
  // by 1 and C notation is used
  // BBMMatrix dSigma (6,1);
  // StateMohrCoulomb MidPoint[8], OldPoint, TrialPoint;

  int methodSteps = 8;
  double methodOrder = 5.0;
  bool errorEstimate = false;

  // double time;
  // time=doRungeutta
  // (A,B,BRes,C,Point,epStrain,absStress,numberIter,methodOrder,methodSteps,false);
  // return time;

  BBMMatrix dSigma(6, 1);
  StateMohrCoulomb MidPoint[8], OldPoint, TrialPoint;

  double DSigma[8][6], DSigmaTemp[7], Result[7];
  double dpZeroStar[8], dpZeroStarTemp, plasticStrainTemp[6],
    plasticStrain[8][6];
  double Error[7], ErrorOther[7], ReUseRes[7], RError, RErrorOther, NewStepSize,
    TotalSize, StepLength, Temp, methodPower;
  double Frequency = 15000 / methodOrder; // how often display info about steps
  double dStressDrift[6], dLambdaDrift;

  bool ReUseStep = false;

  StepLength = 1;
  TotalSize = 0;
  NewStepSize = 1;
  methodPower = pow(2.0, methodOrder) * d_integrationTol;

  for (int i = 0; i < methodOrder; i++) {
    for (int j = 0; j < 6; j++)
      DSigma[i][j] = 0;
    dpZeroStar[i] = 0;
  }

  bool Finished = false, StepAccepted = false;
  double SubstepStrain[7], CurrentStrain[7];
  double MicroStep = 0;
  double StepAccuracycheck;

  for (int i = 0; i < 7; i++) {
    SubstepStrain[i] = epStrain[i];
    absStress[i] = 0;
  }
  int stepNo = 0;
  vector<double> stressInc;
  vector<double> strainInc;
  vector<double> PZeroLambda;
  vector<double>::iterator Iter;

  clock_t StartTime, EndTime;
  StartTime = clock();

  NewStepSize = 0;
  for (int i = 0; i < 6; i++)
    NewStepSize = NewStepSize + fabs(SubstepStrain[i]);
  NewStepSize = 0.01 / (NewStepSize * methodOrder);
  // cout<<"NewStepSize="<<NewStepSize<<"\n";
  if (NewStepSize > 1)
    NewStepSize = 1;

  do {
    StepAccepted = false;
    stepNo++;

    /*if (stepNo>1)
{
	if (NewStepSize>10) NewStepSize=10;
	if (NewStepSize<0.1) NewStepSize=0.1;
}*/ // limiting step increase/decrease does not improve results/ enables faster
    // convergence...

    StepLength = StepLength * NewStepSize; // size of a step
    if ((StepLength + TotalSize) > 1)
      StepLength =
        1 - TotalSize; // check whether the step not exceed the whole increment

    for (int i = 0; i < 7; i++) {
      SubstepStrain[i] =
        StepLength * epStrain[i]; // strain increment in current step
      CurrentStrain[i] = 0;
    }
    RError = 0;

    // cout<<"Step Length="<<StepLength<<"\n";
    // cout<<"Current strain [0]="<<CurrentStrain[0]<<"\n";

    for (int i = 0; i < methodSteps; i++)
      point->Copy(&MidPoint[i]); // point is unchanged in  procedure

    // Below the main R-K loop to calculate the value of intermediate stresses;
    // values stored in DSigmaTemp[][]
    // ReUseStep=false;
    if (ReUseStep) {
      for (int i = 0; i < 6; i++) {
        DSigma[0][i] = ReUseRes[i];
        plasticStrain[0][i] = plasticStrainTemp[i];
      }
      dpZeroStar[0] = ReUseRes[6];
      // add line about plastic strain...
    } else {
      calcPlastic(MidPoint[0], SubstepStrain, &dSigma, plasticStrainTemp,
                  &dpZeroStar[0], computeYieldFunctionNN(Point), dStressDrift,
                  &dLambdaDrift);
      for (int i = 0; i < 6; i++) {
        DSigma[0][i] = dSigma.getElement(i + 1, 1);
        plasticStrain[0][i] = plasticStrainTemp[i];
      }
    }

    for (int rkloop = 1; rkloop < methodSteps; rkloop++) {
      for (int i = 0; i < 6; i++) {
        DSigmaTemp[i] = 0;
        plasticStrainTemp[i] = 0;
      }
      dpZeroStarTemp = 0;
      for (int i = 0; i < 7; i++)
        CurrentStrain[i] =
          C[rkloop] *
          SubstepStrain[i]; // set the beginning point of the procedure
      for (int i = 0; i < rkloop; i++) {
        for (int j = 0; j < 6; j++) {
          DSigmaTemp[j] = DSigmaTemp[j] + A[rkloop][i] * DSigma[i][j];
          plasticStrainTemp[j] =
            plasticStrainTemp[j] + A[rkloop][i] * plasticStrain[i][j];
        }
        dpZeroStarTemp = dpZeroStarTemp + A[rkloop][i] * dpZeroStar[i];
      }
      MidPoint[rkloop].Update(plasticStrainTemp, CurrentStrain, DSigmaTemp,
                              dpZeroStarTemp);
      // double dummy;
      // StateMohrCoulomb TempPoint;
      // BBMMatrix TEMPMATRIX (6,7), TEMPEPSILON(7,1);
      // SubstepStrain[6]=0;
      // for (int i=1; i<8; i++) TEMPEPSILON.PutElement(i,1,SubstepStrain[i-1]);
      // MidPoint[rkloop].Copy(&TempPoint);
      calcPlastic(MidPoint[rkloop], SubstepStrain, &dSigma, plasticStrainTemp,
                  &dpZeroStar[rkloop], computeYieldFunctionNN(Point),
                  dStressDrift, &dLambdaDrift);
      // calculateElastoPlasticTangentMatrix(&TempPoint,&TEMPMATRIX);
      // TEMPMATRIX.Multiply(&TEMPEPSILON,&TEMPMATRIX);
      /*	for (int i=1; i<7; i++)
              {
                      dummy=fabs(dSigma.getElement(i,1)-TEMPMATRIX.getElement(i,1));
                      if (fabs(dSigma.getElement(i,1))>TINY)
         dummy=dummy/fabs(dSigma.getElement(i,1));
                      if (dummy>0.01)
                      {
                              cout<<"Problems with the alternative
         matrix."<<i<<" Change of stress is:"<<dSigma.getElement(i,1)<<"\n";
                              cout<<"calculated change of stress with the DEP
         matrix is:"<<TEMPMATRIX.getElement(i,1)<<"\n";
                              getchar();
                      }
                      else cout<<"*";
              }


      */
      for (int i = 0; i < 6; i++) {
        DSigma[rkloop][i] = dSigma.getElement(i + 1, 1);
        plasticStrain[rkloop][i] = plasticStrainTemp[i];
      }
    }

    // needed: result, error

    for (int i = 0; i < 6; i++) {
      Result[i] = 0;
      plasticStrainTemp[i] = 0;
      for (int j = 0; j < methodSteps; j++) {
        Result[i] = Result[i] + BRes[j] * DSigma[j][i];
        plasticStrainTemp[i] =
          plasticStrainTemp[i] + BRes[j] * plasticStrain[j][i];
      }
    }
    Result[6] = 0;
    for (int j = 0; j < methodSteps; j++)
      Result[6] = Result[6] + BRes[j] * dpZeroStar[j];

    for (int i = 0; i < 7; i++)
      Error[i] = 0;

    for (int i = 0; i < methodSteps; i++) {
      for (int j = 0; j < 6; j++)
        Error[j] = Error[j] + B[i] * DSigma[i][j];
      Error[6] = Error[6] + B[i] * dpZeroStar[i];
    }
    if (!errorEstimate)
      for (int i = 0; i < 7; i++)
        Error[i] = Error[i] - Result[i]; // error estimate calculated in case we
                                         // have lower order solution instead of
                                         // error estimate

    for (int i = 0; i < 7; i++)
      ErrorOther[i] = 0;
    for (int i = 0; i < (methodSteps - 1); i++) {
      for (int j = 0; j < 6; j++)
        ErrorOther[j] = ErrorOther[j] + ErrorCoef[i] * DSigma[i][j];
      ErrorOther[6] = ErrorOther[6] + ErrorCoef[i] * dpZeroStar[i];
    }

    // check the error norm

    switch (int(d_tolMethod)) {
      case 0: {
        // RError=checkNorm (DSigmaTemp, dpZeroStarTemp, point, Error);
        // //returns RError
        RError = checkNorm(Result, Result[6], point, Error); // returns RError
        RErrorOther = checkNorm(Result, Result[6], point, ErrorOther);
        if (RError < RErrorOther)
          RError = RErrorOther;
      } break;
      case 1: {
        // SLOAN NORM
        // RError=checkNormSloan (DSigmaTemp, dpZeroStarTemp, point, Error);
        // //returns RError
        RError = checkNormSloan(Result, Result[6], point, Error); // returns
                                                                  // RError
        RErrorOther = checkNormSloan(Result, Result[6], point,
                                     ErrorOther); // returns RError
        if (RError < RErrorOther)
          RError = RErrorOther;
      } break;
      default: {
        cout << "ERROR !!!! Improper d_tolMethod in increment.dta" << "\n";
        getchar();
      }
    }

    for (int i = 0; i < 7; i++)
      if (!finite(Result[i])) {
        Result[i] = 0;
        if (RError < methodPower)
          RError = methodPower;
      }

    if ((Point->getmeanStress() + (Result[0] + Result[1] + Result[2])) / 3 < 0)
      if (RError < methodPower)
        RError = methodPower;
    if ((Point->p0Star() + Result[6]) < 0)
      if (RError < methodPower)
        RError = methodPower;

    if (RError < d_integrationTol) {
      StepAccepted = true;
    } else {
      StepAccepted = false;
      if (StepLength < 1e-20)
        StepAccepted = true;
    }

    if (d_tolMethod == 0)
      NewStepSize =
        d_betaFactor * pow(d_integrationTol / RError, (1 / (methodOrder - 1.0)));
    else
      NewStepSize =
        d_betaFactor * pow(d_integrationTol / RError, (1 / methodOrder));

    if (!StepAccepted) {
      // here we should take care about correct re - usage of the first
      // evaluation of derivative
      ReUseStep = true;
      for (int i = 0; i < 6; i++) {

        ReUseRes[i] = DSigma[0][i] * NewStepSize;
        plasticStrainTemp[i] = plasticStrain[0][i] * NewStepSize;
      }
      ReUseRes[6] = dpZeroStar[0] * NewStepSize;
    } else {
      // here we need to update all the point data.
      Point->Copy(&OldPoint);
      Point->Update(plasticStrainTemp, SubstepStrain, Result, Result[6]);
      Point->Copy(&TrialPoint);
      // this value is not used in any calculations, it is just to show at the
      // end of the step the p0* increase

      if (d_driftCorrection == 3)
        correctDrift(Point);
      if (d_driftCorrection == 2)
        correctDriftBeg(
          Point,
          &OldPoint); // value of OldPoint copied before updating the point

      // re - evaluate the error in the point:

      for (int i = 0; i < 6; i++) {
        Error[i] = Error[i] + Point->stress[i] - TrialPoint.stress[i];
        ErrorOther[i] = ErrorOther[i] + Point->stress[i] - TrialPoint.stress[i];
      }
      Error[6] = Error[6] + Point->p0Star() - TrialPoint.p0Star();
      ErrorOther[6] = ErrorOther[6] + Point->p0Star() - TrialPoint.p0Star();
      // error vector updated, norm should be re-evaluated:

      switch (int(d_tolMethod)) {
        case 0: {
          // Temp=checkNorm (DSigmaTemp, dpZeroStarTemp, point, Error);
          // //returns RError
          Temp = checkNorm(Result, Result[6], point, Error); // returns RError
          RErrorOther =
            checkNorm(Result, Result[6], point, ErrorOther); // returns RError
          if (Temp > RError)
            RError = Temp;
          if (RErrorOther > RError)
            RError = RErrorOther;
        } break;
        case 1: {
          // SLOAN NORM
          // Temp=checkNormSloan (DSigmaTemp, dpZeroStarTemp, point, Error);
          // //returns RError
          Temp = checkNormSloan(Result, Result[6], point, Error); // returns
                                                                  // RError
          RErrorOther = checkNormSloan(Result, Result[6], point,
                                       ErrorOther); // returns RError
          if (Temp > RError)
            RError = Temp;
          if (RErrorOther > RError)
            RError = RErrorOther;
        } break;
        default: {
          cout << "ERROR !!!! Improper d_tolMethod in increment.dta" << "\n";
          getchar();
        }
      }
      if (!finite(RError))
        if (RError < methodPower)
          RError = methodPower;

      if (d_tolMethod == 0)
        NewStepSize =
          d_betaFactor * pow(d_integrationTol / RError, (1 / (methodOrder - 1.0)));
      else
        NewStepSize =
          d_betaFactor * pow(d_integrationTol / RError, (1 / methodOrder));

      if (RError < d_integrationTol)
        StepAccepted = true;
      else {
        StepAccepted = false;
        if (StepLength < 1e-20)
          StepAccepted = true;
      }

      if (!StepAccepted) {
        ReUseStep = true;
        for (int i = 0; i < 6; i++) {
          ReUseRes[i] = DSigma[0][i] * NewStepSize;
          plasticStrainTemp[i] = plasticStrain[0][i] * NewStepSize;
        }
        ReUseRes[6] = dpZeroStar[0] * NewStepSize;
        OldPoint.Copy(Point);
      } else {
        // this may be done only after successful re-evaluation of the substep
        for (int i = 0; i < 6; i++)
          absStress[i] =
            absStress[i] + fabs(Point->stress[i] - OldPoint.stress[i]);
        absStress[6] =
          absStress[6] + fabs(Point->p0Star() - OldPoint.p0Star());
        ReUseStep = false;
        StepAccuracycheck = TotalSize;
        MicroStep = MicroStep + StepLength;
        TotalSize = TotalSize + MicroStep; // total part of step done updated
        MicroStep = MicroStep - (TotalSize - StepAccuracycheck);
        if (TotalSize >= 1)
          Finished = true;
        Temp = double(stepNo) / Frequency;
        if (modf(Temp, &Temp) == 0) {
          cout << "Step number:" << stepNo << "\n";
          cout << "Total size done is: " << TotalSize
               << " of whole step. Current StepLength=" << StepLength << "\n";
        }
        /*
                //Debug
                strainInc.push_back (StepLength);
                for (int i=0; i<6; i++) stressInc.push_back (Point->stress[i]);
                stressInc.push_back (Point->p0Star());
                //End Debug */
      }
    }

  } while (!Finished);
  EndTime = clock();
  *numberIter = stepNo;
  /*
  //debug
  FILE * ResultsFile;
  ResultsFile = fopen( "ResultsRK.dta", "w" );

  for (Iter=stressInc.begin();Iter!=stressInc.end();)
  {
          for (int i=0; i<7; i++)
                          {
                          fprintf( ResultsFile, "%.15g , ",*Iter);
                          Iter++;
                          }
          fprintf( ResultsFile, "\n");
  }
  fclose (ResultsFile);

  ResultsFile = fopen( "StrainRK.dta", "w" );
  for (Iter=strainInc.begin();Iter!=strainInc.end();)
  {
          fprintf( ResultsFile, "%.15g \n",*Iter);
          Iter++;
  }
  fclose (ResultsFile);
  //end debug */

  return (EndTime - StartTime);
}

double
ShengMohrCoulomb::plasticRKNoExTry(StateMohrCoulomb* point, double* epStrain,
                                   double* absStress, int* numberIter)
/* It is a non - extrapolation procedure !!! */

/*
This procedure calculate stress increment using Runge - Kutta pair as given by
England. The procedure consists of 6 stages and should be less
efficient than the RKErr8544 that uses Bogacki - Shimpine pair.
*/

{
  int methodSteps = 6;

  double methodOrder = 5;
  double time;
  double A[8][8]; // matrix must be this size to be used in the doRungeutta
                  // method
  bool errorEstimate = false; // we give a 2nd order solution
  /*A - matrix with coefficients, B - error estimate, BRes - result
   * coefficients, C - x coefficients. */
  double C[6] = { 0.0, 1.0 / 3.0, 0.4, 1, 2.0 / 3.0, 0.8 };
  double BRes[6] = { 23 / 192.0,   0, 125 / 192.0, 0.0, -81.0 / 192.0,
                     125.0 / 192.0 };
  double B[6] = { -0.5, 1.5, 0, 0, 0, 0 };

  A[0][0] = 0;

  A[1][0] = 1.0 / 3.0;

  A[2][0] = 4.0 / 25.0;
  A[2][1] = 6.0 / 25.0;

  A[3][0] = 0.25;
  A[3][1] = -3.0;
  A[3][2] = 3.75;

  A[4][0] = 6.0 / 81.0;
  A[4][1] = 90 / 81.0;
  A[4][2] = -50 / 81.0;
  A[4][3] = 8.0 / 81.0;

  A[5][0] = 6.0 / 75.0;
  A[5][1] = 36.0 / 75.0;
  A[5][2] = 10.0 / 75.0;
  A[5][3] = 8.0 / 75.0;
  A[5][4] = 0.0;

  time = doRungeutta(A, B, BRes, C, point, epStrain, absStress, numberIter,
                    methodOrder, methodSteps, errorEstimate);
  return time;
}

double
ShengMohrCoulomb::plasticRKDP754(StateMohrCoulomb* point, double* epStrain,
                                 double* absStress, int* numberIter)
/*
This procedure calculate stress increment using Runge - Kutta pair as given by
DORMAND - PRINCE. The
used pair is known also as DOPRI5 or RK5(4)7FM. The procedure consists of 7
stages and should be less
efficient than the RKErr8544 that uses Bogacki - Shimpine pair.
*/

{
  int methodSteps = 7;
  double methodOrder = 5;
  double time;
  double A[8][8];
  bool errorEstimate =
    false; // we give a 4th order solution, not an error estimate

  /*A - matrix with coefficients, B - error estimate, BRes - result
   * coefficients, C - x coefficients. */
  double C[7] = { 0.0, 0.2, 0.3, 0.8, 8.0 / 9.0, 1, 1 };
  double BRes[7] = {
    35.0 / 384.0, 0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0,
    11.0 / 84.0,  0
  };
  double B[7] = { 5179.0 / 57600.0,
                  0,
                  7571.0 / 16695.0,
                  393.0 / 640.0,
                  -92097.0 / 339200.0,
                  187.0 / 2100.0,
                  0.025 };

  A[0][0] = 0;

  A[1][0] = 0.2;

  A[2][0] = 0.075;
  A[2][1] = 0.225;

  A[3][0] = 44.0 / 45.0;
  A[3][1] = -56.0 / 15.0;
  A[3][2] = 32.0 / 9.0;

  A[4][0] = 19372.0 / 6561.0;
  A[4][1] = -25360.0 / 2187.0;
  A[4][2] = 64448.0 / 6561.0;
  A[4][3] = -212.0 / 729.0;

  A[5][0] = 9017.0 / 3168.0;
  A[5][1] = -355.0 / 33.0;
  A[5][2] = 46732.0 / 5247.0;
  A[5][3] = 49.0 / 176.0;
  A[5][4] = -5103.0 / 18656.0;

  A[6][0] = 35.0 / 384.0;
  A[6][1] = 0;
  A[6][2] = 500.0 / 1113.0;
  A[6][3] = 125.0 / 192.0;
  A[6][4] = -2187.0 / 6784.0;
  A[6][5] = 11.0 / 84.0;

  // time=doRungeuttaExtrapol (A, B,BRes, C, point, epStrain, absStress,
  // numberIter, methodOrder, methodSteps, errorEstimate);
  time = doRungeutta(A, B, BRes, C, point, epStrain, absStress, numberIter,
                    methodOrder, methodSteps, errorEstimate);
  return time;
}

double
ShengMohrCoulomb::plasticRKCK654(StateMohrCoulomb* point, double* epStrain,
                                 double* absStress, int* numberIter)
/*
This procedure calculate stress increment using Runge - Kutta pair as given by
Cash - Karp. The coeficients
are given in Numerical Recipes (Cambrige Univ Press) or Cash Karp (1990) ACM
Transactions on Mathematical
software, vol 16, pp 201-222. The procedure consists of 6 stages and should be
less
efficient than the RKErr8544 that uses Bogacki - Shimpine pair.
*/

{
  int methodSteps = 6;

  double methodOrder = 5;
  double time;
  double A[8][8]; // matrix must be this size to be used in the doRungeutta
                  // method
  bool errorEstimate =
    false; // we give a 4th order solution, not an error estimate

  /*A - matrix with coefficients, B - error estimate, BRes - result
   * coefficients, C - x coefficients. */
  double C[6] = { 0.0, 0.2, 0.3, 0.6, 1, 0.875 };
  double BRes[6] = { 37.0 / 378.0,  0, 250.0 / 621.0,
                     125.0 / 594.0, 0, 512.0 / 1771.0 };
  double B[6] = { 2825.0 / 27648.0, 0,   18575.0 / 48384.0, 13525.0 / 55296.0,
                  277.0 / 14336.0,  0.25 };

  A[0][0] = 0;

  A[1][0] = 0.2;

  A[2][0] = 0.075;
  A[2][1] = 0.225;

  A[3][0] = 0.3;
  A[3][1] = -0.9;
  A[3][2] = 1.2;

  A[4][0] = -11.0 / 54.0;
  A[4][1] = 2.5;
  A[4][2] = -70.0 / 27.0;
  A[4][3] = 35.0 / 27.0;

  A[5][0] = 1631.0 / 55296.0;
  A[5][1] = 175.0 / 512.0;
  A[5][2] = 575.0 / 13824.0;
  A[5][3] = 44275.0 / 110592.0;
  A[5][4] = 253.0 / 4096.0;

  time = doRungeutta(A, B, BRes, C, point, epStrain, absStress, numberIter,
                    methodOrder, methodSteps, errorEstimate);
  // time=doRungeuttaExtrapol (A, B,BRes, C, point, epStrain, absStress,
  // numberIter, methodOrder, methodSteps, errorEstimate);
  return time;
}

double
ShengMohrCoulomb::plasticRKEng654(StateMohrCoulomb* point, double* epStrain,
                                  double* absStress, int* numberIter)
/*
This procedure calculate stress increment using Runge - Kutta pair as given by
Sloan (1987). The coeficients
are given in Sloan (1987), Ordinary and Partial Differential Equations routines
in... by Lee & Schiesser, Chapman & Hall 2003
or originally by England (1969) Error Estimates for Runge - Kutta type solutions
to systems of ordinary differential equations,
Computer Journal 12 - 166-170. The procedure consists of 6 stages and should be
the least
efficient from all the R-K pairs presented
*/

{
  int methodSteps = 6;

  double methodOrder = 5;
  double time;
  double A[8][8]; // matrix must be this size to be used in the doRungeutta
                  // method
  bool errorEstimate = true; // we give a 4th order error estimate

  /*A - matrix with coefficients, B - error estimate, BRes - result
   * coefficients, C - x coefficients. */
  double C[6] = { 0.0, 0.5, 0.5, 1, 2.0 / 3.0, 0.2 };
  double BRes[6] = { 14.0 / 336.0, 0, 0, 35.0 / 336.0, 162.0 / 336.0,
                     125.0 / 336.0 };
  double B[6] = { -42.0 / 336.0,  0,
                  -224.0 / 336.0, -21.0 / 336.0,
                  162.0 / 336.0,  125.0 / 336.0 };

  A[0][0] = 0;

  A[1][0] = 0.5;

  A[2][0] = 0.25;
  A[2][1] = 0.25;

  A[3][0] = 0;
  A[3][1] = -1.0;
  A[3][2] = 2.0;

  A[4][0] = 7.0 / 27.0;
  A[4][1] = 10.0 / 27.0;
  A[4][2] = 0;
  A[4][3] = 1.0 / 27.0;

  A[5][0] = 28.0 / 625.0;
  A[5][1] = -125.0 / 625.0;
  A[5][2] = 546.0 / 625.0;
  A[5][3] = 54.0 / 625.0;
  A[5][4] = -378.0 / 625.0;

  time = doRungeutta(A, B, BRes, C, point, epStrain, absStress, numberIter,
                    methodOrder, methodSteps, errorEstimate);
  // time=doRungeuttaExtrapol (A, B,BRes, C, point, epStrain, absStress,
  // numberIter, methodOrder, methodSteps, errorEstimate);
  return time;
}

double
ShengMohrCoulomb::plasticRK543(StateMohrCoulomb* point, double* epStrain,
                               double* absStress, int* numberIter)
/*
This procedure calculate stress increment using 4-3 Runge - Kutta pair. The
coeficients
are given in  Ordinary and Partial Differential Equations routines in... by Lee
& Schiesser, Chapman & Hall 2003
The procedure consists of 5 stages and is 4th order accurate.
*/

{
  int methodSteps = 5;

  double methodOrder = 4;
  double time;
  double A[8][8]; // matrix must be this size to be used in the doRungeutta
                  // method
  bool errorEstimate = true; // we give a 3th order error estimate

  /*A - matrix with coefficients, B - error estimate, BRes - result
   * coefficients, C - x coefficients. */

  /*
  double C[5]={0.0 , 1.0 , 1.0 , 1.5 , 3};
  double LowOrdSolution[5]={0.5 , 0 , 0 , 2 , 0.5}
  double BRes[5]={0.3 , 0 , 0.9 , 1.2 , 0.6};
  double B[5]={ -0.2 , 0 , 0.9 , -0.8 , 0.1};

  A[0][0]=0;

  A[1][0]=1.0;

  A[2][0]=0.5;
  A[2][1]=0.5;

  A[3][0]=0.375;
  A[3][1]=0;
  A[3][2]=1.125;

  A[4][0]=1.5;
  A[4][1]=0;
  A[4][2]= -4.5;
  A[4][3]= 6;

  double TempStrain[7];
  for (int i=0; i<7; i++) TempStrain[i]=epStrain[i]/3; */

  double C[5] = { 0.0, 1.0 / 3.0, 1.0 / 3.0, 0.5, 1 };
  // double LowOrdSolution[5]={1.0/6.0 , 0 , 0 , 2.0/3.0 , 1.0/6.0}
  double BRes[5] = { 0.1, 0, 0.3, 0.4, 0.2 };
  double B[5] = { -1.0 / 15.0, 0, 0.3, -4.0 / 15.0, 1.0 / 30.0 };

  A[0][0] = 0;

  A[1][0] = 1.0 / 3.0;

  A[2][0] = 0.5 / 3.0;
  A[2][1] = 0.5 / 3.0;

  A[3][0] = 0.125;
  A[3][1] = 0;
  A[3][2] = 0.375;

  A[4][0] = 0.5;
  A[4][1] = 0;
  A[4][2] = -1.5;
  A[4][3] = 2.0;

  time = doRungeutta(A, B, BRes, C, point, epStrain, absStress, numberIter,
                    methodOrder, methodSteps, errorEstimate);
  // time=doRungeuttaExtrapol (A, B,BRes, C, point, epStrain, absStress,
  // numberIter, methodOrder, methodSteps, errorEstimate);
  return time;
}

double
ShengMohrCoulomb::plasticRK332(StateMohrCoulomb* point, double* epStrain,
                               double* absStress, int* numberIter)
/*
This procedure calculate stress increment using 3-2 Runge - Kutta pair. The
coeficients
are given in  Ordinary and Partial Differential Equations routines in... by Lee
& Schiesser, Chapman & Hall 2003
The procedure consists of 3 stages and is 3th order accurate.
*/

{
  int methodSteps = 3;

  double methodOrder = 3;
  double time;
  double A[8][8]; // matrix must be this size to be used in the doRungeutta
                  // method
  bool errorEstimate = false; // we give a 2nd order solution

  /*A - matrix with coefficients, B - error estimate, BRes - result
   * coefficients, C - x coefficients. */

  double C[3] = { 0.0, 2.0 / 3.0, 2.0 / 3.0 };
  double BRes[3] = { 0.25, 0.375, 0.375 };
  double B[3] = { 0.25, 0.75, 0 };

  A[0][0] = 0;

  A[1][0] = 2.0 / 3.0;

  A[2][0] = 0;
  A[2][1] = 2.0 / 3.0;

  time = doRungeutta(A, B, BRes, C, point, epStrain, absStress, numberIter,
                    methodOrder, methodSteps, errorEstimate);
  // time=doRungeuttaExtrapol (A, B,BRes, C, point, epStrain, absStress,
  // numberIter, methodOrder, methodSteps, errorEstimate);
  return time;
}

double
ShengMohrCoulomb::plasticRKBog432(StateMohrCoulomb* point, double* epStrain,
                                  double* absStress, int* numberIter)
/*
This procedure calculate stress increment using 3-2 Runge - Kutta pair. The
coeficients
are given in "A 3(2) Pair of Runge-Kutta Formulas" by P. Bogacki and L.F.
Shampine, Appl. Math. Lett., 2, pp. 321-325, 1989.
The procedure consists of 3 stages and is 3th order accurate.
*/

{
  int methodSteps = 4;

  double methodOrder = 3;
  double time;
  double A[8][8]; // matrix must be this size to be used in the doRungeutta
                  // method
  bool errorEstimate = false; // we give a 2nd order solution

  /*A - matrix with coefficients, B - error estimate, BRes - result
   * coefficients, C - x coefficients. */

  double C[4] = { 0.0, 0.5, 0.75, 1.0 };
  double BRes[4] = { 2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0 };
  double B[4] = { 7.0 / 24.0, 0.25, 1 / 3.0, 0.125 };

  A[0][0] = 0;

  A[1][0] = 0.5;

  A[2][0] = 0;
  A[2][1] = 0.75;

  A[3][0] = 2.0 / 9.0;
  A[3][1] = 1.0 / 3.0;
  A[3][2] = 4.0 / 9.0;

  time = doRungeutta(A, B, BRes, C, point, epStrain, absStress, numberIter,
                    methodOrder, methodSteps, errorEstimate);
  // time=doRungeuttaExtrapol (A, B,BRes, C, point, epStrain, absStress,
  // numberIter, methodOrder, methodSteps, errorEstimate);
  return time;
}

double
ShengMohrCoulomb::plasticRKME221(StateMohrCoulomb* point, double* epStrain,
                                 double* absStress, int* numberIter)
/*
This procedure calculate stress increment using Modified Euler method
*/

{
  int methodSteps = 2;
  double methodOrder = 2;
  double time;
  double A[8][8]; // matrix must be this size to be used in the doRungeutta
                  // method
  bool errorEstimate = false; // we give a 2nd order solution

  /*A - matrix with coefficients, B - error estimate, BRes - result
   * coefficients, C - x coefficients. */

  double C[2] = { 0.0, 1.0 };
  double BRes[2] = { 0.5, 0.5 };
  double B[2] = { 1.0, 0 };
  A[0][0] = 0;
  A[1][0] = 1.0;

  // Matrix for midpoint method
  /*
  double C[2]={0.0 , 0.5};
  double BRes[2]={0 , 1.0};
  double B[2]={ 1.0 , 0};
  A[0][0]=0;
  A[1][0]=0.5; */

  *numberIter = 0;

  // time=doRungeuttaEqualStep (A, B,BRes, C, point, epStrain, absStress,
  // &RError, *numberIter, methodOrder, methodSteps, errorEstimate);
  // cout<<"RKME221. Error="<<RError<<"\n";
  // time=doRungeuttaExtrapol (A, B,BRes, C, point, epStrain, absStress,
  // numberIter, methodOrder, methodSteps, errorEstimate);
  // cout<<"Runge Kutta Extrapol Finished. Number of
  // iterations:"<<*numberIter<<"\n";
  time = doRungeutta(A, B, BRes, C, point, epStrain, absStress, numberIter,
                    methodOrder, methodSteps, errorEstimate);
  // time=doRungeuttaExtrapol (A, B,BRes, C, point, epStrain, absStress,
  // numberIter, methodOrder, methodSteps, errorEstimate);

  return time;
}

double
ShengMohrCoulomb::plasticMidpointGallipoli(StateMohrCoulomb* point, double* epStrain,
                                           double* absStress, int* numberIter)
{
  /*this procedure is to calculate the stress increase using the midpoint method
  with given number of iterations numberIter.
  is made mainly for use in the extrapolation procedure. It has no adaptive
  substepping or error control. It just integrates
  the strain to get the stress in given number of substeps using the midpoint
  method
  */
  // clock_t StartTime, EndTime;
  // StartTime=clock();

  StateMohrCoulomb NewPoint, MidPoint;
  BBMMatrix SIGMA(6, 1);
  // vector <double> Result;
  vector<double>::iterator Iter;

  double DSigma[6], CurrentStrain[7], HalfCurrentStrain[7], plasticStrain[6];
  double h, dpZeroStar = 0;
  double dStressDrift[6], dLambdaDrift;

  h = *numberIter;
  h = 1 / h;
  for (int i = 0; i < 7; i++)
    absStress[i] = 0; // 7 components...
  for (int i = 0; i < 7; i++) {
    CurrentStrain[i] = epStrain[i] * h; // strain increment in current step
    HalfCurrentStrain[i] = 0.5 * CurrentStrain[i];
  }
  // cout<<"Step Length="<<StepLength<<"\n";
  // cout<<"Current strain [0]="<<CurrentStrain[0]<<"\n";
  /*
  for (int i=0; i<6; i++) Result.push_back(Point->stress[i]);
  Result.push_back(Point->p0Star());
  Result.push_back(Point->getmeanStress());
  Result.push_back(Point->shearStress());
  */

  Point->Copy(&NewPoint); // point is unchanged in calcPlastic procedure
  Point->Copy(&MidPoint); // point is unchanged in calcPlastic procedure
  calcPlastic(NewPoint, HalfCurrentStrain, &SIGMA, plasticStrain, &dpZeroStar,
              computeYieldFunctionNN(Point), dStressDrift,
              &dLambdaDrift); // calculate the plastic stresses...
  for (int i = 0; i < 6; i++)
    DSigma[i] = SIGMA.getElement(i + 1, 1);
  NewPoint.Update(plasticStrain, HalfCurrentStrain, DSigma, dpZeroStar);
  calcPlastic(NewPoint, HalfCurrentStrain, &SIGMA, plasticStrain, &dpZeroStar,
              computeYieldFunctionNN(Point), dStressDrift,
              &dLambdaDrift); // calculate the plastic stresses...

  for (int loop = 0; loop < 2 * (*numberIter); loop++) {

    MidPoint.Update(plasticStrain, HalfCurrentStrain, DSigma, dpZeroStar);
    MidPoint.Copy(&NewPoint);
    NewPoint.Update(plasticStrain, HalfCurrentStrain, DSigma, dpZeroStar);
    calcPlastic(NewPoint, HalfCurrentStrain, &SIGMA, plasticStrain, &dpZeroStar,
                computeYieldFunctionNN(Point), dStressDrift,
                &dLambdaDrift); // calculate the plastic stresses...
    for (int i = 0; i < 6; i++)
      DSigma[i] = SIGMA.getElement(i + 1, 1);

    /*for (int i=0; i<6; i++) Result.push_back(Point->stress[i]);
    Result.push_back(Point->p0Star());
    Result.push_back(Point->getmeanStress());
    Result.push_back(Point->shearStress());*/
  }
  MidPoint.Copy(Point);
  // EndTime=clock();
  // return (EndTime-StartTime);
  /*FILE *File;
  File=fopen ("ConstStep.dta", "a+");

  fprintf(File, "\n");
  fprintf(File, "Data for %d substeps. \n", *numberIter);
  for (Iter=Result.begin() ; Iter!=Result.end(); )
                  {
                          for (int i=0; i<9; i++)
                          {
                          fprintf( File, "%.20f , ",*Iter);
                          Iter++;
                          }
                  fprintf(File, "\n");
                  }
  fprintf(File, "\n");
  fclose (File); */
  return 0;
}

double
ShengMohrCoulomb::plasticMidpoint(StateMohrCoulomb* point, double* epStrain,
                                  double* absStress, int* numberIter)
{
  /*this procedure is to calculate the stress increase using the midpoint method
  with given number of iterations numberIter.
  is made mainly for use in the extrapolation procedure. It has no adaptive
  substepping or error control. It just integrates
  the strain to get the stress in given number of substeps using the midpoint
  method
  */
  // clock_t StartTime, EndTime;
  // StartTime=clock();

  StateMohrCoulomb NewPoint, MidPoint;
  BBMMatrix SIGMA(6, 1);
  // vector <double> Result;
  vector<double>::iterator Iter;

  double DSigma[6], CurrentStrain[7], HalfCurrentStrain[7], plasticStrain[6];
  double h, dpZeroStar = 0;
  double dStressDrift[6], dLambdaDrift;

  h = *numberIter;
  h = 1 / h;
  for (int i = 0; i < 7; i++)
    absStress[i] = 0; // 7 components...
  for (int i = 0; i < 7; i++) {
    CurrentStrain[i] = epStrain[i] * h; // strain increment in current step
    HalfCurrentStrain[i] = 0.5 * CurrentStrain[i];
  }
  // cout<<"Step Length="<<StepLength<<"\n";
  // cout<<"Current strain [0]="<<CurrentStrain[0]<<"\n";
  /*
  for (int i=0; i<6; i++) Result.push_back(Point->stress[i]);
  Result.push_back(Point->p0Star());
  Result.push_back(Point->getmeanStress());
  Result.push_back(Point->shearStress());
  */
  for (int loop = 0; loop < *numberIter; loop++) {
    Point->Copy(&NewPoint); // point is unchanged in calcPlastic procedure
    Point->Copy(&MidPoint);
    calcPlastic(NewPoint, HalfCurrentStrain, &SIGMA, plasticStrain, &dpZeroStar,
                computeYieldFunctionNN(Point), dStressDrift,
                &dLambdaDrift); // calculate the plastic stresses...
    for (int i = 0; i < 6; i++)
      DSigma[i] = SIGMA.getElement(i + 1, 1);
    MidPoint.Update(plasticStrain, HalfCurrentStrain, DSigma, dpZeroStar);
    calcPlastic(MidPoint, HalfCurrentStrain, &SIGMA, plasticStrain, &dpZeroStar,
                computeYieldFunctionNN(Point), dStressDrift,
                &dLambdaDrift); // calculate the plastic stresses...
    for (int i = 0; i < 6; i++) {
      DSigma[i] = 2 * SIGMA.getElement(i + 1, 1);
      plasticStrain[i] = 2 * plasticStrain[i];
    }
    dpZeroStar = 2 * dpZeroStar;

    for (int i = 0; i < 6; i++)
      absStress[i] = absStress[i] + fabs(DSigma[i]);
    absStress[6] = absStress[6] + fabs(dpZeroStar);
    Point->Update(plasticStrain, CurrentStrain, DSigma, dpZeroStar);
    /*for (int i=0; i<6; i++) Result.push_back(Point->stress[i]);
    Result.push_back(Point->p0Star());
    Result.push_back(Point->getmeanStress());
    Result.push_back(Point->shearStress());*/
  }
  // EndTime=clock();
  // return (EndTime-StartTime);
  /*FILE *File;
  File=fopen ("ConstStep.dta", "a+");

  fprintf(File, "\n");
  fprintf(File, "Data for %d substeps. \n", *numberIter);
  for (Iter=Result.begin() ; Iter!=Result.end(); )
                  {
                          for (int i=0; i<9; i++)
                          {
                          fprintf( File, "%.20f , ",*Iter);
                          Iter++;
                          }
                  fprintf(File, "\n");
                  }
  fprintf(File, "\n");
  fclose (File); */
  return 0;
}

double
ShengMohrCoulomb::plasticExtrapol(StateMohrCoulomb* point, double* epStrain,
                                  double* absStress, int* numberIter)
{
  // Here the d_driftCorrection parameter is used. 0=no correction, 1 -
  // correction at the beginning (A),
  // 2 - correction at the end (B), 3 - Zero Drift Algorithm.
  /*
  switch (int(d_driftCorrection))
  {
  case 1 :
          {
                  cout<<"Procedure proceed with standard algorithm and no drift
  correction."<<"\n";
                  break;
          }
  case 2 :
          {
                  cout<<"Procedure proceed with standard algorithm and drift
  correction at the beginning (point A)."<<"\n";
                  break;
          }

  default :
          {
                  cout<<"Unknown d_driftCorrection parameter. Parameter read
  is:"<<d_driftCorrection<<"\n";
                  break;
          }
  }
  */

  // no drift correction as may worsen the results. at least possibly.

  clock_t StartTime, EndTime;
  StartTime = clock();

  int STEPMAX = 15;
  // int
  // DivisionsInt[15]={2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768};
  // int DivisionsInt[15]={2,4,6,8,10,12,14,16,18,20,22,24,26,28,30};
  // int DivisionsInt[15]={2,4,6,8,12,16,24,32,48,64,96,128,192,256,384};
  // int DivisionsInt[15]={2,4,6,8,12,16,24,32,48,64,96,128,192,256,384};
  // int
  // DivisionsInt[15]={12,16,24,32,48,64,96,128,160,192,256,320,384,448,512};
  // //fastests so far
  int DivisionsInt[15] = {
    32, 48, 64, 96, 128, 160, 192, 256, 320, 384, 448, 512, 608, 736, 992
  }; // fastests so far
  // int
  // DivisionsInt[15]={2,4,6,8,10,14,20,28,38,52,72,100,138,190,262,362,500,690,952};
  // //N+N-4
  // int
  // DivisionsInt[15]={2,4,6,8,10,12,16,22,30,40,52,68,90,120,170,222,290,380,500,670,892};
  // //N+N-5
  // int
  // DivisionsInt[15]={28,32,40,52,64,78,94,120,154,200,240,290,330,380,440};
  // int
  // DivisionsInt[15]={32,40,52,64,78,94,120,154,200,240,290,330,380,440,520};
  // int DivisionsInt[15]={32,36,40,46,52,60,70,82,96,112,130,150,176,220,380};
  // int
  // DivisionsInt[15]={32,36,40,52,68,92,114,154,200,240,290,330,380,440,520};
  // int DivisionsInt[15]={32,36,40,44,50,56,62,68,76,84,92,102,112,124,136};
  // //n=n-1*1.1, doesn't converge too often
  // int
  // DivisionsInt[15]={32,38,46,56,66,80,96,114,138,166,198,238,286,344,412};
  // //n=n-1*1.2
  double Zeros[6] = { 0, 0, 0, 0, 0, 0 };

  double Divisions[15];
  for (int i = 0; i < 15; i++)
    Divisions[i] = DivisionsInt[i];

  double ApproximationTable[16][20]; // six stresses, p0*, six absolute
                                     // stresses, absolute p0*, 6 plastic
                                     // strains
  double ApproximationTableOld[16][20];
  double hSquareTable[16][15];
  double DError[15], RError;
  double InitialStress[7], InitialplasticStrain[6], plasticStrain[6];
  double DSigma[6];
  /*
  double TestTable[5]={0.375,0.37109375,0.36945588,0.36879683,0.36829712};
  double TestApproxTable[5]; double TestApproxTableOld[5];
  for (int i=0; i<5; i++)
  {
          TestApproxTable[i]=0;
          TestApproxTableOld[i]=0;
  }
  */

  StateMohrCoulomb NewPoint, CopyPoint, OldPoint;
  Point->Copy(&CopyPoint);
  Point->Copy(&OldPoint);

  for (int i = 0; i < 16; i++) {
    for (int j = 0; j < 20; j++) {
      ApproximationTable[i][j] = 0;
      ApproximationTableOld[i][j] = 0;
    }
    for (int j = 0; j < 15; j++)
      hSquareTable[i][j] = 0;
  }

  for (int i = 1; i < 16; i++) {
    for (int j = 0; j < i; j++) {
      hSquareTable[i][j] =
        (Divisions[i] / Divisions[j]) * (Divisions[i] / Divisions[j]);
      // cout<<"Divisions["<<i<<"]="<<Divisions[i]<<"
      // Divisions["<<j<<"]="<<Divisions[j]<<"\n";
      // cout<<"hSquareTable["<<i<<"]["<<j<<"]="<<hSquareTable[i][j]<<"\n";
    }
  }
  for (int i = 0; i < 6; i++) {
    InitialStress[i] = Point->stress[i];
    InitialplasticStrain[i] = Point->plastic_strain[i];
  }
  InitialStress[6] = Point->p0Star();

  int loop = 0;
  for (; loop < STEPMAX; loop++) {
    // calculating stress increment using MidPoint rule

    CopyPoint.Copy(&NewPoint);
    plasticMidpoint(&NewPoint, epStrain, absStress,
                    &DivisionsInt[loop]); // calculation of stress increment
                                          // using the midpoint procedure
    // TestApproxTable[0]=TestTable[loop];

    // saving data using MidPoint rule
    for (int i = 0; i < 6; i++) {
      ApproximationTable[0][i] = NewPoint.stress[i] - InitialStress[i];
      ApproximationTable[0][i + 7] = absStress[i];
      ApproximationTable[0][i + 14] =
        NewPoint.plastic_strain[i] - InitialplasticStrain[i];
    }
    ApproximationTable[0][6] = NewPoint.p0Star() - InitialStress[6];
    ApproximationTable[0][13] = absStress[6];

    // all data from the Midpoint Rule saved in the ResultsTable.

    for (int i = 0; i < loop; i++) {

      for (int j = 0; j < 20; j++) {
        ApproximationTable[i + 1][j] =
          ApproximationTable[i][j] +
          (ApproximationTable[i][j] - ApproximationTableOld[i][j]) /
            (hSquareTable[loop][loop - i - 1] - 1);
      }
    }

    // for (i=0; i<loop+1; i++) cout<<ApproximationTable[i][0]<<"  ";
    // cout<<"\n"<<"\n";

    // approximations are calculated
    // two possibilities of error control. In literature rather the second one;
    // this is the more stringent, but first one should be enough
    //(1) use ApproximationTable[loop][i] - ApproximationTable[loop-1][i]
    //(2) use ApproximationTable[loop][i] - ApproximationTableOld [loop-1][i]
    // still using relative error definition...
    // OldPoint (used for calculating the norm in q) is set accordingly

    // FIRST ONE OK, SEE Deuflhard Bornemann 2002 Springer p.206

    RError = 0;
    if (loop > 0) {
      for (int i = 0; i < 6; i++)
        DSigma[i] = ApproximationTable[loop][i];
      for (int i = 0; i < 7; i++) {
        // DError[i]=ApproximationTable[loop][i] -
        // ApproximationTableOld[loop-1][i];
        DError[i] =
          ApproximationTable[loop][i] - ApproximationTable[loop - 1][i];
      }

      switch (int(d_tolMethod)) {
        case 0: {
          RError = checkNorm(DSigma, ApproximationTable[loop][6], &CopyPoint,
                             DError); // returns RError
        } break;
        case 1: {
          // SLOAN NORM
          // cout<<"SLOAN NORM"<<"\n";
          RError =
            checkNormSloan(DSigma, ApproximationTable[loop][6], &CopyPoint,
                           DError); // returns RError
        } break;
        default: {
          cout << "ERROR !!!! Improper d_tolMethod in increment.dta" << "\n";
          getchar();
        }
      } // end switch
    } else
      RError = 2 * d_integrationTol;

    // cout<<"Relative error after iteration "<<loop<<" is equal
    // to:"<<RError<<"\n";
    for (int i = 0; i < loop + 1; i++)
      for (int j = 0; j < 20; j++)
        ApproximationTableOld[i][j] = ApproximationTable[i][j];
    // for (int i=0; i<loop+1; i++) cout<<ApproximationTable[i][0]<<"  ";
    // cout<<"\n";
    if (RError < d_integrationTol) {
      // error less then requested...check for drift and add the error
      CopyPoint.Copy(&OldPoint);
      OldPoint.Update(Zeros, epStrain, DSigma, ApproximationTable[loop][6]);
      OldPoint.Copy(&NewPoint);
      if (d_driftCorrection == 3)
        correctDrift(&NewPoint);
      if (d_driftCorrection == 2)
        correctDriftBeg(
          &NewPoint,
          &CopyPoint); // value of OldPoint copied before updating the point
      for (int i = 0; i < 6; i++)
        DError[i] = DError[i] + (NewPoint.stress[i] - OldPoint.stress[i]);
      DError[6] = DError[6] + NewPoint.p0Star() - OldPoint.p0Star();
      switch (int(d_tolMethod)) {
        case 0: {
          RError = checkNorm(DSigma, ApproximationTable[loop][6], &CopyPoint,
                             DError); // returns RError
        } break;
        case 1: {
          // SLOAN NORM
          // cout<<"SLOAN NORM"<<"\n";
          RError =
            checkNormSloan(DSigma, ApproximationTable[loop][6], &CopyPoint,
                           DError); // returns RError
        } break;
        default: {
          cout << "ERROR !!!! Improper d_tolMethod in increment.dta" << "\n";
          getchar();
        }
      } // end switch
      if (RError < d_integrationTol) {
        loop = loop + 100; // no more interations after - finished
      } else {
        // one more loop, do not need to change anything, as point anyway copied
        // back later on.
      }
    } else { // NewPoint.Copy(&OldPoint);
    }
  }

  if (loop > 100) {
    loop = loop - 101;
    // done - time to update everything and sent time and no of iterations...

    for (int i = 0; i < 6; i++) {
      absStress[i] = ApproximationTable[loop][i + 7];
      plasticStrain[i] = ApproximationTable[loop][i + 14];
      // cout<<"DSigma="<<DSigma[i]<<"   absStress="<<absStress[i]<<"\n";
    }
    absStress[6] = ApproximationTable[loop][13];
    // cout<<absStress[6]<<"\n";
    Point->Update(plasticStrain, epStrain, DSigma, ApproximationTable[loop][6]);
    EndTime = clock();
    *numberIter = 0;
    for (int i = 0; i < loop + 1; i++)
      *numberIter = *numberIter + Divisions[i];
    // cout<<"Procedure has coverged after"<<*numberIter<<" iterations."<<"\n";
    // getchar();
  } else {
    loop--;
    double DSigma[6];
    for (int i = 0; i < 6; i++) {
      DSigma[i] = ApproximationTable[loop][i];
      absStress[i] = ApproximationTable[loop][i + 7];
      plasticStrain[i] = ApproximationTable[loop][i + 14];
      // cout<<DSigma[i]<<"\n";
    }
    absStress[6] = ApproximationTable[loop][13];
    Point->Update(plasticStrain, epStrain, DSigma, ApproximationTable[loop][6]);
    EndTime = clock();
    *numberIter = 0;
    for (int i = 0; i < loop + 1; i++)
      *numberIter = *numberIter + Divisions[i];
    // cout<<"Procedure has NOT CONVERGED after"<<*numberIter<<"
    // iterations."<<"\n";
    // getchar();
  }
  // cout<<"calculation took:"<<double(EndTime-StartTime)/CLOCKS_PER_SEC<<"
  // s."<<"\n";

  /*
  switch (int(d_driftCorrection))
  {case 1 :{cout<<"Values calculated with standard algorithm and no drift
  correction."<<"\n";break;}
  case 2 :{cout<<"Values calculated with standard algorithm and drift correction
  at the beginning (point A)."<<"\n";break;}
  case 3 :{cout<<"Values calculated with standard algorithm and drift correction
  at the end (point B)."<<"\n";break;}
  case 4 :{cout<<"Values calculated with zero drift algorithm."<<"\n";break;}
  default :{cout<<"Unknown d_driftCorrection parameter. Parameter read
  is:"<<d_driftCorrection<<"\n";break;}
  }

  cout<<"The d_integrationTol parameter is equal to:"<<d_integrationTol<<"\n";
  cout<<"Total number of steps done:"<<stepNo<<"\n";
  cout<<"Total parameter lambda change: "<<LambdaTot<<"\n";
  cout<<"Total PZeroStar change: "<<PZeroStarTot<<"\n";
  checkYield (Point->state, point->stress, point->strain[6],&fValue);
  cout<<"Yield Function value is:"<<fValue<<"\n";
  cout<<"Number of drift correction performed:"<<correctdriftCount<<"\n";
  cout<<"Over the whole step change of stress is:"<<"\n";
  for (int i=0; i<6; i++)
  {
          cout<<"s["<<i<<"]="<<Point->stress[i]-OldPoint.stress[i]<<"\n";
  }

  //double DEpsV=0, dNu, SpecificVolume;

  //SpecificVolume=OldPoint.getSpecVol();
  //LambdaS=LambdaZero*(r+(1-r)*exp(-Beta*OldPoint.suction()));



  cout<<"Initial specific volume="<<OldPoint.getSpecVol()<<"\n";
  cout<<"Final specific volume="<<Point->getSpecVol()<<"\n";
  cout<<"Change of specific volume is equal
  to="<<(OldPoint.getSpecVol()-Point->getSpecVol())<<"\n";
  cout<<"Change of mean stress p is equal
  to:"<<Point->getmeanStress()-OldPoint.getmeanStress()<<"\n";
  cout<<"Change of shear stress q is equal
  to:"<<Point->shearStress()-OldPoint.shearStress()<<"\n";



  FILE * StressFile;
  FILE * StrainFile;
  FILE * PZeroFile;
  StressFile = fopen( "stressInc.dta", "w" );

  for (Iter=stressInc.begin() ; Iter!=stressInc.end(); )
  {
  for (int i=0; i<6; i++)
          {
          fprintf( StressFile, "%.20f , ",*Iter);
          Iter++;
          }
  fprintf( StressFile, "\n");
  }
  fclose (StressFile);

  StrainFile = fopen( "strainInc.dta", "w" );

  for (Iter=strainInc.begin() ; Iter!=strainInc.end(); )
  {
  for (int i=0; i<7; i++)
          {
          fprintf( StrainFile, "%.20f , ",*Iter);
          Iter++;
          }
  fprintf( StrainFile, "\n");
  }
  fclose (StrainFile);


  PZeroFile = fopen( "PZeroInc.dta", "w" );
  for (Iter=PZeroLambda.begin() ; Iter!=PZeroLambda.end(); )
  {
  for (int i=0; i<2; i++)
          {
          fprintf( PZeroFile, "%.25f , ",*Iter);
          Iter++;
          }
  fprintf( PZeroFile, "\n");
  }
  fclose (PZeroFile);
  */

  // cout<<"Press any key..."<<"\n";
  // getchar();
  //*numberIter=stepNo;
  /*
  FILE *File;
  File=fopen ("ConstStep.dta", "a+");

  fprintf(File, "\n");
  fprintf(File, "Extrapolation method Results.\n", *numberIter);
  for (int i=0; i<7;i++ )
                  {
                          fprintf( File, "%.20f ,
  ",ApproximationTable[loop][i]);
                  }
  fprintf( File, "%.20f ,
  ",(ApproximationTable[loop][0]+ApproximationTable[loop][1]+ApproximationTable[loop][2])/3);
  fprintf( File, "%.20f ,
  ",ApproximationTable[loop][0]-ApproximationTable[loop][1]);

  fprintf(File, "\n");
  fclose (File); */
  return (EndTime - StartTime);
}

double
ShengMohrCoulomb::doRKExtrapolation(double A[][8], double* B, double* BRes,
                                  double* C, StateMohrCoulomb* point, double* epStrain,
                                  double* absStress, StateMohrCoulomb* OldPoint,
                                  double* RelError, int* numberIter,
                                  double methodOrder, int methodSteps,
                                  bool errorEstimate)
{
  // Here the d_driftCorrection parameter is used. 0=no correction, 1 -
  // correction at the beginning (A),
  // 2 - correction at the end (B), 3 - Zero Drift Algorithm.
  /*
  switch (int(d_driftCorrection))
  {
  case 1 :
          {
                  cout<<"Procedure proceed with standard algorithm and no drift
  correction."<<"\n";
                  break;
          }
  case 2 :
          {
                  cout<<"Procedure proceed with standard algorithm and drift
  correction at the beginning (point A)."<<"\n";
                  break;
          }

  default :
          {
                  cout<<"Unknown d_driftCorrection parameter. Parameter read
  is:"<<d_driftCorrection<<"\n";
                  break;
          }
  }
  */

  // no drift correction as may worsen the results. at least possibly.

  clock_t StartTime, EndTime;
  StartTime = clock();

  int STEPMAX = 15;
  // int
  // DivisionsInt[15]={2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768};
  // int DivisionsInt[15]={2,4,6,8,10,12,14,16,18,20,22,24,26,28,30};
  int DivisionsInt[15] = { 2,  4,  6,  8,   12,  16,  24, 32,
                           48, 64, 96, 128, 192, 256, 384 };
  // int
  // DivisionsInt[15]={32,48,64,96,128,160,192,256,320,384,448,512,608,736,992};
  // //fastests so far
  // int
  // DivisionsInt[15]={2,4,6,8,10,14,20,28,38,52,72,100,138,190,262,362,500,690,952};
  // //N+N-4
  // int
  // DivisionsInt[15]={2,4,6,8,10,12,16,22,30,40,52,68,90,120,170,222,290,380,500,670,892};
  // //N+N-5
  // int
  // DivisionsInt[15]={28,32,40,52,64,78,94,120,154,200,240,290,330,380,440};
  // int
  // DivisionsInt[15]={32,40,52,64,78,94,120,154,200,240,290,330,380,440,520};
  // int DivisionsInt[15]={32,36,40,46,52,60,70,82,96,112,130,150,176,220,380};
  // int
  // DivisionsInt[15]={32,36,40,52,68,92,114,154,200,240,290,330,380,440,520};
  // int DivisionsInt[15]={32,36,40,44,50,56,62,68,76,84,92,102,112,124,136};
  // //n=n-1*1.1, doesn't converge too often
  // int
  // DivisionsInt[15]={32,38,46,56,66,80,96,114,138,166,198,238,286,344,412};
  // //n=n-1*1.2

  double Divisions[15];
  for (int i = 0; i < 15; i++)
    Divisions[i] = DivisionsInt[i];

  double ApproximationTable[16][20]; // 6 stresses, p0*, 6 absolute stresses,
                                     // absolute p0*, 6 plastic strains
  double ApproximationTableOld[16][20];
  double hSquareTable[16][15];
  double DError[15], RError;
  double InitialStress[7], InitialplasticStrain[6], plasticStrain[6];
  double DSigma[6];
  bool StepAccepted = false;
  /*
  double TestTable[5]={0.375,0.37109375,0.36945588,0.36879683,0.36829712};
  double TestApproxTable[5]; double TestApproxTableOld[5];
  for (int i=0; i<5; i++)
  {
          TestApproxTable[i]=0;
          TestApproxTableOld[i]=0;
  }
  */

  StateMohrCoulomb NewPoint, CopyPoint;
  Point->Copy(&CopyPoint);
  // Point->Copy(&OldPoint);

  for (int i = 0; i < 16; i++) {
    for (int j = 0; j < 20; j++) {
      ApproximationTable[i][j] = 0;
      ApproximationTableOld[i][j] = 0;
    }
    for (int j = 0; j < 15; j++)
      hSquareTable[i][j] = 0;
  }

  for (int i = 1; i < 16; i++) {
    for (int j = 0; j < i; j++) {
      hSquareTable[i][j] =
        (Divisions[i] / Divisions[j]) * (Divisions[i] / Divisions[j]);
      // cout<<"Divisions["<<i<<"]="<<Divisions[i]<<"
      // Divisions["<<j<<"]="<<Divisions[j]<<"\n";
      // cout<<"hSquareTable["<<i<<"]["<<j<<"]="<<hSquareTable[i][j]<<"\n";
    }
  }
  for (int i = 0; i < 6; i++) {
    InitialStress[i] = Point->stress[i];
    InitialplasticStrain[i] = Point->plastic_strain[i];
  }
  InitialStress[6] = Point->p0Star();

  for (int i = 0; i < 6; i++) {
    ApproximationTable[0][i] = OldPoint->stress[i] - InitialStress[i];
    ApproximationTable[0][i + 7] = absStress[i];
    ApproximationTable[0][i + 14] =
      OldPoint->plastic_strain[i] - InitialplasticStrain[i];
  }
  ApproximationTable[0][6] = OldPoint->p0Star() - InitialStress[6];
  ApproximationTable[0][13] = absStress[6];

  int loop = 1;

  // for (int i=0; i<loop+1; i++) cout<<ApproximationTable[i][0]<<"  ";
  // cout<<"\n";

  for (int i = 0; i < loop; i++) {

    for (int j = 0; j < 20; j++) {
      ApproximationTable[i + 1][j] =
        ApproximationTable[i][j] +
        (ApproximationTable[i][j] - ApproximationTableOld[i][j]) /
          (hSquareTable[loop][loop - i - 1] - 1);
    }
  }
  for (int i = 0; i < loop + 1; i++)
    for (int j = 0; j < 20; j++)
      ApproximationTableOld[i][j] = ApproximationTable[i][j];

  // cout<<"Initial"<<"\n";
  // for (int i=0; i<loop; i++) cout<<ApproximationTable[i][0]<<"  ";
  // cout<<"\n"<<"\n";
  // loop=0;

  for (; loop < STEPMAX; loop++) {
    // calculating stress increment using MidPoint rule

    CopyPoint.Copy(&NewPoint);
    doRungeuttaEqualStep(A, B, BRes, C, &NewPoint, epStrain, absStress, &RError,
                        DivisionsInt[loop], methodOrder, methodSteps,
                        errorEstimate); // calculation of stress increment using
                                        // the RK procedure
    // TestApproxTable[0]=TestTable[loop];

    // saving data using MidPoint rule
    for (int i = 0; i < 6; i++) {
      ApproximationTable[0][i] = NewPoint.stress[i] - InitialStress[i];
      ApproximationTable[0][i + 7] = absStress[i];
      ApproximationTable[0][i + 14] =
        NewPoint.plastic_strain[i] - InitialplasticStrain[i];
    }
    ApproximationTable[0][6] = NewPoint.p0Star() - InitialStress[6];
    ApproximationTable[0][13] = absStress[6];

    // all data from the Midpoint Rule saved in the ResultsTable.

    for (int i = 0; i < loop; i++) {

      for (int j = 0; j < 20; j++) {
        ApproximationTable[i + 1][j] =
          ApproximationTable[i][j] +
          (ApproximationTable[i][j] - ApproximationTableOld[i][j]) /
            (hSquareTable[loop][loop - i - 1] - 1);
      }
    }

    // for (i=0; i<loop+1; i++) cout<<ApproximationTable[i][0]<<"  ";
    // cout<<"\n";
    // getchar();

    // approximations are calculated
    // two possibilities of error control. In literature rather the second one;
    // this is the more stringent, but first one should be enough
    //(1) use ApproximationTable[loop][i] - ApproximationTable[loop-1][i]
    //(2) use ApproximationTable[loop][i] - ApproximationTableOld [loop-1][i]
    // still using relative error definition...
    // OldPoint (used for calculating the norm in q) is set accordingly

    StepAccepted = false;
    if (RError < d_integrationTol)
      StepAccepted = true;

    for (int i = 0; i < 6; i++)
      DSigma[i] = ApproximationTable[loop][i];
    for (int i = 0; i < 7; i++) {
      DError[i] =
        ApproximationTable[loop][i] - ApproximationTableOld[loop - 1][i];
    }

    switch (int(d_tolMethod)) {
      case 0: {
        RError = checkNorm(DSigma, ApproximationTable[loop][6], OldPoint,
                           DError); // returns RError
      } break;
      case 1: {
        // SLOAN NORM
        // cout<<"SLOAN NORM"<<"\n";
        RError = checkNormSloan(DSigma, ApproximationTable[loop][6], OldPoint,
                                DError); // returns RError
      } break;
      default: {
        cout << "ERROR !!!! Improper d_tolMethod in increment.dta" << "\n";
        getchar();
      }
    } // end switch

    // cout<<"Relative error after iteration "<<loop<<" is equal
    // to:"<<RError<<"\n";
    for (int i = 0; i < loop + 1; i++)
      for (int j = 0; j < 20; j++)
        ApproximationTableOld[i][j] = ApproximationTable[i][j];
    // for (int i=0; i<loop+1; i++) cout<<ApproximationTable[i][0]<<"  ";
    // cout<<"\n";
    if (RError < d_integrationTol)
      StepAccepted = true;
    if (StepAccepted)
      loop = loop + 100; // no more interations after - finished
    else
      NewPoint.Copy(OldPoint);
  }

  if (loop > 100) {
    loop = loop - 101;
    // done - time to update everything and sent time and no of iterations...

    for (int i = 0; i < 6; i++) {
      absStress[i] = ApproximationTable[loop][i + 7];
      plasticStrain[i] = ApproximationTable[loop][i + 14];
      // cout<<"DSigma="<<DSigma[i]<<"   absStress="<<absStress[i]<<"\n";
    }
    absStress[6] = ApproximationTable[loop][13];
    // cout<<absStress[6]<<"\n";
    Point->Update(plasticStrain, epStrain, DSigma, ApproximationTable[loop][6]);
    EndTime = clock();
    *numberIter = 0;
    for (int i = 0; i < loop + 1; i++)
      *numberIter = *numberIter + Divisions[i];
    // cout<<"Procedure has coverged after"<<*numberIter<<" iterations."<<"\n";
    // getchar();
    return (EndTime - StartTime);
  } else {
    loop--;
    double DSigma[6];
    for (int i = 0; i < 6; i++) {
      DSigma[i] = ApproximationTable[loop][i];
      absStress[i] = ApproximationTable[loop][i + 7];
      plasticStrain[i] = ApproximationTable[loop][i + 14];
      // cout<<DSigma[i]<<"\n";
    }
    absStress[6] = ApproximationTable[loop][13];
    Point->Update(plasticStrain, epStrain, DSigma, ApproximationTable[loop][6]);
    EndTime = clock();
    *numberIter = 0;
    for (int i = 0; i < loop + 1; i++)
      *numberIter = *numberIter + Divisions[i];
    cout << "Procedure has NOT CONVERGED after" << *numberIter << " iterations."
         << "\n";
    // getchar();
  }
  return (EndTime - StartTime);
}

double ShengMohrCoulomb::checkNorm(double* DSigma, double dpZeroStar,
                                   StateMohrCoulomb* initialState,
                                   double* DError) // returns RError
{
  /*
  Procedure returns value of the relative error RError
  The value is calculated as maximum of all the values of RErrors
  RErrors used:
          - each stress component
          - error in p0*
          - error in p
          - error in q

  INPUT: DSigma [6] - values of stress increment (min 6, but may be more; values
  over six ignored)
                  dpZeroStar - value of the p0* increment
                  initialState - point used to calculate the values of the step
                  DError [7] - error estimates for all values of DSigma
  DError[i]  - error estimate for DSigma[i];
  */

  double DRError[9];
  double RError = 0;

  // standard norm:
  for (int i = 0; i < 6; i++) {
    if (fabs(DSigma[i]) > TINY)
      DRError[i] = DError[i] / DSigma[i];
    else
      DRError[i] = DError[i] / TINY;
  }
  if (fabs(dpZeroStar) > TINY)
    DRError[6] = DError[6] / dpZeroStar;
  else
    DRError[6] = DError[6] / TINY;

  // norm in mean stress p:

  double P, ErrorP;

  P = fabs(DSigma[0] + DSigma[1] + DSigma[2]) / 3;
  ErrorP = fabs(DError[0] + DError[1] + DError[2]) / 3;
  if (P > TINY)
    DRError[7] = ErrorP / P;
  else
    DRError[7] = ErrorP / TINY;

  // norm of q...

  double InitialShear, FinalShearMin, FinalShearMax, FinalShear;
  double ShearError, DShear;
  double SigmaTemp[6], InitialSigma[6];

  InitialShear = initialState->shearStress();

  for (int i = 0; i < 6; i++) {
    InitialSigma[i] = initialState->stress[i];
    SigmaTemp[i] = InitialSigma[i] + DSigma[i];
  }

  // computing shear stress like this is much more effective than calling the
  // point to give the shear stress
  // as the point would have to be updated and copied. It is not time effective,
  // as there are more
  // things to update than just the stresses; Even if zeros are given, they are
  // still updated what is time consuming
  // This is less elegant but more efficient

  FinalShear = (SigmaTemp[0] - SigmaTemp[1]) * (SigmaTemp[0] - SigmaTemp[1]) +
               (SigmaTemp[0] - SigmaTemp[2]) * (SigmaTemp[0] - SigmaTemp[2]) +
               (SigmaTemp[1] - SigmaTemp[2]) * (SigmaTemp[1] - SigmaTemp[2]);
  FinalShear = FinalShear +
               6 * (SigmaTemp[3] * SigmaTemp[3] + SigmaTemp[4] * SigmaTemp[4] +
                    SigmaTemp[5] * SigmaTemp[5]);
  FinalShear = sqrt(0.5 * FinalShear);
  DShear = fabs(FinalShear - InitialShear);

  for (int i = 0; i < 6; i++)
    SigmaTemp[i] = SigmaTemp[i] + DError[i];
  FinalShearMax =
    (SigmaTemp[0] - SigmaTemp[1]) * (SigmaTemp[0] - SigmaTemp[1]) +
    (SigmaTemp[0] - SigmaTemp[2]) * (SigmaTemp[0] - SigmaTemp[2]) +
    (SigmaTemp[1] - SigmaTemp[2]) * (SigmaTemp[1] - SigmaTemp[2]);
  FinalShearMax =
    FinalShearMax +
    6 * (SigmaTemp[3] * SigmaTemp[3] + SigmaTemp[4] * SigmaTemp[4] +
         SigmaTemp[5] * SigmaTemp[5]);
  FinalShearMax = sqrt(0.5 * FinalShearMax);

  for (int i = 0; i < 6; i++)
    SigmaTemp[i] = SigmaTemp[i] - 2 * DError[i];
  FinalShearMin =
    (SigmaTemp[0] - SigmaTemp[1]) * (SigmaTemp[0] - SigmaTemp[1]) +
    (SigmaTemp[0] - SigmaTemp[2]) * (SigmaTemp[0] - SigmaTemp[2]) +
    (SigmaTemp[1] - SigmaTemp[2]) * (SigmaTemp[1] - SigmaTemp[2]);
  FinalShearMin =
    FinalShearMin +
    6 * (SigmaTemp[3] * SigmaTemp[3] + SigmaTemp[4] * SigmaTemp[4] +
         SigmaTemp[5] * SigmaTemp[5]);
  FinalShearMin = sqrt(0.5 * FinalShearMin);

  ShearError = fabs(FinalShearMax - FinalShear);
  if (fabs(FinalShearMin - FinalShear) > ShearError)
    ShearError = fabs(FinalShearMin - FinalShear);
  if (DShear > TINY)
    DRError[8] = ShearError / DShear;
  else
    DRError[8] = ShearError / TINY;

  // final check: max norm
  for (int i = 0; i < 9; i++) {
    DRError[i] = fabs(DRError[i]);
    if (DRError[i] > RError)
      RError = DRError[i];
  }
  return RError;
}

double ShengMohrCoulomb::checkNormSloan(double* DSigma, double dpZeroStar,
                                        StateMohrCoulomb* initialState,
                                        double* DError) // returns RError
{
  /*
  Procedure returns value of the relative error RError
  The value is calculated as maximum of all the values of RErrors
  RErrors used:
          - each stress component
          - error in p0*
          - error in p
          - error in q

  INPUT: DSigma [6] - values of stress increment (min 6, but may be more; values
  over six ignored)
                  dpZeroStar - value of the p0* increment
                  initialState - point used to calculate the values of the step
                  DError [7] - error estimates for all values of DSigma
  DError[i]  - error estimate for DSigma[i];
  */

  double DRError[9];
  double RError = 0;
  double InitialSigma[6], SigmaTemp[6], dpZeroStarEnd;

  for (int i = 0; i < 6; i++) {
    InitialSigma[i] = initialState->stress[i];
    SigmaTemp[i] = InitialSigma[i] + DSigma[i];
  }
  dpZeroStarEnd = initialState->p0Star() + dpZeroStar;

  // standard norm:
  for (int i = 0; i < 6; i++) {
    if (fabs(SigmaTemp[i]) > TINY)
      DRError[i] = DError[i] / SigmaTemp[i];
    else
      DRError[i] = DError[i] / TINY;
  }
  if (fabs(dpZeroStarEnd) > TINY)
    DRError[6] = DError[6] / dpZeroStarEnd;
  else
    DRError[6] = DError[6] / TINY;

  // norm in mean stress p:

  double P, ErrorP;

  P = fabs(InitialSigma[0] + InitialSigma[1] + InitialSigma[2] + DSigma[0] +
           DSigma[1] + DSigma[2]) /
      3;
  ErrorP = fabs(DError[0] + DError[1] + DError[2]) / 3;
  if (P > TINY)
    DRError[7] = ErrorP / P;
  else
    DRError[7] = ErrorP / TINY;

  // norm of q...

  double FinalShearMin, FinalShearMax, FinalShear; // InitialShear;
  double ShearError;

  // InitialShear=initialState->shearStress();

  // computing shear stress like this is much more effective than calling the
  // point to give the shear stress
  // as the point would have to be updated and copied. It is not time effective,
  // as there are more
  // things to update than just the stresses; Even if zeros are given, they are
  // still updated what is time consuming
  // This is less elegant but more efficient

  FinalShear = (SigmaTemp[0] - SigmaTemp[1]) * (SigmaTemp[0] - SigmaTemp[1]) +
               (SigmaTemp[0] - SigmaTemp[2]) * (SigmaTemp[0] - SigmaTemp[2]) +
               (SigmaTemp[1] - SigmaTemp[2]) * (SigmaTemp[1] - SigmaTemp[2]);
  FinalShear = FinalShear +
               6 * (SigmaTemp[3] * SigmaTemp[3] + SigmaTemp[4] * SigmaTemp[4] +
                    SigmaTemp[5] * SigmaTemp[5]);
  FinalShear = sqrt(0.5 * FinalShear);

  for (int i = 0; i < 6; i++)
    SigmaTemp[i] = SigmaTemp[i] + DError[i];
  FinalShearMax =
    (SigmaTemp[0] - SigmaTemp[1]) * (SigmaTemp[0] - SigmaTemp[1]) +
    (SigmaTemp[0] - SigmaTemp[2]) * (SigmaTemp[0] - SigmaTemp[2]) +
    (SigmaTemp[1] - SigmaTemp[2]) * (SigmaTemp[1] - SigmaTemp[2]);
  FinalShearMax =
    FinalShearMax +
    6 * (SigmaTemp[3] * SigmaTemp[3] + SigmaTemp[4] * SigmaTemp[4] +
         SigmaTemp[5] * SigmaTemp[5]);
  FinalShearMax = sqrt(0.5 * FinalShearMax);

  for (int i = 0; i < 6; i++)
    SigmaTemp[i] = SigmaTemp[i] - 2 * DError[i];
  FinalShearMin =
    (SigmaTemp[0] - SigmaTemp[1]) * (SigmaTemp[0] - SigmaTemp[1]) +
    (SigmaTemp[0] - SigmaTemp[2]) * (SigmaTemp[0] - SigmaTemp[2]) +
    (SigmaTemp[1] - SigmaTemp[2]) * (SigmaTemp[1] - SigmaTemp[2]);
  FinalShearMin =
    FinalShearMin +
    6 * (SigmaTemp[3] * SigmaTemp[3] + SigmaTemp[4] * SigmaTemp[4] +
         SigmaTemp[5] * SigmaTemp[5]);
  FinalShearMin = sqrt(0.5 * FinalShearMin);

  ShearError = fabs(FinalShearMax - FinalShear);
  if (fabs(FinalShearMin - FinalShear) > ShearError)
    ShearError = fabs(FinalShearMin - FinalShear);
  if (FinalShear > TINY)
    DRError[8] = ShearError / FinalShear;
  else
    DRError[8] = ShearError / TINY;

  // final check: max norm
  for (int i = 0; i < 9; i++) {
    DRError[i] = fabs(DRError[i]);
    if (DRError[i] > RError)
      RError = DRError[i];
  }
  return RError;
}

void
ShengMohrCoulomb::computeG1(StateMohrCoulomb* initialState, int retentionModel,
                            double* retentionParameters, double* G1)
{
  /*
  See Word document for details

  Retention Model:
  1 - State surface
  2 - Van Genuchten
  3 - Gallipoli Wheeler & Karstunen

  Retention Parameters:
  1: [a, b]
  2: irrelevant (matrix ==0)
  3: [fi, psi, n, m]

  G1- double[6];
  */
  switch (retentionModel) {
    case 1: // state surface model.
    {
      double a, b, s;
      a = retentionParameters[0];
      b = retentionParameters[1];
      s = initialState->suction();
      G1[0] = (1 - exp(b * s)) * a / 3;
      G1[1] = G1[0];
      G1[2] = G1[0];
      G1[3] = 0;
      G1[4] = 0;
      G1[5] = 0;
    } break;
    case 2: // Van Genuchten Model
    {
      for (int i = 0; i < 6; i++)
        G1[i] = 0;
    } break;
    case 3: // Gallipoli, Wheeler, Karstunen
    {
      cout << "calculation of G1 Matrix for Gallipoli Wheeler Karstunen model "
              "not implemented yet..."
           << "\n";
      for (int i = 0; i < 6; i++)
        G1[i] = 0;
    } break;
    default: {
      cout << "Procedure Compute G1. Unknown Retention Model... No matrix G1 "
              "calculated. Please press any key to continue..."
           << "\n";
      getchar();
    }
  }
}
void
ShengMohrCoulomb::computeG2(StateMohrCoulomb* initialState, int retentionModel,
                            double* retentionParameters, double* G2)
{
  /*
  See Word document for details

  Retention Model:
  1 - State surface
  2 - Van Genuchten
  3 - Gallipoli Wheeler & Karstunen

  Retention Parameters:
  1: [a, b]
  2: [Ew, Fw, Ssn, Sir]
  3: [fi, psi, n, m]

  G2- double;
  */

  switch (retentionModel) {
    case 1: // state surface model.
    {
      double a, b, s;
      a = retentionParameters[0];
      b = retentionParameters[1];
      s = initialState->suction();
      *G2 = b * exp(b * s) * (0.5 - a * initialState->getmeanStress());
    } break;
    case 2: // Van Genuchten Model
    {
      double Ew, Fw, Ssn, Sir, s, Numerator, Denominator;
      Ew = retentionParameters[0];
      Fw = retentionParameters[1];
      Ssn = retentionParameters[2];
      Sir = retentionParameters[3];
      s = initialState->suction();
      Numerator = (Ssn - Sir) * (1 - Fw) * Ew * pow(Ew * s, Fw - 1);
      Denominator = 1 + pow(Ew * s, Fw);
      Denominator = Fw * pow(Denominator, 1 - 1 / Fw);
      *G2 = Numerator / Denominator;
    } break;
    case 3: // Gallipoli, Wheeler, Karstunen
    {
      double Fi, psi, n, m, nu, s, Numerator, Denominator;
      Fi = retentionParameters[0];
      psi = retentionParameters[1];
      n = retentionParameters[2];
      m = retentionParameters[3];
      s = initialState->suction();
      nu = initialState->getSpecVol();
      Numerator = Fi * (1 - nu) * psi * s;
      Numerator = Fi * (1 - nu) * psi * pow(Numerator, n - 1);
      Denominator = Numerator * s + 1;
      Denominator = pow(Denominator, m + 1);
      *G2 = Numerator / Denominator;
    } break;
    default: {
      cout << "Procedure Compute G1. Unknown Retention Model... No matrix G1 "
              "calculated. Please press any key to continue..."
           << "\n";
      getchar();
    }
  }
}

double
ShengMohrCoulomb::getLambda(double* deriv, double stresspq[3],
                            double strainpq[3])
{
  double numerator[3], denominator[3];
  // get K from calcstresselast, hmmmm, later
  double K = 1;
  numerator[0] = K * strainpq[0];
  numerator[1] = 3 * G * strainpq[1];
  numerator[2] = strainpq[2];
  double DeltaPZero = 0;

  for (int i = 0; i < 3; i++)
    numerator[i] = numerator[i] * deriv[i]; // now it needed to be added to get
                                            // the first part of the numerator
  for (int i = 0; i < 3; i++)
    numerator[i] =
      numerator[i] + deriv[i] * DeltaPZero; // here it does nothing at the
                                            // moment, because DeltaPZero==0

  numerator[0] =
    numerator[0] + numerator[1] +
    numerator[2]; // numerator finished; it's a number stored in [0]

  denominator[0] = K * deriv[0];
  denominator[1] = 3 * G * deriv[1];
  denominator[2] = deriv[2];

  for (int i = 0; i < 3; i++)
    denominator[i] = denominator[i] * deriv[i];
  denominator[0] =
    denominator[0] + denominator[1] +
    denominator[2]; // denominator finished; it's a number stored in [0]

  return numerator[0] / denominator[0]; // numerator/denominator=lambda
}

double
ShengMohrCoulomb::getk()
{
  return 0; // k;
}
double
ShengMohrCoulomb::getLambdaZero()
{
  return 0; // LambdaZero;
}
double
ShengMohrCoulomb::getr()
{
  return 0; // r;
}
double
ShengMohrCoulomb::getBeta()
{
  return 0; // Beta;
}
double
ShengMohrCoulomb::getKappaP()
{
  return 0; // KappaP;
}

double
ShengMohrCoulomb::getpc()
{
  return 0; // pc;
}
