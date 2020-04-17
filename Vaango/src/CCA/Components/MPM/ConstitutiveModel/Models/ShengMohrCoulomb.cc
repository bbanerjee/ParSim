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
#include <chrono>

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
  // for (int i=0; i<3; i++) dbg << strainIncrement[i]<<"  "<<"\n";
  // for (int i=0; i<3; i++) dbg << finalState.stress[i]<<"  "<<"\n";

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
  double cosinus = findGradient(initialState->stress, stressInc, df, 0, 0);

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
  Matrix66 DEP_66 = calculateElasticTangentMatrix(K, G);

  Matrix67 DEP = Matrix67::Zero();
  DEP.block<6,6>(0,0) = DEP_66;
  return DEP;
}

Matrix66
ShengMohrCoulomb::calculateElasticTangentMatrix(double K, double G) const
{
  double K43G = K + 4.0 * G / 3.0;
  double K23G = K - 2.0 * G / 3.0;

  Matrix66 DEP = Matrix66::Zero();
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
ShengMohrCoulomb::findGradient(const Vector6& s, const Vector6& ds, 
                               Vector6s& df_dSigma, 
                               double /*suction*/, double /*dsuction*/) const
{
  // compute gradient
  Vector6 dg_dSigma;
  std::tie(df_dSigma, dg_dSigma) = computeDfDsigma(s);

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
    cosin += df_dSigma(i) * ds(i) / (dslength * dFlength); 
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
 * Actually compute the gradient of the yield finction wrt stress
 */
std::tuple<Vector6, Vector6>
ShengMohrCoulomb::computeDfDsigma(const Vector6& stress) const
{
  double I1 = firstInvariant(stress);
  double I2 = secondInvariant(stress);
  double J3 = thirdDevInvariant(stress);
  double J2, shearStress;
  std::tie(J2, shearStress) = vonMisesStress(stress);

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
  double M    = (3 - d_yield.d_sin_phi) / (6 * d_yield.d_alpha * d_yield.d_sin_phi);
  double Mpsi = (3 - d_yield.d_sin_psi) / (6 * d_yield.d_alpha * d_yield.d_sin_psi);

  Vector6 dJ2_dSig = Vector6::Zero();
  dJ2_dSig(0) =  (2 * s(0) - s(1) - s(2)) / 3.0;
  dJ2_dSig(1) =  (2 * s(1) - s(0) - s(2)) / 3.0;
  dJ2_dSig(2) =  (2 * s(2) - s(0) - s(1)) / 3.0;
  dJ2_dSig(3) =  2 * s(3);
  dJ2_dSig(4) =  2 * s(4);
  dJ2_dSig(5) =  2 * s(5);

  Vector6 df_dSigma = dJ2_dSig;
  Vector6 dg_dSigma = dJ2_dSig;
  df_dSigma *= (0.21022410391343 * M    * factor025 / shearStress);
  dg_dSigma *= (0.21022410391343 * Mpsi * factor025 / shearStress);

  for (int i = 0; i < 6; i++) {
    if (!std::isfinite(df_dSigma(i))) {
      std::ostringstream err;
      err << "*ERROR** df/dSigma not finite at df_dsigma(" << i
           << ")." << df_dSigma << "\n";
      throw InvalidValue(err.str(), __FILE__, __LINE__);
    }
  }

  if (dbg.active()) {
    dbg << "dF/dJ2 part\n";
    dbg << setprecision(15) << df_dSigma << "\n";
  }

  Vector6 dq_dSig = dJ2_dSig;
  dq_dSig *= std::sqrt(3.0) * 0.5 / std::sqrt(J2);

  Vector6 dq_dSig_g = dq_dSig;
  dq_dSig *= (-2.838025401481287 * (d_yield.alpha4 - 1) * d_yield.cohesion * M * J3 * factor075 /
               (shearStress * shearStress * shearStress * shearStress) +
              4.257038102221929 * (d_yield.alpha4 - 1) * M * std::sqrt(J2) * J3 * factor075 /
               (std::sqrt(3.0) * shearStress * shearStress * shearStress * shearStress));
  dq_dSig_g *= (-2.838025401481287 * (d_yield.alpha4 - 1) * d_yield.cohesion * Mpsi * J3 * factor075 /
                 (shearStress * shearStress * shearStress * shearStress) +
                4.257038102221929 * (d_yield.alpha4 - 1) * Mpsi * std::sqrt(J2) * J3 * factor075 /
                 (std::sqrt(3.0) * shearStress * shearStress * shearStress * shearStress));


  if (dbg.active()) {
    dbg << "dF/dq part\n";
    dbg << setprecision(15) << dq_dSig << "\n";
  }

  df_dSigma += dq_dSig;
  dg_dSigma += dq_dSig_g;

  Vector6 dI1_dSig = Vector6::Zero();
  dI1_dSig(0) =  1.0;
  dI1_dSig(1) =  1.0;
  dI1_dSig(2) =  1.0; //{1,1,1,0,0,0}

  Vector6 dJ3_dSig = dI1_dSig;
  dJ3_dSig *= (2.0 / 9.0 * I1 * I1 - I2 / 3.0);

  if (dbg.active()) {
    dbg << setprecision(15) << dJ3_dSig << "\n";
  }

  Vector6 dI2_dSig = Vector6::Zero();
  dI2_dSig(0) =  s(1) + s(2);
  dI2_dSig(1) =  s(0) + s(2);
  dI2_dSig(2) =  s(0) + s(1);
  dI2_dSig(3) =  -2 * s(3);
  dI2_dSig(4) =  -2 * s(4);
  dI2_dSig(5) =  -2 * s(5);

  if (dbg.active()) {
    dbg << setprecision(15) << dI2_dSig << "\n";
  }

  Vector6 dI3_dSig = Vector6::Zero(); 
  dI3_dSig(0) =  s(1) * s(2) - s(5) * s(5);
  dI3_dSig(1) =  s(0) * s(2) - s(4) * s(4);
  dI3_dSig(2) =  s(0) * s(1) - s(3) * s(3);
  dI3_dSig(3) =  2 * s(5) * s(4) - 2 * s(2) * s(3);
  dI3_dSig(4) =  2 * s(3) * s(5) - 2 * s(1) * s(4);
  dI3_dSig(5) =  2 * s(3) * s(4) - 2 * s(0) * s(5);

  if (dbg.active()) {
    dbg << setprecision(15) << dI3_dSig << "\n";
  }

  dJ3_dSig += (dI2_dSig * (-I1 / 3.0) + dI3_dSig);

  if (dbg.active()) {
    dbg << setprecision(15) << dJ3_dSig << "\n";
  }

  Vector6 dJ3_dSig_g = dJ3_dSig;
  dJ3_dSig *= (-1.419012700740643 * (d_yield.alpha4 - 1) * M * sqrt(J2) * factor075 /
                (std::sqrt(3.0) * shearStress * shearStress * shearStress) +
               0.94600846716043 * (d_yield.alpha4 - 1) * d_yield.cohesion * M * factor075 /
                (shearStress * shearStress * shearStress));
  dJ3_dSig_g *= (-1.419012700740643 * (d_yield.alpha4 - 1) * Mpsi * sqrt(J2) * factor075 /
                  (std::sqrt(3.0) * shearStress * shearStress * shearStress) +
                 0.94600846716043 * (d_yield.alpha4 - 1) * d_yield.cohesion * Mpsi * factor075 /
                  (shearStress * shearStress * shearStress));

  if (dbg.active()) {
    dbg << "dF/dJ3 part\n";
    dbg << setprecision(15) << dJ3_dSig << "\n";
  }

  df_dSigma += dJ3_dSig;
  dg_dSigma += dJ3_dSig_g;

  for (int i = 0; i < 6; i++) {
    if (!std::isfinite(df_dSigma(i))) {
      std::ostringstream err;
      err << "*ERROR** df/dSigma not finite at df_dsigma(" << i
           << ")." << df_dSigma << "\n";
      throw InvalidValue(err.str(), __FILE__, __LINE__);
    }
  }

  df_dsigma += dI1_dSig * (-1.0 / 3.0);
  dg_dsigma += dI1_dSig * (-1.0 / 3.0);

  for (int i = 0; i < 6; i++) {
    if (!std::isfinite(df_dSigma(i))) {
      std::ostringstream err;
      err << "*ERROR** df/dSigma not finite at df_dsigma(" << i
           << ")." << df_dSigma << "\n";
      throw InvalidValue(err.str(), __FILE__, __LINE__);
    }
  }

  if (dbg.active()) {
    dbg << "df_dSigma = " << df_dSigma << "\n";
  }

  return std::make_tuple(df_dsigma, dg_dsigma);
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

/**
  Calculate stress increment during plastic deformation
 */
double
ShengMohrCoulomb::calculatePlastic(const Vector7& purelyPlasticStrain, 
                                   const StateMohrCoulomb& state) const
{
  double time; 
  int    numIter;
  double stressIncrAbs[7];

  switch (d_int.d_solutionAlgorithm) {
    case SolutionAlgorithm::RUNGE_KUTTA_SECOND_ORDER_MODIFIED_EULER: {
      // calculate using Modified Euler scheme
      std::tie(time, numIter) = plasticRKME221(state, purelyPlasticStrain);
      break;
    }
    case SolutionAlgorithm::RUNGE_KUTTA_THIRD_ORDER_NYSTROM: {
      // calculate using 3rd order RK scheme (Nystrom)
      std::tie(time, numIter) = plasticRK332(state, purelyPlasticStrain);
      break;
    }
    case SolutionAlgorithm::RUNGE_KUTTA_THIRD_ORDER_BOGACKI: {
      // calculate using 3rd order RK scheme (Bogacki - Shampine)
      std::tie(time, numIter) = plasticRKBog432(state, purelyPlasticStrain);
      break;
    }
    case SolutionAlgorithm::RUNGE_KUTTA_FOURTH_ORDER: {
      // calculate using 4th order RK scheme
      std::tie(time, numIter) = plasticRK543(state, purelyPlasticStrain);
      break;
    }
    case SolutionAlgorithm::RUNGE_KUTTA_FIFTH_ORDER_ENGLAND: {
      // calculate using 5th order RK scheme (England)
      std::tie(time, numIter) = plasticRKEng654(state, purelyPlasticStrain);
      break;
    }
    case SolutionAlgorithm::RUNGE_KUTTA_FIFTH_ORDER_CASH: {
      // calculate using 5th order RK scheme (Cash - Karp)
      std::tie(time, numIter) = plasticRKCK654(state, purelyPlasticStrain);
      break;
    }
    case SolutionAlgorithm::RUNGE_KUTTA_FIFTH_ORDER_DORMAND: {
      // calculate using 5th order RK scheme (Dormand - Prince)
      std::tie(time, numIter) = plasticRKDP754(state, purelyPlasticStrain);
      break;
    }
    case SolutionAlgorithm::RUNGE_KUTTA_FIFTH_ORDER_BOGACKI: {
      // calculate using 5th order RK scheme (Bogacki - Shampine)
      std::tie(time, numIter) = plasticRKErr8544(state, purelyPlasticStrain);
      break;
    }
    case SolutionAlgorithm::EXTRAPOLATION_BULIRSCH: {
      // calculate using extrapolation method (Bulirsch - Stoer)
      std::tie(time, numIter) = plasticExtrapol(state, purelyPlasticStrain);
      break;
    }
    default: {
      std::ostringstream err;
      err << "**ERROR** Unknown Solution Algorithm. Value of d_solutionAlgorithm variable "
             "is set to:" << d_solutionAlgorithm << "\n";
      err << "Acceptable values are ints from 1 to 9. Procedure calculate "
             "Plastic exits without any calculations done.\n";
      throw InvalidValue(err.str(), __FILE__, __LINE__);
    }
  }

  return time;
}

/**
  This procedure calculate stress increment using Modified Euler method
  A - matrix with coefficients, 
  C - x coefficients. 
  BRes - result coefficients, 
  B - error estimate, 
 */
std::tuple<double, int>
ShengMohrCoulomb::plasticRKME221(StateMohrCoulomb& state, 
                                 const Vector7& epStrain)

{
  constexpr int Order = 2;
  constexpr int Steps = 2;

  using RKMatrix = Eigen::Matrix<double, Steps, Steps>;
  using RKVector = Eigen::Matrix<double, Steps, 1>;

  RKMatrix A = RKMatrix::Zero();
  A(0, 0) = 0.0;
  A(1, 0) = 1.0;

  RKVector C = RKVector::Zero();
  C(0) = 0.0;
  C(1) = 1.0;

  RKVector BRes = RKVector::Zero();
  BRes(0) = 0.5;
  BRes(1) = 0.5;

  RKVector B = RKVector::Zero();
  B(0) = 1.0;
  B(1) = 0.0;

  bool errorEstimate = false;

  double time = 0;
  int numIter = 0;

  auto result = doRungeKutta<Order, Steps>(A, B, BRes, C, 
                                           state, epStrain,
                                           errorEstimate);

  return result;
}

/**
  This procedure calculate any Runge-Kutta pair, given the coefficient of
  stress in A for each x used to calculate values in C, where B gives the
  coefficient to calculate error estimate (or lower order solution) and 
  BRes to calculate result. The procedure will integrate whole step of strain
  re-using the initial derivative for rejected steps.  

  Order = order of the method (required for substep prediction)
  Steps = the number of stages to get the result
  errorEstimate = if true - the B table contains error estimate. 
                  if false, it contains 4th order solution
*/
template<int Order, int Steps>
std::tuple<double, int>
ShengMohrCoulomb::doRungeKutta(const Eigen::Matrix<double, Steps, Steps>& AA, 
                               const Eigen::Matrix<double, Steps, 1>&     BB, 
                               const Eigen::Matrix<double, Steps, 1>&     BRes, 
                               const Eigen::Matrix<double, Steps, 1>&     CC,
                               StateMohrCoulomb&       state, 
                               const Vector7&          epStrain,
                               bool                    errorEstimate)
{
  // *WARNING** prevents algorithm to have more than 1e4 steps; cost: accuracy
  // reduction in some unusual cases, but the whole thing will keep on going
  constexpr double criticalStepSize = 1.0e-4; 
  double newStepSize = 0;
  for (int i = 0; i < 6; i++) {
    newStepSize += std::abs(substepStrain[i]);
  }
  newStepSize = 0.01 / (newStepSize * Order);
  if (newStepSize > 1) {
    newStepSize = 1;
  }
  if (dbg.active()) {
    dbg << "In doRungeKutta::newStepSize = " << newStepSize << "\n";
  }
  double stepLength        = 1;
  double totalSize         = 0;
  bool   reuseStep         = false;
  double microStep         = 0;
  double stepAccuracyCheck = 0;
  double frequency         = 100000.0 / Order;    // how often display info about steps
  double methodPower       = std::pow(2.0, Order) * d_int.d_integrationTol;

  std::vector<StateMohrCoulomb> midStates(Steps);

  Eigen::Matrix<double, 6, Steps> dSigma        = Eigen::Matrix<double, 6, Steps>::Zero();
  Eigen::Matrix<double, 6, Steps> plasticStrain = Eigen::Matrix<double, 6, Steps>::Zero();
  Eigen::Matrix<double, 1, Steps> dP0Star       = Eigen::Matrix<double, 1, Steps>::Zero();

  Vector7 reuseRes         = Vector7::Zero();
  Vector6 stressInc        = Vector6::Zero();
  Vector6 plasticStrainInc = Vector6::Zero();
  double p0StarInc         = 0.0;

  std::chrono::time_point<std::chrono::system_clock> startTime, endTime;
  startTime = std::chrono::system_clock::now();

  int    stepNo = 0;
  bool finished = false;
  do {
    bool stepAccepted = false;
    ++stepNo;

    // compute size of a step
    stepLength *= newStepSize; 
    if (stepLength < criticalStepSize) {
      stepLength = criticalStepSize; // HACK here, it is a fairly large value
    }
    // check whether the step exceeds the whole increment
    if ((stepLength + totalSize) > 1) {
      stepLength = 1 - totalSize; 
    }

    // strain increment in current step
    Vector7 substepStrain = epStrain * stepLength;
    Vector7 currentStrain = Vector7::Zero();

    double rError = 0;

    if (dbg.active()) {
      dbg << "Step Length = "<< stepLength 
          << " Current strain (0) = " << currentStrain(0) << "\n";
    }

    // Make copies of the initial state
    for (auto& midState : midStates) {
      midState = state;
    }

    // Below the main R-K loop to calculate the value of intermediate stresses;
    // values stored in dSigma_temp
    if (reuseStep) {

      dSigma.col(0)        = reuseRes.block<6,1>(0,0);
      plasticStrain.col(0) = plasticStrainInc.block<6,1>(0,0);
      dP0Star.col(0)       = reuseRes(6);

    } else {

      stressInc        = Vector6::Zero();
      plasticStrainInc = Vector6::Zero();
      calcPlastic(midStates[0], substepStrain, 
                  stressInc, plasticStrainInc, p0StarInc);

      // **TODO** Check if plastic strain increment should not be computed
      // plasticStrainInc = Vector6::Zero();
      dSigma.col(0)        = stressInc;
      plasticStrain.col(0) = plasticStrainInc;
      dP0Star.col(0)       = p0StarInc;

    }

    for (int rkloop = 1; rkloop < Steps; rkloop++) {

      stressInc        = Vector6::Zero();
      plasticStrainInc = Vector6::Zero();
      p0StarInc        = 0;

      currentStrain = substepStrain * CC[rkloop];

      for (int i = 0; i < rkloop; i++) {

        stressInc        += (dSigma.col(i) * AA(rkloop, i));
        plasticStrainInc += (plasticStrain.col(i) * AA(rkloop, i));
        p0StarInc        += (dP0Star.col(i) * AA(rkloop, i);

      }

      midStates[rkloop].update(plasticStrainInc, currentStrain, stressInc,
                               p0StarInc);

      calcPlastic(midStates[rkloop], substepStrain, 
                  stressInc, plasticStrainInc, p0StarInc);

      // **TODO** Check if plastic strain increment should not be computed
      // plasticStrainInc = Vector6::Zero();
      dSigma.col(rkloop)        = stressInc;
      plasticStrain.col(rkloop) = plasticStrainInc;
      dP0Star.col(rkloop)       = p0StarInc;
    }

    Vector6 sigBRes  = dSigma        * BRes;
    Vector6 epsPBres = plasticStrain * BRes;
    double  p0Bres   = dP0Star       * BRes;
    Vector6 sigB     = dSigma        * BB;
    double  p0B      = dP0Star       * BB;
    
    Vector7 result, error;
    result.block<6,1>(0,0) = sigBRes;
    error.block<6,1>(0,0)  = sigB;
    result(6)              = p0BRes;
    error(6)               = p0B;

    // error estimate calculated in case we
    // have lower order solution instead of
    // error estimate
    if (!errorEstimate) {
      error -= result;
    }

    // check the error norm
    switch (d_int.d_tolMethod) {
      case EPUS_RELATIVE_ERROR:
        rError = checkNorm(result, result(6), point, error); 
        break;
      case SLOAN:
        rError = checkNormSloan(result, result(6), point, error);
        break;
      default: 
        std::ostringstream err;
        err << "**ERROR** Improper d_tolMethod in increment.dta" << "\n";
        throw InvalidValue(err.str(), __FILE__, __LINE__);
    }

    for (int i = 0; i < 7; i++) {
      if (!std::isfinite(result(i))) {
        std::ostringstream err;
        err << "**ERROR** Results not a number. Try correcting issue, results may be "
               "incorrect."
            << " Results = \n" << result << "\n";
        throw InvalidValue(err.str(), __FILE__, __LINE__);
      }
    }

    if (rError < d_int.d_integrationTol) {
      stepAccepted = true;
    } else {
      stepAccepted = false;
      if (stepLength <= criticalStepSize) {
        stepAccepted = true;
      }
    }

    if (rError < TINY) {
      rError = TINY;
    }
    
    if (d_int.d_tolMethod == EPUS_RELATIVE_ERROR) {
      newStepSize = d_int.d_betaFactor * 
        std::pow(d_int.d_integrationTol / rError, (1 / (Order - 1.0)));
    } else {
      newStepSize = d_int.d_betaFactor * 
        std::pow(d_int.d_integrationTol / rError, (1 / Order));
    }

    if (!stepAccepted) {

      // here we should take care about correct re - usage of the first
      // evaluation of derivative
      reuseStep = true;
      reuseRes.block<6,1>(0,0)         = dSigma.col(0)        * newStepSize;
      reuseRes(6)                      = dP0Star.col(0)       * newStepSize;
      plasticStrainInc.block<6,1>(0,0) = plasticStrain.col(0) * newStepSize;

    } else {

      StateMohrCoulomb oldState = state; 

      //******
      // here we update all the state data.
      state.update(plasticStrainInc, substepStrain, result, result(6));
      //******

      // Keeping a copy to check values
      StateMohrCoulomb trialState = state;

      switch (d_int.d_driftCorrection) {
        case NO_CORRECTION:
          break;
        case CORRECTION_AT_BEGIN:
          correctDriftBeg(state, &oldState);
          break;
        case CORRECTION_AT_END:
          correctDrift(state);
          break;
        default:
          break;
      }

      // reevaluate the error in the point:
      error.block<6,1>(0,0) += (state.stress - trialState.stress);
      error(6) += state->p0Star() - trialState.p0Star();

      // error vector updated, norm should be re-evaluated:
      switch (d_int.d_tolMethod) {
        case EPUS_RELATIVE_ERROR:
          rError = checkNorm(result, result(6), point, error); 
          break;
        case SLOAN:
          rError = checkNormSloan(result, result(6), point, error); 
          break;
        default: 
          std::ostringstream err;
          err << "**ERROR** Improper d_tolMethod in increment.dta" << "\n";
          throw InvalidValue(err.str(), __FILE__, __LINE__);
      }

      if (!std::isfinite(rError)) {
        if (rError < methodPower) {
          rError = methodPower;
        }
      }
    
      if (d_int.d_tolMethod == EPUS_RELATIVE_ERROR) {
        newStepSize = d_int.d_betaFactor * 
          std::pow(d_int.d_integrationTol / rError, (1 / (Order - 1.0)));
      } else {
        newStepSize = d_int.d_betaFactor * 
          std::pow(d_int.d_integrationTol / rError, (1 / Order));
      }

      if (rError < d_int.d_integrationTol) {
        stepAccepted = true;
      } else {
        stepAccepted = false;
        if (stepLength <= criticalStepSize) {
          stepAccepted = true;
        }
      }

      if (!stepAccepted) {

        reuseStep = true;
        reuseRes.block<6,1>(0,0)         = dSigma.col(0)        * newStepSize;
        reuseRes(6)                      = dP0Star.col(0)       * newStepSize;
        plasticStrainInc.block<6,1>(0,0) = plasticStrain.col(0) * newStepSize;

        oldState = state;
        double temp = double(stepNo) / frequency;
        if (std::modf(temp, &temp) == 0) {
          std::cout << "Step number : " << stepNo << "\n";
          std::cout << "Total size done is : " << totalSize
                    << " of whole step. Current stepLength = " << stepLength << "\n";
          std::cout << "Stress state :\n" << state.stress << "\n";
        }

      } else {

        // this may be done only after successful re-evaluation of the substep
        reuseStep = false;
        stepAccuracyCheck = totalSize;
        microStep += stepLength;
        totalSize += totalSize + microStep; // total part of step done updated
        microStep -= (totalSize - stepAccuracyCheck);
        if (totalSize >= 1) {
          finished = true;
        }
        double temp = double(stepNo) / frequency;
        if (std::modf(temp, &temp) == 0) {
          std::cout << "Step number : " << stepNo << "\n";
          std::cout << "Total size done is : " << totalSize
                    << " of whole step. Current stepLength = " << stepLength << "\n";
          std::cout << "Stress state :\n" << state.stress << "\n";
        }
      }
    }

  } while (!finished);

  endTime = std::chrono::system_clock::now();

  double time = std::chrono::duration<double>(endTime - startTime).count();
  return std::make_tuple(time, stepNo);
}

/**
 * Compute the stress and p0star increment, 
 */
int
ShengMohrCoulomb::calcPlastic(const StateMohrCoulomb& state, 
                              const Vector7&          epStrainInc,
                              Vector6&                dSigma, 
                              Vector6&                dEps_p, 
                              double&                 dP0Star)
{
  if (!state.checkIfFinite()) {
    std::ostringstream err;
    err << "**Error** in the calcPlastic Procedure. "
        << "Point internal values are not finite.\n";.
    throw InvalidValue(err.str(), __FILE__, __LINE__);
  }

  auto stress = state.stress;

  Vector6 df_dsigma = Vector6::Zero();
  Vector6 dg_dsigma = Vector6::Zero();
  std::tie(df_dsigma, dg_dsigma) = computeDfDsigma(stress);

  double K = d_elastic.K;
  double G = d_elastic.G;
  Matrix66 elasticTangent = calculateElasticTangentMatrix(K, G);

  Vector6 dEps = epStrainInc.block<6,1>(0,0);
  auto df_dsigma_x_D = df_dsigma.transpose() * elasticTangent;
  double numerator   = df_dsigma_x_D * dEps;
  double denominator = df_dsigma_x_D * dg_dsigma;
  if (denominator < TINY) {
    std::cout << "**WARNING** Denominator of plastic multiplier is very small."
                 " Some error may arise and results may be incorrect" << "\n";
  }

  double ratio = numerator/denominator;
  dEps_p = dg_dsigma * ratio;
  auto dSigma = elasticTangent * (dEps - dEps_p);

  for (int i = 0; i < 6; i++) {
    if (!std::isfinite(dSigma(i))) {
      std::ostringstream err;
      err << "calcPlastic: dSigma not finite. dSigma = \n"
          << dSigma << "\n";
      throw InvalidValue(err.str(), __FILE__, __LINE__);
           << "\n";
    }
  }

  dP0Star = 1.0;
  return 0;
}

/**
  Procedure returns value of the relative error rError
  The value is calculated as maximum of all the values of rErrors
  rErrors used:
          - each stress component
          - error in p0*
          - error in p
          - error in q

  INPUT: dSigma [6] - values of stress increment (min 6, but may be more; values
  over six ignored)
                  dP0Star - value of the p0* increment
                  initialState - point used to calculate the values of the step
                  dError [7] - error estimates for all values of dSigma
  dError[i]  - error estimate for dSigma[i];
  returns rError
 */
double 
ShengMohrCoulomb::checkNorm(const Vector7& dSigma, double dP0Star,
                            const StateMohrCoulomb& initialState,
                            const Vector7& dError) 
{
  using Vector9 = Eigen::Matrix<double, 9, 1>;

  Vector9 drError = Vector9::Zero();
  double rError   = 0;

  // standard norm:
  for (int i = 0; i < 6; i++) {
    if (std::abs(dSigma(i)) > TINY) {
      drError(i) = dError(i) / dSigma(i);
    } else {
      drError(i) = dError(i) / TINY;
    }
  }
  if (std::abs(dP0Star) > TINY) {
    drError(6) = dError(6) / dP0Star;
  } else {
    drError(6) = dError(6) / TINY;
  }

  // norm in mean stress p:
  double p = std::abs(dSigma(0) + dSigma(1) + dSigma(2)) / 3;
  double errorP = std::abs(dError(0) + dError(1) + dError(2)) / 3;
  if (p > TINY) {
    drError(7) = errorP / P;
  } else {
    drError(7) = errorP / TINY;
  }

  // norm of q...
  Vector6 initialSigma  = initialState.stress;
  Vector6 finalSigma    = initialSigma + dSigma.block<6,1>(0,0);
  double  initialShear  = initialState.shearStress();
  double  finalShear    = 
    (finalSigma(0) - finalSigma(1)) * (finalSigma(0) - finalSigma(1)) +
    (finalSigma(0) - finalSigma(2)) * (finalSigma(0) - finalSigma(2)) +
    (finalSigma(1) - finalSigma(2)) * (finalSigma(1) - finalSigma(2)) +
    6 * (finalSigma(3) * finalSigma(3) + finalSigma(4) * finalSigma(4) +
         finalSigma(5) * finalSigma(5));
  finalShear = std::sqrt(0.5 * finalShear);

  finalSigma += dError.block<6,1>(0,0);
  double finalShearMax =
    (finalSigma(0) - finalSigma(1)) * (finalSigma(0) - finalSigma(1)) +
    (finalSigma(0) - finalSigma(2)) * (finalSigma(0) - finalSigma(2)) +
    (finalSigma(1) - finalSigma(2)) * (finalSigma(1) - finalSigma(2)) +
    6 * (finalSigma(3) * finalSigma(3) + finalSigma(4) * finalSigma(4) +
         finalSigma(5) * finalSigma(5));
  finalShearMax = std::sqrt(0.5 * finalShearMax);

  finalSigma -= dError.block<6,1>(0,0) * 2.0;
  double finalShearMin =
    (finalSigma(0) - finalSigma(1)) * (finalSigma(0) - finalSigma(1)) +
    (finalSigma(0) - finalSigma(2)) * (finalSigma(0) - finalSigma(2)) +
    (finalSigma(1) - finalSigma(2)) * (finalSigma(1) - finalSigma(2)) +
    6 * (finalSigma(3) * finalSigma(3) + finalSigma(4) * finalSigma(4) +
         finalSigma(5) * finalSigma(5));
  finalShearMin = std::sqrt(0.5 * finalShearMin);

  double shearError = std::abs(finalShearMax - finalShear);
  if (std::abs(finalShearMin - finalShear) > shearError) {
    shearError = std::abs(finalShearMin - finalShear);
  }

  double dShear = std::abs(finalShear - initialShear);
  if (dShear > TINY) {
    drError(8) = shearError / DShear;
  } else {
    drError(8) = shearError / TINY;
  }

  // final check: max norm
  for (int i = 0; i < 9; i++) {
    drError(i) = std::abs(drError(i));
    if (drError(i) > rError) {
      rError = drError(i);
    }
  }
  return rError;
}

double 
ShengMohrCoulomb::checkNormSloan(const Vector7& dSigma, double dP0Star,
                                 const StateMohrCoulomb& initialState,
                                 const Vector7& dError) 
{
  using Vector9 = Eigen::Matrix<double, 9, 1>;

  Vector9 drError = Vector9::Zero();
  double rError   = 0;

  double InitialSigma[6], finalSigma[6], dP0StarEnd;

  Vector6 initialSigma  = initialState.stress;
  Vector6 finalSigma    = initialSigma + dSigma.block<6,1>(0,0);
  double  dP0StarEnd    = initialState.p0Star() + dP0Star;

  // standard norm:
  for (int i = 0; i < 6; i++) {
    if (std::abs(finalSigma(i)) > TINY) {
      drError(i) = dError(i) / finalSigma(i);
    } else {
      drError(i) = dError(i) / TINY;
    }
  }
  if (std::abs(dP0StarEnd) > TINY) {
    drError(6) = dError(6) / dP0StarEnd;
  } else {
    drError(6) = dError(6) / TINY;
  }

  // norm in mean stress p:
  double p = std::abs(finalSigma(0) + finalSigma(1) + finalSigma(2)) / 3;
  double errorP = std::abs(dError(0) + dError(1) + dError(2)) / 3;
  if (p > TINY) {
    drError(7) = errorP / P;
  } else {
    drError(7) = errorP / TINY;
  }

  // norm of q...
  double  initialShear  = initialState.shearStress();
  double  finalShear    = 
    (finalSigma(0) - finalSigma(1)) * (finalSigma(0) - finalSigma(1)) +
    (finalSigma(0) - finalSigma(2)) * (finalSigma(0) - finalSigma(2)) +
    (finalSigma(1) - finalSigma(2)) * (finalSigma(1) - finalSigma(2)) +
    6 * (finalSigma(3) * finalSigma(3) + finalSigma(4) * finalSigma(4) +
         finalSigma(5) * finalSigma(5));
  finalShear = std::sqrt(0.5 * finalShear);

  finalSigma += dError.block<6,1>(0,0);
  double finalShearMax =
    (finalSigma(0) - finalSigma(1)) * (finalSigma(0) - finalSigma(1)) +
    (finalSigma(0) - finalSigma(2)) * (finalSigma(0) - finalSigma(2)) +
    (finalSigma(1) - finalSigma(2)) * (finalSigma(1) - finalSigma(2)) +
    6 * (finalSigma(3) * finalSigma(3) + finalSigma(4) * finalSigma(4) +
         finalSigma(5) * finalSigma(5));
  finalShearMax = std::sqrt(0.5 * finalShearMax);

  finalSigma -= (dError.block<6,1>(0,0) * 2.0);
  double finalShearMin =
    (finalSigma(0) - finalSigma(1)) * (finalSigma(0) - finalSigma(1)) +
    (finalSigma(0) - finalSigma(2)) * (finalSigma(0) - finalSigma(2)) +
    (finalSigma(1) - finalSigma(2)) * (finalSigma(1) - finalSigma(2)) +
    6 * (finalSigma(3) * finalSigma(3) + finalSigma(4) * finalSigma(4) +
         finalSigma(5) * finalSigma(5));
  finalShearMin = std::sqrt(0.5 * finalShearMin);

  double shearError = std::abs(finalShearMax - finalShear);
  if (std::abs(finalShearMin - finalShear) > shearError) {
    shearError = std::abs(finalShearMin - finalShear);
  }
  if (finalShear > TINY) {
    drError(8) = shearError / finalShear;
  } else {
    drError(8) = shearError / TINY;
  }

  // final check: max norm
  for (int i = 0; i < 9; i++) {
    drError(i) = std::abs(drError(i));
    if (drError(i) > rError) {
      rError = drError(i);
    }
  }
  return rError;
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

  // for (int i=0; i<6; i++) dbg << "Delta stress i="<<ds[i]<<"\n";
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
      dbg << "check Yield: Mean Stress="<<meanStress;
          dbg << " Shear Stress="<<shearStress;
          dbg << " cohesion="<<cohesion;
          dbg << " M="<<M;
          dbg << " Yield Function="<<Value<<"\n";
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
    // dbg << "Elasto-Plastic Matrix used"<<"\n";
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
  BBMMatrix denominator (1,1);

  BBMMatrix dg_dsigma (6,1); //dG/dsigma, in case of associated flow
  (nonAssociated==false) same as A
  BBMMatrix MM (6,1); // will be vector(1,1,1,0,0,0)T
  BBMMatrix elasticTangent (6,6); //D elastic matrix...

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
                       dbg << "Factor in calculate Elasto-Plastic Tangent Matrix
  is not finite. set to 1"<<"\n";
                       dbg << "alpha4="<<alpha4<<" J3="<<J3<<"
  shearStress="<<shearStress<<"\n";
                  }
  Factor=1+alpha4-(1-alpha4)*Factor;
  double Factor025=pow(Factor,0.25);
  double Factor075=pow(Factor,-0.75);
  double M=(3-sin_phi)/(6*alpha*sin_phi);
  double Mpsi=(3-sin_psi)/(6*alpha*sin_psi);



  BBMMatrix dJ2_dSig (6,1), dJ3_dSig (6,1), dq_dSig (6,1), dI1_dSig (6,1), dI2_dSig
  (6,1), dI3_dSig (6,1), TEMP(6,1), numerator(6,1);

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

  dJ2_dSig.Copy(&dg_dsigma);
  dg_dsigma.Multiply (0.21022410391343*Mpsi*Factor025/shearStress,&dg_dsigma);
  dJ3_dSig.Copy(&TEMP);
  TEMP.Multiply
  (1.419012700740643*(alpha4-1)*Mpsi*sqrt(J2)*Factor075/(sqrt(3.0)*shearStress*shearStress*shearStress)-
                             0.94600846716043*(alpha4-1)*cohesion*Mpsi*Factor075/(shearStress*shearStress*shearStress),&TEMP);
  dg_dsigma.Add (&TEMP,&dg_dsigma);
  dq_dSig.Copy(&TEMP);
  TEMP.Multiply
  (2.838025401481287*(alpha4-1)*cohesion*Mpsi*J3*Factor075/(shearStress*shearStress*shearStress*shearStress)-
                             4.257038102221929*(alpha4-1)*Mpsi*sqrt(J2)*J3*Factor075/(sqrt(3.0)*shearStress*shearStress*shearStress*shearStress),&TEMP);
  dg_dsigma.Add (&TEMP,&dg_dsigma);
  dI1_dSig.Copy(&TEMP);
  TEMP.Multiply(-1/3.0, &TEMP);  //First correction: sign here, but most likely
  not enough!!!
  dg_dsigma.Add (&TEMP,&dg_dsigma);

  //FINISHED dQ/dSigma

  double K43G=K+4.0*G/3.0;
  double K23G=K-2.0*G/3.0;

  elasticTangent.PutElement (1,1,K43G);
  elasticTangent.PutElement (1,2,K23G);
  elasticTangent.PutElement (1,3,K23G); //rest of the line are zeros and rightly so
  elasticTangent.PutElement (2,1,K23G);
  elasticTangent.PutElement (2,2,K43G);
  elasticTangent.PutElement (2,3,K23G); //yes, the matrix is symmetrical, but it is faster
  to put this 3 additional elements
  elasticTangent.PutElement (3,1,K23G); //than just mirror all, including zeros, which are
  there already
  elasticTangent.PutElement (3,2,K23G);
  elasticTangent.PutElement (3,3,K43G);
  elasticTangent.PutElement (4,4,2.0*G);
  elasticTangent.PutElement (5,5,2.0*G);
  elasticTangent.PutElement (6,6,2.0*G); //rest of the matrix is filled with zeros...

  //getting lambda and Dep

  A.Transpose(&numerator);
  numerator.Multiply(&elasticTangent,&numerator); //numerator=aT*Del -->Numerator of Lambda
  without multiplication by dEpsilon
  numerator.Multiply(&dg_dsigma,&denominator);

  elasticTangent.Multiply(&dg_dsigma,&TEMP);
  TEMP.Multiply(&numerator,&TEMP);
  TEMP.Multiply(1/denominator.getElement(1,1),&TEMP);
  elasticTangent.Substract(&TEMP,&elasticTangent);
  for (int i=1;i<7; i++) for (int j=1; j<7; j++)
  DEP->PutElement(i,j,elasticTangent.getElement(i,j));
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
          //dbg << "1:"<<deriv[2]<<"\n";
          deriv[2]=deriv[2]*-1*numerator/(denominator*denominator);
          //dbg << "2:"<<deriv[2]<<"\n";
          deriv[2]=deriv[2]*pc*log(P0Star/pc);
          //dbg << "3:"<<deriv[2]<<"\n";
          deriv[2]=deriv[2]*pow((P0Star/pc),numerator/denominator); //it's
     finished with dp0/ds
          //dbg << "4:"<<deriv[2]<<"\n";
          deriv[2]=M*M*(meanStress+SuctionPressure)*deriv[2];
          //dbg << "5:"<<deriv[2]<<"\n";
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
                  dbg << "yield of the suction surface occurred"<<"\n";
                  if (dsuction>0)
                  {
                          cout <<"Unable to find gradient; Whole step
  plastic"<<"\n";
                          return -2;
                  //this check may be unused...
                  }
          }
  //check of the standard yield surface

          //dbg << "meanStress="<<meanStress<<"\n";
          //dbg << "shearStress="<<shearStress<<"\n";

          PZero=LambdaZero*((1-r)*exp(-1*Beta*Suction)+r);
      PZero=(LambdaZero-KappaP)/(PZero-KappaP);
          PZero=pc*pow((P0Star/pc),PZero);
          SuctionPressure=k*Suction;
  //	dbg << "suction="<<suction<<"\n";
  //	dbg << "PZero="<<PZero<<"\n";
          fValue=shearStress*shearStress-M*M*(meanStress+SuctionPressure)*(PZero-meanStress);
          fValue=fValue/((PZero+SuctionPressure)*(PZero+SuctionPressure));
  //value of Yield function calculated and normalized
          if (fValue<-d_yieldTol)
          {
                  dbg << "!!!there is no yield at the beginning !!!"<<"\n";
                  dbg << "fValue is="<<fValue<<"\n";
                  dbg << "!!!find gradient procedure
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
          //dbg << "1:"<<dF[2]<<"\n";
          dF[2]=dF[2]*-1*numerator/(denominator*denominator);
          //dbg << "2:"<<dF[2]<<"\n";
          dF[2]=dF[2]*pc*log(P0Star/pc);
          //dbg << "3:"<<dF[2]<<"\n";
          dF[2]=dF[2]*pow((P0Star/pc),numerator/denominator); //it's finished
  with dp0/ds
          //dbg << "4:"<<dF[2]<<"\n";
          dF[2]=M*M*(meanStress+SuctionPressure)*dF[2];
          //dbg << "5:"<<dF[2]<<"\n";
          dF[2]=-1*M*M*k*(PZero-meanStress)-dF[2];	//final result
          //dbg << dF[2]<<"\n";
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
          dbg << "Mean Stress="<<meanStress<<"  Shear Stress="<<shearStress<< "
  Suction="<<Suction<<"\n";
          dbg << "dMean Stress="<<dmeanStress<<"  dShear
  Stress="<<dshearStress<<" dSuction="<<dsuction<<"\n";
          for (int i=0; i<3; i++) dbg << "value of F["<<i<<"] is:"<<dF[i]<<"\n";
          dbg << "df length is:"<<dFlength<<"\n";
          dbg << "cosinus is:"<<cosin<<"\n";
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

  // dbg << "Correct Drift Procedure entered!"<<"\n";

  BBMMatrix A(6, 1); // dF/dsigma
  BBMMatrix denominator(1, 1);

  BBMMatrix dg_dsigma(6, 1);  // dG/dsigma, in case of associated flow
                       // (nonAssociated==false) same as A
  BBMMatrix MM(6, 1);  // will be vector(1,1,1,0,0,0)T
  BBMMatrix elasticTangent(6, 6); // D elastic matrix...
  BBMMatrix dEps(6, 1);
  BBMMatrix dSigma(6, 1);

  double dSigma[6], epStrain[6], zeros[7];
  for (int i = 0; i < 7; i++)
    zeros[i] = 0;

  for (int i = 1; i < 6; i++)
    dEps.PutElement(i, 1, epStrain[i - 1]); // increase of epsilon copied
  // we do not have g, as a == g in associated flow rule.

  bool correct;

  double I1 = Point->firstInvariant();
  double I2 = Point->secondInvariant();
  double J2 = Point->secondDevInvariant();
  if (std::abs(J2) < TINY)
    J2 = TINY;
  double J3 = Point->thirdDevInvariant();
  double shearStress = Point->shearStress();
  if (std::abs(shearStress) < TINY)
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
      // dbg << "Drift Correction, Iteration="<<numberIter<<" Function
      // Value="<<fValue<<"\n";
      // CORRECT FOR DRIFT
      // HERE THE DRIFT WILL BE CORRECTED BY USING THE D MATRIX FROM THE
      // FORBIDDEN SPACE
      // ALTHOUGH BECAUSE THE DRIFT WILL BE CHECKED AGAIN, IT SHOULDN'T POSE
      // MUCH PROBLEM.

      BBMMatrix dJ2_dSig(6, 1), dJ3_dSig(6, 1), dq_dSig(6, 1), dI1_dSig(6, 1),
        dI2_dSig(6, 1), dI3_dSig(6, 1), TEMP(6, 1), numerator(6, 1);

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

      // dbg << "dF/dJ2 part"<<"\n";
      // A.PrintPrecise();
      dJ3_dSig.Copy(&TEMP);
      TEMP.Multiply(-1.419012700740643 * (alpha4 - 1) * M * sqrt(J2) *
                        Factor075 /
                        (sqrt(3.0) * shearStress * shearStress * shearStress) +
                      0.94600846716043 * (alpha4 - 1) * cohesion * M *
                        Factor075 / (shearStress * shearStress * shearStress),
                    &TEMP);
      // dbg << "dF/dJ3 part"<<"\n";
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
      // dbg << "dF/dq part"<<"\n";
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
      dbg << "dSigma="<<dSs<<"\n";
      dbg << "dF/dSigma11="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[0]-=dSs;

      CopyPoint.stress[1]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      dbg << "dF/dSigma22="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[1]-=dSs;

      CopyPoint.stress[2]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      dbg << "dF/dSigma33="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[2]-=dSs;

      dSs=10*dSs;
      CopyPoint.stress[0]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      dbg << "dSigma="<<dSs<<"\n";
      dbg << "dF/dSigma11="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[0]-=dSs;

      CopyPoint.stress[1]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      dbg << "dF/dSigma22="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[1]-=dSs;

      CopyPoint.stress[2]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      dbg << "dF/dSigma33="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[2]-=dSs;

      dSs=10*dSs;
      CopyPoint.stress[0]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      dbg << "dSigma="<<dSs<<"\n";
      dbg << "dF/dSigma11="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[0]-=dSs;

      CopyPoint.stress[1]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      dbg << "dF/dSigma22="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[1]-=dSs;

      CopyPoint.stress[2]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      dbg << "dF/dSigma33="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[2]-=dSs;

      dSs=10*dSs;
      CopyPoint.stress[0]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      dbg << "dSigma="<<dSs<<"\n";
      dbg << "dF/dSigma11="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[0]-=dSs;

      CopyPoint.stress[1]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      dbg << "dF/dSigma22="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[1]-=dSs;

      CopyPoint.stress[2]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      dbg << "dF/dSigma33="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[2]-=dSs;

      dSs=10*dSs;
      CopyPoint.stress[0]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      dbg << "dSigma="<<dSs<<"\n";
      dbg << "dF/dSigma11="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[0]-=dSs;

      CopyPoint.stress[1]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      dbg << "dF/dSigma22="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[1]-=dSs;

      CopyPoint.stress[2]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      dbg << "dF/dSigma33="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[2]-=dSs;

      dSs=10*dSs;
      CopyPoint.stress[0]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      dbg << "dSigma="<<dSs<<"\n";
      dbg << "dF/dSigma11="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[0]-=dSs;

      CopyPoint.stress[1]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      dbg << "dF/dSigma22="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[1]-=dSs;

      CopyPoint.stress[2]+=dSs;
      Yield2=computeYieldFunctionNN(&CopyPoint);
      dbg << "dF/dSigma33="<<(Yield2-Yield1)/dSs<<"\n";
      CopyPoint.stress[2]-=dSs;
      */
      dJ2_dSig.Copy(&dg_dsigma);
      dg_dsigma.Multiply(0.21022410391343 * Mpsi * Factor025 / shearStress, &dg_dsigma);
      dJ3_dSig.Copy(&TEMP);
      TEMP.Multiply(-1.419012700740643 * (alpha4 - 1) * Mpsi * sqrt(J2) *
                        Factor075 /
                        (sqrt(3.0) * shearStress * shearStress * shearStress) +
                      0.94600846716043 * (alpha4 - 1) * cohesion * Mpsi *
                        Factor075 / (shearStress * shearStress * shearStress),
                    &TEMP);
      dg_dsigma.Add(&TEMP, &dg_dsigma);
      dq_dSig.Copy(&TEMP);
      TEMP.Multiply(
        -2.838025401481287 * (alpha4 - 1) * cohesion * Mpsi * J3 * Factor075 /
            (shearStress * shearStress * shearStress * shearStress) +
          4.257038102221929 * (alpha4 - 1) * Mpsi * sqrt(J2) * J3 * Factor075 /
            (sqrt(3.0) * shearStress * shearStress * shearStress * shearStress),
        &TEMP);
      dg_dsigma.Add(&TEMP, &dg_dsigma);
      dI1_dSig.Copy(&TEMP);
      TEMP.Multiply(
        -1 / 3.0,
        &TEMP); // First correction: sign here, but most likely not enough!!!
      dg_dsigma.Add(&TEMP, &dg_dsigma);

      // FINISHED dQ/dSigma

      if (meanStress < -cohesion * M) {
        // tension cut-off plane
        dg_dsigma(0, 0) =  1.0 / 3.0);
        dg_dsigma(1, 0) =  1.0 / 3.0);
        dg_dsigma(2, 0) =  1.0 / 3.0);
        dg_dsigma(3, 0) =  0.0);
        dg_dsigma(4, 0) =  0.0);
        dg_dsigma(5, 0) =  0.0);
      }

      double K43G = K + 4.0 * G / 3.0;
      double K23G = K - 2.0 * G / 3.0;

      elasticTangent(0, 0) =  K43G);
      elasticTangent(0, 1) =  K23G);
      elasticTangent(0, 2) =  K23G); // rest of the line are zeros and rightly so
      elasticTangent(1, 0) =  K23G);
      elasticTangent(1, 1) =  K43G);
      elasticTangent(1, 2) =  K23G); // yes, the matrix is symmetrical, but it is
                                  // faster to put this 3 additional elements
      elasticTangent(2, 0) =  K23G); // than just mirror all, including zeros,
                                  // which are there already
      elasticTangent.PutElement(3, 2, K23G);
      elasticTangent.PutElement(3, 3, K43G);
      elasticTangent.PutElement(4, 4, 2.0 * G);
      elasticTangent.PutElement(5, 5, 2.0 * G);
      elasticTangent.PutElement(6, 6,
                     2.0 * G); // rest of the matrix is filled with zeros...

      // getting lambda and Dep

      A.Transpose(&numerator);
      numerator.Multiply(&elasticTangent, &numerator); // numerator=aT*Del -->Numerator of
                                            // Lambda without multiplication by
                                            // dEpsilon
      numerator.Multiply(&dg_dsigma, &denominator);

      double Lambda = fValue / denominator.getElement(1, 1);

      A.Multiply(Lambda, &TEMP); // delta epsilon plastic= -delta epsilon
                                 // elastic

      // dbg << "Delta Epsilon Plastic:"<<"\n";
      // TEMP.Print();

      for (int i = 1; i < 7; i++)
        epStrain[i - 1] = TEMP.getElement(i, 1); // epsilon pl change

      elasticTangent.Multiply(&dg_dsigma, &TEMP);
      TEMP.Multiply(-Lambda, &dSigma);
      // final result for stress change, with the negative sign, so the stress
      // should be ADDED
      // be ADDED to get the right result of corrected stresses.
      for (int i = 0; i < 6; i++)
        dSigma[i] = dSigma.getElement(i + 1, 1);
      Point->Update(epStrain, zeros, dSigma, 0);
    }
    if (numberIter > 10) {
      // dbg << "Drift Correction Procedure failed"<<"\n";
      correct = FALSE;
    }
  } while (correct == TRUE);

  /*


bool correct;
double fValue, dP0Star, Lambda, dSigma[6], epStrain[6], zeros[7];
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
				BBMMatrix dg_dsigma(6,1);	//dG/dsigma
				BBMMatrix PEP (1,1); //dp0* /depsvpl
				BBMMatrix MM (6,1); // will be vector(1,1,1,0,0,0)T
				BBMMatrix elasticTangent (6,6); //D elastic matrix...
				BBMMatrix DPZEROSTAR (1,1);
				BBMMatrix TEMP (1,1);
				BBMMatrix dSigma (6,1);
				MM.PutElement (1,1,1);
				MM.PutElement (2,1,1);
				MM.PutElement (3,1,1); //rest is zero as initialized


				SpecificVolume=PointCopy.getSpecVol();  //specific volume need to be used the right one
				//dbg << "Specific Volume:"<<SpecificVolume<<"\n";
				PZeroStar=PointCopy.p0Star ();

				//convention - matrices names are made from CAPITALIZED letters


				meanStress=PointCopy.getmeanStress();
				shearStress=PointCopy.shearStress();
				Suction=PointCopy.suction();

				LambdaS=(1-r)*exp(-Beta*Suction)+r;
				LambdaS=LambdaS*LambdaZero;
				Fraction=(LambdaZero-KappaP)/(LambdaS- KappaP);
				PZero=pc*pow(PZeroStar/pc,Fraction);  //get.calculated.pzero;

				//dbg << "PZero = "<<PZero<<"\n";
				//dbg << "PZeroStar = "<<PZeroStar<<"\n";
				//dbg << "p = "<<meanStress<<"\n";
				//dbg << "q = "<<shearStress<<"\n";
				//dbg << "s = "<<Suction<<"\n";

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
				//dbg << "A:"<<"\n"; A.Print();
				//dF/dsigma - inserted into A

				if (nonAssociated)
				{
					temp=alfa*(2*Point->stress[0]-Point->stress[1]-Point->stress[2])+M*M/3*(2*meanStress+k*Suction-PZero);
					dg_dsigma.PutElement(1,1,temp);
					temp=alfa*(2*Point->stress[1]-Point->stress[0]-Point->stress[2])+M*M/3*(2*meanStress+k*Suction-PZero);
					dg_dsigma.PutElement(2,1,temp);
					temp=alfa*(2*Point->stress[2]-Point->stress[0]-Point->stress[1])+M*M/3*(2*meanStress+k*Suction-PZero);
					dg_dsigma.PutElement(3,1,temp);
					temp=6*alfa*Point->stress[3];
					dg_dsigma.PutElement(4,1,temp);
					temp=6*alfa*Point->stress[4];
					dg_dsigma.PutElement(5,1,temp);
					temp=6*alfa*Point->stress[5];
					dg_dsigma.PutElement(6,1,temp);
				}
				else */ /*A.Copy (&dg_dsigma);


                                //d

                                temp=0;
                                temp=-M*M*(meanStress+k*Suction)*Fraction*pow((PZeroStar/pc),
                                Fraction-1);

                                P.PutElement (1,1,temp);

                                //dbg << "P:"<<"\n"; P.Print();

                                temp=PZeroStar*SpecificVolume/(LambdaZero-KappaP);
                                PEP.PutElement (1,1,temp); //dP0* /depsplv
                                //dbg << "dP0Star/Depsvpl:"<<temp<<"\n";
                                //elasticTangent... elastic matrix... values of K. Here we
                                need values of K at the point...
                                We need to have the K - bulk modulus of the soil
                                calculated, and then it is
                                possible to fill into the elasticTangent Matrix...
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
                                          dbg << "WARNING !!! Mean stress too
  low. Mean stress is adjusted to d_minMeanStress value!!!"<<"\n";
                                          meanStress=d_minMeanStress;
                                  }

                                  K=meanStress*SpecificVolume/KappaP;  //tangent
  bulk modulus K=p*specificVol/KappaP, from eq dSpecVol=-KappaP*ln(p/pzero);
                                  //dbg << "K="<<K<<"\n";

                                  // ****************************** K need
  correcting, but with const. suction seems to be ok
  ******************************

                                  // Tiny set in stdafx to 1e-12
                                  // calculate helpful variables:
                                  double K43G=K+4.0*G/3.0;
                                  double K23G=K-2.0*G/3.0;


                                  // Fill in the matrix:

                                  elasticTangent.PutElement (1,1,K43G);
                                  elasticTangent.PutElement (1,2,K23G);
                                  elasticTangent.PutElement (1,3,K23G); //rest of the line
  are zeros and rightly so
                                  elasticTangent.PutElement (2,1,K23G);
                                  elasticTangent.PutElement (2,2,K43G);
                                  elasticTangent.PutElement (2,3,K23G); //yes, the matrix
  is symmetrical, but it is faster to put this 3 additional elements
                                  elasticTangent.PutElement (3,1,K23G); //than just mirror
  all, including zeros, which are there already
                                  elasticTangent.PutElement (3,2,K23G);
                                  elasticTangent.PutElement (3,3,K43G);
                                  elasticTangent.PutElement (4,4,2*G);
                                  elasticTangent.PutElement (5,5,2*G);
                                  elasticTangent.PutElement (6,6,2*G); //rest of the matrix
  is filled with zeros...

                                  A.Transpose (&TEMP);
                                  TEMP.Multiply (&elasticTangent,&TEMP);
                                  TEMP.Multiply (&dg_dsigma, &TEMP);

                                  temp=TEMP.getElement (1,1);  //first part of
  the denominator
                                  //dbg << "First Part of Denominator done...
  temp="<<temp<<"\n";


                                  P.Multiply (&PEP,&TEMP);
                                  MM.Transpose (&MM);
                                  TEMP.Multiply (&MM,&TEMP);	//MM is
  transposed
                                  TEMP.Multiply (&dg_dsigma,&TEMP);

                                  temp=temp+TEMP.getElement (1,1); //'end of the
  denominator

                                  //dbg << "Denominator="<<temp<<"\n";

                                  Lambda=fValue*(PZero+Suction*k)*(PZero+Suction*k)/temp;
  //because we need the value, not the reduced value...

                                  //dbg << "Lambda="<<Lambda<<"\n";

                                  A.Multiply (Lambda, &TEMP); //delta epsilon
  plastic= -delta epsilon elastic

                                  //dbg << "Delta Epsilon Plastic:"<<"\n";
                                  //TEMP.Print();

                                  for (int i=1; i<7; i++) epStrain
  [i-1]=TEMP.getElement (i,1); //epsilon pl change
                                  temp=epStrain[0]+epStrain[1]+epStrain[2];
  //DepsilonV
                                  dP0Star=PEP.getElement(1,1)*temp;
                                  //dbg << "dP0Star="<<*dP0Star<<"\n";

                                  elasticTangent.Multiply (&dg_dsigma, &TEMP);
                                  TEMP.Multiply (-Lambda, &dSigma);
                          //final result for stress change, with the negative
  sign, so the stress should be ADDED
                          //be ADDED to get the right result of corrected
  stresses.

                          //dbg << "Delta Sigma="<<"\n";
                          //dSigma->Print();


                          //dbg << "Press any key (end of Correct Drift
  Procedure)"<<"\n";
                          //getchar();
                          for (int i=0; i<6; i++) dSigma[i]=dSigma.getElement
  (i+1,1);
                          PointCopy.Update (epStrain, zeros, dSigma,
  dP0Star);

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

  // dbg << "Correct Drift Procedure entered!"<<"\n";

  BBMMatrix A(6, 1); // dF/dsigma
  BBMMatrix denominator(1, 1);

  BBMMatrix dg_dsigma(6, 1);  // dG/dsigma, in case of associated flow
                       // (nonAssociated==false) same as A
  BBMMatrix MM(6, 1);  // will be vector(1,1,1,0,0,0)T
  BBMMatrix elasticTangent(6, 6); // D elastic matrix...
  BBMMatrix dEps(6, 1);
  BBMMatrix dSigma(6, 1);

  double dSigma[6], epStrain[6], zeros[7];
  for (int i = 0; i < 7; i++)
    zeros[i] = 0;

  for (int i = 1; i < 6; i++)
    dEps.PutElement(i, 1, epStrain[i - 1]); // increase of epsilon copied
  // we do not have g, as a == g in associated flow rule.

  bool correct;

  double I1 = PointOld->firstInvariant();
  double I2 = PointOld->secondInvariant();
  double J2 = PointOld->secondDevInvariant();
  if (std::abs(J2) < TINY)
    J2 = TINY;
  double J3 = PointOld->thirdDevInvariant();
  double shearStress = PointOld->shearStress();
  if (std::abs(shearStress) < TINY)
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
    dI2_dSig(6, 1), dI3_dSig(6, 1), TEMP(6, 1), numerator(6, 1);

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

  dJ2_dSig.Copy(&dg_dsigma);
  dg_dsigma.Multiply(0.21022410391343 * Mpsi * Factor025 / shearStress, &dg_dsigma);
  dJ3_dSig.Copy(&TEMP);
  TEMP.Multiply(-1.419012700740643 * (alpha4 - 1) * Mpsi * sqrt(J2) *
                    Factor075 /
                    (sqrt(3.0) * shearStress * shearStress * shearStress) +
                  0.94600846716043 * (alpha4 - 1) * cohesion * Mpsi *
                    Factor075 / (shearStress * shearStress * shearStress),
                &TEMP);
  dg_dsigma.Add(&TEMP, &dg_dsigma);
  dq_dSig.Copy(&TEMP);
  TEMP.Multiply(
    -2.838025401481287 * (alpha4 - 1) * cohesion * Mpsi * J3 * Factor075 /
        (shearStress * shearStress * shearStress * shearStress) +
      4.257038102221929 * (alpha4 - 1) * Mpsi * sqrt(J2) * J3 * Factor075 /
        (sqrt(3.0) * shearStress * shearStress * shearStress * shearStress),
    &TEMP);
  dg_dsigma.Add(&TEMP, &dg_dsigma);
  dI1_dSig.Copy(&TEMP);
  TEMP.Multiply(
    -1 / 3.0,
    &TEMP); // First correction: sign here, but most likely not enough!!!
  dg_dsigma.Add(&TEMP, &dg_dsigma);

  // FINISHED dQ/dSigma

  double K43G = K + 4.0 * G / 3.0;
  double K23G = K - 2.0 * G / 3.0;

  elasticTangent(0, 0) =  K43G);
  elasticTangent(0, 1) =  K23G);
  elasticTangent(0, 2) =  K23G); // rest of the line are zeros and rightly so
  elasticTangent(1, 0) =  K23G);
  elasticTangent(1, 1) =  K43G);
  elasticTangent(1, 2) =  K23G); // yes, the matrix is symmetrical, but it is
                              // faster to put this 3 additional elements
  elasticTangent.PutElement(
    3, 1,
    K23G); // than just mirror all, including zeros, which are there already
  elasticTangent.PutElement(3, 2, K23G);
  elasticTangent.PutElement(3, 3, K43G);
  elasticTangent.PutElement(4, 4, 2.0 * G);
  elasticTangent.PutElement(5, 5, 2.0 * G);
  elasticTangent.PutElement(6, 6, 2.0 * G); // rest of the matrix is filled with zeros...

  // getting lambda and Dep

  A.Transpose(&numerator);
  numerator.Multiply(&elasticTangent, &numerator); // numerator=aT*Del -->Numerator of
                                        // Lambda without multiplication by
                                        // dEpsilon
  numerator.Multiply(&dg_dsigma, &denominator);

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

      double Lambda = fValue / denominator.getElement(1, 1);
      A.Multiply(Lambda, &TEMP); // delta epsilon plastic= -delta epsilon
                                 // elastic

      // dbg << "Delta Epsilon Plastic:"<<"\n";
      // TEMP.Print();

      for (int i = 1; i < 7; i++)
        epStrain[i - 1] = TEMP.getElement(i, 1); // epsilon pl change
      elasticTangent.Multiply(&dg_dsigma, &TEMP);
      TEMP.Multiply(-Lambda, &dSigma);
      // final result for stress change, with the negative sign, so the stress
      // should be ADDED
      // be ADDED to get the right result of corrected stresses.
      for (int i = 0; i < 6; i++)
        dSigma[i] = dSigma.getElement(i + 1, 1);
      Point->Update(epStrain, zeros, dSigma, 0);
    }
    if (numberIter > 10) {
      // dbg << "Drift Correction Procedure failed"<<"\n";
      correct = FALSE;
    }
  } while (correct == TRUE);

  /*
StateMohrCoulomb pointCopy, pointEnd;
PointOld->Copy(&PointCopy);
Point->Copy(&PointEnd);

bool correct;
double fValue, dP0Star, Lambda, dSigma[6], epStrain[6], zeros[7];
for (int i=0; i<7; i++) zeros[i]=0;

double temp, PZero, meanStress, shearStress, Suction, PZeroStar, SpecificVolume, LambdaS, Fraction;
// we do not have g, as a == g in associated flow rule.
	do
	{
	checkYield (PointEnd.state, pointEnd.stress, pointEnd.suction(), &fValue); //20 Feb 2006, preparations for the drift correction

	if (std::abs(fValue)>d_yieldTol) correct=TRUE; else correct=FALSE;
		if (correct==TRUE)
		{
				// CORRECT FOR DRIFT
				//HERE THE DRIFT WILL BE CORRECTED BY USING THE D MATRIX FROM THE FORBIDDEN SPACE
				//ALTHOUGH BECAUSE THE DRIFT WILL BE CHECKED AGAIN, IT SHOULDN'T POSE MUCH PROBLEM.
				BBMMatrix A (6,1); //dF/dsigma
				BBMMatrix P (1,1); //dF/dP0*
				BBMMatrix dg_dsigma (6,1);
				BBMMatrix PEP (1,1); //dp0* /depsvpl
				BBMMatrix MM (6,1); // will be vector(1,1,1,0,0,0)T
				BBMMatrix elasticTangent (6,6); //D elastic matrix...
				BBMMatrix DPZEROSTAR (1,1);
				BBMMatrix TEMP (1,1);
				BBMMatrix dSigma (6,1);
				MM.PutElement (1,1,1);
				MM.PutElement (2,1,1);
				MM.PutElement (3,1,1); //rest is zero as initialized



				SpecificVolume=PointCopy.getSpecVol();  //specific volume need to be used the right one
				//dbg << "Specific Volume:"<<SpecificVolume<<"\n";
				PZeroStar=PointCopy.p0Star ();

				//convention - matrices names are made from CAPITALIZED letters


				meanStress=PointCopy.getmeanStress();
				shearStress=PointCopy.shearStress();
				Suction=PointCopy.suction();

				LambdaS=(1-r)*exp(-Beta*Suction)+r;
				LambdaS=LambdaS*LambdaZero;
				Fraction=(LambdaZero-KappaP)/(LambdaS- KappaP);
				PZero=pc*pow(PZeroStar/pc,Fraction);  //get.calculated.pzero;

				//dbg << "PZero = "<<PZero<<"\n";
				//dbg << "PZeroStar = "<<PZeroStar<<"\n";
				//dbg << "p = "<<meanStress<<"\n";
				//dbg << "q = "<<shearStress<<"\n";
				//dbg << "s = "<<Suction<<"\n";

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
				//dbg << "A:"<<"\n"; A.Print();
				//dF/dsigma - inserted into A

				if (nonAssociated)
				{
					temp=alfa*(2*Point->stress[0]-Point->stress[1]-Point->stress[2])+M*M/3*(2*meanStress+k*Suction-PZero);
					dg_dsigma.PutElement(1,1,temp);
					temp=alfa*(2*Point->stress[1]-Point->stress[0]-Point->stress[2])+M*M/3*(2*meanStress+k*Suction-PZero);
					dg_dsigma.PutElement(2,1,temp);
					temp=alfa*(2*Point->stress[2]-Point->stress[0]-Point->stress[1])+M*M/3*(2*meanStress+k*Suction-PZero);
					dg_dsigma.PutElement(3,1,temp);
					temp=6*alfa*Point->stress[3];
					dg_dsigma.PutElement(4,1,temp);
					temp=6*alfa*Point->stress[4];
					dg_dsigma.PutElement(5,1,temp);
					temp=6*alfa*Point->stress[5];
					dg_dsigma.PutElement(6,1,temp);
				}
				else */ /* A.Copy (&dg_dsigma);

                                //d

                                temp=0;
                                temp=-M*M*(meanStress+k*Suction)*Fraction*pow((PZeroStar/pc),
                                Fraction-1);

                                P.PutElement (1,1,temp);

                                //dbg << "P:"<<"\n"; P.Print();

                                temp=PZeroStar*SpecificVolume/(LambdaZero-KappaP);
                                PEP.PutElement (1,1,temp); //dP0* /depsplv
                                //dbg << "dP0Star/Depsvpl:"<<temp<<"\n";
                                //elasticTangent... elastic matrix... values of K. Here we
                                need values of K at the point...
                                We need to have the K - bulk modulus of the soil
                                calculated, and then it is
                                possible to fill into the elasticTangent Matrix...
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
                                          dbg << "WARNING !!! Mean stress too
  low. Mean stress is adjusted to d_minMeanStress value!!!"<<"\n";
                                          meanStress=d_minMeanStress;
                                  }

                                  K=meanStress*SpecificVolume/KappaP;  //tangent
  bulk modulus K=p*specificVol/KappaP, from eq dSpecVol=-KappaP*ln(p/pzero);
                                  //dbg << "K="<<K<<"\n";

                                  // ****************************** K need
  correcting, but with const. suction seems to be ok
  ******************************

                                  // Tiny set in stdafx to 1e-12
                                  // calculate helpful variables:
                                  double K43G=K+4.0*G/3.0;
                                  double K23G=K-2.0*G/3.0;


                                  // Fill in the matrix:

                                  elasticTangent.PutElement (1,1,K43G);
                                  elasticTangent.PutElement (1,2,K23G);
                                  elasticTangent.PutElement (1,3,K23G); //rest of the line
  are zeros and rightly so
                                  elasticTangent.PutElement (2,1,K23G);
                                  elasticTangent.PutElement (2,2,K43G);
                                  elasticTangent.PutElement (2,3,K23G); //yes, the matrix
  is symmetrical, but it is faster to put this 3 additional elements
                                  elasticTangent.PutElement (3,1,K23G); //than just mirror
  all, including zeros, which are there already
                                  elasticTangent.PutElement (3,2,K23G);
                                  elasticTangent.PutElement (3,3,K43G);
                                  elasticTangent.PutElement (4,4,2*G);
                                  elasticTangent.PutElement (5,5,2*G);
                                  elasticTangent.PutElement (6,6,2*G); //rest of the matrix
  is filled with zeros...

                                  A.Transpose (&TEMP);
                                  TEMP.Multiply (&elasticTangent,&TEMP);
                                  TEMP.Multiply (&dg_dsigma, &TEMP);

                                  temp=TEMP.getElement (1,1);  //first part of
  the denominator
                                  //dbg << "First Part of Denominator done...
  temp="<<temp<<"\n";


                                  P.Multiply (&PEP,&TEMP);
                                  MM.Transpose (&MM);
                                  TEMP.Multiply (&MM,&TEMP);	//MM is
  transposed
                                  TEMP.Multiply (&dg_dsigma,&TEMP);

                                  temp=temp+TEMP.getElement (1,1); //'end of the
  denominator

                                  //dbg << "Denominator="<<temp<<"\n";

                                  Lambda=fValue*(PZero+Suction*k)*(PZero+Suction*k)/temp;
  //because we need the value, not the reduced value...

                                  //dbg << "Lambda="<<Lambda<<"\n";

                                  A.Multiply (Lambda, &TEMP); //delta epsilon
  plastic= -delta epsilon elastic

                                  //dbg << "Delta Epsilon Plastic:"<<"\n";
                                  //TEMP.Print();

                                  for (int i=1; i<7; i++) epStrain
  [i-1]=TEMP.getElement (i,1); //epsilon pl change
                                  temp=epStrain[0]+epStrain[1]+epStrain[2];
  //DepsilonV
                                  dP0Star=PEP.getElement(1,1)*temp;
                                  //dbg << "dP0Star="<<*dP0Star<<"\n";

                                  elasticTangent.Multiply (&dg_dsigma, &TEMP);
                                  TEMP.Multiply (-Lambda, &dSigma);
                          //final result for stress change, with the negative
  sign, so the stress should be ADDED
                          //be ADDED to get the right result of corrected
  stresses.

                          //dbg << "Delta Sigma="<<"\n";
                          //dSigma->Print();


                          //dbg << "Press any key (end of Correct Drift
  Procedure)"<<"\n";
                          //getchar();
                          for (int i=0; i<6; i++) dSigma[i]=dSigma.getElement
  (i+1,1);
                          PointEnd.Update (epStrain, zeros, dSigma, dP0Star);

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
  // dbg << purelyElastic<<"\n";

  if (purelyElastic) {
    // finalState.Copy(initialState);
    // update with what we have - all updated after the procedure
    // dbg << "Purely Elastic Step"<<"\n";
  } else // we have elasto-plastic step
  {
    if (onTheYieldLocus) {
      unloading =
        checkGradient(initialState, &finalState); // checking if there is any
                                                  // unloading part; in
                                                  // finalState purely elastic
                                                  // values; we care only for
                                                  // correct change of suction
      // dbg << "\n"<<"\n"<<"unloading="<<unloading<<"\n";
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
  // for (int i=0; i<3; i++) dbg << strainIncrement[i]<<"  "<<"\n";
  // for (int i=0; i<3; i++) dbg << finalState.stress[i]<<"  "<<"\n";
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
    // dbg << "In iteration "<<iter<<" alfa="<<alfa<<" and F="<<fAlfa<<"\n";

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
  // dbg << "Error in procedure findYieldOriginal"<<"\n";
  // dbg << "After "<<d_maxIter<<" iterations crossing point not found"<<"\n";
  // dbg << "This is likely to cause incorrect results... Results obtained should
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
                  dbg << "Non associated flow rule used. Value of
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
          dbg << "Model Parameters:"<<"\n";
          dbg << "Shear Modulus G="<<G<<"\n";
          dbg << "Kappa p="<<KappaP<<"\n";
          dbg << "Kappa s="<<KappaS<<"\n";
          dbg << "Pc="<<pc<<"\n";
          dbg << "k="<<k<<"\n";
          dbg << "r="<<r<<"\n";
          dbg << "Beta="<<Beta<<"\n";
          dbg << "Lambda (0)="<<LambdaZero<<"\n";
          dbg << "N (0)="<<NZero<<"\n";
          dbg << "M="<<M<<" deg"<<"\n";

          dbg << "Computed parameters:"<<"\n";
          dbg << "K="<<K<<"\n";
          dbg << "Ks="<<Ks<<"\n";
          dbg << "Kp="<<Kp<<"\n";*/
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
  //dbg << Max<<"\n"; // in first line we put how many points we have...
  fprintf( stream, "%d\n", Max );
  minP=-k*suction;
  PZero=LambdaZero*((1-r)*exp(-1*Beta*suction)+r);
  PZero=(LambdaZero-KappaP)/(PZero-KappaP);
  PZero=pc*pow((state[0]/pc),PZero);
  maxP=PZero;
  //dbg << minP<<"\n";	//minimum
  //dbg << maxP<<"\n";	//maximum, used to set the view...
  fprintf( stream, "%f\n",minP );
  fprintf( stream, "%f\n",maxP );
  difference=maxP-minP;	//this is the max difference...
  //dbg << difference<<"\n";
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
  int Steps = 6;

  double Order = 5;
  double A[8][8]; // matrix must be this size to be used in the doRungeKutta
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

  time = doRungeKuttaEqualStep(A, B, BRes, C, point, purelyPlasticStrain,
                             stressIncrAbs, &RelativeError, numberIter,
                             Order, Steps, errorEstimate);
  return time;
}


double
ShengMohrCoulomb::plasticEuler(StateMohrCoulomb* point, double* epStrain,
                               double* absStress, int numberIterations)
{

  BBMMatrix dSigma(6, 1);
  double dP0Star;
  double dSigma[6];
  double currentStrain[7], plasticStrain[6];
  double PZeroStarTot = 0;
  double dSigmaDrift[6], dLambdaDrift;
  StateMohrCoulomb oldState;
  Point->Copy(&oldState);

  vector<double> stressInc;
  vector<double> strainInc;
  vector<double> p0Lambda;
  vector<double>::iterator Iter;

  // dbg << "Euler Procedure Entered"<<"\n";
  // dbg << "Number of iterations in Euler Algorithm:"<<numberIterations<<"\n";

  clock_t startTime, endTime;
  startTime = clock();

  for (int i = 0; i < 7; i++)
    absStress[i] = 0;
  for (int i = 0; i < 7; i++)
    currentStrain[i] = epStrain[i] / numberIterations;
  for (int loop = 0; loop < numberIterations; loop++) {
    double fValue = 0;
    if (USE_NICE_SCHEME > 0)
      fValue = computeYieldFunction(Point);
    calcPlastic(*Point, currentStrain, &dSigma, plasticStrain, &dP0Star,
                fValue, dSigmaDrift, &dLambdaDrift);
    for (int i = 0; i < 6; i++)
      dSigma[i] = dSigma.getElement(i + 1, 1);
    Point->Update(plasticStrain, currentStrain, dSigma, dP0Star);

    p0Lambda.push_back(dP0Star);
    p0Lambda.push_back(plasticStrain[0]);
    PZeroStarTot = PZeroStarTot + dP0Star;

    for (int i = 0; i < 6; i++)
      absStress[i] = absStress[i] + std::abs(dSigma[i]);
    absStress[6] = absStress[6] + std::abs(dP0Star);

    /*for (int i=0; i<6; i++)
            {
            stressInc.push_back (dSigma[i]);  //is ok, as updated total stress?
            strainInc.push_back (currentStrain [i]);
            }
    strainInc.push_back (currentStrain [6]);*/
  }

  endTime = clock();

  /*
  dbg << "calculation took:"<<double(endTime-startTime)/CLOCKS_PER_SEC<<"
  s."<<"\n";
  dbg << "The d_integrationTol parameter is equal to:"<<d_integrationTol<<"\n";
  dbg << "Total number of steps done:"<<numberIterations<<"\n";

  dbg << "Total PZeroStar change: "<<PZeroStarTot<<"\n";
  checkYield (Point->state, point->stress, point->strain[6],&fValue);
  dbg << "Yield Function value is:"<<fValue<<"\n";

  dbg << "Over the whole step change of stress is:"<<"\n";
  for (int i=0; i<6; i++)
  {
          dbg << "s["<<i<<"]="<<Point->stress[i]-oldState.stress[i]<<"\n";
  }

  dbg << "Initial specific volume="<<oldState.getSpecVol()<<"\n";
  dbg << "Final specific volume="<<Point->getSpecVol()<<"\n";
  dbg << "Change of specific volume is equal
  to="<<(oldState.getSpecVol()-Point->getSpecVol())<<"\n";
  dbg << "Change of mean stress p is equal
  to:"<<Point->getmeanStress()-oldState.getmeanStress()<<"\n";
  dbg << "Change of shear stress q is equal
  to:"<<Point->shearStress()-oldState.shearStress()<<"\n";



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
  for (Iter=p0Lambda.begin() ; Iter!=p0Lambda.end(); )
  {
  for (int i=0; i<2; i++)
          {
          fprintf( PZeroFile, "%.25f , ",*Iter);
          Iter++;
          }
  fprintf( PZeroFile, "\n");
  }
  fclose (PZeroFile);



  dbg << "Press any key..."<<"\n";
  getchar();
  */

  return (endTime - startTime);
}


double
ShengMohrCoulomb::doRungeKuttaEqualStep(double A[][8], double* B, double* BRes,
                                      double* C, StateMohrCoulomb* point,
                                      double* epStrain, double* absStress,
                                      double* RelError, int numberIter,
                                      double Order, int Steps,
                                      bool errorEstimate)
/*
This procedure calculate any Runge - Kutta pair, given the coefficient of
stress in A for each x used to calculate values in C, where B gives the
coefficient to calculate
error estimate (or lower order solution) and BRes to calculate result. The
procedure will integrate whole step of strain
re-using the initial derivative for rejected steps.
Order contain order of the method (required for substep prediction)
Steps contain the number of stages to get the result
errorEstimate if true - the B table contains error estimate. If false, it
contains 4th order solution
*/

{
  BBMMatrix dSigma(6, 1);
  StateMohrCoulomb midState[8], oldState, trialState;

  double dSigma[8][6], stressInc[7], Result[7];
  double dP0Star[8], p0StarInc, plasticStrain[8][6],
    plasticStrainInc[6];
  double Error[7], rError, TotrError, totalSize, stepLength, Temp, methodPower;
  double frequency = 15000 / Order; // how often display info about steps

  for (int i = 0; i < Order; i++) {
    for (int j = 0; j < 6; j++)
      dSigma[i][j] = 0;
    dP0Star[i] = 0;
  }

  // bool finished=false, stepAccepted=false;
  double substepStrain[7], currentStrain[7];
  double microStep = 0;
  double stepAccuracyCheck;
  double dStressDrift[6], dLambdaDrift;

  int stepNo = 0;
  vector<double> stressInc;
  vector<double> strainInc;
  vector<double> p0Lambda;
  vector<double>::iterator Iter;

  stepLength = 1.0 / numberIter;

  for (int i = 0; i < 7; i++) {
    currentStrain[i] = 0;
    substepStrain[i] =
      stepLength * epStrain[i]; // strain increment in all steps (equally sized)
    absStress[i] = 0;
  }
  totalSize = 0; // newStepSize=1;
  methodPower = pow(2.0, Order) * d_integrationTol;
  // stepAccepted=true;
  rError = 0;
  TotrError = 0;

  for (int loop = 0; loop < numberIter; loop++) {

    stepNo++;

    /*if (stepNo>1)
{
	if (newStepSize>10) newStepSize=10;
	if (newStepSize<0.1) newStepSize=0.1;
}*/ // limiting step increase/decrease does not improve results/ enables faster
    // convergence...

    // dbg << "Step Length="<<stepLength<<"\n";
    // dbg << "Current strain [0]="<<currentStrain[0]<<"\n";

    for (int i = 0; i < Steps; i++)
      point->Copy(&midState[i]); // point is unchanged in calcPlastic procedure

    // Below the main R-K loop to calculate the value of intermediate stresses;
    // values stored in stressInc[][]
    // reuseStep=false;

    for (int rkloop = 0; rkloop < Steps; rkloop++) {
      for (int i = 0; i < 6; i++) {
        stressInc[i] = 0;
        plasticStrainInc[i] = 0;
      }
      p0StarInc = 0;
      for (int i = 0; i < 7; i++)
        currentStrain[i] =
          C[rkloop] *
          substepStrain[i]; // set the beginning point of the procedure
      for (int i = 0; i < rkloop; i++) {
        for (int j = 0; j < 6; j++) {
          stressInc[j] = stressInc[j] + A[rkloop][i] * dSigma[i][j];
          plasticStrainInc[j] =
            plasticStrainInc[j] + A[rkloop][i] * plasticStrain[i][j];
        }
        p0StarInc = p0StarInc + A[rkloop][i] * dP0Star[i];
      }
      midState[rkloop].Update(plasticStrainInc, currentStrain, stressInc,
                              p0StarInc);
      calcPlastic(midState[rkloop], substepStrain, &dSigma, plasticStrainInc,
                  &dP0Star[rkloop], computeYieldFunctionNN(Point),
                  dStressDrift, &dLambdaDrift);
      for (int i = 0; i < 6; i++) {
        dSigma[rkloop][i] = dSigma.getElement(i + 1, 1);
        plasticStrain[rkloop][i] = plasticStrainInc[i];
      }
    }

    // needed: result, error

    for (int i = 0; i < 6; i++) {
      Result[i] = 0;
      plasticStrainInc[i] = 0;
      for (int j = 0; j < Steps; j++) {
        Result[i] = Result[i] + BRes[j] * dSigma[j][i];
        plasticStrainInc[i] =
          plasticStrainInc[i] + BRes[j] * plasticStrain[j][i];
      }
    }
    Result[6] = 0;
    for (int j = 0; j < Steps; j++)
      Result[6] = Result[6] + BRes[j] * dP0Star[j];

    for (int i = 0; i < 7; i++)
      Error[i] = 0;

    for (int i = 0; i < Steps; i++) {
      for (int j = 0; j < 6; j++)
        Error[j] = Error[j] + B[i] * dSigma[i][j];
      Error[6] = Error[6] + B[i] * dP0Star[i];
    }
    if (!errorEstimate)
      for (int i = 0; i < 7; i++)
        Error[i] = Error[i] - Result[i]; // error estimate calculated in case we
                                         // have lower order solution instead of
                                         // error estimate

    // check the error norm

    switch (int(d_tolMethod)) {
      case 0: {
        rError = checkNorm(Result, Result[6], point, Error); // returns rError
      } break;
      case 1: {
        // SLOAN NORM
        rError = checkNormSloan(Result, Result[6], point, Error); // returns
                                                                  // rError
        // dbg << "Runge - Kutta constant size procedure. rError="<<rError<<"\n";
        // dbg << "Error[0]="<<Error[0]<<"  Result[0]="<<Result[0]<<"\n";
      } break;
      default: {
        cout << "ERROR !!!! Improper d_tolMethod in increment.dta" << "\n";
        getchar();
      }
    }

    // dbg << "Procedure R-K constant. rError="<<rError<<"\n";

    for (int i = 0; i < 7; i++)
      if (!finite(Result[i])) {
        Result[i] = 0;
        if (rError < methodPower)
          rError = methodPower;
      }

    if ((Point->getmeanStress() + (Result[0] + Result[1] + Result[2])) / 3 <
        0) {
      // dbg << "Mean Stress less then 0!!!
      // Result:"<<((Result[0]+Result[1]+Result[2])/3)<<"  Mean
      // stress:"<<Point->getmeanStress()<<"\n";

      if (rError < methodPower)
        rError = methodPower;
    }
    if ((Point->p0Star() + Result[6]) < 0) {
      // dbg << "P Zero Star less then 0!!!"<<"\n";
      if (rError < methodPower)
        rError = methodPower;
    }
    // here we need to update all the point data.
    Point->Copy(&oldState);
    Point->Update(plasticStrainInc, substepStrain, Result, Result[6]);
    Point->Copy(&trialState);
    // this value is not used in any calculations, it is just to show at the end
    // of the step the p0* increase

    // dbg << "Procedure R-K constant. rError="<<rError<<"\n";

    if (d_driftCorrection == 3)
      correctDrift(Point);
    if (d_driftCorrection == 2)
      correctDriftBeg(
        Point,
        &oldState); // value of oldState copied before updating the point

    // re - evaluate the error in the point:

    for (int i = 0; i < 6; i++) {
      Error[i] = Error[i] + Point->stress[i] - trialState.stress[i];
    }
    Error[6] = Error[6] + Point->p0Star() - trialState.p0Star();
    // error vector updated, norm should be re-evaluated:

    switch (int(d_tolMethod)) {
      case 0: {
        Temp = checkNorm(Result, Result[6], point, Error); // returns rError
        if (Temp > rError)
          rError = Temp;
      } break;
      case 1: {
        // SLOAN NORM
        Temp =
          checkNormSloan(Result, Result[6], point, Error); // returns rError
        if (Temp > rError)
          rError = Temp;
      } break;
      default: {
        cout << "ERROR !!!! Improper d_tolMethod in increment.dta" << "\n";
        getchar();
      }
    }
    if (!finite(rError))
      if (rError < methodPower)
        rError = methodPower;

    // dbg << "Procedure R-K constant. rError="<<rError<<"\n";
    TotrError = TotrError + rError;
    // dbg << "Procedure R-K constant. Total rError="<<TotrError<<"\n";

    // this is done anyway
    for (int i = 0; i < 6; i++)
      absStress[i] = absStress[i] + std::abs(Point->stress[i] - oldState.stress[i]);
    absStress[6] = absStress[6] + std::abs(Point->p0Star() - oldState.p0Star());
    stepAccuracyCheck = totalSize;
    microStep = microStep + stepLength;
    totalSize = totalSize + microStep; // total part of step done updated
    microStep = microStep - (totalSize - stepAccuracyCheck);
    // if (totalSize>=1) finished=true;
    Temp = double(stepNo) / frequency;
    /*if (modf(Temp,&Temp)==0)
            {
                    dbg << "Step number:"<<stepNo<<"\n";
                    dbg << "Total size done is: "<<totalSize<<" of whole step.
       Current stepLength="<<stepLength<<"\n";
             }*/
    /*
            //Debug
            strainInc.push_back (stepLength);
            for (int i=0; i<6; i++) stressInc.push_back (Point->stress[i]);
            stressInc.push_back (Point->p0Star());
            //End Debug */
  }

  *RelError = TotrError / numberIter;

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
ShengMohrCoulomb::doRungeKuttaExtrapol(double A[][8], double* B, double* BRes,
                                     double* C, StateMohrCoulomb* point,
                                     double* epStrain, double* absStress,
                                     int* numberIter, double Order,
                                     int Steps, bool errorEstimate)
{
  /*
  this procedure tend to use Runge Kutta scheme. In case the step is not
  accepted, and the substep size is just a bit smaller
  than current one , instead of cancelling the step, it goes into extrapolation
  procedure, just to save the substep done.
  */

  BBMMatrix dSigma(6, 1);
  StateMohrCoulomb trialState, oldState;
  double WORTH_EXTRAPOL = 3;

  // double dSigma[8][6];
  // double dP0Star[8];
  double rError, newStepSize, totalSize, stepLength, Temp; // methodPower;
  double frequency = 15000 / Order; // how often display info about steps
  // bool reuseStep=false;

  stepLength = 1;
  totalSize = 0;
  newStepSize = 1;
  // methodPower=pow(2.0,Order)*d_integrationTol;

  // for (int i=0; i<Order; i++) {
  // for (int j=0; j<6; j++) dSigma[i][j]=0;
  // dP0Star[i]=0;	}

  bool finished = false, stepAccepted = false;
  double substepStrain[7]; // currentStrain[7];
  double microStep = 0;
  double stepAccuracyCheck;

  for (int i = 0; i < 7; i++) {
    substepStrain[i] = epStrain[i];
    absStress[i] = 0;
  }
  int stepNo = 0;
  vector<double> stressInc;
  vector<double> strainInc;
  vector<double> p0Lambda;
  vector<double>::iterator Iter;

  clock_t startTime, endTime;
  startTime = clock();

  // assuming that maximum tolerable step is 0.5% for order 2 Initial step size
  // is being calculated:
  newStepSize = 0;
  for (int i = 0; i < 6; i++)
    newStepSize = newStepSize + std::abs(substepStrain[i]);
  newStepSize = 0.01 / (newStepSize * Order);
  // dbg << "newStepSize="<<newStepSize<<"\n";
  if (newStepSize > 1)
    newStepSize = 1;
  // getchar(); */

  do {
    stepAccepted = false;

    /*if (stepNo>1)
{
	if (newStepSize>10) newStepSize=10;
	if (newStepSize<0.1) newStepSize=0.1;
}*/ // limiting step increase/decrease does not improve results/ enables faster
    // convergence...

    stepLength = stepLength * newStepSize; // size of a step
    if ((stepLength + totalSize) > 1)
      stepLength =
        1 - totalSize; // check whether the step not exceed the whole increment

    for (int i = 0; i < 7; i++) {
      substepStrain[i] =
        stepLength * epStrain[i]; // strain increment in current step
                                  // currentStrain[i]=0;
    }
    rError = 0;

    Point->Copy(&trialState);
    Point->Copy(&oldState);
    doRungeKuttaEqualStep(A, B, BRes, C, &trialState, substepStrain, absStress,
                        &rError, 2, Order, Steps, errorEstimate);

    stepNo = stepNo + 2;
    // dbg << "Step Length="<<stepLength<<"\n";
    // dbg << "Current strain [0]="<<substepStrain[0]<<"\n";
    // dbg << "rError="<<rError<<"\n";
    // getchar();

    if (rError < d_integrationTol) {
      stepAccepted = true;
    } else {
      stepAccepted = false;
      if (stepLength < 1e-20)
        stepAccepted = true;
    }

    if (d_tolMethod == 0)
      newStepSize =
        d_betaFactor * pow(d_integrationTol / rError, (1 / (Order - 1.0)));
    else
      newStepSize =
        d_betaFactor * pow(d_integrationTol / rError, (1 / Order));

    if (!stepAccepted) {
      // check the extrapolation and go into it, if it looks sensible. otherwise
      // reject the step...
      if (newStepSize < WORTH_EXTRAPOL) {
        // double Result[7];
        int TempNumber = 0;
        // for (int i=0; i<6; i++)
        // Result[i]=trialState.stress[i]-oldState.stress[i];
        // Result[6]=trialState.p0Star()-oldState.p0Star();
        trialState.Copy(&oldState);
        Point->Copy(&trialState);
        doRKExtrapolation(
          A, B, BRes, C, &trialState, substepStrain, absStress, &oldState,
          &rError, &TempNumber, Order, Steps,
          errorEstimate); // Extrapolate and finally accept the step
        stepNo = stepNo + TempNumber;
        stepAccepted = true;
        Point->Copy(&oldState);
      } else
        ;
      /*	//here we should take care about correct re - usage of the first evaluation of derivative
	reuseStep=true;
	for (int i=0; i<6; i++) reuseRes[i]=dSigma[0][i]*newStepSize;
	reuseRes[6]=dP0Star[0]*newStepSize; */ // No step re-using so far
      // reject the step.
    }

    if (stepAccepted) {
      // here we need to update all the point data. and nothing else...
      for (int i = 0; i < 6; i++)
        absStress[i] =
          absStress[i] + std::abs(oldState.stress[i] - trialState.stress[i]);
      absStress[6] =
        absStress[6] + std::abs(Point->p0Star() - trialState.p0Star());
      trialState.Copy(Point);
      // Point->Update(0,substepStrain,stressInc,oldState.p0Star()-trialState.p0Star());
      // drift is already corrected
      // reuseStep=false;
      stepAccuracyCheck = totalSize;
      microStep = microStep + stepLength;
      totalSize = totalSize + microStep; // total part of step done updated
      microStep = microStep - (totalSize - stepAccuracyCheck);
      if (totalSize >= 1)
        finished = true;
      Temp = double(stepNo) / frequency;
      if (modf(Temp, &Temp) == 0) {
        cout << "Step number:" << stepNo << "\n";
        cout << "Total size done is: " << totalSize
             << " of whole step. Current stepLength=" << stepLength << "\n";
      }
      /*
              //Debug
              strainInc.push_back (stepLength);
              for (int i=0; i<6; i++) stressInc.push_back (Point->stress[i]);
              stressInc.push_back (Point->p0Star());
              //End Debug */
    }

  } while (!finished);
  endTime = clock();
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
  return (endTime - startTime);
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
  // StateMohrCoulomb midState[8], oldState, trialState;

  int Steps = 8;
  double Order = 5.0;
  bool errorEstimate = false;

  // double time;
  // time=doRungeKutta
  // (A,B,BRes,C,Point,epStrain,absStress,numberIter,Order,Steps,false);
  // return time;

  BBMMatrix dSigma(6, 1);
  StateMohrCoulomb midState[8], oldState, trialState;

  double dSigma[8][6], stressInc[7], Result[7];
  double dP0Star[8], p0StarInc, plasticStrainInc[6],
    plasticStrain[8][6];
  double Error[7], ErrorOther[7], reuseRes[7], rError, rErrorOther, newStepSize,
    totalSize, stepLength, Temp, methodPower;
  double frequency = 15000 / Order; // how often display info about steps
  double dStressDrift[6], dLambdaDrift;

  bool reuseStep = false;

  stepLength = 1;
  totalSize = 0;
  newStepSize = 1;
  methodPower = pow(2.0, Order) * d_integrationTol;

  for (int i = 0; i < Order; i++) {
    for (int j = 0; j < 6; j++)
      dSigma[i][j] = 0;
    dP0Star[i] = 0;
  }

  bool finished = false, stepAccepted = false;
  double substepStrain[7], currentStrain[7];
  double microStep = 0;
  double stepAccuracyCheck;

  for (int i = 0; i < 7; i++) {
    substepStrain[i] = epStrain[i];
    absStress[i] = 0;
  }
  int stepNo = 0;
  vector<double> stressInc;
  vector<double> strainInc;
  vector<double> p0Lambda;
  vector<double>::iterator Iter;

  clock_t startTime, endTime;
  startTime = clock();

  newStepSize = 0;
  for (int i = 0; i < 6; i++)
    newStepSize = newStepSize + std::abs(substepStrain[i]);
  newStepSize = 0.01 / (newStepSize * Order);
  // dbg << "newStepSize="<<newStepSize<<"\n";
  if (newStepSize > 1)
    newStepSize = 1;

  do {
    stepAccepted = false;
    stepNo++;

    /*if (stepNo>1)
{
	if (newStepSize>10) newStepSize=10;
	if (newStepSize<0.1) newStepSize=0.1;
}*/ // limiting step increase/decrease does not improve results/ enables faster
    // convergence...

    stepLength = stepLength * newStepSize; // size of a step
    if ((stepLength + totalSize) > 1)
      stepLength =
        1 - totalSize; // check whether the step not exceed the whole increment

    for (int i = 0; i < 7; i++) {
      substepStrain[i] =
        stepLength * epStrain[i]; // strain increment in current step
      currentStrain[i] = 0;
    }
    rError = 0;

    // dbg << "Step Length="<<stepLength<<"\n";
    // dbg << "Current strain [0]="<<currentStrain[0]<<"\n";

    for (int i = 0; i < Steps; i++)
      point->Copy(&midState[i]); // point is unchanged in  procedure

    // Below the main R-K loop to calculate the value of intermediate stresses;
    // values stored in stressInc[][]
    // reuseStep=false;
    if (reuseStep) {
      for (int i = 0; i < 6; i++) {
        dSigma[0][i] = reuseRes[i];
        plasticStrain[0][i] = plasticStrainInc[i];
      }
      dP0Star[0] = reuseRes[6];
      // add line about plastic strain...
    } else {
      calcPlastic(midState[0], substepStrain, &dSigma, plasticStrainInc,
                  &dP0Star[0], computeYieldFunctionNN(Point), dStressDrift,
                  &dLambdaDrift);
      for (int i = 0; i < 6; i++) {
        dSigma[0][i] = dSigma.getElement(i + 1, 1);
        plasticStrain[0][i] = plasticStrainInc[i];
      }
    }

    for (int rkloop = 1; rkloop < Steps; rkloop++) {
      for (int i = 0; i < 6; i++) {
        stressInc[i] = 0;
        plasticStrainInc[i] = 0;
      }
      p0StarInc = 0;
      for (int i = 0; i < 7; i++)
        currentStrain[i] =
          C[rkloop] *
          substepStrain[i]; // set the beginning point of the procedure
      for (int i = 0; i < rkloop; i++) {
        for (int j = 0; j < 6; j++) {
          stressInc[j] = stressInc[j] + A[rkloop][i] * dSigma[i][j];
          plasticStrainInc[j] =
            plasticStrainInc[j] + A[rkloop][i] * plasticStrain[i][j];
        }
        p0StarInc = p0StarInc + A[rkloop][i] * dP0Star[i];
      }
      midState[rkloop].Update(plasticStrainInc, currentStrain, stressInc,
                              p0StarInc);
      // double dummy;
      // StateMohrCoulomb TempPoint;
      // BBMMatrix TEMPMATRIX (6,7), TEMPEPSILON(7,1);
      // substepStrain[6]=0;
      // for (int i=1; i<8; i++) TEMPEPSILON.PutElement(i,1,substepStrain[i-1]);
      // midState[rkloop].Copy(&TempPoint);
      calcPlastic(midState[rkloop], substepStrain, &dSigma, plasticStrainInc,
                  &dP0Star[rkloop], computeYieldFunctionNN(Point),
                  dStressDrift, &dLambdaDrift);
      // calculateElastoPlasticTangentMatrix(&TempPoint,&TEMPMATRIX);
      // TEMPMATRIX.Multiply(&TEMPEPSILON,&TEMPMATRIX);
      /*	for (int i=1; i<7; i++)
              {
                      dummy=std::abs(dSigma.getElement(i,1)-TEMPMATRIX.getElement(i,1));
                      if (std::abs(dSigma.getElement(i,1))>TINY)
         dummy=dummy/std::abs(dSigma.getElement(i,1));
                      if (dummy>0.01)
                      {
                              dbg << "Problems with the alternative
         matrix."<<i<<" Change of stress is:"<<dSigma.getElement(i,1)<<"\n";
                              dbg << "calculated change of stress with the DEP
         matrix is:"<<TEMPMATRIX.getElement(i,1)<<"\n";
                              getchar();
                      }
                      else dbg << "*";
              }


      */
      for (int i = 0; i < 6; i++) {
        dSigma[rkloop][i] = dSigma.getElement(i + 1, 1);
        plasticStrain[rkloop][i] = plasticStrainInc[i];
      }
    }

    // needed: result, error

    for (int i = 0; i < 6; i++) {
      Result[i] = 0;
      plasticStrainInc[i] = 0;
      for (int j = 0; j < Steps; j++) {
        Result[i] = Result[i] + BRes[j] * dSigma[j][i];
        plasticStrainInc[i] =
          plasticStrainInc[i] + BRes[j] * plasticStrain[j][i];
      }
    }
    Result[6] = 0;
    for (int j = 0; j < Steps; j++)
      Result[6] = Result[6] + BRes[j] * dP0Star[j];

    for (int i = 0; i < 7; i++)
      Error[i] = 0;

    for (int i = 0; i < Steps; i++) {
      for (int j = 0; j < 6; j++)
        Error[j] = Error[j] + B[i] * dSigma[i][j];
      Error[6] = Error[6] + B[i] * dP0Star[i];
    }
    if (!errorEstimate)
      for (int i = 0; i < 7; i++)
        Error[i] = Error[i] - Result[i]; // error estimate calculated in case we
                                         // have lower order solution instead of
                                         // error estimate

    for (int i = 0; i < 7; i++)
      ErrorOther[i] = 0;
    for (int i = 0; i < (Steps - 1); i++) {
      for (int j = 0; j < 6; j++)
        ErrorOther[j] = ErrorOther[j] + ErrorCoef[i] * dSigma[i][j];
      ErrorOther[6] = ErrorOther[6] + ErrorCoef[i] * dP0Star[i];
    }

    // check the error norm

    switch (int(d_tolMethod)) {
      case 0: {
        // rError=checkNorm (stressInc, p0StarInc, point, Error);
        // //returns rError
        rError = checkNorm(Result, Result[6], point, Error); // returns rError
        rErrorOther = checkNorm(Result, Result[6], point, ErrorOther);
        if (rError < rErrorOther)
          rError = rErrorOther;
      } break;
      case 1: {
        // SLOAN NORM
        // rError=checkNormSloan (stressInc, p0StarInc, point, Error);
        // //returns rError
        rError = checkNormSloan(Result, Result[6], point, Error); // returns
                                                                  // rError
        rErrorOther = checkNormSloan(Result, Result[6], point,
                                     ErrorOther); // returns rError
        if (rError < rErrorOther)
          rError = rErrorOther;
      } break;
      default: {
        cout << "ERROR !!!! Improper d_tolMethod in increment.dta" << "\n";
        getchar();
      }
    }

    for (int i = 0; i < 7; i++)
      if (!finite(Result[i])) {
        Result[i] = 0;
        if (rError < methodPower)
          rError = methodPower;
      }

    if ((Point->getmeanStress() + (Result[0] + Result[1] + Result[2])) / 3 < 0)
      if (rError < methodPower)
        rError = methodPower;
    if ((Point->p0Star() + Result[6]) < 0)
      if (rError < methodPower)
        rError = methodPower;

    if (rError < d_integrationTol) {
      stepAccepted = true;
    } else {
      stepAccepted = false;
      if (stepLength < 1e-20)
        stepAccepted = true;
    }

    if (d_tolMethod == 0)
      newStepSize =
        d_betaFactor * pow(d_integrationTol / rError, (1 / (Order - 1.0)));
    else
      newStepSize =
        d_betaFactor * pow(d_integrationTol / rError, (1 / Order));

    if (!stepAccepted) {
      // here we should take care about correct re - usage of the first
      // evaluation of derivative
      reuseStep = true;
      for (int i = 0; i < 6; i++) {

        reuseRes[i] = dSigma[0][i] * newStepSize;
        plasticStrainInc[i] = plasticStrain[0][i] * newStepSize;
      }
      reuseRes[6] = dP0Star[0] * newStepSize;
    } else {
      // here we need to update all the point data.
      Point->Copy(&oldState);
      Point->Update(plasticStrainInc, substepStrain, Result, Result[6]);
      Point->Copy(&trialState);
      // this value is not used in any calculations, it is just to show at the
      // end of the step the p0* increase

      if (d_driftCorrection == 3)
        correctDrift(Point);
      if (d_driftCorrection == 2)
        correctDriftBeg(
          Point,
          &oldState); // value of oldState copied before updating the point

      // re - evaluate the error in the point:

      for (int i = 0; i < 6; i++) {
        Error[i] = Error[i] + Point->stress[i] - trialState.stress[i];
        ErrorOther[i] = ErrorOther[i] + Point->stress[i] - trialState.stress[i];
      }
      Error[6] = Error[6] + Point->p0Star() - trialState.p0Star();
      ErrorOther[6] = ErrorOther[6] + Point->p0Star() - trialState.p0Star();
      // error vector updated, norm should be re-evaluated:

      switch (int(d_tolMethod)) {
        case 0: {
          // Temp=checkNorm (stressInc, p0StarInc, point, Error);
          // //returns rError
          Temp = checkNorm(Result, Result[6], point, Error); // returns rError
          rErrorOther =
            checkNorm(Result, Result[6], point, ErrorOther); // returns rError
          if (Temp > rError)
            rError = Temp;
          if (rErrorOther > rError)
            rError = rErrorOther;
        } break;
        case 1: {
          // SLOAN NORM
          // Temp=checkNormSloan (stressInc, p0StarInc, point, Error);
          // //returns rError
          Temp = checkNormSloan(Result, Result[6], point, Error); // returns
                                                                  // rError
          rErrorOther = checkNormSloan(Result, Result[6], point,
                                       ErrorOther); // returns rError
          if (Temp > rError)
            rError = Temp;
          if (rErrorOther > rError)
            rError = rErrorOther;
        } break;
        default: {
          cout << "ERROR !!!! Improper d_tolMethod in increment.dta" << "\n";
          getchar();
        }
      }
      if (!finite(rError))
        if (rError < methodPower)
          rError = methodPower;

      if (d_tolMethod == 0)
        newStepSize =
          d_betaFactor * pow(d_integrationTol / rError, (1 / (Order - 1.0)));
      else
        newStepSize =
          d_betaFactor * pow(d_integrationTol / rError, (1 / Order));

      if (rError < d_integrationTol)
        stepAccepted = true;
      else {
        stepAccepted = false;
        if (stepLength < 1e-20)
          stepAccepted = true;
      }

      if (!stepAccepted) {
        reuseStep = true;
        for (int i = 0; i < 6; i++) {
          reuseRes[i] = dSigma[0][i] * newStepSize;
          plasticStrainInc[i] = plasticStrain[0][i] * newStepSize;
        }
        reuseRes[6] = dP0Star[0] * newStepSize;
        oldState.Copy(Point);
      } else {
        // this may be done only after successful re-evaluation of the substep
        for (int i = 0; i < 6; i++)
          absStress[i] =
            absStress[i] + std::abs(Point->stress[i] - oldState.stress[i]);
        absStress[6] =
          absStress[6] + std::abs(Point->p0Star() - oldState.p0Star());
        reuseStep = false;
        stepAccuracyCheck = totalSize;
        microStep = microStep + stepLength;
        totalSize = totalSize + microStep; // total part of step done updated
        microStep = microStep - (totalSize - stepAccuracyCheck);
        if (totalSize >= 1)
          finished = true;
        Temp = double(stepNo) / frequency;
        if (modf(Temp, &Temp) == 0) {
          cout << "Step number:" << stepNo << "\n";
          cout << "Total size done is: " << totalSize
               << " of whole step. Current stepLength=" << stepLength << "\n";
        }
        /*
                //Debug
                strainInc.push_back (stepLength);
                for (int i=0; i<6; i++) stressInc.push_back (Point->stress[i]);
                stressInc.push_back (Point->p0Star());
                //End Debug */
      }
    }

  } while (!finished);
  endTime = clock();
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

  return (endTime - startTime);
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
  int Steps = 6;

  double Order = 5;
  double time;
  double A[8][8]; // matrix must be this size to be used in the doRungeKutta
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

  time = doRungeKutta(A, B, BRes, C, point, epStrain, absStress, numberIter,
                    Order, Steps, errorEstimate);
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
  int Steps = 7;
  double Order = 5;
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

  // time=doRungeKuttaExtrapol (A, B,BRes, C, point, epStrain, absStress,
  // numberIter, Order, Steps, errorEstimate);
  time = doRungeKutta(A, B, BRes, C, point, epStrain, absStress, numberIter,
                    Order, Steps, errorEstimate);
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
  int Steps = 6;

  double Order = 5;
  double time;
  double A[8][8]; // matrix must be this size to be used in the doRungeKutta
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

  time = doRungeKutta(A, B, BRes, C, point, epStrain, absStress, numberIter,
                    Order, Steps, errorEstimate);
  // time=doRungeKuttaExtrapol (A, B,BRes, C, point, epStrain, absStress,
  // numberIter, Order, Steps, errorEstimate);
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
  int Steps = 6;

  double Order = 5;
  double time;
  double A[8][8]; // matrix must be this size to be used in the doRungeKutta
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

  time = doRungeKutta(A, B, BRes, C, point, epStrain, absStress, numberIter,
                    Order, Steps, errorEstimate);
  // time=doRungeKuttaExtrapol (A, B,BRes, C, point, epStrain, absStress,
  // numberIter, Order, Steps, errorEstimate);
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
  int Steps = 5;

  double Order = 4;
  double time;
  double A[8][8]; // matrix must be this size to be used in the doRungeKutta
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

  time = doRungeKutta(A, B, BRes, C, point, epStrain, absStress, numberIter,
                    Order, Steps, errorEstimate);
  // time=doRungeKuttaExtrapol (A, B,BRes, C, point, epStrain, absStress,
  // numberIter, Order, Steps, errorEstimate);
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
  int Steps = 3;

  double Order = 3;
  double time;
  double A[8][8]; // matrix must be this size to be used in the doRungeKutta
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

  time = doRungeKutta(A, B, BRes, C, point, epStrain, absStress, numberIter,
                    Order, Steps, errorEstimate);
  // time=doRungeKuttaExtrapol (A, B,BRes, C, point, epStrain, absStress,
  // numberIter, Order, Steps, errorEstimate);
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
  int Steps = 4;

  double Order = 3;
  double time;
  double A[8][8]; // matrix must be this size to be used in the doRungeKutta
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

  time = doRungeKutta(A, B, BRes, C, point, epStrain, absStress, numberIter,
                    Order, Steps, errorEstimate);
  // time=doRungeKuttaExtrapol (A, B,BRes, C, point, epStrain, absStress,
  // numberIter, Order, Steps, errorEstimate);
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
  // clock_t startTime, endTime;
  // startTime=clock();

  StateMohrCoulomb newState, midState;
  BBMMatrix SIGMA(6, 1);
  // vector <double> Result;
  vector<double>::iterator Iter;

  double dSigma[6], currentStrain[7], HalfcurrentStrain[7], plasticStrain[6];
  double h, dP0Star = 0;
  double dStressDrift[6], dLambdaDrift;

  h = *numberIter;
  h = 1 / h;
  for (int i = 0; i < 7; i++)
    absStress[i] = 0; // 7 components...
  for (int i = 0; i < 7; i++) {
    currentStrain[i] = epStrain[i] * h; // strain increment in current step
    HalfcurrentStrain[i] = 0.5 * currentStrain[i];
  }
  // dbg << "Step Length="<<stepLength<<"\n";
  // dbg << "Current strain [0]="<<currentStrain[0]<<"\n";
  /*
  for (int i=0; i<6; i++) Result.push_back(Point->stress[i]);
  Result.push_back(Point->p0Star());
  Result.push_back(Point->getmeanStress());
  Result.push_back(Point->shearStress());
  */

  Point->Copy(&newState); // point is unchanged in calcPlastic procedure
  Point->Copy(&midState); // point is unchanged in calcPlastic procedure
  calcPlastic(newState, HalfcurrentStrain, &SIGMA, plasticStrain, &dP0Star,
              computeYieldFunctionNN(Point), dStressDrift,
              &dLambdaDrift); // calculate the plastic stresses...
  for (int i = 0; i < 6; i++)
    dSigma[i] = SIGMA.getElement(i + 1, 1);
  newState.Update(plasticStrain, HalfcurrentStrain, dSigma, dP0Star);
  calcPlastic(newState, HalfcurrentStrain, &SIGMA, plasticStrain, &dP0Star,
              computeYieldFunctionNN(Point), dStressDrift,
              &dLambdaDrift); // calculate the plastic stresses...

  for (int loop = 0; loop < 2 * (*numberIter); loop++) {

    midState.Update(plasticStrain, HalfcurrentStrain, dSigma, dP0Star);
    midState.Copy(&newState);
    newState.Update(plasticStrain, HalfcurrentStrain, dSigma, dP0Star);
    calcPlastic(newState, HalfcurrentStrain, &SIGMA, plasticStrain, &dP0Star,
                computeYieldFunctionNN(Point), dStressDrift,
                &dLambdaDrift); // calculate the plastic stresses...
    for (int i = 0; i < 6; i++)
      dSigma[i] = SIGMA.getElement(i + 1, 1);

    /*for (int i=0; i<6; i++) Result.push_back(Point->stress[i]);
    Result.push_back(Point->p0Star());
    Result.push_back(Point->getmeanStress());
    Result.push_back(Point->shearStress());*/
  }
  midState.Copy(Point);
  // endTime=clock();
  // return (endTime-startTime);
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
  // clock_t startTime, endTime;
  // startTime=clock();

  StateMohrCoulomb newState, midState;
  BBMMatrix SIGMA(6, 1);
  // vector <double> Result;
  vector<double>::iterator Iter;

  double dSigma[6], currentStrain[7], HalfcurrentStrain[7], plasticStrain[6];
  double h, dP0Star = 0;
  double dStressDrift[6], dLambdaDrift;

  h = *numberIter;
  h = 1 / h;
  for (int i = 0; i < 7; i++)
    absStress[i] = 0; // 7 components...
  for (int i = 0; i < 7; i++) {
    currentStrain[i] = epStrain[i] * h; // strain increment in current step
    HalfcurrentStrain[i] = 0.5 * currentStrain[i];
  }
  // dbg << "Step Length="<<stepLength<<"\n";
  // dbg << "Current strain [0]="<<currentStrain[0]<<"\n";
  /*
  for (int i=0; i<6; i++) Result.push_back(Point->stress[i]);
  Result.push_back(Point->p0Star());
  Result.push_back(Point->getmeanStress());
  Result.push_back(Point->shearStress());
  */
  for (int loop = 0; loop < *numberIter; loop++) {
    Point->Copy(&newState); // point is unchanged in calcPlastic procedure
    Point->Copy(&midState);
    calcPlastic(newState, HalfcurrentStrain, &SIGMA, plasticStrain, &dP0Star,
                computeYieldFunctionNN(Point), dStressDrift,
                &dLambdaDrift); // calculate the plastic stresses...
    for (int i = 0; i < 6; i++)
      dSigma[i] = SIGMA.getElement(i + 1, 1);
    midState.Update(plasticStrain, HalfcurrentStrain, dSigma, dP0Star);
    calcPlastic(midState, HalfcurrentStrain, &SIGMA, plasticStrain, &dP0Star,
                computeYieldFunctionNN(Point), dStressDrift,
                &dLambdaDrift); // calculate the plastic stresses...
    for (int i = 0; i < 6; i++) {
      dSigma[i] = 2 * SIGMA.getElement(i + 1, 1);
      plasticStrain[i] = 2 * plasticStrain[i];
    }
    dP0Star = 2 * dP0Star;

    for (int i = 0; i < 6; i++)
      absStress[i] = absStress[i] + std::abs(dSigma[i]);
    absStress[6] = absStress[6] + std::abs(dP0Star);
    Point->Update(plasticStrain, currentStrain, dSigma, dP0Star);
    /*for (int i=0; i<6; i++) Result.push_back(Point->stress[i]);
    Result.push_back(Point->p0Star());
    Result.push_back(Point->getmeanStress());
    Result.push_back(Point->shearStress());*/
  }
  // endTime=clock();
  // return (endTime-startTime);
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
                  dbg << "Procedure proceed with standard algorithm and no drift
  correction."<<"\n";
                  break;
          }
  case 2 :
          {
                  dbg << "Procedure proceed with standard algorithm and drift
  correction at the beginning (point A)."<<"\n";
                  break;
          }

  default :
          {
                  dbg << "Unknown d_driftCorrection parameter. Parameter read
  is:"<<d_driftCorrection<<"\n";
                  break;
          }
  }
  */

  // no drift correction as may worsen the results. at least possibly.

  clock_t startTime, endTime;
  startTime = clock();

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
  double Zero[6] = { 0, 0, 0, 0, 0, 0 };

  double Divisions[15];
  for (int i = 0; i < 15; i++)
    Divisions[i] = DivisionsInt[i];

  double approximationTable[16][20]; // six stresses, p0*, six absolute
                                     // stresses, absolute p0*, 6 plastic
                                     // strains
  double approximationTableOld[16][20];
  double hSquareTable[16][15];
  double dError[15], rError;
  double InitialStress[7], InitialplasticStrain[6], plasticStrain[6];
  double dSigma[6];
  /*
  double TestTable[5]={0.375,0.37109375,0.36945588,0.36879683,0.36829712};
  double TestApproxTable[5]; double TestApproxTableOld[5];
  for (int i=0; i<5; i++)
  {
          TestApproxTable[i]=0;
          TestApproxTableOld[i]=0;
  }
  */

  StateMohrCoulomb newState, CopyPoint, oldState;
  Point->Copy(&CopyPoint);
  Point->Copy(&oldState);

  for (int i = 0; i < 16; i++) {
    for (int j = 0; j < 20; j++) {
      approximationTable[i][j] = 0;
      approximationTableOld[i][j] = 0;
    }
    for (int j = 0; j < 15; j++)
      hSquareTable[i][j] = 0;
  }

  for (int i = 1; i < 16; i++) {
    for (int j = 0; j < i; j++) {
      hSquareTable[i][j] =
        (Divisions[i] / Divisions[j]) * (Divisions[i] / Divisions[j]);
      // dbg << "Divisions["<<i<<"]="<<Divisions[i]<<"
      // Divisions["<<j<<"]="<<Divisions[j]<<"\n";
      // dbg << "hSquareTable["<<i<<"]["<<j<<"]="<<hSquareTable[i][j]<<"\n";
    }
  }
  for (int i = 0; i < 6; i++) {
    InitialStress[i] = Point->stress[i];
    InitialplasticStrain[i] = Point->plastic_strain[i];
  }
  InitialStress[6] = Point->p0Star();

  int loop = 0;
  for (; loop < STEPMAX; loop++) {
    // calculating stress increment using midState rule

    CopyPoint.Copy(&newState);
    plasticMidpoint(&newState, epStrain, absStress,
                    &DivisionsInt[loop]); // calculation of stress increment
                                          // using the midpoint procedure
    // TestApproxTable[0]=TestTable[loop];

    // saving data using midState rule
    for (int i = 0; i < 6; i++) {
      approximationTable[0][i] = newState.stress[i] - InitialStress[i];
      approximationTable[0][i + 7] = absStress[i];
      approximationTable[0][i + 14] =
        newState.plastic_strain[i] - InitialplasticStrain[i];
    }
    approximationTable[0][6] = newState.p0Star() - InitialStress[6];
    approximationTable[0][13] = absStress[6];

    // all data from the Midpoint Rule saved in the ResultsTable.

    for (int i = 0; i < loop; i++) {

      for (int j = 0; j < 20; j++) {
        approximationTable[i + 1][j] =
          approximationTable[i][j] +
          (approximationTable[i][j] - approximationTableOld[i][j]) /
            (hSquareTable[loop][loop - i - 1] - 1);
      }
    }

    // for (i=0; i<loop+1; i++) dbg << approximationTable[i][0]<<"  ";
    // dbg << "\n"<<"\n";

    // approximations are calculated
    // two possibilities of error control. In literature rather the second one;
    // this is the more stringent, but first one should be enough
    //(1) use approximationTable[loop][i] - approximationTable[loop-1][i]
    //(2) use approximationTable[loop][i] - approximationTableOld [loop-1][i]
    // still using relative error definition...
    // oldState (used for calculating the norm in q) is set accordingly

    // FIRST ONE OK, SEE Deuflhard Bornemann 2002 Springer p.206

    rError = 0;
    if (loop > 0) {
      for (int i = 0; i < 6; i++)
        dSigma[i] = approximationTable[loop][i];
      for (int i = 0; i < 7; i++) {
        // dError[i]=approximationTable[loop][i] -
        // approximationTableOld[loop-1][i];
        dError[i] =
          approximationTable[loop][i] - approximationTable[loop - 1][i];
      }

      switch (int(d_tolMethod)) {
        case 0: {
          rError = checkNorm(dSigma, approximationTable[loop][6], &CopyPoint,
                             dError); // returns rError
        } break;
        case 1: {
          // SLOAN NORM
          // dbg << "SLOAN NORM"<<"\n";
          rError =
            checkNormSloan(dSigma, approximationTable[loop][6], &CopyPoint,
                           dError); // returns rError
        } break;
        default: {
          cout << "ERROR !!!! Improper d_tolMethod in increment.dta" << "\n";
          getchar();
        }
      } // end switch
    } else
      rError = 2 * d_integrationTol;

    // dbg << "Relative error after iteration "<<loop<<" is equal
    // to:"<<rError<<"\n";
    for (int i = 0; i < loop + 1; i++)
      for (int j = 0; j < 20; j++)
        approximationTableOld[i][j] = approximationTable[i][j];
    // for (int i=0; i<loop+1; i++) dbg << approximationTable[i][0]<<"  ";
    // dbg << "\n";
    if (rError < d_integrationTol) {
      // error less then requested...check for drift and add the error
      CopyPoint.Copy(&oldState);
      oldState.Update(Zero, epStrain, dSigma, approximationTable[loop][6]);
      oldState.Copy(&newState);
      if (d_driftCorrection == 3)
        correctDrift(&newState);
      if (d_driftCorrection == 2)
        correctDriftBeg(
          &newState,
          &CopyPoint); // value of oldState copied before updating the point
      for (int i = 0; i < 6; i++)
        dError[i] = dError[i] + (newState.stress[i] - oldState.stress[i]);
      dError[6] = dError[6] + newState.p0Star() - oldState.p0Star();
      switch (int(d_tolMethod)) {
        case 0: {
          rError = checkNorm(dSigma, approximationTable[loop][6], &CopyPoint,
                             dError); // returns rError
        } break;
        case 1: {
          // SLOAN NORM
          // dbg << "SLOAN NORM"<<"\n";
          rError =
            checkNormSloan(dSigma, approximationTable[loop][6], &CopyPoint,
                           dError); // returns rError
        } break;
        default: {
          cout << "ERROR !!!! Improper d_tolMethod in increment.dta" << "\n";
          getchar();
        }
      } // end switch
      if (rError < d_integrationTol) {
        loop = loop + 100; // no more interations after - finished
      } else {
        // one more loop, do not need to change anything, as point anyway copied
        // back later on.
      }
    } else { // newState.Copy(&oldState);
    }
  }

  if (loop > 100) {
    loop = loop - 101;
    // done - time to update everything and sent time and no of iterations...

    for (int i = 0; i < 6; i++) {
      absStress[i] = approximationTable[loop][i + 7];
      plasticStrain[i] = approximationTable[loop][i + 14];
      // dbg << "dSigma="<<dSigma[i]<<"   absStress="<<absStress[i]<<"\n";
    }
    absStress[6] = approximationTable[loop][13];
    // dbg << absStress[6]<<"\n";
    Point->Update(plasticStrain, epStrain, dSigma, approximationTable[loop][6]);
    endTime = clock();
    *numberIter = 0;
    for (int i = 0; i < loop + 1; i++)
      *numberIter = *numberIter + Divisions[i];
    // dbg << "Procedure has coverged after"<<*numberIter<<" iterations."<<"\n";
    // getchar();
  } else {
    loop--;
    double dSigma[6];
    for (int i = 0; i < 6; i++) {
      dSigma[i] = approximationTable[loop][i];
      absStress[i] = approximationTable[loop][i + 7];
      plasticStrain[i] = approximationTable[loop][i + 14];
      // dbg << dSigma[i]<<"\n";
    }
    absStress[6] = approximationTable[loop][13];
    Point->Update(plasticStrain, epStrain, dSigma, approximationTable[loop][6]);
    endTime = clock();
    *numberIter = 0;
    for (int i = 0; i < loop + 1; i++)
      *numberIter = *numberIter + Divisions[i];
    // dbg << "Procedure has NOT CONVERGED after"<<*numberIter<<"
    // iterations."<<"\n";
    // getchar();
  }
  // dbg << "calculation took:"<<double(endTime-startTime)/CLOCKS_PER_SEC<<"
  // s."<<"\n";

  /*
  switch (int(d_driftCorrection))
  {case 1 :{dbg << "Values calculated with standard algorithm and no drift
  correction."<<"\n";break;}
  case 2 :{dbg << "Values calculated with standard algorithm and drift correction
  at the beginning (point A)."<<"\n";break;}
  case 3 :{dbg << "Values calculated with standard algorithm and drift correction
  at the end (point B)."<<"\n";break;}
  case 4 :{dbg << "Values calculated with zero drift algorithm."<<"\n";break;}
  default :{dbg << "Unknown d_driftCorrection parameter. Parameter read
  is:"<<d_driftCorrection<<"\n";break;}
  }

  dbg << "The d_integrationTol parameter is equal to:"<<d_integrationTol<<"\n";
  dbg << "Total number of steps done:"<<stepNo<<"\n";
  dbg << "Total parameter lambda change: "<<LambdaTot<<"\n";
  dbg << "Total PZeroStar change: "<<PZeroStarTot<<"\n";
  checkYield (Point->state, point->stress, point->strain[6],&fValue);
  dbg << "Yield Function value is:"<<fValue<<"\n";
  dbg << "Number of drift correction performed:"<<correctdriftCount<<"\n";
  dbg << "Over the whole step change of stress is:"<<"\n";
  for (int i=0; i<6; i++)
  {
          dbg << "s["<<i<<"]="<<Point->stress[i]-oldState.stress[i]<<"\n";
  }

  //double DEpsV=0, dNu, SpecificVolume;

  //SpecificVolume=oldState.getSpecVol();
  //LambdaS=LambdaZero*(r+(1-r)*exp(-Beta*oldState.suction()));



  dbg << "Initial specific volume="<<oldState.getSpecVol()<<"\n";
  dbg << "Final specific volume="<<Point->getSpecVol()<<"\n";
  dbg << "Change of specific volume is equal
  to="<<(oldState.getSpecVol()-Point->getSpecVol())<<"\n";
  dbg << "Change of mean stress p is equal
  to:"<<Point->getmeanStress()-oldState.getmeanStress()<<"\n";
  dbg << "Change of shear stress q is equal
  to:"<<Point->shearStress()-oldState.shearStress()<<"\n";



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
  for (Iter=p0Lambda.begin() ; Iter!=p0Lambda.end(); )
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

  // dbg << "Press any key..."<<"\n";
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
  ",approximationTable[loop][i]);
                  }
  fprintf( File, "%.20f ,
  ",(approximationTable[loop][0]+approximationTable[loop][1]+approximationTable[loop][2])/3);
  fprintf( File, "%.20f ,
  ",approximationTable[loop][0]-approximationTable[loop][1]);

  fprintf(File, "\n");
  fclose (File); */
  return (endTime - startTime);
}

double
ShengMohrCoulomb::doRKExtrapolation(double A[][8], double* B, double* BRes,
                                  double* C, StateMohrCoulomb* point, double* epStrain,
                                  double* absStress, StateMohrCoulomb* oldState,
                                  double* RelError, int* numberIter,
                                  double Order, int Steps,
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
                  dbg << "Procedure proceed with standard algorithm and no drift
  correction."<<"\n";
                  break;
          }
  case 2 :
          {
                  dbg << "Procedure proceed with standard algorithm and drift
  correction at the beginning (point A)."<<"\n";
                  break;
          }

  default :
          {
                  dbg << "Unknown d_driftCorrection parameter. Parameter read
  is:"<<d_driftCorrection<<"\n";
                  break;
          }
  }
  */

  // no drift correction as may worsen the results. at least possibly.

  clock_t startTime, endTime;
  startTime = clock();

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

  double approximationTable[16][20]; // 6 stresses, p0*, 6 absolute stresses,
                                     // absolute p0*, 6 plastic strains
  double approximationTableOld[16][20];
  double hSquareTable[16][15];
  double dError[15], rError;
  double InitialStress[7], InitialplasticStrain[6], plasticStrain[6];
  double dSigma[6];
  bool stepAccepted = false;
  /*
  double TestTable[5]={0.375,0.37109375,0.36945588,0.36879683,0.36829712};
  double TestApproxTable[5]; double TestApproxTableOld[5];
  for (int i=0; i<5; i++)
  {
          TestApproxTable[i]=0;
          TestApproxTableOld[i]=0;
  }
  */

  StateMohrCoulomb newState, CopyPoint;
  Point->Copy(&CopyPoint);
  // Point->Copy(&oldState);

  for (int i = 0; i < 16; i++) {
    for (int j = 0; j < 20; j++) {
      approximationTable[i][j] = 0;
      approximationTableOld[i][j] = 0;
    }
    for (int j = 0; j < 15; j++)
      hSquareTable[i][j] = 0;
  }

  for (int i = 1; i < 16; i++) {
    for (int j = 0; j < i; j++) {
      hSquareTable[i][j] =
        (Divisions[i] / Divisions[j]) * (Divisions[i] / Divisions[j]);
      // dbg << "Divisions["<<i<<"]="<<Divisions[i]<<"
      // Divisions["<<j<<"]="<<Divisions[j]<<"\n";
      // dbg << "hSquareTable["<<i<<"]["<<j<<"]="<<hSquareTable[i][j]<<"\n";
    }
  }
  for (int i = 0; i < 6; i++) {
    InitialStress[i] = Point->stress[i];
    InitialplasticStrain[i] = Point->plastic_strain[i];
  }
  InitialStress[6] = Point->p0Star();

  for (int i = 0; i < 6; i++) {
    approximationTable[0][i] = oldState->stress[i] - InitialStress[i];
    approximationTable[0][i + 7] = absStress[i];
    approximationTable[0][i + 14] =
      oldState->plastic_strain[i] - InitialplasticStrain[i];
  }
  approximationTable[0][6] = oldState->p0Star() - InitialStress[6];
  approximationTable[0][13] = absStress[6];

  int loop = 1;

  // for (int i=0; i<loop+1; i++) dbg << approximationTable[i][0]<<"  ";
  // dbg << "\n";

  for (int i = 0; i < loop; i++) {

    for (int j = 0; j < 20; j++) {
      approximationTable[i + 1][j] =
        approximationTable[i][j] +
        (approximationTable[i][j] - approximationTableOld[i][j]) /
          (hSquareTable[loop][loop - i - 1] - 1);
    }
  }
  for (int i = 0; i < loop + 1; i++)
    for (int j = 0; j < 20; j++)
      approximationTableOld[i][j] = approximationTable[i][j];

  // dbg << "Initial"<<"\n";
  // for (int i=0; i<loop; i++) dbg << approximationTable[i][0]<<"  ";
  // dbg << "\n"<<"\n";
  // loop=0;

  for (; loop < STEPMAX; loop++) {
    // calculating stress increment using midState rule

    CopyPoint.Copy(&newState);
    doRungeKuttaEqualStep(A, B, BRes, C, &newState, epStrain, absStress, &rError,
                        DivisionsInt[loop], Order, Steps,
                        errorEstimate); // calculation of stress increment using
                                        // the RK procedure
    // TestApproxTable[0]=TestTable[loop];

    // saving data using midPoint rule
    for (int i = 0; i < 6; i++) {
      approximationTable[0][i] = newState.stress[i] - InitialStress[i];
      approximationTable[0][i + 7] = absStress[i];
      approximationTable[0][i + 14] =
        newState.plastic_strain[i] - InitialplasticStrain[i];
    }
    approximationTable[0][6] = newState.p0Star() - InitialStress[6];
    approximationTable[0][13] = absStress[6];

    // all data from the Midpoint Rule saved in the ResultsTable.

    for (int i = 0; i < loop; i++) {

      for (int j = 0; j < 20; j++) {
        approximationTable[i + 1][j] =
          approximationTable[i][j] +
          (approximationTable[i][j] - approximationTableOld[i][j]) /
            (hSquareTable[loop][loop - i - 1] - 1);
      }
    }

    // for (i=0; i<loop+1; i++) dbg << approximationTable[i][0]<<"  ";
    // dbg << "\n";
    // getchar();

    // approximations are calculated
    // two possibilities of error control. In literature rather the second one;
    // this is the more stringent, but first one should be enough
    //(1) use approximationTable[loop][i] - approximationTable[loop-1][i]
    //(2) use approximationTable[loop][i] - approximationTableOld [loop-1][i]
    // still using relative error definition...
    // oldState (used for calculating the norm in q) is set accordingly

    stepAccepted = false;
    if (rError < d_integrationTol)
      stepAccepted = true;

    for (int i = 0; i < 6; i++)
      dSigma[i] = approximationTable[loop][i];
    for (int i = 0; i < 7; i++) {
      dError[i] =
        approximationTable[loop][i] - approximationTableOld[loop - 1][i];
    }

    switch (int(d_tolMethod)) {
      case 0: {
        rError = checkNorm(dSigma, approximationTable[loop][6], oldState,
                           dError); // returns rError
      } break;
      case 1: {
        // SLOAN NORM
        // dbg << "SLOAN NORM"<<"\n";
        rError = checkNormSloan(dSigma, approximationTable[loop][6], oldState,
                                dError); // returns rError
      } break;
      default: {
        cout << "ERROR !!!! Improper d_tolMethod in increment.dta" << "\n";
        getchar();
      }
    } // end switch

    // dbg << "Relative error after iteration "<<loop<<" is equal
    // to:"<<rError<<"\n";
    for (int i = 0; i < loop + 1; i++)
      for (int j = 0; j < 20; j++)
        approximationTableOld[i][j] = approximationTable[i][j];
    // for (int i=0; i<loop+1; i++) dbg << approximationTable[i][0]<<"  ";
    // dbg << "\n";
    if (rError < d_integrationTol)
      stepAccepted = true;
    if (stepAccepted)
      loop = loop + 100; // no more interations after - finished
    else
      newState.Copy(oldState);
  }

  if (loop > 100) {
    loop = loop - 101;
    // done - time to update everything and sent time and no of iterations...

    for (int i = 0; i < 6; i++) {
      absStress[i] = approximationTable[loop][i + 7];
      plasticStrain[i] = approximationTable[loop][i + 14];
      // dbg << "dSigma="<<dSigma[i]<<"   absStress="<<absStress[i]<<"\n";
    }
    absStress[6] = approximationTable[loop][13];
    // dbg << absStress[6]<<"\n";
    Point->Update(plasticStrain, epStrain, dSigma, approximationTable[loop][6]);
    endTime = clock();
    *numberIter = 0;
    for (int i = 0; i < loop + 1; i++)
      *numberIter = *numberIter + Divisions[i];
    // dbg << "Procedure has coverged after"<<*numberIter<<" iterations."<<"\n";
    // getchar();
    return (endTime - startTime);
  } else {
    loop--;
    double dSigma[6];
    for (int i = 0; i < 6; i++) {
      dSigma[i] = approximationTable[loop][i];
      absStress[i] = approximationTable[loop][i + 7];
      plasticStrain[i] = approximationTable[loop][i + 14];
      // dbg << dSigma[i]<<"\n";
    }
    absStress[6] = approximationTable[loop][13];
    Point->Update(plasticStrain, epStrain, dSigma, approximationTable[loop][6]);
    endTime = clock();
    *numberIter = 0;
    for (int i = 0; i < loop + 1; i++)
      *numberIter = *numberIter + Divisions[i];
    cout << "Procedure has NOT CONVERGED after" << *numberIter << " iterations."
         << "\n";
    // getchar();
  }
  return (endTime - startTime);
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
