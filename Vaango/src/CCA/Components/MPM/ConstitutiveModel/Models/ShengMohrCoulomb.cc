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

#include <CCA/Components/MPM/ConstitutiveModel/Models/ShengMohrCoulomb.h>

#include <Core/Exceptions/InvalidValue.h>
#include <Core/Util/DebugStream.h>

#include <chrono>
#include <cmath>
#include <iomanip>

using namespace Uintah;

static DebugStream dbg("ShengMC", false);
static DebugStream dbg_unloading("ShengMC_unloading", false);

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
  d_int.setDefaults(d_yield);
}

ShengMohrCoulomb::ShengMohrCoulomb(double G, double K, double cohesion,
                                   double phi, double psi)
{
  setModelParameters(G, K, cohesion, phi, psi);
  d_int.setDefaults(d_yield);
}

/**
 * Initialization methods
 */
void
ShengMohrCoulomb::setModelParameters(double G, double K, double cohesion,
                                     double phi, double psi)
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
ShengMohrCoulomb::setIntegrationParameters(
  int maxIterPegasus, double integrationTolerance, double betaFactor,
  double yieldLocTolerance, SolutionAlgorithm solutionAlgorithm,
  ToleranceMethod toleranceMethod, DriftCorrection driftCorrection)
{
  d_int.d_maxIter = maxIterPegasus;
  d_int.d_integrationTol = integrationTolerance;
  d_int.d_yieldTol = yieldLocTolerance;
  d_int.d_betaFactor = betaFactor;

  d_int.d_solutionAlgorithm = solutionAlgorithm;
  d_int.d_tolMethod = toleranceMethod;
  d_int.d_driftCorrection = driftCorrection;
}

/**
 * Integration
 */
StateMohrCoulomb
ShengMohrCoulomb::integrate(const Vector7& strainIncrement,
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

  // checkYieldNormalized returns true  if in the plastic part of the YL. If on
  // the YL
  // (within Tolerance INT_TOL returns false)
  bool elasticPlastic = checkYieldNormalized(finalState);

  if (dbg.active()) {
    dbg << "p = " << finalState.meanStress()
        << " q = " << finalState.shearStress() << "  Yield? " << std::boolalpha
        << elasticPlastic << "\n";
  }

  bool unloading = false;
  if (elasticPlastic) {
    if (onTheYieldLocus) {

      if (dbg_unloading.active()) {
        // checking if there is any unloading part; in finalState purely elastic
        // values;
        // we care only for correct change of suction
        unloading = checkGradient(initialState, finalState);
      }

      if (unloading) {

        if (dbg_unloading.active()) {
          dbg_unloading << "\n\n Elasto-Plastic unloading=" << unloading
                        << "\n";
        }

        // find elastic part
        findIntersectionUnloading(strainIncrement, initialState,
                                  purelyElasticStrain, purelyPlasticStrain);

        // calculating elastic part, updated in the same point.
        calcElastic(purelyElasticStrain, initialState, finalState);

        // Update specific volume
        double spVol =
          computeNu(finalState.stress, finalState.state, finalState.suction());
        finalState.specificVolume(spVol);

        calculatePlastic(purelyPlasticStrain, finalState);

      } else {

        calculatePlastic(strainIncrement, finalState);
      }

    } else {
      // not on the yield locus, finding intersection and subsequently
      // calculate elastic and plastic stress increment
      findIntersection(strainIncrement, initialState, purelyElasticStrain,
                       purelyPlasticStrain);

      calcElastic(purelyElasticStrain, initialState, finalState);

      double spVol =
        computeNu(finalState.stress, finalState.state, finalState.suction());
      finalState.specificVolume(spVol);

      calculatePlastic(purelyPlasticStrain, finalState);
    }
  }

  if (dbg.active()) {
    dbg << "strain = " << strainIncrement.transpose() << "\n"
        << "stress = " << finalState.stress << "\n";
  }

  double spVol =
    computeNu(finalState.stress, finalState.state, finalState.suction());
  finalState.specificVolume(spVol);

  return finalState;
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
    double factor =
      -27.0 * J3 / (2.0 * shearStress * shearStress * shearStress);
    if (!std::isfinite(factor)) {
      factor = 1.0;
      std::cout
        << "**WARNING** Factor in checkYield Normalised not finite. set to 1\n"
           "  alpha4 = "
        << d_yield.d_alpha4 << " J3 = " << J3
        << " shearStress = " << shearStress << "\n";
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
  double yieldFnValue = shearStress / M - 2.0 * cohesion / M - meanStress;

  // normalisation
  yieldFnValue /= (std::abs(meanStress) + 2.0 * cohesion);

  if (dbg.active()) {
    std::cout << "Check Yield: Mean Stress = " << meanStress
              << " Shear Stress=" << shearStress << " cohesion = " << cohesion
              << " M = " << M << " Yield Function = " << yieldFnValue << "\n";
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

  Vector6 stressInc =
    calcStressIncElast(initialState.specificVolume(), initialState.stress,
                       initialState.strain, strainInc);

  std::cout << "Stress increment is: " << stressInc << "\n";

  Vector6 plasticStrainInc = Vector6::Zero();
  double p0StarInc = 0.0;
  finalState.update(plasticStrainInc, strainInc, stressInc, p0StarInc);
}

Vector6
ShengMohrCoulomb::calcStressIncElast(double nu0, const Vector6& s0,
                                     const Vector7& eps0,
                                     const Vector7& deps) const
{
  double K = d_elastic.d_K;
  double G = d_elastic.d_G;

  double K43G = K + 4.0 * G / 3.0;
  double K23G = K - 2.0 * G / 3.0;

  Vector6 ds;
  ds(0) = K43G * deps(0) + K23G * (deps(1) + deps(2));
  ds(1) = K43G * deps(1) + K23G * (deps(0) + deps(2));
  ds(2) = K43G * deps(2) + K23G * (deps(0) + deps(1));
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
  Vector7 strainInc = finalState.strain - initialState.strain;

  double max = 0.0;
  for (int i = 0; i < 6; i++) {
    if (std::abs(strainInc(i)) > max) {
      max = std::abs(strainInc(i));
    }
  }
  for (int i = 0; i < 6; i++) {
    strainInc(i) /= (max * 10E-10);
  }
  strainInc(6) = finalState.suction() - initialState.suction();

  // The normalisation is important to catch the unloading due to shear stress
  // 12 13 23. As the q is always positive, large unloading - loading values
  // of 12 13 or 23 component lead to larger q, the shear stress change is
  // taken as positive and there is no unloading. This fix, though not the
  // most elegant, should work.
  Matrix67 tangentElastic = calculateElasticTangentMatrix(initialState);

  Vector6 stressInc = tangentElastic * strainInc;

  if (dbg.active()) {
    dbg << "Strain inc = " << strainInc << "\n";
    dbg << "TangentElastic: " << tangentElastic << "\n";
    dbg << "Stress inc = " << stressInc << "\n";
    dbg << "Suction Increment = " << strainInc(7) << "\n";
  }

  Vector6 df; // not initialised; will contain values of derivatives of the
              // yield locus function
  double cosinus = findGradient(initialState.stress, stressInc, df, 0, 0);

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
ShengMohrCoulomb::calculateElasticTangentMatrix(
  const StateMohrCoulomb& state) const
{
  double K = d_elastic.d_K;
  double G = d_elastic.d_G;
  Matrix66 DEP_66 = calculateElasticTangentMatrix(K, G);

  Matrix67 DEP = Matrix67::Zero();
  DEP.block<6, 6>(0, 0) = DEP_66;
  return DEP;
}

Matrix66
ShengMohrCoulomb::calculateElasticTangentMatrix(double K, double G) const
{
  double K43G = K + 4.0 * G / 3.0;
  double K23G = K - 2.0 * G / 3.0;

  Matrix66 DEP = Matrix66::Zero();
  DEP(0, 0) = K43G;
  DEP(0, 1) = K23G;
  DEP(0, 2) = K23G;
  DEP(1, 0) = K23G;
  DEP(1, 1) = K43G;
  DEP(1, 2) = K23G;
  DEP(2, 0) = K23G;
  DEP(2, 1) = K23G;
  DEP(2, 2) = K43G;
  DEP(3, 3) = 2.0 * G;
  DEP(4, 4) = 2.0 * G;
  DEP(5, 5) = 2.0 * G;
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
                               Vector6& df_dSigma, double /*suction*/,
                               double /*dsuction*/) const
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
                 "  alpha4 = "
              << d_yield.d_alpha4 << " J3 = " << J3
              << " shearStress = " << shearStress << "\n";
  } else if (factor > 1) {
    factor = 1;
  } else if (factor < -1) {
    factor = -1;
  }
  factor = 1 + d_yield.d_alpha4 - (1 - d_yield.d_alpha4) * factor;

  double factor025 = std::pow(factor, 0.25);
  double factor075 = std::pow(factor, -0.75);
  double M =
    (3 - d_yield.d_sin_phi) / (6 * d_yield.d_alpha * d_yield.d_sin_phi);
  double Mpsi =
    (3 - d_potential.d_sin_psi) / (6 * d_yield.d_alpha * d_potential.d_sin_psi);

  Vector6 dJ2_dSig = Vector6::Zero();
  dJ2_dSig(0) = (2 * stress(0) - stress(1) - stress(2)) / 3.0;
  dJ2_dSig(1) = (2 * stress(1) - stress(0) - stress(2)) / 3.0;
  dJ2_dSig(2) = (2 * stress(2) - stress(0) - stress(1)) / 3.0;
  dJ2_dSig(3) = 2 * stress(3);
  dJ2_dSig(4) = 2 * stress(4);
  dJ2_dSig(5) = 2 * stress(5);

  Vector6 df_dSigma = dJ2_dSig;
  Vector6 dg_dSigma = dJ2_dSig;
  df_dSigma *= (0.21022410391343 * M * factor025 / shearStress);
  dg_dSigma *= (0.21022410391343 * Mpsi * factor025 / shearStress);

  for (int i = 0; i < 6; i++) {
    if (!std::isfinite(df_dSigma(i))) {
      std::ostringstream err;
      err << "*ERROR** df/dSigma not finite at df_dsigma(" << i << ")."
          << df_dSigma << "\n";
      throw InvalidValue(err.str(), __FILE__, __LINE__);
    }
  }

  if (dbg.active()) {
    dbg << "dF/dJ2 part\n";
    dbg << std::setprecision(15) << df_dSigma << "\n";
  }

  Vector6 dq_dSig = dJ2_dSig;
  dq_dSig *= std::sqrt(3.0) * 0.5 / std::sqrt(J2);

  Vector6 dq_dSig_g = dq_dSig;
  dq_dSig *=
    (-2.838025401481287 * (d_yield.d_alpha4 - 1) * d_yield.d_cohesion * M * J3 *
       factor075 / (shearStress * shearStress * shearStress * shearStress) +
     4.257038102221929 * (d_yield.d_alpha4 - 1) * M * std::sqrt(J2) * J3 *
       factor075 / (std::sqrt(3.0) * shearStress * shearStress * shearStress *
                    shearStress));
  dq_dSig_g *=
    (-2.838025401481287 * (d_yield.d_alpha4 - 1) * d_yield.d_cohesion * Mpsi *
       J3 * factor075 /
       (shearStress * shearStress * shearStress * shearStress) +
     4.257038102221929 * (d_yield.d_alpha4 - 1) * Mpsi * std::sqrt(J2) * J3 *
       factor075 / (std::sqrt(3.0) * shearStress * shearStress * shearStress *
                    shearStress));

  if (dbg.active()) {
    dbg << "dF/dq part\n";
    dbg << std::setprecision(15) << dq_dSig << "\n";
  }

  df_dSigma += dq_dSig;
  dg_dSigma += dq_dSig_g;

  Vector6 dI1_dSig = Vector6::Zero();
  dI1_dSig(0) = 1.0;
  dI1_dSig(1) = 1.0;
  dI1_dSig(2) = 1.0; //{1,1,1,0,0,0}

  Vector6 dJ3_dSig = dI1_dSig;
  dJ3_dSig *= (2.0 / 9.0 * I1 * I1 - I2 / 3.0);

  if (dbg.active()) {
    dbg << std::setprecision(15) << dJ3_dSig << "\n";
  }

  Vector6 dI2_dSig = Vector6::Zero();
  dI2_dSig(0) = stress(1) + stress(2);
  dI2_dSig(1) = stress(0) + stress(2);
  dI2_dSig(2) = stress(0) + stress(1);
  dI2_dSig(3) = -2 * stress(3);
  dI2_dSig(4) = -2 * stress(4);
  dI2_dSig(5) = -2 * stress(5);

  if (dbg.active()) {
    dbg << std::setprecision(15) << dI2_dSig << "\n";
  }

  Vector6 dI3_dSig = Vector6::Zero();
  dI3_dSig(0) = stress(1) * stress(2) - stress(5) * stress(5);
  dI3_dSig(1) = stress(0) * stress(2) - stress(4) * stress(4);
  dI3_dSig(2) = stress(0) * stress(1) - stress(3) * stress(3);
  dI3_dSig(3) = 2 * stress(5) * stress(4) - 2 * stress(2) * stress(3);
  dI3_dSig(4) = 2 * stress(3) * stress(5) - 2 * stress(1) * stress(4);
  dI3_dSig(5) = 2 * stress(3) * stress(4) - 2 * stress(0) * stress(5);

  if (dbg.active()) {
    dbg << std::setprecision(15) << dI3_dSig << "\n";
  }

  dJ3_dSig += (dI2_dSig * (-I1 / 3.0) + dI3_dSig);

  if (dbg.active()) {
    dbg << std::setprecision(15) << dJ3_dSig << "\n";
  }

  Vector6 dJ3_dSig_g = dJ3_dSig;
  dJ3_dSig *=
    (-1.419012700740643 * (d_yield.d_alpha4 - 1) * M * sqrt(J2) * factor075 /
       (std::sqrt(3.0) * shearStress * shearStress * shearStress) +
     0.94600846716043 * (d_yield.d_alpha4 - 1) * d_yield.d_cohesion * M *
       factor075 / (shearStress * shearStress * shearStress));
  dJ3_dSig_g *=
    (-1.419012700740643 * (d_yield.d_alpha4 - 1) * Mpsi * sqrt(J2) * factor075 /
       (std::sqrt(3.0) * shearStress * shearStress * shearStress) +
     0.94600846716043 * (d_yield.d_alpha4 - 1) * d_yield.d_cohesion * Mpsi *
       factor075 / (shearStress * shearStress * shearStress));

  if (dbg.active()) {
    dbg << "dF/dJ3 part\n";
    dbg << std::setprecision(15) << dJ3_dSig << "\n";
  }

  df_dSigma += dJ3_dSig;
  dg_dSigma += dJ3_dSig_g;

  for (int i = 0; i < 6; i++) {
    if (!std::isfinite(df_dSigma(i))) {
      std::ostringstream err;
      err << "*ERROR** df/dSigma not finite at df_dsigma(" << i << ")."
          << df_dSigma << "\n";
      throw InvalidValue(err.str(), __FILE__, __LINE__);
    }
  }

  df_dSigma += dI1_dSig * (-1.0 / 3.0);
  dg_dSigma += dI1_dSig * (-1.0 / 3.0);

  for (int i = 0; i < 6; i++) {
    if (!std::isfinite(df_dSigma(i))) {
      std::ostringstream err;
      err << "*ERROR** df/dSigma not finite at df_dsigma(" << i << ")."
          << df_dSigma << "\n";
      throw InvalidValue(err.str(), __FILE__, __LINE__);
    }
  }

  // tension cut-off plane
  double meanStress = I1 / 3.0;
  if (meanStress < -d_yield.d_cohesion * M) {
    df_dSigma(0) = 1.0 / 3.0;
    dg_dSigma(0) = 1.0 / 3.0;
    df_dSigma(1) = 1.0 / 3.0;
    dg_dSigma(1) = 1.0 / 3.0;
    df_dSigma(2) = 1.0 / 3.0;
    dg_dSigma(2) = 1.0 / 3.0;
    df_dSigma(3) = 0.0;
    dg_dSigma(3) = 0.0;
    df_dSigma(4) = 0.0;
    dg_dSigma(4) = 0.0;
    df_dSigma(5) = 0.0;
    dg_dSigma(5) = 0.0;
  }

  if (dbg.active()) {
    dbg << "df_dSigma = \n"
        << df_dSigma.transpose() << "dg_dsigma = \n"
        << dg_dSigma.transpose() << "\n";
  }

  return std::make_tuple(df_dSigma, dg_dSigma);
}

/**
 Find intersection with yield surface during unloading
 */
void
ShengMohrCoulomb::findIntersectionUnloading(
  const Vector7& strainIncrement, const StateMohrCoulomb& initialState,
  Vector7& elasticStrainInc, Vector7& plasticStrainInc)
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
                                   Vector7& plasticStrainInc) const
{
  double alpha = findYieldModified(initialState.state, initialState.stress,
                                   initialState.strain, strainIncrement);
  elasticStrainInc = strainIncrement * alpha;
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
ShengMohrCoulomb::findYieldAlpha(const Vector3& state, const Vector6& s0,
                                 const Vector7& eps0, const Vector7& dep)
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
    delta = (f0 + d_int.d_yieldTol) / 2;
    d_int.d_yieldTol = 0.9 * (d_int.d_yieldTol - f0) / 2;
  }

  if (dbg.active()) {
    dbg << "Procedure Find Yield Yield ... s = " << s0 << "\n";
    dbg << "Procedure Find Yield Yield ... dEps = " << dep << "\n";
    dbg << "Initial normalised Yield Locus... f0 =" << f0 << "\n";
    dbg << "Procedure Find Yield Yield ... delta =" << delta << "\n";
    dbg
      << "Call of this procedure is generally rare and most likely improper\n";
    dbg << "Parameters are set to:\n";
    dbg << "Maximum number of iteration d_maxIter: " << d_int.d_maxIter << "\n";
    ;
    dbg << "Value of d_alfaCheck - when the additional iter. is to "
           "enter: "
        << d_int.d_alfaCheck << "\n";
    dbg << "Value of change of alfa in the step: " << d_int.d_alfaChange
        << "\n";
    ;
    dbg << "alfa old/alfa ratio: " << d_int.d_alfaRatio << "\n";
  }

  double alfa0 = 0;
  double alfa1 = 1;
  // Vector7 epsIni = dep * alfa0;
  Vector7 epsFin = dep * alfa1;

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
    dbg << "f0 should lie on the yield locus within tolerance: " << yieldTol_old
        << "\n";
  }

  double alphaYield = 0.0;
  double alfa = 0.5;
  double alfa_old = 0;
  for (int iter = 0; iter < d_int.d_maxIter; iter++) {
    bool problems = false;
    alfa_old = alfa;

    alfa = f0 / (f0 - f1);
    Vector7 epsAlfa = dep * (alfa0 + alfa * (alfa1 - alfa0));

    // calculate stress increment for current alfa
    Vector6 sAlfa = calcStressIncElast(state(2), s0, eps0, epsAlfa);

    // update stress
    sAlfa += s0;

    // calculate yield function for current alfa
    double fAlfa = computeYieldNormalized(sAlfa);
    fAlfa -= delta;

    if (dbg.active()) {
      dbg << "In iteration " << iter << " alfa = " << alfa
          << " and f = " << fAlfa << "\n";
    }

    // if fAlfa within tolerance, we have the solution
    if (std::abs(fAlfa) < d_int.d_yieldTol) {

      alphaYield = alfa0 + alfa * (alfa1 - alfa0);

      if (dbg.active()) {
        dbg << "Solution in findYieldAlpha procedure was found after " << iter
            << " iterations."
            << "\n";
        dbg << "Value of the yield function is equal to: " << fAlfa << "\n";
        dbg << "Modified value of tolerance is: " << d_int.d_yieldTol << "\n";
        dbg << "Value of delta is: " << delta << "\n";
        dbg << "Alfa is equal to: " << alphaYield << "\n";
      }

      if (iter > 50) {
        std::cout << "**WARNING** Large number of iterations!!! Solution is "
                     "however correct... "
                  << "\n";
      }

      d_int.d_yieldTol = yieldTol_old;

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
          dbg << "Problematic iteration entered !!!"
              << "\n";
        }

        alfa = alfa1 - d_int.d_alfaChange * (alfa1 - alfa0);
        epsAlfa = dep * (alfa0 + alfa * (alfa1 - alfa0));

        sAlfa = calcStressIncElast(state(2), s0, eps0, epsAlfa);
        sAlfa += s0;

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

      if ((alfa < d_int.d_alfaCheck) && (alfa_old / alfa) < d_int.d_alfaRatio) {
        problems = true;
      }
      alfa0 = alfa0 + alfa * (alfa1 - alfa0);
      f0 = fAlfa;

      if (problems) {
        if (dbg.active()) {
          dbg << "Problematic iteration entered !!!"
              << "\n";
        }

        alfa = alfa0 + d_int.d_alfaChange * (alfa1 - alfa0);
        epsAlfa = dep * (alfa0 + alfa * (alfa1 - alfa0));

        sAlfa = calcStressIncElast(state(2), s0, eps0, epsAlfa);
        sAlfa += s0;

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

  d_int.d_yieldTol = yieldTol_old;

  // if we are here, we must have perforemed to many iterations...
  std::ostringstream err;
  err << "**ERROR** in procedure findYieldAlpha"
      << "\n";
  err << "  After " << d_int.d_maxIter << " iterations crossing point not found"
      << "\n";
  err << "  This is likely to cause incorrect results... Results obtained "
         "  should not be taken too seriously..."
      << "\n";
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
ShengMohrCoulomb::findYieldModified(const Vector3& state, const Vector6& s0,
                                    const Vector7& eps0,
                                    const Vector7& deps) const
{

  if (dbg.active()) {
    dbg << "Parameters are set to:\n";
    dbg << "Maximum number of iteration d_maxIter: " << d_int.d_maxIter << "\n";
    ;
    dbg << "Value of d_alfaCheck - when the additional iter. is to "
           "enter: "
        << d_int.d_alfaCheck << "\n";
    dbg << "Value of change of alfa in the step: " << d_int.d_alfaChange
        << "\n";
    ;
    dbg << "alfa old/alfa ratio: " << d_int.d_alfaRatio << "\n";
  }

  double alfa0 = 0;
  double alfa1 = 1;
  // Vector7 epsIni = deps * alfa0;
  Vector7 epsFin = deps * alfa1;

  Vector6 sIni = s0;
  Vector6 sFin = calcStressIncElast(state(2), s0, eps0, epsFin);
  sFin += s0;

  double f0 = computeYieldNormalized(s0);
  double f1 = computeYieldNormalized(sFin);

  if (dbg.active()) {
    dbg << "Procedure findYieldModified. Initial values of f0 = " << f0
        << " f1 = " << f1 << "\n";
    dbg << "Value of f0 should be negative, and value of f1 should be "
           "positive.\n";
    dbg << "Values should be larger than tolerance for yield: "
        << d_int.d_yieldTol << "\n";
  }

  double alfa = 0.5;
  double alfa_old = 0;

  for (int iter = 0; iter < d_int.d_maxIter; iter++) {

    bool problems = false;
    alfa_old = alfa;
    alfa = f0 / (f0 - f1);

    Vector7 epsAlfa = deps * (alfa0 + alfa * (alfa1 - alfa0));
    Vector6 sAlfa = calcStressIncElast(state(2), s0, eps0, epsAlfa);
    sAlfa += s0;

    double fAlfa = computeYieldNormalized(sAlfa);

    if (dbg.active()) {
      dbg << "In iteration " << iter << " alfa = " << alfa
          << " and f = " << fAlfa;
      dbg << " alfa0 = " << alfa0 << "\n";
    }

    // if the difference is below numerical the accuracy, we have the solution
    if ((alfa1 - alfa0) < TINY) {
      fAlfa = 0.0;
    }

    if (std::abs(fAlfa) < d_int.d_yieldTol) {

      // if fAlfa within tolerance, we have the solution
      double alphaYield = alfa0 + alfa * (alfa1 - alfa0);

      if (iter > 50) {
        std::cout << "**WARNING** Large number of iterations in "
                     "findYieldModified procedure!!! "
                     "Solution is however correct..."
                  << "\n";
      }
      if (dbg.active()) {
        dbg << "Solution in findYieldModified procedure was found after "
            << iter << " iterations."
            << "\n";
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
          dbg << "Problematic iteration entered !!!"
              << "\n";
        }
        alfa = alfa1 - d_int.d_alfaChange * (alfa1 - alfa0);
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
          dbg << "Problematic iteration entered !!!"
              << "\n";
        }

        alfa = alfa0 + d_int.d_alfaChange * (alfa1 - alfa0);
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
  err << "Error in procedure findYieldModified"
      << "\n";
  err << "After " << d_int.d_maxIter << " iterations crossing point not found"
      << "\n";
  err << "alphamin = " << alfa0 << " alphamax = " << alfa1
      << " dalpha = " << alfa1 - alfa0;
  err << "Yield Function value Min=" << f0 << " Max=" << f1 << "\n";
  err << "Stress:" << s0 << "\n";
  err << "Strain:" << deps << "\n";
  err << "G: " << d_elastic.d_G << " K: " << d_elastic.d_K
      << " cohesion: " << d_yield.d_cohesion << " phi: " << d_yield.d_phi
      << "\n";
  err << "This is likely to cause incorrect results... Results obtained "
         "should not be taken too seriously..."
      << "\n";
  throw InvalidValue(err.str(), __FILE__, __LINE__);

  return -999.0;
}

double
ShengMohrCoulomb::computeNu(const Vector6& s, const Vector3& state,
                            double suction) const
{
  // does nothing for SMC
  return 1.0;
}

/**
  Calculate stress increment during plastic deformation
 */
double
ShengMohrCoulomb::calculatePlastic(const Vector7& purelyPlasticStrain,
                                   StateMohrCoulomb& state) const
{
  double time;
  int numIter;

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
      err << "**ERROR** Unknown Solution Algorithm. Value of "
             "d_solutionAlgorithm variable "
             "is set to:"
          << d_int.d_solutionAlgorithm << "\n";
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
                                 const Vector7& epStrain) const

{
  constexpr int Order = 2;
  constexpr int Steps = 2;

  using RKMatrix = Eigen::Matrix<double, Steps, Steps>;
  using RKVector = Eigen::Matrix<double, Steps, 1>;

  RKMatrix A = RKMatrix::Zero();
  RKVector C = RKVector::Zero();
  RKVector BRes = RKVector::Zero();
  RKVector B = RKVector::Zero();

  getParamRKME221(A, B, BRes, C);

  bool errorEstimate = false;
  auto result =
    doRungeKutta<Order, Steps>(A, B, BRes, C, state, epStrain, errorEstimate);
  return result;
}

void
ShengMohrCoulomb::getParamRKME221(Eigen::Matrix<double, 2, 2>& A,
                                  Eigen::Matrix<double, 2, 1>& B,
                                  Eigen::Matrix<double, 2, 1>& BRes,
                                  Eigen::Matrix<double, 2, 1>& C) const
{
  A(0, 0) = 0.0;
  A(1, 0) = 1.0;

  C(0) = 0.0;
  C(1) = 1.0;

  BRes(0) = 0.5;
  BRes(1) = 0.5;

  B(0) = 1.0;
  B(1) = 0.0;
}

/**
  This procedure calculates stress increment using 3-2 Runge - Kutta pair. The
  coeficients are given in  Ordinary and Partial Differential Equations
  routines in... by Lee & Schiesser, Chapman & Hall 2003
  The procedure consists of 3 stages and is 3th order accurate.
  A - matrix with coefficients, B - error estimate, BRes - result
  coefficients, C - x coefficients.
 */
std::tuple<double, int>
ShengMohrCoulomb::plasticRK332(StateMohrCoulomb& state,
                               const Vector7& epStrain) const

{
  constexpr int Steps = 3;
  constexpr int Order = 3;

  using RKMatrix = Eigen::Matrix<double, Steps, Steps>;
  using RKVector = Eigen::Matrix<double, Steps, 1>;

  RKMatrix A = RKMatrix::Zero();
  RKVector C = RKVector::Zero();
  RKVector BRes = RKVector::Zero();
  RKVector B = RKVector::Zero();

  getParamRK332(A, B, BRes, C);

  bool errorEstimate = false;
  auto result =
    doRungeKutta<Order, Steps>(A, B, BRes, C, state, epStrain, errorEstimate);
  return result;
}

void
ShengMohrCoulomb::getParamRK332(Eigen::Matrix<double, 3, 3>& A,
                                Eigen::Matrix<double, 3, 1>& B,
                                Eigen::Matrix<double, 3, 1>& BRes,
                                Eigen::Matrix<double, 3, 1>& C) const
{
  A(0, 0) = 0;
  A(1, 0) = 2.0 / 3.0;
  A(2, 0) = 0;
  A(2, 1) = 2.0 / 3.0;

  C << 0.0, 2.0 / 3.0, 2.0 / 3.0;

  BRes << 0.25, 0.375, 0.375;

  B << 0.25, 0.75, 0;
}

/**
  This procedure calculate stress increment using 3-2 Runge - Kutta pair. The
  coeficients are given in "A 3(2) Pair of Runge-Kutta Formulas"
  by P. Bogacki and L.F.  Shampine, Appl. Math. Lett., 2, pp. 321-325, 1989.
  The procedure consists of 3 stages and is 3th order accurate.

  A - matrix with coefficients, B - error estimate, BRes - result
  coefficients, C - x coefficients.
 */
std::tuple<double, int>
ShengMohrCoulomb::plasticRKBog432(StateMohrCoulomb& state,
                                  const Vector7& epStrain) const
{
  constexpr int Steps = 4;
  constexpr int Order = 3;

  using RKMatrix = Eigen::Matrix<double, Steps, Steps>;
  using RKVector = Eigen::Matrix<double, Steps, 1>;

  RKMatrix A = RKMatrix::Zero();
  RKVector C = RKVector::Zero();
  RKVector BRes = RKVector::Zero();
  RKVector B = RKVector::Zero();

  getParamRKBog432(A, B, BRes, C);

  bool errorEstimate = false;
  auto result =
    doRungeKutta<Order, Steps>(A, B, BRes, C, state, epStrain, errorEstimate);
  return result;
}

void
ShengMohrCoulomb::getParamRKBog432(Eigen::Matrix<double, 4, 4>& A,
                                   Eigen::Matrix<double, 4, 1>& B,
                                   Eigen::Matrix<double, 4, 1>& BRes,
                                   Eigen::Matrix<double, 4, 1>& C) const
{
  A(0, 0) = 0;
  A(1, 0) = 0.5;
  A(2, 0) = 0;
  A(2, 1) = 0.75;
  A(3, 0) = 2.0 / 9.0;
  A(3, 1) = 1.0 / 3.0;
  A(3, 2) = 4.0 / 9.0;

  C << 0.0, 0.5, 0.75, 1.0;

  BRes << 2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0;

  B << 7.0 / 24.0, 0.25, 1 / 3.0, 0.125;
}

/**
 This procedure calculate stress increment using 4-3 Runge - Kutta pair. The
 coeficients are given in  Ordinary and Partial Differential Equations routines
 in... by Lee & Schiesser, Chapman & Hall 2003
 The procedure consists of 5 stages and is 4th order accurate.
*/
std::tuple<double, int>
ShengMohrCoulomb::plasticRK543(StateMohrCoulomb& state,
                               const Vector7& epStrain) const
{
  constexpr int Steps = 5;
  constexpr int Order = 4;

  using RKMatrix = Eigen::Matrix<double, Steps, Steps>;
  using RKVector = Eigen::Matrix<double, Steps, 1>;

  RKMatrix A = RKMatrix::Zero();
  RKVector C = RKVector::Zero();
  RKVector BRes = RKVector::Zero();
  RKVector B = RKVector::Zero();

  getParamRK543(A, B, BRes, C);

  bool errorEstimate = false;
  auto result =
    doRungeKutta<Order, Steps>(A, B, BRes, C, state, epStrain, errorEstimate);
  return result;
}

void
ShengMohrCoulomb::getParamRK543(Eigen::Matrix<double, 5, 5>& A,
                                Eigen::Matrix<double, 5, 1>& B,
                                Eigen::Matrix<double, 5, 1>& BRes,
                                Eigen::Matrix<double, 5, 1>& C) const
{
  A(0, 0) = 0;
  A(1, 0) = 1.0 / 3.0;
  A(2, 0) = 0.5 / 3.0;
  A(2, 1) = 0.5 / 3.0;
  A(3, 0) = 0.125;
  A(3, 1) = 0;
  A(3, 2) = 0.375;
  A(4, 0) = 0.5;
  A(4, 1) = 0;
  A(4, 2) = -1.5;
  A(4, 3) = 2.0;

  C << 0.0, 1.0 / 3.0, 1.0 / 3.0, 0.5, 1;

  BRes << 0.1, 0, 0.3, 0.4, 0.2;

  B << -1.0 / 15.0, 0, 0.3, -4.0 / 15.0, 1.0 / 30.0;
}

/**
  This procedure calculate stress increment using Runge - Kutta pair as given by
  Sloan (1987). The coeficients are given in Sloan (1987), Ordinary and Partial
  Differential Equations routines in... by Lee & Schiesser, Chapman & Hall 2003
  or originally by England (1969) Error Estimates for Runge - Kutta type
  solutions
  to systems of ordinary differential equations, Computer Journal 12 - 166-170.
  The procedure consists of 6 stages and should be the least
  efficient from all the R-K pairs presented
 */
std::tuple<double, int>
ShengMohrCoulomb::plasticRKEng654(StateMohrCoulomb& state,
                                  const Vector7& epStrain) const
{
  constexpr int Steps = 6;
  constexpr int Order = 5;

  using RKMatrix = Eigen::Matrix<double, Steps, Steps>;
  using RKVector = Eigen::Matrix<double, Steps, 1>;

  RKMatrix A = RKMatrix::Zero();
  RKVector C = RKVector::Zero();
  RKVector BRes = RKVector::Zero();
  RKVector B = RKVector::Zero();

  getParamRKEng654(A, B, BRes, C);

  bool errorEstimate = false;
  auto result =
    doRungeKutta<Order, Steps>(A, B, BRes, C, state, epStrain, errorEstimate);
  return result;
}

void
ShengMohrCoulomb::getParamRKEng654(Eigen::Matrix<double, 6, 6>& A,
                                   Eigen::Matrix<double, 6, 1>& B,
                                   Eigen::Matrix<double, 6, 1>& BRes,
                                   Eigen::Matrix<double, 6, 1>& C) const
{
  A(0, 0) = 0;
  A(1, 0) = 0.5;
  A(2, 0) = 0.25;
  A(2, 1) = 0.25;
  A(3, 0) = 0;
  A(3, 1) = -1.0;
  A(3, 2) = 2.0;
  A(4, 0) = 7.0 / 27.0;
  A(4, 1) = 10.0 / 27.0;
  A(4, 2) = 0;
  A(4, 3) = 1.0 / 27.0;
  A(5, 0) = 28.0 / 625.0;
  A(5, 1) = -125.0 / 625.0;
  A(5, 2) = 546.0 / 625.0;
  A(5, 3) = 54.0 / 625.0;
  A(5, 4) = -378.0 / 625.0;

  C << 0.0, 0.5, 0.5, 1, 2.0 / 3.0, 0.2;

  BRes << 14.0 / 336.0, 0, 0, 35.0 / 336.0, 162.0 / 336.0, 125.0 / 336.0;

  B << -42.0 / 336.0, 0, -224.0 / 336.0, -21.0 / 336.0, 162.0 / 336.0,
    125.0 / 336.0;
}

/**
  This procedure calculate stress increment using Runge - Kutta pair as given by
  Cash - Karp. The coeficients are given in Numerical Recipes (Cambrige Univ
  Press)
  or Cash Karp (1990) ACM Transactions on Mathematical software, vol 16, pp
  201-222.
  The procedure consists of 6 stages and should be
  less efficient than the RKErr8544 that uses Bogacki - Shimpine pair.
 */
std::tuple<double, int>
ShengMohrCoulomb::plasticRKCK654(StateMohrCoulomb& state,
                                 const Vector7& epStrain) const
{
  constexpr int Steps = 6;
  constexpr int Order = 5;

  using RKMatrix = Eigen::Matrix<double, Steps, Steps>;
  using RKVector = Eigen::Matrix<double, Steps, 1>;

  RKMatrix A = RKMatrix::Zero();
  RKVector C = RKVector::Zero();
  RKVector BRes = RKVector::Zero();
  RKVector B = RKVector::Zero();

  getParamRKCK654(A, B, BRes, C);

  bool errorEstimate = false;
  auto result =
    doRungeKutta<Order, Steps>(A, B, BRes, C, state, epStrain, errorEstimate);
  return result;
}

void
ShengMohrCoulomb::getParamRKCK654(Eigen::Matrix<double, 6, 6>& A,
                                  Eigen::Matrix<double, 6, 1>& B,
                                  Eigen::Matrix<double, 6, 1>& BRes,
                                  Eigen::Matrix<double, 6, 1>& C) const
{
  A(0, 0) = 0;
  A(1, 0) = 0.2;
  A(2, 0) = 0.075;
  A(2, 1) = 0.225;
  A(3, 0) = 0.3;
  A(3, 1) = -0.9;
  A(3, 2) = 1.2;
  A(4, 0) = -11.0 / 54.0;
  A(4, 1) = 2.5;
  A(4, 2) = -70.0 / 27.0;
  A(4, 3) = 35.0 / 27.0;
  A(5, 0) = 1631.0 / 55296.0;
  A(5, 1) = 175.0 / 512.0;
  A(5, 2) = 575.0 / 13824.0;
  A(5, 3) = 44275.0 / 110592.0;
  A(5, 4) = 253.0 / 4096.0;

  C << 0.0, 0.2, 0.3, 0.6, 1, 0.875;

  BRes << 37.0 / 378.0, 0, 250.0 / 621.0, 125.0 / 594.0, 0, 512.0 / 1771.0;

  B << 2825.0 / 27648.0, 0, 18575.0 / 48384.0, 13525.0 / 55296.0,
    277.0 / 14336.0, 0.25;
}

/**
  This procedure calculate stress increment using Runge - Kutta pair as given by
  DORMAND - PRINCE. The used pair is known also as DOPRI5 or RK5(4)7FM. The
  procedure consists of 7 stages and should be less
  efficient than the RKErr8544 that uses Bogacki - Shimpine pair.
 */
std::tuple<double, int>
ShengMohrCoulomb::plasticRKDP754(StateMohrCoulomb& state,
                                 const Vector7& epStrain) const
{
  constexpr int Steps = 7;
  constexpr int Order = 5;

  using RKMatrix = Eigen::Matrix<double, Steps, Steps>;
  using RKVector = Eigen::Matrix<double, Steps, 1>;

  RKMatrix A = RKMatrix::Zero();
  RKVector C = RKVector::Zero();
  RKVector BRes = RKVector::Zero();
  RKVector B = RKVector::Zero();

  getParamRKDP754(A, B, BRes, C);

  bool errorEstimate = false;
  auto result =
    doRungeKutta<Order, Steps>(A, B, BRes, C, state, epStrain, errorEstimate);
  return result;
}

void
ShengMohrCoulomb::getParamRKDP754(Eigen::Matrix<double, 7, 7>& A,
                                  Eigen::Matrix<double, 7, 1>& B,
                                  Eigen::Matrix<double, 7, 1>& BRes,
                                  Eigen::Matrix<double, 7, 1>& C) const
{
  A(0, 0) = 0;
  A(1, 0) = 0.2;
  A(2, 0) = 0.075;
  A(2, 1) = 0.225;
  A(3, 0) = 44.0 / 45.0;
  A(3, 1) = -56.0 / 15.0;
  A(3, 2) = 32.0 / 9.0;
  A(4, 0) = 19372.0 / 6561.0;
  A(4, 1) = -25360.0 / 2187.0;
  A(4, 2) = 64448.0 / 6561.0;
  A(4, 3) = -212.0 / 729.0;
  A(5, 0) = 9017.0 / 3168.0;
  A(5, 1) = -355.0 / 33.0;
  A(5, 2) = 46732.0 / 5247.0;
  A(5, 3) = 49.0 / 176.0;
  A(5, 4) = -5103.0 / 18656.0;
  A(6, 0) = 35.0 / 384.0;
  A(6, 1) = 0;
  A(6, 2) = 500.0 / 1113.0;
  A(6, 3) = 125.0 / 192.0;
  A(6, 4) = -2187.0 / 6784.0;
  A(6, 5) = 11.0 / 84.0;

  C << 0.0, 0.2, 0.3, 0.8, 8.0 / 9.0, 1, 1;

  BRes << 35.0 / 384.0, 0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0,
    11.0 / 84.0, 0;

  B << 5179.0 / 57600.0, 0, 7571.0 / 16695.0, 393.0 / 640.0,
    -92097.0 / 339200.0, 187.0 / 2100.0, 0.025;
}

/**
  The procedure uses the embedded Runge - Kutta integration scheme with
  Adaptive Stepsize Control the constants are as proposed by Bogacki and
  Shampine (1996), An efficient R-K (4,5) pair, Computers Math Applic,
  Vol 32 No 6 pp 15-28 with FSAL feauture the method allows for getting the
  error estimate and calculating value in one go
  It is arguably better than any other 5(4) RK pair; It has double 4th order
  error estimate

  The coefficients E(*) refer to an estimate of the local error based on
  the first formula of order 4.  It is the difference of the fifth order
  result, here located in A(8,*), and the fourth order result.  By
  construction both ErrorCoef(1) and ErrorCoef(6) are zero.
 */
std::tuple<double, int>
ShengMohrCoulomb::plasticRKErr8544(StateMohrCoulomb& state,
                                   const Vector7& epStrain) const
{
  constexpr int Steps = 8;
  constexpr int Order = 5;

  using RKMatrix = Eigen::Matrix<double, Steps, Steps>;
  using RKVector = Eigen::Matrix<double, Steps, 1>;
  using RKVectorErr = Eigen::Matrix<double, Steps - 1, 1>;

  RKMatrix A = RKMatrix::Zero();
  RKVector C = RKVector::Zero();
  RKVector BRes = RKVector::Zero();
  RKVector B = RKVector::Zero();
  RKVectorErr ErrorCoef = RKVectorErr::Zero();

  getParamRKErr8544(A, B, BRes, C, ErrorCoef);

  bool errorEstimate = false;
  auto result = doRungeKuttaErr<Order, Steps>(A, B, BRes, C, ErrorCoef, state,
                                              epStrain, errorEstimate);
  return result;
}

void
ShengMohrCoulomb::getParamRKErr8544(Eigen::Matrix<double, 8, 8>& A,
                                    Eigen::Matrix<double, 8, 1>& B,
                                    Eigen::Matrix<double, 8, 1>& BRes,
                                    Eigen::Matrix<double, 8, 1>& C,
                                    Eigen::Matrix<double, 7, 1>& ErrorCoef) const
{
  A(0, 0) = 0;
  A(1, 0) = 1.0 / 6.0;
  A(2, 0) = 2.0 / 27.0;
  A(2, 1) = 4.0 / 27.0;
  A(3, 0) = 183.0 / 1372.0;
  A(3, 1) = -162.0 / 343.0;
  A(3, 2) = 1053.0 / 1372.0;
  A(4, 0) = 68.0 / 297.0;
  A(4, 1) = -4.0 / 11.0;
  A(4, 2) = 42.0 / 143.0;
  A(4, 3) = 1960.0 / 3861.0;
  A(5, 0) = 597.0 / 22528.0;
  A(5, 1) = 81.0 / 352.0;
  A(5, 2) = 63099.0 / 585728.0;
  A(5, 3) = 58653.0 / 366080.0;
  A(5, 4) = 4617.0 / 20480.0;
  A(6, 0) = 174197.0 / 959244.0;
  A(6, 1) = -30942.0 / 79937.0;
  A(6, 2) = 8152137.0 / 19744439.0;
  A(6, 3) = 666106.0 / 1039181.0;
  A(6, 4) = -29421.0 / 29068.0;
  A(6, 5) = 482048.0 / 414219.0;
  A(7, 0) = 587.0 / 8064.0;
  A(7, 1) = 0.0;
  A(7, 2) = 4440339.0 / 15491840.0;
  A(7, 3) = 24353.0 / 124800.0;
  A(7, 4) = 387.0 / 44800.0;
  A(7, 5) = 2152.0 / 5985.0;
  A(7, 6) = 7267.0 / 94080.0;

  C(0) = 0.0;
  C(1) = 1.0 / 6.0;
  C(2) = 2.0 / 9.0;
  C(3) = 3.0 / 7.0;
  C(4) = 2.0 / 3.0;
  C(5) = 3.0 / 4.0;
  C(6) = 1.0;
  C(7) = 1.0;

  BRes = A.row(7);
  BRes(7) = 0; // as this scheme is the FSAL scheme

  B(0) = 2479.0 / 34992.0;
  B(1) = 0.0;
  B(2) = 123.0 / 416.0;
  B(3) = 612941.0 / 3411720.0;
  B(4) = 43.0 / 1440.0;
  B(5) = 2272.0 / 6561.0;
  B(6) = 79937.0 / 1113912.0;
  B(7) = 3293.0 / 556956.0;

  ErrorCoef(0) = -3.0 / 1280.0;
  ErrorCoef(1) = 0.0;
  ErrorCoef(2) = 6561.0 / 632320.0;
  ErrorCoef(3) = -343.0 / 20800.0;
  ErrorCoef(4) = 243.0 / 12800.0;
  ErrorCoef(5) = -1.0 / 95.0;
  ErrorCoef(6) = 0.0;
}


/**
 * Use table to integrate
 */
std::tuple<double, int>
ShengMohrCoulomb::plasticExtrapol(StateMohrCoulomb& state,
                                  const Vector7& epStrain) const
{
  // six stresses, p0*, six absolute stresses, absolute p0*,
  // 6 plastic strains
  constexpr int TableRows = 20;
  constexpr int StepMax = 15;
  using ApproxTableMatrix = Eigen::Matrix<double, TableRows, StepMax + 1>;
  using HTableMatrix = Eigen::Matrix<double, StepMax, StepMax + 1>;

  Vector6 Zero = Vector6::Zero();

  std::chrono::time_point<std::chrono::system_clock> startTime, endTime;
  startTime = std::chrono::system_clock::now();

  HTableMatrix hSquareTable = HTableMatrix::Zero();
  for (int i = 1; i < StepMax + 1; i++) {
    for (int j = 0; j < i; j++) {
      hSquareTable(i, j) = (DivInt[i] / DivInt[j]) * (DivInt[i] / DivInt[j]);
      if (dbg.active()) {
        dbg << "DivInt[" << i << "] = " << DivInt[i] << " DivInt[" << j
            << "] = " << DivInt[j] << "\n";
        dbg << "hSquareTable[" << i << ", " << j << "]=" << hSquareTable(i, j)
            << "\n";
      }
    }
  }

  StateMohrCoulomb state_new = state;
  StateMohrCoulomb state_copy = state;
  StateMohrCoulomb state_old = state;

  Vector6 initPlasticStrain = state.plasticStrain;
  Vector7 initStress = Vector7::Zero();
  initStress.block<6, 1>(0, 0) = state.stress;
  initStress(6) = state.p0Star();

  ApproxTableMatrix approxTable = ApproxTableMatrix::Zero();
  ApproxTableMatrix approxTableOld = ApproxTableMatrix::Zero();

  double dP0Star = 0.0;
  Vector6 plasticStrain = Vector6::Zero();
  Vector7 dSigma = Vector7::Zero();
  Vector7 absStress = Vector7::Zero();
  Vector7 dError = Vector7::Zero();

  // calculating stress increment using midpoint rule
  int loop = 0;
  for (; loop < StepMax; loop++) {

    // calculation of stress increment using the midpoint procedure
    state_copy = state_new;
    plasticMidpoint(state_new, epStrain, absStress, DivInt[loop]);

    // saving data using midpoint rule
    auto curCol = approxTable.col(0);
    curCol.block<6, 1>(0, 0) = state_new.stress - initStress.block<6, 1>(0, 0);
    curCol(6) = state_new.p0Star() - initStress(6);
    curCol.block<7, 1>(7, 0) = absStress;
    curCol.block<6, 1>(14, 0) = state_new.plasticStrain - initPlasticStrain;
    approxTable.col(0) = curCol;

    for (int i = 0; i < loop; i++) {
      for (int j = 0; j < TableRows; j++) {
        approxTable(j, i + 1) = approxTable(j, i) +
                                (approxTable(j, i) - approxTableOld(j, i)) /
                                  (hSquareTable(loop - i - 1, loop) - 1);
      }
    }

    if (dbg.active()) {
      for (int i = 0; i < loop + 1; i++)
        dbg << approxTable.col(i) << "\n";
    }

    // approximations are calculated
    // two possibilities of error control. In literature rather the second one;
    // this is the more stringent, but first one should be enough
    // (1) use approxTable[loop][i] - approxTable[loop-1][i]
    // (2) use approxTable[loop][i] - approxTableOld [loop-1][i]
    // still using relative error definition...
    // state_old (used for calculating the norm in q) is set accordingly

    // FIRST ONE OK, SEE Deuflhard Bornemann 2002 Springer p.206
    double rError = 0;
    if (loop > 0) {
      auto prevData = approxTable.col(loop - 1);
      auto currData = approxTable.col(loop);
      dSigma = currData.block<7, 1>(0, 0);
      dP0Star = currData(6);
      dError = currData.block<7, 1>(0, 0) - prevData.block<7, 1>(0, 0);

      switch (d_int.d_tolMethod) {
        case ToleranceMethod::EPUS_RELATIVE_ERROR:
          rError = checkNorm(dSigma, dP0Star, state_copy, dError);
          break;
        case ToleranceMethod::SLOAN:
          rError = checkNormSloan(dSigma, dP0Star, state_copy, dError);
          break;
        default:
          std::ostringstream err;
          err << "**ERROR** Improper d_tolMethod in increment.dta"
              << "\n";
          throw InvalidValue(err.str(), __FILE__, __LINE__);
      } // end switch

    } else {
      rError = 2 * d_int.d_integrationTol;
    }

    if (dbg.active()) {
      dbg << "Relative error after iteration " << loop
          << " is equal to: " << rError << "\n";
    }

    approxTableOld = approxTable;

    // error less then requested...check for drift and add the error
    if (rError < d_int.d_integrationTol) {
      auto currData = approxTable.col(loop);
      dSigma = currData.block<7, 1>(0, 0);
      dP0Star = currData(6);

      state_copy = state_old;
      state_old.update(Zero, epStrain, dSigma.block<6, 1>(0, 0), dP0Star);
      state_old = state_new;

      switch (d_int.d_driftCorrection) {
        case DriftCorrection::NO_CORRECTION:
          break;
        case DriftCorrection::CORRECTION_AT_BEGIN:
          correctDriftBeg(state_new, state_copy);
          break;
        case DriftCorrection::CORRECTION_AT_END:
          correctDriftEnd(state_new);
          break;
        default:
          break;
      }

      dError.block<6, 1>(0, 0) = state_new.stress - state_old.stress;
      dError(6) = state_new.p0Star() - state_old.p0Star();
      switch (d_int.d_tolMethod) {
        case ToleranceMethod::EPUS_RELATIVE_ERROR:
          rError = checkNorm(dSigma, dP0Star, state_copy, dError);
          break;
        case ToleranceMethod::SLOAN:
          rError = checkNormSloan(dSigma, dP0Star, state_copy, dError);
          break;
        default:
          std::ostringstream err;
          err << "**ERROR** Improper d_tolMethod in increment.dta"
              << "\n";
          throw InvalidValue(err.str(), __FILE__, __LINE__);
      } // end switch

      if (rError < d_int.d_integrationTol) {
        loop = loop + 100; // no more interations after - finished
      }
    }
  }

  // done - time to update everything and sent time and no of iterations...
  int numIter = 0.0;
  if (loop > 100) {
    loop = loop - 101;

    auto currData = approxTable.col(loop);
    absStress = currData.block<7, 1>(7, 0);
    plasticStrain = currData.block<6, 1>(14, 0);
    dP0Star = currData(6);
    state.update(plasticStrain, epStrain, dSigma.block<6, 1>(0, 0), dP0Star);

    endTime = std::chrono::system_clock::now();

    for (int i = 0; i < loop + 1; i++)
      numIter += numIter + DivInt[i];
    if (dbg.active()) {
      dbg << "Procedure has coverged after " << numIter << " iterations.\n";
    }
  } else {
    loop--;

    auto currData = approxTable.col(loop);
    dSigma = currData.block<7, 1>(0, 0);
    absStress = currData.block<7, 1>(7, 0);
    plasticStrain = currData.block<6, 1>(14, 0);
    dP0Star = currData(6);

    state.update(plasticStrain, epStrain, dSigma.block<6, 1>(0, 0), dP0Star);

    endTime = std::chrono::system_clock::now();

    for (int i = 0; i < loop + 1; i++) {
      numIter += DivInt[i];
      if (dbg.active()) {
        dbg << "Procedure has NOT CONVERGED after " << numIter
            << "iterations.\n";
      }
    }
  }

  double time = std::chrono::duration<double>(endTime - startTime).count();
  if (dbg.active()) {
    dbg << "Calculation took: " << time << "s."
        << "\n";
  }

  return std::make_tuple(time, numIter);
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
template <int Order, int Steps>
std::tuple<double, int>
ShengMohrCoulomb::doRungeKutta(const Eigen::Matrix<double, Steps, Steps>& AA,
                               const Eigen::Matrix<double, Steps, 1>& BB,
                               const Eigen::Matrix<double, Steps, 1>& BRes,
                               const Eigen::Matrix<double, Steps, 1>& CC,
                               StateMohrCoulomb& state, const Vector7& epStrain,
                               bool errorEstimate) const
{
  Vector7 substepStrain = epStrain;
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
  double stepLength = 1;
  double totalSize = 0;
  bool reuseStep = false;
  double microStep = 0;
  double stepAccuracyCheck = 0;
  double frequency = 100000.0 / Order; // how often display info about steps
  double methodPower = std::pow(2.0, Order) * d_int.d_integrationTol;

  std::vector<StateMohrCoulomb> midStates(Steps);

  Eigen::Matrix<double, 6, Steps> dSigma =
    Eigen::Matrix<double, 6, Steps>::Zero();
  Eigen::Matrix<double, 6, Steps> plasticStrain =
    Eigen::Matrix<double, 6, Steps>::Zero();
  Eigen::Matrix<double, 1, Steps> dP0Star =
    Eigen::Matrix<double, 1, Steps>::Zero();

  Vector7 reuseRes = Vector7::Zero();
  Vector6 stressInc = Vector6::Zero();
  Vector6 plasticStrainInc = Vector6::Zero();
  double p0StarInc = 0.0;

  std::chrono::time_point<std::chrono::system_clock> startTime, endTime;
  startTime = std::chrono::system_clock::now();

  int stepNo = 0;
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
      dbg << "Step Length = " << stepLength
          << " Current strain (0) = " << currentStrain(0) << "\n";
    }

    // Make copies of the initial state
    for (auto& midState : midStates) {
      midState = state;
    }

    // Below the main R-K loop to calculate the value of intermediate stresses;
    // values stored in dSigma_temp
    if (reuseStep) {

      dSigma.col(0) = reuseRes.block<6, 1>(0, 0);
      plasticStrain.col(0) = plasticStrainInc.block<6, 1>(0, 0);
      dP0Star(0, 0) = reuseRes(6);

    } else {

      stressInc = Vector6::Zero();
      plasticStrainInc = Vector6::Zero();
      calcPlastic(midStates[0], substepStrain, stressInc, plasticStrainInc,
                  p0StarInc);

      // **TODO** Check if plastic strain increment should not be computed
      // plasticStrainInc = Vector6::Zero();
      dSigma.col(0) = stressInc;
      plasticStrain.col(0) = plasticStrainInc;
      dP0Star(0, 0) = p0StarInc;
    }

    for (int rkloop = 1; rkloop < Steps; rkloop++) {

      stressInc = Vector6::Zero();
      plasticStrainInc = Vector6::Zero();
      p0StarInc = 0;

      currentStrain = substepStrain * CC[rkloop];

      for (int i = 0; i < rkloop; i++) {

        stressInc += (dSigma.col(i) * AA(rkloop, i));
        plasticStrainInc += (plasticStrain.col(i) * AA(rkloop, i));
        p0StarInc += (dP0Star(0, i) * AA(rkloop, i));
      }

      midStates[rkloop].update(plasticStrainInc, currentStrain, stressInc,
                               p0StarInc);

      calcPlastic(midStates[rkloop], substepStrain, stressInc, plasticStrainInc,
                  p0StarInc);

      // **TODO** Check if plastic strain increment should not be computed
      // plasticStrainInc = Vector6::Zero();
      dSigma.col(rkloop) = stressInc;
      plasticStrain.col(rkloop) = plasticStrainInc;
      dP0Star(0, rkloop) = p0StarInc;
    }

    Vector6 sigBRes = dSigma * BRes;
    // Vector6 epsPBRes = plasticStrain * BRes;
    double p0BRes = dP0Star * BRes;
    Vector6 sigB = dSigma * BB;
    double p0B = dP0Star * BB;

    Vector7 result, error;
    result.block<6, 1>(0, 0) = sigBRes;
    error.block<6, 1>(0, 0) = sigB;
    result(6) = p0BRes;
    error(6) = p0B;

    // error estimate calculated in case we
    // have lower order solution instead of
    // error estimate
    if (!errorEstimate) {
      error -= result;
    }

    // check the error norm
    switch (d_int.d_tolMethod) {
      case ToleranceMethod::EPUS_RELATIVE_ERROR:
        rError = checkNorm(result, result(6), state, error);
        break;
      case ToleranceMethod::SLOAN:
        rError = checkNormSloan(result, result(6), state, error);
        break;
      default:
        std::ostringstream err;
        err << "**ERROR** Improper d_tolMethod in increment.dta"
            << "\n";
        throw InvalidValue(err.str(), __FILE__, __LINE__);
    }

    for (int i = 0; i < 7; i++) {
      if (!std::isfinite(result(i))) {
        std::ostringstream err;
        err << "**ERROR** Results not a number. Try correcting issue, results "
               "may be "
               "incorrect."
            << " Results = \n"
            << result << "\n";
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

    if (d_int.d_tolMethod == ToleranceMethod::EPUS_RELATIVE_ERROR) {
      newStepSize =
        d_int.d_betaFactor *
        std::pow(d_int.d_integrationTol / rError, (1 / (Order - 1.0)));
    } else {
      newStepSize = d_int.d_betaFactor *
                    std::pow(d_int.d_integrationTol / rError, (1 / Order));
    }

    if (!stepAccepted) {

      // here we should take care about correct re - usage of the first
      // evaluation of derivative
      reuseStep = true;
      reuseRes.block<6, 1>(0, 0) = dSigma.col(0) * newStepSize;
      reuseRes(6) = dP0Star(0, 0) * newStepSize;
      plasticStrainInc.block<6, 1>(0, 0) = plasticStrain.col(0) * newStepSize;

    } else {

      StateMohrCoulomb state_old = state;

      //******
      //******
      // here we update all the state data.
      state.update(plasticStrainInc, substepStrain, result.block<6, 1>(0, 0),
                   result(6));
      //******
      //******

      // Keeping a copy to check values
      StateMohrCoulomb trialState = state;

      switch (d_int.d_driftCorrection) {
        case DriftCorrection::CORRECTION_AT_BEGIN:
          correctDriftBeg(state, state_old);
          break;
        case DriftCorrection::CORRECTION_AT_END:
          correctDriftEnd(state);
          break;
        case DriftCorrection::NO_CORRECTION:
        default:
          break;
      }

      // reevaluate the error in the point:
      error.block<6, 1>(0, 0) += (state.stress - trialState.stress);
      error(6) += state.p0Star() - trialState.p0Star();

      // error vector updated, norm should be re-evaluated:
      switch (d_int.d_tolMethod) {
        case ToleranceMethod::EPUS_RELATIVE_ERROR:
          rError = checkNorm(result, result(6), state, error);
          break;
        case ToleranceMethod::SLOAN:
          rError = checkNormSloan(result, result(6), state, error);
          break;
        default:
          std::ostringstream err;
          err << "**ERROR** Improper d_tolMethod in increment.dta"
              << "\n";
          throw InvalidValue(err.str(), __FILE__, __LINE__);
      }

      if (!std::isfinite(rError)) {
        if (rError < methodPower) {
          rError = methodPower;
        }
      }

      if (d_int.d_tolMethod == ToleranceMethod::EPUS_RELATIVE_ERROR) {
        newStepSize =
          d_int.d_betaFactor *
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
        reuseRes.block<6, 1>(0, 0) = dSigma.col(0) * newStepSize;
        reuseRes(6) = dP0Star(0, 0) * newStepSize;
        plasticStrainInc.block<6, 1>(0, 0) = plasticStrain.col(0) * newStepSize;

        state_old = state;
        double temp = double(stepNo) / frequency;
        if (std::modf(temp, &temp) == 0) {
          std::cout << "Step number : " << stepNo << "\n";
          std::cout << "Total size done is : " << totalSize
                    << " of whole step. Current stepLength = " << stepLength
                    << "\n";
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
                    << " of whole step. Current stepLength = " << stepLength
                    << "\n";
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
 * Slight variation of above with different error estimate and adaptive step
 * size
 */
template <int Order, int Steps>
std::tuple<double, int>
ShengMohrCoulomb::doRungeKuttaErr(
  const Eigen::Matrix<double, Steps, Steps>& AA,
  const Eigen::Matrix<double, Steps, 1>& BB,
  const Eigen::Matrix<double, Steps, 1>& BRes,
  const Eigen::Matrix<double, Steps, 1>& CC,
  const Eigen::Matrix<double, Steps - 1, 1>& ErrCoef, StateMohrCoulomb& state,
  const Vector7& epStrain, bool errorEstimate) const
{
  constexpr double criticalStepSize = 1.0e-4;
  Vector7 substepStrain = epStrain;
  double newStepSize = 0;
  for (int i = 0; i < 6; i++) {
    newStepSize += std::abs(substepStrain[i]);
  }
  newStepSize = 0.01 / (newStepSize * Order);
  if (newStepSize > 1) {
    newStepSize = 1;
  }
  if (dbg.active()) {
    dbg << "In doRungeKuttaErr::newStepSize = " << newStepSize << "\n";
  }

  double stepLength = 1;
  double totalSize = 0;
  bool reuseStep = false;
  double microStep = 0;
  double stepAccuracyCheck = 0;
  double frequency = 100000.0 / Order; // how often display info about steps
  double methodPower = std::pow(2.0, Order) * d_int.d_integrationTol;

  std::vector<StateMohrCoulomb> midStates(Steps);

  Eigen::Matrix<double, 6, Steps> dSigma =
    Eigen::Matrix<double, 6, Steps>::Zero();
  Eigen::Matrix<double, 6, Steps> plasticStrain =
    Eigen::Matrix<double, 6, Steps>::Zero();
  Eigen::Matrix<double, 1, Steps> dP0Star =
    Eigen::Matrix<double, 1, Steps>::Zero();

  Vector7 reuseRes = Vector7::Zero();
  Vector6 stressInc = Vector6::Zero();
  Vector6 plasticStrainInc = Vector6::Zero();
  double p0StarInc = 0.0;

  std::chrono::time_point<std::chrono::system_clock> startTime, endTime;
  startTime = std::chrono::system_clock::now();

  int stepNo = 0;
  bool finished = false;
  do {
    bool stepAccepted = false;
    ++stepNo;

    // compute size of a step
    stepLength *= newStepSize;

    // check whether the step exceeds the whole increment
    if ((stepLength + totalSize) > 1) {
      stepLength = 1 - totalSize;
    }

    // strain increment in current step
    Vector7 substepStrain = epStrain * stepLength;
    Vector7 currentStrain = Vector7::Zero();

    double rError = 0;
    double rErrorOther = 0;

    if (dbg.active()) {
      dbg << "In doRungeKuttaErr:: Step Length = " << stepLength
          << " Current strain (0) = " << currentStrain(0) << "\n";
    }

    // Make copies of the initial state
    for (auto& midState : midStates) {
      midState = state;
    }

    // Below the main R-K loop to calculate the value of intermediate stresses;
    // values stored in dSigma_temp
    if (reuseStep) {

      dSigma.col(0) = reuseRes.block<6, 1>(0, 0);
      plasticStrain.col(0) = plasticStrainInc.block<6, 1>(0, 0);
      dP0Star(0, 0) = reuseRes(6);

    } else {

      stressInc = Vector6::Zero();
      plasticStrainInc = Vector6::Zero();
      calcPlastic(midStates[0], substepStrain, stressInc, plasticStrainInc,
                  p0StarInc);

      // **TODO** Check if plastic strain increment should not be computed
      // plasticStrainInc = Vector6::Zero();
      dSigma.col(0) = stressInc;
      plasticStrain.col(0) = plasticStrainInc;
      dP0Star(0, 0) = p0StarInc;
    }

    for (int rkloop = 1; rkloop < Steps; rkloop++) {

      stressInc = Vector6::Zero();
      plasticStrainInc = Vector6::Zero();
      p0StarInc = 0;

      currentStrain = substepStrain * CC[rkloop];

      for (int i = 0; i < rkloop; i++) {

        stressInc += (dSigma.col(i) * AA(rkloop, i));
        plasticStrainInc += (plasticStrain.col(i) * AA(rkloop, i));
        p0StarInc += (dP0Star(0, i) * AA(rkloop, i));
      }

      midStates[rkloop].update(plasticStrainInc, currentStrain, stressInc,
                               p0StarInc);

      calcPlastic(midStates[rkloop], substepStrain, stressInc, plasticStrainInc,
                  p0StarInc);

      // **TODO** Check if plastic strain increment should not be computed
      // plasticStrainInc = Vector6::Zero();
      dSigma.col(rkloop) = stressInc;
      plasticStrain.col(rkloop) = plasticStrainInc;
      dP0Star(0, rkloop) = p0StarInc;
    }

    Vector6 sigBRes = dSigma * BRes;
    // Vector6 epsPBRes = plasticStrain * BRes;
    double p0BRes = dP0Star * BRes;
    Vector6 sigB = dSigma * BB;
    double p0B = dP0Star * BB;

    Vector7 result, error;
    result.block<6, 1>(0, 0) = sigBRes;
    error.block<6, 1>(0, 0) = sigB;
    result(6) = p0BRes;
    error(6) = p0B;

    // error estimate calculated in case we
    // have lower order solution instead of
    // error estimate
    if (!errorEstimate) {
      error -= result;
    }

    // For the a-priori error estimate
    Vector7 errorOther = Vector7::Zero();
    for (int i = 0; i < (Steps - 1); i++) {
      for (int j = 0; j < 6; j++) {
        errorOther[j] += dSigma(j, i) * ErrCoef(i);
        errorOther[6] += dP0Star(i) * ErrCoef(i);
      }
    }

    // check the error norm
    switch (d_int.d_tolMethod) {
      case ToleranceMethod::EPUS_RELATIVE_ERROR:
        rError = checkNorm(result, result(6), state, error);
        rErrorOther = checkNorm(result, result(6), state, errorOther);
        break;
      case ToleranceMethod::SLOAN:
        rError = checkNormSloan(result, result(6), state, error);
        rErrorOther = checkNormSloan(result, result(6), state, errorOther);
        break;
      default:
        std::ostringstream err;
        err << "**ERROR** Improper d_tolMethod in increment.dta"
            << "\n";
        throw InvalidValue(err.str(), __FILE__, __LINE__);
    }

    for (int i = 0; i < 7; i++) {
      if (!std::isfinite(result(i))) {
        std::cout << "**ERROR** Results not a number. Try correcting issue, "
                     "results may be "
                     "incorrect. Results = \n"
                  << result << "\n";
        result(i) = 0;
        if (rError < methodPower) {
          rError = methodPower;
        }
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

    if (d_int.d_tolMethod == ToleranceMethod::EPUS_RELATIVE_ERROR) {
      newStepSize =
        d_int.d_betaFactor *
        std::pow(d_int.d_integrationTol / rError, (1 / (Order - 1.0)));
    } else {
      newStepSize = d_int.d_betaFactor *
                    std::pow(d_int.d_integrationTol / rError, (1 / Order));
    }

    if (!stepAccepted) {

      // here we should take care about correct re - usage of the first
      // evaluation of derivative
      reuseStep = true;
      reuseRes.block<6, 1>(0, 0) = dSigma.col(0) * newStepSize;
      reuseRes(6) = dP0Star(0, 0) * newStepSize;
      plasticStrainInc.block<6, 1>(0, 0) = plasticStrain.col(0) * newStepSize;

    } else {

      StateMohrCoulomb state_old = state;

      //******
      //******
      // here we update all the state data.
      state.update(plasticStrainInc, substepStrain, result.block<6, 1>(0, 0),
                   result(6));
      //******
      //******

      // Keeping a copy to check values
      StateMohrCoulomb trialState = state;

      switch (d_int.d_driftCorrection) {
        case DriftCorrection::CORRECTION_AT_BEGIN:
          correctDriftBeg(state, state_old);
          break;
        case DriftCorrection::CORRECTION_AT_END:
          correctDriftEnd(state);
          break;
        case DriftCorrection::NO_CORRECTION:
        default:
          break;
      }

      // reevaluate the error in the point:
      error.block<6, 1>(0, 0) += (state.stress - trialState.stress);
      error(6) += state.p0Star() - trialState.p0Star();
      errorOther.block<6, 1>(0, 0) += (state.stress - trialState.stress);
      errorOther(6) += state.p0Star() - trialState.p0Star();

      // error vector updated, norm should be re-evaluated:
      double temp = 0;
      switch (d_int.d_tolMethod) {
        case ToleranceMethod::EPUS_RELATIVE_ERROR:
          temp = checkNorm(result, result(6), state, error);
          rErrorOther = checkNorm(result, result(6), state, errorOther);
          if (temp > rError) {
            rError = temp;
          }
          if (rErrorOther > rError) {
            rError = rErrorOther;
          }
          break;
        case ToleranceMethod::SLOAN:
          temp = checkNormSloan(result, result(6), state, error);
          rErrorOther = checkNormSloan(result, result(6), state, errorOther);
          if (temp > rError) {
            rError = temp;
          }
          if (rErrorOther > rError) {
            rError = rErrorOther;
          }
          break;
        default:
          std::ostringstream err;
          err << "**ERROR** Improper d_tolMethod in increment.dta"
              << "\n";
          throw InvalidValue(err.str(), __FILE__, __LINE__);
      }

      if (!std::isfinite(rError)) {
        if (rError < methodPower) {
          rError = methodPower;
        }
      }

      if (d_int.d_tolMethod == ToleranceMethod::EPUS_RELATIVE_ERROR) {
        newStepSize =
          d_int.d_betaFactor *
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
        reuseRes.block<6, 1>(0, 0) = dSigma.col(0) * newStepSize;
        reuseRes(6) = dP0Star(0, 0) * newStepSize;
        plasticStrainInc.block<6, 1>(0, 0) = plasticStrain.col(0) * newStepSize;

        state_old = state;
        double temp = double(stepNo) / frequency;
        if (std::modf(temp, &temp) == 0) {
          std::cout << "Step number : " << stepNo << "\n";
          std::cout << "Total size done is : " << totalSize
                    << " of whole step. Current stepLength = " << stepLength
                    << "\n";
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
                    << " of whole step. Current stepLength = " << stepLength
                    << "\n";
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
                              const Vector7& epStrainInc, Vector6& dSigma,
                              Vector6& dEps_p, double& dP0Star) const
{
  if (!state.checkIfFinite()) {
    std::ostringstream err;
    err << "**Error** in the calcPlastic Procedure. "
        << "Point internal values are not finite.\n";
    throw InvalidValue(err.str(), __FILE__, __LINE__);
  }

  auto stress = state.stress;

  Vector6 df_dsigma = Vector6::Zero();
  Vector6 dg_dsigma = Vector6::Zero();
  std::tie(df_dsigma, dg_dsigma) = computeDfDsigma(stress);

  double K = d_elastic.d_K;
  double G = d_elastic.d_G;
  Matrix66 elasticTangent = calculateElasticTangentMatrix(K, G);

  Vector6 dEps = epStrainInc.block<6, 1>(0, 0);
  auto numerator = df_dsigma.transpose() * elasticTangent;
  double denominator = numerator * dg_dsigma;
  if (denominator < TINY) {
    std::cout << "**WARNING** denominator of plastic multiplier is very small."
                 " Some error may arise and results may be incorrect"
              << "\n";
  }

  dEps_p = (dg_dsigma * numerator * dEps) / denominator;
  dSigma = elasticTangent * (dEps - dEps_p);

  for (int i = 0; i < 6; i++) {
    if (!std::isfinite(dSigma(i))) {
      std::ostringstream err;
      err << "calcPlastic: dSigma not finite. dSigma = \n" << dSigma << "\n";
      throw InvalidValue(err.str(), __FILE__, __LINE__);
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
                            const Vector7& dError) const
{
  using Vector9 = Eigen::Matrix<double, 9, 1>;

  Vector9 drError = Vector9::Zero();
  double rError = 0;

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
    drError(7) = errorP / p;
  } else {
    drError(7) = errorP / TINY;
  }

  // norm of q...
  Vector6 initialSigma = initialState.stress;
  Vector6 finalSigma = initialSigma + dSigma.block<6, 1>(0, 0);
  double initialShear = initialState.shearStress();
  double finalShear =
    (finalSigma(0) - finalSigma(1)) * (finalSigma(0) - finalSigma(1)) +
    (finalSigma(0) - finalSigma(2)) * (finalSigma(0) - finalSigma(2)) +
    (finalSigma(1) - finalSigma(2)) * (finalSigma(1) - finalSigma(2)) +
    6 * (finalSigma(3) * finalSigma(3) + finalSigma(4) * finalSigma(4) +
         finalSigma(5) * finalSigma(5));
  finalShear = std::sqrt(0.5 * finalShear);

  finalSigma += dError.block<6, 1>(0, 0);
  double finalShearMax =
    (finalSigma(0) - finalSigma(1)) * (finalSigma(0) - finalSigma(1)) +
    (finalSigma(0) - finalSigma(2)) * (finalSigma(0) - finalSigma(2)) +
    (finalSigma(1) - finalSigma(2)) * (finalSigma(1) - finalSigma(2)) +
    6 * (finalSigma(3) * finalSigma(3) + finalSigma(4) * finalSigma(4) +
         finalSigma(5) * finalSigma(5));
  finalShearMax = std::sqrt(0.5 * finalShearMax);

  finalSigma -= dError.block<6, 1>(0, 0) * 2.0;
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
    drError(8) = shearError / dShear;
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
                                 const Vector7& dError) const
{
  using Vector9 = Eigen::Matrix<double, 9, 1>;

  Vector9 drError = Vector9::Zero();
  double rError = 0;

  Vector6 initialSigma = initialState.stress;
  Vector6 finalSigma = initialSigma + dSigma.block<6, 1>(0, 0);
  double dP0StarEnd = initialState.p0Star() + dP0Star;

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
    drError(7) = errorP / p;
  } else {
    drError(7) = errorP / TINY;
  }

  // norm of q...
  // double  initialShear  = initialState.shearStress();
  double finalShear =
    (finalSigma(0) - finalSigma(1)) * (finalSigma(0) - finalSigma(1)) +
    (finalSigma(0) - finalSigma(2)) * (finalSigma(0) - finalSigma(2)) +
    (finalSigma(1) - finalSigma(2)) * (finalSigma(1) - finalSigma(2)) +
    6 * (finalSigma(3) * finalSigma(3) + finalSigma(4) * finalSigma(4) +
         finalSigma(5) * finalSigma(5));
  finalShear = std::sqrt(0.5 * finalShear);

  finalSigma += dError.block<6, 1>(0, 0);
  double finalShearMax =
    (finalSigma(0) - finalSigma(1)) * (finalSigma(0) - finalSigma(1)) +
    (finalSigma(0) - finalSigma(2)) * (finalSigma(0) - finalSigma(2)) +
    (finalSigma(1) - finalSigma(2)) * (finalSigma(1) - finalSigma(2)) +
    6 * (finalSigma(3) * finalSigma(3) + finalSigma(4) * finalSigma(4) +
         finalSigma(5) * finalSigma(5));
  finalShearMax = std::sqrt(0.5 * finalShearMax);

  finalSigma -= (dError.block<6, 1>(0, 0) * 2.0);
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

/**
  This procedure corrects the drift. The method of drift correction is based on
  the Potts&Zdravkovic book; It is however slightly changes so to use it in
  unsaturated soil model. It is noted however, that in current version the
  algorithm
  does calculate all the derivatives in the forbidden space outside the yield
  locus which may be cause of concerns.

  Input: state to correct
  Output: Updated state, corrected using values of stress at the beginning

  Need dF/Ds, D, dF/dP0star, dP0star/DEpsPl, m, dF/dS
  at the beginning we need to calculate the derivatives a, b, c, d, g, p...)
 */
void
ShengMohrCoulomb::correctDriftBeg(StateMohrCoulomb& state,
                                  const StateMohrCoulomb& stateOld) const
{
  if (dbg.active()) {
    dbg << "Correct Drift Procedure entered!"
        << "\n";
  }

  Vector6 df_dsigma = Vector6::Zero();
  Vector6 dg_dsigma = Vector6::Zero();
  std::tie(df_dsigma, dg_dsigma) = computeDfDsigma(stateOld.stress);

  double K = d_elastic.d_K;
  double G = d_elastic.d_G;
  Matrix66 elasticTangent = calculateElasticTangentMatrix(K, G);

  auto numerator = df_dsigma.transpose() * elasticTangent;
  double denominator = numerator * dg_dsigma;

  bool correctDrift = true;
  int numIter = 0;
  do {
    ++numIter;
    double fValue = computeYieldNormalized(state.stress);
    if (std::abs(fValue) > d_int.d_yieldTol) {
      correctDrift = true;
    } else {
      correctDrift = false;
    }

    if (correctDrift) {
      double lambda = fValue / denominator;
      Vector6 dEps_p = df_dsigma * lambda;

      if (dbg.active()) {
        dbg << "Delta Epsilon Plastic:\n" << dEps_p << "\n";
      }

      auto dSigma = (elasticTangent * dg_dsigma) * (-lambda);

      Vector7 zeros = Vector7::Zero();
      state.update(dEps_p, zeros, dSigma, 0);
    }

    if (numIter > 10) {
      if (dbg.active()) {
        dbg << "**WARNING** Drift Correction Procedure failed."
            << "\n";
      }
      correctDrift = false;
    }
  } while (correctDrift);
}

/**
  Same as above except current stress value is used for correction.
  Input: state to correct
  Output: Updated state, corrected using values of stress at the end
 */
void
ShengMohrCoulomb::correctDriftEnd(StateMohrCoulomb& state) const
{
  if (dbg.active()) {
    dbg << "Correct Drift at End Procedure entered!"
        << "\n";
  }

  Vector6 df_dsigma = Vector6::Zero();
  Vector6 dg_dsigma = Vector6::Zero();
  double K = d_elastic.d_K;
  double G = d_elastic.d_G;

  bool correctDrift = true;
  int numIter = 0;
  do {
    double fValue = computeYieldNormalized(state.stress);

    if (std::abs(fValue) > d_int.d_yieldTol) {
      correctDrift = true;
    } else {
      correctDrift = false;
    }
    if (correctDrift) {
      // Correct for drift
      // here the drift will be corrected by using the D matrix from the
      // forbidden space.  Although because the drift will be checked again,
      // it shouldn't pose much problem.
      ++numIter;

      if (dbg.active()) {
        dbg << "Drift Correction, Iteration = " << numIter
            << " Function Value = " << fValue << "\n";
      }

      std::tie(df_dsigma, dg_dsigma) = computeDfDsigma(state.stress);

      Matrix66 elasticTangent = calculateElasticTangentMatrix(K, G);

      auto numerator = df_dsigma.transpose() * elasticTangent;
      double denominator = numerator * dg_dsigma;

      double lambda = fValue / denominator;

      Vector6 dEps_p = df_dsigma * lambda;

      if (dbg.active()) {
        dbg << "Delta Epsilon Plastic:\n" << dEps_p << "\n";
      }

      auto dSigma = (elasticTangent * dg_dsigma) * (-lambda);

      Vector7 zeros = Vector7::Zero();
      state.update(dEps_p, zeros, dSigma, 0);
    }
    if (numIter > 10) {
      if (dbg.active()) {
        dbg << "**WARNING** Drift Correction Procedure failed"
            << "\n";
      }
      correctDrift = false;
    }
  } while (correctDrift);
}

/**
  This procedure is to calculate the stress increase using the midpoint method
  with given number of iterations numIter.
  is made mainly for use in the extrapolation procedure. It has no adaptive
  substepping or error control. It just integrates
  the strain to get the stress in given number of substeps using the midpoint
  method
  */
double
ShengMohrCoulomb::plasticMidpoint(StateMohrCoulomb& state,
                                  const Vector7& epStrain, Vector7& absStress,
                                  int numIter) const
{
  double h = 1.0 / (double)numIter;

  absStress = Vector7::Zero();
  Vector7 currentStrain = epStrain * h;
  Vector7 halfCurrentStrain = currentStrain * 0.5;

  if (dbg.active()) {
    dbg << "Step Length = " << h << "\n";
    dbg << "Current strain [0] = " << currentStrain[0] << "\n";
  }

  StateMohrCoulomb state_new, state_mid;
  Vector6 dSigma = Vector6::Zero();
  Vector6 plasticStrain = Vector6::Zero();
  double dP0Star = 0.0;
  for (int loop = 0; loop < numIter; loop++) {
    state_new = state;
    state_mid = state;

    calcPlastic(state_new, halfCurrentStrain, dSigma, plasticStrain, dP0Star);
    state_mid.update(plasticStrain, halfCurrentStrain, dSigma, dP0Star);

    calcPlastic(state_mid, halfCurrentStrain, dSigma, plasticStrain, dP0Star);
    dSigma *= 2.0;
    plasticStrain *= 2.0;
    dP0Star *= 2.0;
    state.update(plasticStrain, currentStrain, dSigma, dP0Star);

    absStress.block<6, 1>(0, 0) += dSigma.cwiseAbs();
    absStress(6) += std::abs(dP0Star);
  }
  return 0;
}

/**
  This procedure returns tangent matrix. If given point lies on the yield locus,
  the matrix is elastoplastic.
  When the point lies inside the yield locus, the matrix is elastic.
  Point is supposed to have all the required data about describing its state,
  like stress state, specific volume nu and the hardening parameter p0*.

  Input: point, Matrix (7x6 - six 6 rows, seven 7 columns) that will be filled
  with coefficients corresponding to the tangent matrix
  Initialise BBMMatrix (6,7);
 */
Matrix67
ShengMohrCoulomb::getTangentMatrix(const StateMohrCoulomb& state) const
{

  Matrix67 DEP;

  if (checkYieldNormalized(state)) {
    DEP = calculateElastoPlasticTangentMatrix(state);
  } else {
    DEP = calculateElasticTangentMatrix(state);
  }

  return DEP;
}

/**
  Here the ElastoPlastic tangent matrix is calculated.
  For more theoretical explanation what is being done below please see
  'Evaluation of the tangent elastic matrix Del & elastoplastic matrix Dep'
  (evaluation of the tangent matrix D.doc)

  This procedure is quite similar to procedure calcPlastic where the stress
  increment using the tangent elasto-plastic matrix is calculated. The
  matrix itself, however, is not explicitly calculated there.
  */
Matrix67
ShengMohrCoulomb::calculateElastoPlasticTangentMatrix(
  const StateMohrCoulomb& state) const
{
  Vector6 df_dsigma, dg_dsigma;
  std::tie(df_dsigma, dg_dsigma) = computeDfDsigma(state.stress);

  Matrix66 tangent =
    calculateElasticTangentMatrix(d_elastic.d_K, d_elastic.d_G);

  auto numerator = df_dsigma.transpose() * tangent;
  double denominator = numerator * dg_dsigma;

  tangent -= (tangent * dg_dsigma * numerator) / denominator;

  Matrix67 DEP = Matrix67::Zero();
  DEP.block<6, 6>(0, 0) = tangent;

  return DEP;
}

/**
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
Vector7
ShengMohrCoulomb::computeG1(
  const StateMohrCoulomb& initialState, RetentionModel retentionModel,
  const std::vector<double>& retentionParameters) const
{
  Vector7 G1 = Vector7::Zero();

  switch (retentionModel) {
    case RetentionModel::STATE_SURFACE:
      double a, b, s;
      a = retentionParameters[0];
      b = retentionParameters[1];
      s = initialState.suction();
      G1(0) = (1 - std::exp(b * s)) * a / 3;
      G1(1) = G1(0);
      G1(2) = G1(0);
      G1(3) = 0;
      G1(4) = 0;
      G1(5) = 0;
      break;
    case RetentionModel::VAN_GENUCHTEN:
    case RetentionModel::GALLIPOLI:
    default:
      break;
  }

  return G1;
}

/**
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
double
ShengMohrCoulomb::computeG2(
  const StateMohrCoulomb& initialState, RetentionModel retentionModel,
  const std::vector<double>& retentionParameters) const
{
  double G2 = 0.0;

  switch (retentionModel) {
    case RetentionModel::STATE_SURFACE: {
      double a = retentionParameters[0];
      double b = retentionParameters[1];
      double s = initialState.suction();
      G2 = b * std::exp(b * s) * (0.5 - a * initialState.meanStress());
      break;
    }
    case RetentionModel::VAN_GENUCHTEN: {
      double Ew = retentionParameters[0];
      double Fw = retentionParameters[1];
      double Ssn = retentionParameters[2];
      double Sir = retentionParameters[3];
      double s = initialState.suction();
      double numerator = (Ssn - Sir) * (1 - Fw) * Ew * std::pow(Ew * s, Fw - 1);
      double denominator = 1 + std::pow(Ew * s, Fw);
      denominator = Fw * std::pow(denominator, 1 - 1 / Fw);
      G2 = numerator / denominator;
      break;
    }
    case RetentionModel::GALLIPOLI: {
      double Fi = retentionParameters[0];
      double psi = retentionParameters[1];
      double n = retentionParameters[2];
      double m = retentionParameters[3];
      double s = initialState.suction();
      double nu = initialState.specificVolume();
      double numerator = Fi * (1 - nu) * psi * s;
      numerator = Fi * (1 - nu) * psi * std::pow(numerator, n - 1);
      double denominator = numerator * s + 1;
      denominator = std::pow(denominator, m + 1);
      G2 = numerator / denominator;
      break;
    }
    default:
      std::cout << "Procedure Compute G2. Unknown Retention Model... No G2 "
                   "calculated.\n";
      break;
  }

  return G2;
}
