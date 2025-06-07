/*
 * The MIT License
 *
 * Copyright (c) 1997-2019 Center for the Simulation of Accidental Fires and
 * Explosions (CSAFE), and  Scientific Computing and Imaging Institute (SCI),
 * University of Utah.
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#include <CCA/Components/MPM/ConstitutiveModel/RockSoilModels/MohrCoulombSheng.h>

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
MohrCoulombSheng::MohrCoulombSheng()
  : MohrCoulombBase()
{
}

MohrCoulombSheng::MohrCoulombSheng(double G, double K, double cohesion,
                                   double phi, double psi, double pMin)
  : MohrCoulombBase(G, K, cohesion, phi, psi, pMin)
{
}

MohrCoulombSheng::MohrCoulombSheng(const MohrCoulombSheng* cm)
  : MohrCoulombBase(cm)
{
}

/**
 * Integration
 */
MohrCoulombState
MohrCoulombSheng::integrate(const Vector7& strainIncrement,
                            const MohrCoulombState& initialState)
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
  MohrCoulombState finalState;
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

      Vector7 plasticStrainInc;
      if (unloading) {

        if (dbg_unloading.active()) {
          dbg_unloading << "\n\n Elasto-Plastic unloading=" << unloading
                        << "\n";
        }

        // find elastic part
        Vector7 elasticStrainInc;
        std::tie(elasticStrainInc, plasticStrainInc) = 
          findIntersectionUnloading(strainIncrement, initialState);

        // calculating elastic part, updated in the same point.
        calcElastic(elasticStrainInc, initialState, finalState);

        // Update specific volume
        double spVol =
          computeNu(finalState.stress, finalState.state, finalState.suction());
        finalState.specificVolume(spVol);

      } else {
        
        plasticStrainInc = strainIncrement;

      }

      calculatePlastic(plasticStrainInc, finalState);

    } else {

      // not on the yield locus, finding intersection and subsequently
      // calculate elastic and plastic stress increment
      Vector7 elasticStrainInc, plasticStrainInc;
      std::tie(elasticStrainInc, plasticStrainInc) = 
        findIntersection(strainIncrement, initialState);

      calcElastic(elasticStrainInc, initialState, finalState);

      double spVol =
        computeNu(finalState.stress, finalState.state, finalState.suction());
      finalState.specificVolume(spVol);

      calculatePlastic(plasticStrainInc, finalState);
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

MohrCoulombState 
MohrCoulombSheng::integrate(const Vector7& strainIncrement,
                            const MohrCoulombState& initialState, 
                            RegionType& region)
{
  return initialState;
}

/**
 * Check yield
 *  check of the standard yield surface
 *  Note: value of function is normalised by the mean stress+2*cohesion!
 */
bool
MohrCoulombSheng::checkYieldNormalized(const MohrCoulombState& state) const
{
  double yieldFnValue = computeYieldNormalized(state.stress);
  if (yieldFnValue > (-d_int.d_yieldTol)) {
    return true;
  }
  return false;
}

double
MohrCoulombSheng::computeYieldNormalized(const Vector6& stress) const
{
  double p = firstInvariant(stress) / 3.0;
  double q = std::sqrt(3.0 * secondDevInvariant(stress));
  double J3 = thirdDevInvariant(stress);
  double M = (6.0 * d_yield.d_sin_phi) / (3.0 - d_yield.d_sin_phi);
  if (q > TINY) {
    double sin3thetabar = -27.0 * J3 / (2.0 * q * q * q);
    if (!std::isfinite(sin3thetabar)) {
      sin3thetabar = 1.0;
      std::cout
        << "**WARNING** sin3thetabar in checkYield Normalised not finite. set to 1\n"
           "  alpha4 = "
        << d_yield.d_alpha4 << " J3 = " << J3
        << " q = " << q << "\n";
    } else if (sin3thetabar > 1) {
      sin3thetabar = 1;
    } else if (sin3thetabar < -1) {
      sin3thetabar = -1;
    } 
    double factor = 1 + d_yield.d_alpha4 - (1 - d_yield.d_alpha4) * sin3thetabar;
    double scaledAlpha = d_yield.d_alpha / std::pow(0.5 * factor, 0.25);
    if (scaledAlpha > -1.0 && scaledAlpha < 1.0) {
      M *= scaledAlpha;
    }
  } else {
    M *= d_yield.d_alpha / std::pow(0.5 * (1 + d_yield.d_alpha4), 0.25);
  }

  double cohesion = d_yield.d_cohesion;
  double yieldFnValue = q / M - 2.0 * cohesion / M - p;

  // normalisation
  yieldFnValue /= (std::abs(p) + 2.0 * cohesion);

  if (dbg.active()) {
    std::cout << "Check Yield: Mean Stress = " << p
              << " Shear Stress=" << q << " cohesion = " << cohesion
              << " M = " << M << " Yield Function = " << yieldFnValue << "\n";
  }
  return yieldFnValue;
}

/**
 * Actually compute the gradient of the yield finction wrt stress
 */
std::tuple<Vector6, Vector6>
MohrCoulombSheng::computeDfDsigma(const Vector6& stress) const
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

  df_dSigma /= df_dSigma.norm();
  dg_dSigma /= dg_dSigma.norm();

  if (dbg.active()) {
    dbg << "df_dSigma = \n"
        << df_dSigma.transpose() << "dg_dsigma = \n"
        << dg_dSigma.transpose() << "\n";
  }

  return std::make_tuple(df_dSigma, dg_dSigma);
}

/**
  This procedure calculate stress increment using Modified Euler method
  A - matrix with coefficients,
  C - x coefficients.
  BRes - result coefficients,
  B - error estimate,
 */
std::tuple<double, int>
MohrCoulombSheng::plasticRKME221(MohrCoulombState& state,
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

/**
  This procedure calculates stress increment using 3-2 Runge - Kutta pair. The
  coeficients are given in  Ordinary and Partial Differential Equations
  routines in... by Lee & Schiesser, Chapman & Hall 2003
  The procedure consists of 3 stages and is 3th order accurate.
  A - matrix with coefficients, B - error estimate, BRes - result
  coefficients, C - x coefficients.
 */
std::tuple<double, int>
MohrCoulombSheng::plasticRK332(MohrCoulombState& state,
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

/**
  This procedure calculate stress increment using 3-2 Runge - Kutta pair. The
  coeficients are given in "A 3(2) Pair of Runge-Kutta Formulas"
  by P. Bogacki and L.F.  Shampine, Appl. Math. Lett., 2, pp. 321-325, 1989.
  The procedure consists of 3 stages and is 3th order accurate.

  A - matrix with coefficients, B - error estimate, BRes - result
  coefficients, C - x coefficients.
 */
std::tuple<double, int>
MohrCoulombSheng::plasticRKBog432(MohrCoulombState& state,
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

/**
 This procedure calculate stress increment using 4-3 Runge - Kutta pair. The
 coeficients are given in  Ordinary and Partial Differential Equations routines
 in... by Lee & Schiesser, Chapman & Hall 2003
 The procedure consists of 5 stages and is 4th order accurate.
*/
std::tuple<double, int>
MohrCoulombSheng::plasticRK543(MohrCoulombState& state,
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
MohrCoulombSheng::plasticRKEng654(MohrCoulombState& state,
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
MohrCoulombSheng::plasticRKCK654(MohrCoulombState& state,
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

/**
  This procedure calculate stress increment using Runge - Kutta pair as given by
  DORMAND - PRINCE. The used pair is known also as DOPRI5 or RK5(4)7FM. The
  procedure consists of 7 stages and should be less
  efficient than the RKErr8544 that uses Bogacki - Shimpine pair.
 */
std::tuple<double, int>
MohrCoulombSheng::plasticRKDP754(MohrCoulombState& state,
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


/**
  The procedure uses the embedded Runge - Kutta integration scheme with
  Adaptive StepSize Control the constants are as proposed by Bogacki and
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
MohrCoulombSheng::plasticRKErr8544(MohrCoulombState& state,
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


/**
 * Use table to integrate
 */
std::tuple<double, int>
MohrCoulombSheng::plasticExtrapol(MohrCoulombState& state,
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

  MohrCoulombState state_new = state;
  MohrCoulombState state_copy = state;
  MohrCoulombState state_old = state;

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
MohrCoulombSheng::doRungeKutta(const Eigen::Matrix<double, Steps, Steps>& AA,
                               const Eigen::Matrix<double, Steps, 1>& BB,
                               const Eigen::Matrix<double, Steps, 1>& BRes,
                               const Eigen::Matrix<double, Steps, 1>& CC,
                               MohrCoulombState& state, const Vector7& epStrain,
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

  std::vector<MohrCoulombState> midStates(Steps);

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
      std::tie(stressInc, p0StarInc) = calcPlastic(midStates[0], substepStrain);

      dSigma.col(0) = stressInc;
      plasticStrain.col(0) = substepStrain.block<6,1>(0,0);
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

      std::tie(stressInc, p0StarInc) = calcPlastic(midStates[rkloop], substepStrain);

      dSigma.col(rkloop) = stressInc;
      plasticStrain.col(rkloop) = substepStrain.block<6,1>(0,0);
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

      MohrCoulombState state_old = state;

      //******
      //******
      // here we update all the state data.
      state.update(plasticStrainInc, substepStrain, result.block<6, 1>(0, 0),
                   result(6));
      //******
      //******

      // Keeping a copy to check values
      MohrCoulombState trialState = state;

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
MohrCoulombSheng::doRungeKuttaErr(
  const Eigen::Matrix<double, Steps, Steps>& AA,
  const Eigen::Matrix<double, Steps, 1>& BB,
  const Eigen::Matrix<double, Steps, 1>& BRes,
  const Eigen::Matrix<double, Steps, 1>& CC,
  const Eigen::Matrix<double, Steps - 1, 1>& ErrCoef, MohrCoulombState& state,
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

  std::vector<MohrCoulombState> midStates(Steps);

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
      std::tie(stressInc, p0StarInc) = calcPlastic(midStates[0], substepStrain);

      dSigma.col(0) = stressInc;
      plasticStrain.col(0) = substepStrain.block<6,1>(0,0);
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

      std::tie(stressInc, p0StarInc) = calcPlastic(midStates[rkloop], substepStrain);

      dSigma.col(rkloop) = stressInc;
      plasticStrain.col(rkloop) = substepStrain.block<6,1>(0,0);
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

      MohrCoulombState state_old = state;

      //******
      //******
      // here we update all the state data.
      state.update(plasticStrainInc, substepStrain, result.block<6, 1>(0, 0),
                   result(6));
      //******
      //******

      // Keeping a copy to check values
      MohrCoulombState trialState = state;

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
  This procedure is to calculate the stress increase using the midpoint method
  with given number of iterations numIter.
  is made mainly for use in the extrapolation procedure. It has no adaptive
  substepping or error control. It just integrates
  the strain to get the stress in given number of substeps using the midpoint
  method
  */
double
MohrCoulombSheng::plasticMidpoint(MohrCoulombState& state,
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

  MohrCoulombState state_new, state_mid;
  Vector6 dSigma = Vector6::Zero();
  Vector6 plasticStrain = Vector6::Zero();
  double dP0Star = 0.0;
  for (int loop = 0; loop < numIter; loop++) {
    state_new = state;
    state_mid = state;

    std::tie(dSigma, dP0Star) = calcPlastic(state_new, halfCurrentStrain);
    state_mid.update(plasticStrain, halfCurrentStrain, dSigma, dP0Star);

    std::tie(dSigma, dP0Star) = calcPlastic(state_mid, halfCurrentStrain);
    dSigma *= 2.0;
    plasticStrain *= 2.0;
    dP0Star *= 2.0;
    state.update(plasticStrain, currentStrain, dSigma, dP0Star);

    absStress.block<6, 1>(0, 0) += dSigma.cwiseAbs();
    absStress(6) += std::abs(dP0Star);
  }
  return 0;
}

