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

#include "ClassicMohrCoulomb.h"

#include <Core/Exceptions/InvalidValue.h>
#include <Core/Util/DebugStream.h>

#include <chrono>
#include <cmath>
#include <iomanip>

using namespace Uintah;

static DebugStream dbg("ClassicMC", false);
static DebugStream dbg_unloading("ClassicMC_unloading", false);

/**
 * Constructors
 */
ClassicMohrCoulomb::ClassicMohrCoulomb()
  : ShengMohrCoulomb()
{
}

ClassicMohrCoulomb::ClassicMohrCoulomb(double G, double K, double cohesion,
                                   double phi, double psi)
  : ShengMohrCoulomb(G, K, cohesion, phi, psi)
{
}

/**
 Integrate
 */
StateMohrCoulomb 
ClassicMohrCoulomb::integrate(const Vector7& strainIncrement,
                              const StateMohrCoulomb& initialState)
{
  Vector7 purelyElasticStrain, purelyplasticStrain;

  if (dbg.active()) {
    dbg << "Strain increment:" << strainIncrement << "\n";
  }

  bool onTheYieldLocus = checkYieldNormalized(initialState);

  StateMohrCoulomb finalState;
  calcElastic(strainIncrement, initialState, finalState); 

  if (dbg.active()) {
    dbg << "Stress" << finalState.stress.transpose() << "\n";
  }

  bool elasticPlastic = checkYieldNormalized(finalState); 

  if (dbg.active()) {
    dbg << "p = " << finalState.meanStress()
        << " q = " << finalState.shearStress() << "  Yield? " << std::boolalpha
        << elasticPlastic << "\n";
  }

  bool unloading = false;
  if (elasticPlastic) {
    if (onTheYieldLocus) {
      
      unloading = checkGradient(initialState, finalState);

      if (unloading) {

        if (dbg_unloading.active()) {
          dbg_unloading << "\n\n Elasto-Plastic unloading=" << unloading
                        << "\n";
        }

        findIntersectionUnloading(strainIncrement, initialState,
                                  purelyElasticStrain, purelyplasticStrain);

        calcElastic(purelyElasticStrain, initialState, finalState);

        double spVol =
          computeNu(finalState.stress, finalState.state, finalState.suction());
        finalState.specificVolume(spVol);

        calculatePlastic(purelyplasticStrain, finalState);

      } else {

        calculatePlastic(strainIncrement, finalState);
      }

    } else {

      findIntersection(strainIncrement, initialState, purelyElasticStrain,
                       purelyplasticStrain);

      calcElastic(purelyElasticStrain, initialState, finalState);

      double spVol =
        computeNu(finalState.stress, finalState.state, finalState.suction());
      finalState.specificVolume(spVol);

      calculatePlastic(purelyplasticStrain, finalState);
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
 *  Note: value of function is normalised by the sum of eigenvalues
 */
bool
ClassicMohrCoulomb::checkYieldNormalized(const StateMohrCoulomb& state) const
{
  double yieldFnValue = computeYieldNormalized(state.stress);
  if (yieldFnValue > d_int.d_yieldTol) {
    return true;
  }
  return false;
}

double
ClassicMohrCoulomb::computeYieldNormalized(const Vector6& stress) const
{
  Vector eigenval;
  Matrix33 dummy;

  std::tie(eigenval, dummy) = getEigen(stress);
  double eig_min = eigenval(0);
  double eig_max = eigenval(2);

  double yieldFnValue = (eig_max - eig_min) - 
                        (eig_max + eig_min) * d_yield.sin_phi -
                        2.0 * d_yield.d_cohesion * d_yield.cos_phi;

  // Normalisation
  yieldFnValue /= (std::abs(eig_max) + std::abs(eig_min) + 2.0 * d_yield.d_cohesion);
  
  if (dbg.active()) {
    std::cout << "Check Yield: Mean Stress = " << state.meanStress()
              << " Shear Stress=" << state.shearStress 
              << " cohesion = " << d_yield.d_cohesion
              << " Yield Function = " << yieldFnValue << "\n";
  }
  return yieldFnValue;
}

/**
 * Eigenvalues are sorted from smallest to largest
 */
std::tuple<Vector3, Matrix33>
ClassicMohrCoulomb::getEigen(const Vector6& stress) const
{
  Matrix33 stressMat;
  stressMat(0, 0) = stress(0); 
  stressMat(0, 1) = stress(3); 
  stressMat(0, 2) = stress(4);
  stressMat(1, 0) = stressMat(0, 1); 
  stressMat(1, 1) = stress(1);
  stressMat(1, 2) = stress(5);
  stressMat(2, 0) = stressMat(0, 2); 
  stressMat(2, 1) = stressMat(1, 2); 
  stressMat(2, 2) = stress(2);
  
  Eigen::SelfAdjointEigenSolver<Matrix33> solver(stressMat);
  Vector3 eigenval = solver.eigenvalues();
  Matrix33 eigenvec = solver.eigenvectors();

  return std::make_tuple(eigenval, eigenvec);
}

Vector6
ClassicMohrCoulomb::rotateToOrigin(const Vector6& vec,
                                   const Matrix33& eigenVecs) 
{
  Matrix33 mat = toMatrix33(vec);
  Matrix33 rotMat = eigenVecs * mat * eigenVecs.transpose();
  Vector6 rotVec = toVector6(rotMat);

  return rotVec;
}

Vector6
ClassicMohrCoulomb::rotateToEigen(const Vector6& vec,
                                  const Matrix33& eigenVecs) 
{
  Matrix33 mat = toMatrix33(vec);
  Matrix33 rotMat = eigenVecs.transpose() * mat * eigenVecs;
  Vector6 rotVec = toVector6(rotMat);

  return rotVec;
}

Matrix33
ClassicMohrCoulomb::toMatrix33(const Vector6& vec) 
{
  Matrix33 mat = toMatrix33(vec);;
  mat(0, 0) = vec(0); 
  mat(0, 1) = vec(3); 
  mat(0, 2) = vec(4);
  mat(1, 0) = mat(0, 1); 
  mat(1, 1) = vec(1);
  mat(1, 2) = vec(5);
  mat(2, 0) = mat(0, 2); 
  mat(2, 1) = mat(1, 2); 
  mat(2, 2) = vec(2);
  return mat;
}

Vector6
ClassicMohrCoulomb::toVector6(const Matrix33& mat) 
{
  Vector6 vec;
  vec(0) = mat(0, 0);
  vec(1) = mat(1, 1);
  vec(2) = mat(2, 2);
  vec(3) = mat(0, 1);
  vec(4) = mat(0, 2);
  vec(5) = mat(1, 2);
  return vec;
}

/**
 * Actually compute the gradient of the yield finction wrt stress
 */
std::tuple<Vector6, Vector6>
ClassicMohrCoulomb::computeDfDsigma(const Vector6& stress) const
{
  Vector6 df_dsigma = Vector6::Zero();
  Vector6 dg_dsigma = Vector6::Zero();

  double meanStress = firstInvariant(stress)/3.0;
  if (meanStress > -d_yield.d_pMin) {

    double df_dsigma1 = 1.0 - d_yield.d_sin_phi;
    double df_dsigma3 = -1.0 - d_yield.d_sin_phi;
    double dg_dsigma1 = 1.0 - d_yield.d_sin_psi;
    double dg_dsigma3 = -1.0 - d_yield.d_sin_psi;

    df_dsigma(0) = df_dsigma1;
    df_dsigma(2) = df_dsigma3;

    dg_dsigma(0) = dg_dsigma1;
    dg_dsigma(2) = dg_dsigma3;

  } else {

    double one_third = 1.0 / 3.0;
    df_dsigma(0) = - one_third;
    df_dsigma(1) = - one_third;
    df_dsigma(2) = - one_third;

    dg_dsigma(0) = - one_third;
    dg_dsigma(1) = - one_third;
    dg_dsigma(2) = - one_third;
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
ClassicMohrCoulomb::plasticRKME221(StateMohrCoulomb& state,
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
    doRungeKuttaEig<Order, Steps>(A, B, BRes, C, state, epStrain, errorEstimate);
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
ClassicMohrCoulomb::plasticRK332(StateMohrCoulomb& state,
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
    doRungeKuttaEig<Order, Steps>(A, B, BRes, C, state, epStrain, errorEstimate);
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
ClassicMohrCoulomb::plasticRKBog432(StateMohrCoulomb& state,
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
    doRungeKuttaEig<Order, Steps>(A, B, BRes, C, state, epStrain, errorEstimate);
  return result;
}

/**
 This procedure calculate stress increment using 4-3 Runge - Kutta pair. The
 coeficients are given in  Ordinary and Partial Differential Equations routines
 in... by Lee & Schiesser, Chapman & Hall 2003
 The procedure consists of 5 stages and is 4th order accurate.
*/
std::tuple<double, int>
ClassicMohrCoulomb::plasticRK543(StateMohrCoulomb& state,
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
    doRungeKuttaEig<Order, Steps>(A, B, BRes, C, state, epStrain, errorEstimate);
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
ClassicMohrCoulomb::plasticRKEng654(StateMohrCoulomb& state,
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
    doRungeKuttaEig<Order, Steps>(A, B, BRes, C, state, epStrain, errorEstimate);
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
ClassicMohrCoulomb::plasticRKCK654(StateMohrCoulomb& state,
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
    doRungeKuttaEig<Order, Steps>(A, B, BRes, C, state, epStrain, errorEstimate);
  return result;
}

/**
  This procedure calculate stress increment using Runge - Kutta pair as given by
  DORMAND - PRINCE. The used pair is known also as DOPRI5 or RK5(4)7FM. The
  procedure consists of 7 stages and should be less
  efficient than the RKErr8544 that uses Bogacki - Shimpine pair.
 */
std::tuple<double, int>
ClassicMohrCoulomb::plasticRKDP754(StateMohrCoulomb& state,
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
    doRungeKuttaEig<Order, Steps>(A, B, BRes, C, state, epStrain, errorEstimate);
  return result;
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
ClassicMohrCoulomb::plasticRKErr8544(StateMohrCoulomb& state,
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
  auto result = doRungeKuttaEigErr<Order, Steps>(A, B, BRes, C, ErrorCoef, state,
                                                 epStrain, errorEstimate);
  return result;
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
        case DriftCorrection::CORRECTION_AT_END:
          Vector3 eigenVal;
          Matrix33 eigenVec;
          std::tie(eigenVal, eigenVec) = getEigen(state_new.stress);
          state_new.setStressEigen(eigenVal);
          state_new.strain.block<6,1>(0,0) =
            rotateToEigen(state_new.strain.block<6,1>(0,0), eigenVec);
          state_new.plasticStrain =
            rotateToEigen(midStates[0].plasticStrain, eigenVec);

          correctDriftEnd(state_new);

          state_new.stress = rotateToOrigin(state_new.stress, eigenVec);
          state_new.strain.block<6,1>(0,0) =
            rotateToOrigin(state_new.strain.block<6,1>(0,0), eigenVec);
          state_new.plasticStrain =
            rotateToOrigin(state_new.plasticStrain, eigenVec);

          break;
        case DriftCorrection::NO_CORRECTION:
        case DriftCorrection::CORRECTION_AT_BEGIN:
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
ClassicMohrCoulomb::doRungeKuttaEig(const Eigen::Matrix<double, Steps, Steps>& AA,
                                  const Eigen::Matrix<double, Steps, 1>& BB,
                                  const Eigen::Matrix<double, Steps, 1>& BRes,
                                  const Eigen::Matrix<double, Steps, 1>& CC,
                                  StateMohrCoulomb& state, const Vector7& epStrain,
                                  bool errorEstimate) const
{
  // Check eigenvalues and rotation are being computed correctly
  std::vector<StateMohrCoulomb> midStates(Steps);
  midStates[0] = state;

  Vector3 eigenValues;
  Matrix33 eigenVectors;
  std::tie(eigenValues, eigenVectors) = getEigen(midStates[0].stress);
  midStates[0].setStressEigen(eigenValues);
  midStates[0].stress = rotateToOrigin(midStates[0].stress, eigenVectors);

  if ((midStates[0].stress(0) - state.stress(0)) > 0.1 ||
      (midStates[0].stress(1) - state.stress(1)) > 0.1 ||
      (midStates[0].stress(2) - state.stress(2)) > 0.1) {
    std::cout << "**WARNING** "
                 "Problem with initial test rotation procedure in the Runge-Kutta "
                 "procedure. Stress tensor rotated back to origin is not equal to "
                 "the initial stress. Please investigate. \n "
  }

  Eigen::Matrix<double, 6, Steps> dSigma =
    Eigen::Matrix<double, 6, Steps>::Zero();
  Eigen::Matrix<double, 6, Steps> plasticStrain =
    Eigen::Matrix<double, 6, Steps>::Zero();
  Eigen::Matrix<double, 1, Steps> dP0Star =
    Eigen::Matrix<double, 1, Steps>::Zero();

  Vector7 absStress = Vector7::Zero();
  Vector7 reuseRes = Vector7::Zero();

  double newStepSize = 1.0;
  double stepLength = 1.0;
  double microState = 0.0;
  double totalSize = 0.0;
  double frequency = 100000.0 / Order; // how often display info about steps
  double methodPower = std::pow(2.0, Order) * d_int.d_integrationTol;
  double stepAccuracyCheck = 0.0;
  bool reuseStep = false;

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
    Vector7 currentStrain = substepStrain;

    double rError = 0;

    if (dbg.active()) {
      dbg << "Step Length = " << stepLength
          << " Current strain (0) = " << currentStrain(0) << "\n";
    }

    // Make copies of the initial state
    for (auto& midState : midStates) {
      midState = state;
    }

    if (reuseStep) {

      dSigma.col(0) = reuseRes.block<6, 1>(0, 0);
      plasticStrain.col(0) = plasticStrainInc.block<6, 1>(0, 0);
      dP0Star(0, 0) = reuseRes(6);

    } else {

      Vector3 eigenVal;
      Matrix33 eigenVec;
      std::tie(eigenVal, eigenVec) = getEigen(midStates[0].stress);

      midStates[0].setStressEigen(eigenVal);
      midStates[0].strain.block<6,1>(0,0) = 
        rotateToEigen(midStates[0].strain.block<6,1>(0,0), eigenVec);
      midStates[0].plasticStrain = 
        rotateToEigen(midStates[0].plasticStrain, eigenVec);
      substepStrain = rotateToEigen(substepStrain, eigenVec);
     
      Vector6 stressInc = Vector6::Zero();
      Vector6 plasticStrainInc = Vector6::Zero();
      calcPlastic(midStates[0], substepStrain, stressInc, plasticStrainInc,
                  p0StarInc);
  
      // 4:Rotate back the computed stress
      midStates[0].stress = 
        rotateToOrigin(midStates[0].stress, eigenVec);
      midStates[0].strain.block<6,1>(0,0) = 
        rotateToOrigin(midStates[0].strain.block<6,1>(0,0), eigenVec);
      midStates[0].plasticStrain = 
        rotateToOrigin(midStates[0].plasticStrain, eigenVec);
      substepStrain = rotateToOrigin(substepStrain, eigenVec);
      plasticStrainInc = rotateToOrigin(plasticStrainInc, eigenVec);
      stressInc = rotateToOrigin(stressInc, eigenVec);

      if (dbg.active()) {
        dbg << "Stress midState[" << 0 << "] = " << midStates[0].stress << "\n";
        dbg << "Stress increment = " << stressInc << "\n";
      }

      // **TODO** Check if plastic strain increment should not be computed
      // plasticStrainInc = Vector6::Zero();
      dSigma.col(0) = stressInc;
      plasticStrain.col(0) = plasticStrainInc;
      dP0Star(0, 0) = p0StarInc;
    }

    for (int rkloop = 1; rkloop < Steps; rkloop++) {

      Vector6 stressInc = Vector6::Zero();
      Vector6 plasticStrainInc = Vector6::Zero();
      double p0StarInc = 0;

      currentStrain = substepStrain * CC[rkloop];
     
      for (int i = 0; i < rkloop; i++) {
        stressInc += (dSigma.col(i) * AA(rkloop, i));
        plasticStrainInc += (plasticStrain.col(i) * AA(rkloop, i));
        p0StarInc += (dP0Star(0, i) * AA(rkloop, i));
      }

      midStates[rkloop].update(plasticStrainInc, currentStrain, stressInc,
                               p0StarInc);

      Vector3 eigenVal;
      Matrix33 eigenVec;
      std::tie(eigenVal, eigenVec) = getEigen(midStates[rkloop].stress);

      midStates[rkloop].setStressEigen(eigenVal);
      midStates[rkloop].strain.block<6,1>(0,0) = 
        rotateToEigen(midStates[0].strain.block<6,1>(0,0), eigenVec);
      midStates[rkloop].plasticStrain = 
        rotateToEigen(midStates[0].plasticStrain, eigenVec);
      substepStrain = rotateToEigen(substepStrain, eigenVec);
     
      if (dbg.active()) {
        dbg << "Stress midState[" << rkloop << "] = " << midStates[rkloop].stress << "\n";
        dbg << "Strain increment = " << substepStrain << "\n";
      }

      Vector6 stressInc = Vector6::Zero();
      Vector6 plasticStrainInc = Vector6::Zero();
      calcPlastic(midStates[rkloop], substepStrain, stressInc, plasticStrainInc,
                  p0StarInc);
  
      midStates[rkloop].stress = 
        rotateToOrigin(midStates[rkloop].stress, eigenVec);
      midStates[rkloop].strain.block<6,1>(0,0) = 
        rotateToOrigin(midStates[rkloop].strain.block<6,1>(0,0), eigenVec);
      midStates[rkloop].plasticStrain = 
        rotateToOrigin(midStates[rkloop].plasticStrain, eigenVec);
      substepStrain = rotateToOrigin(substepStrain, eigenVec);
      plasticStrainInc = rotateToOrigin(plasticStrainInc, eigenVec);
      stressInc = rotateToOrigin(stressInc, eigenVec);

      if (dbg.active()) {
        dbg << "Stress increment = " << stressInc << "\n";
      }

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

    if (!errorEstimate) {
      error -= result;
    }

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
        case DriftCorrection::CORRECTION_AT_END:
          {
            Vector3 eigenVal;
            Matrix33 eigenVec;
            std::tie(eigenVal, eigenVec) = getEigen(state.stress);
            state.setStressEigen(eigenVal);
            state.strain.block<6,1>(0,0) = 
              rotateToEigen(state.strain.block<6,1>(0,0), eigenVec);
            state.plasticStrain = 
              rotateToEigen(midStates[0].plasticStrain, eigenVec);

            correctDriftEnd(state);

            state.stress = rotateToOrigin(state.stress, eigenVec);
            state.strain.block<6,1>(0,0) = 
              rotateToOrigin(state.strain.block<6,1>(0,0), eigenVec);
            state.plasticStrain = 
              rotateToOrigin(state.plasticStrain, eigenVec);
          }
          break;
        case DriftCorrection::CORRECTION_AT_BEGIN:
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
          std::cout << "Strain inc :\n" << substepStrain << "\n";
        }

      } else {

        absStress.block<6,1>(0,0) = (state.stress - state_old.stress).cwiseAbs();
        absStress(6) = std::abs(state.p0Star - state_old.p0Star);

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
          std::cout << "Strain inc :\n" << substepStrain << "\n";
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
ClassicMohrCoulomb::doRungeKuttaEigErr(
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

      Vector3 eigenVal;
      Matrix33 eigenVec;
      std::tie(eigenVal, eigenVec) = getEigen(midStates[0].stress);

      midStates[0].setStressEigen(eigenVal);
      midStates[0].strain.block<6,1>(0,0) =
        rotateToEigen(midStates[0].strain.block<6,1>(0,0), eigenVec);
      midStates[0].plasticStrain =
        rotateToEigen(midStates[0].plasticStrain, eigenVec);
      substepStrain = rotateToEigen(substepStrain, eigenVec);

      stressInc = Vector6::Zero();
      plasticStrainInc = Vector6::Zero();
      calcPlastic(midStates[0], substepStrain, stressInc, plasticStrainInc,
                  p0StarInc);

      midStates[0].stress =
        rotateToOrigin(midStates[0].stress, eigenVec);
      midStates[0].strain.block<6,1>(0,0) =
        rotateToOrigin(midStates[0].strain.block<6,1>(0,0), eigenVec);
      midStates[0].plasticStrain =
        rotateToOrigin(midStates[0].plasticStrain, eigenVec);
      substepStrain = rotateToOrigin(substepStrain, eigenVec);
      plasticStrainInc = rotateToOrigin(plasticStrainInc, eigenVec);
      stressInc = rotateToOrigin(stressInc, eigenVec);


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

      Vector3 eigenVal;
      Matrix33 eigenVec;
      std::tie(eigenVal, eigenVec) = getEigen(midStates[rkloop].stress);

      midStates[rkloop].setStressEigen(eigenVal);
      midStates[rkloop].strain.block<6,1>(0,0) =
        rotateToEigen(midStates[0].strain.block<6,1>(0,0), eigenVec);
      midStates[rkloop].plasticStrain =
        rotateToEigen(midStates[0].plasticStrain, eigenVec);
      substepStrain = rotateToEigen(substepStrain, eigenVec);

      calcPlastic(midStates[rkloop], substepStrain, stressInc, plasticStrainInc,
                  p0StarInc);

      midStates[rkloop].stress =
        rotateToOrigin(midStates[rkloop].stress, eigenVec);
      midStates[rkloop].strain.block<6,1>(0,0) =
        rotateToOrigin(midStates[rkloop].strain.block<6,1>(0,0), eigenVec);
      midStates[rkloop].plasticStrain =
        rotateToOrigin(midStates[rkloop].plasticStrain, eigenVec);
      substepStrain = rotateToOrigin(substepStrain, eigenVec);
      plasticStrainInc = rotateToOrigin(plasticStrainInc, eigenVec);
      stressInc = rotateToOrigin(stressInc, eigenVec);

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
        case DriftCorrection::CORRECTION_AT_END:
          {
            Vector3 eigenVal;
            Matrix33 eigenVec;
            std::tie(eigenVal, eigenVec) = getEigen(state.stress);
            state.setStressEigen(eigenVal);
            state.strain.block<6,1>(0,0) =
              rotateToEigen(state.strain.block<6,1>(0,0), eigenVec);
            state.plasticStrain =
              rotateToEigen(midStates[0].plasticStrain, eigenVec);

            correctDriftEnd(state);
 
            state.stress = rotateToOrigin(state.stress, eigenVec);
            state.strain.block<6,1>(0,0) =
              rotateToOrigin(state.strain.block<6,1>(0,0), eigenVec);
            state.plasticStrain =
              rotateToOrigin(state.plasticStrain, eigenVec);
          }
          break;
        case DriftCorrection::CORRECTION_AT_BEGIN:
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

        absStress.block<6,1>(0,0) = (state.stress - state_old.stress).cwiseAbs();
        absStress(6) = std::abs(state.p0Star - state_old.p0Star);

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

double
ClassicMohrCoulomb::plasticMidpoint(StateMohrCoulomb& state,
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

    Vector3 eigenVal;
    Matrix33 eigenVec;
    std::tie(eigenVal, eigenVec) = getEigen(state_new.stress);
    state_new.setStressEigen(eigenVal);
    state_mid.setStressEigen(eigenVal);

    state_new.strain.block<6,1>(0,0) =
      rotateToEigen(state_new.strain.block<6,1>(0,0), eigenVec);
    state_mid.strain.block<6,1>(0,0) =
      rotateToEigen(state_new.strain.block<6,1>(0,0), eigenVec);
    state_new.plasticStrain =
      rotateToEigen(state_new.plasticStrain, eigenVec);
    state_mid.plasticStrain =
      rotateToEigen(state_new.plasticStrain, eigenVec);

    halfCurrentStrain = rotateToEigen(halfCurrentStrain, eigenVec);
    
    calcPlastic(state_new, halfCurrentStrain, dSigma, plasticStrain, dP0Star);

    state_mid.update(plasticStrain, halfCurrentStrain, dSigma, dP0Star);

    calcPlastic(state_mid, halfCurrentStrain, dSigma, plasticStrain, dP0Star);

    plasticStrain = rotateToOrigin(plasticStrain, eigenVec);
    dSigma = rotateToOrigin(dSigma, eigenVec);

    dSigma *= 2.0;
    plasticStrain *= 2.0;
    dP0Star *= 2.0;
    state.update(plasticStrain, currentStrain, dSigma, dP0Star);

    absStress.block<6, 1>(0, 0) += dSigma.cwiseAbs();
    absStress(6) += std::abs(dP0Star);
  }
  return 0;
}

void ClassicMohrCoulomb::IntegrateMCIClassic (double* StrainIncrement, BBMPoint* InitialPoint, int* Region)
{
// Modified implicit, should be accurate and equal to the explicit
// take into account that different surfaces may be active during integration

/*written in mind that InitialPoint may be equal to the FinalPoint
Basically the procedure made before should be re-introduced, with some changes
due to the different approach to parameters setting. It should also incorporate the
correct value of plastic strain, instead of useless lambda. This should have been done
long time ago, but well, better late then never.
The plastic/elastic part of the strain is also linked to the
*/


bool PurelyElastic;
//bool OnTheYieldLocus;
//bool Unloading;
//double PurelyElasticStrain[7],PurelyPlasticStrain[7], IntersectionStrain[7], NewIntersectionStrain[7];

BBMMatrix EigenVectors(3,3);
//double EigenValues [3];

//cout<<"Strains:"<<StrainIncrement[0]<<' '<<StrainIncrement[1]<<' '<<StrainIncrement[2]<<' '<<StrainIncrement[3]<<' '<<StrainIncrement[4]<<' '<<StrainIncrement[5]<<endl;



//InitialPoint->GetEigen(EigenValues,&EigenVectors);
//InitialPoint->SetStressEigen(EigenValues);
//EigenVectors.RotateToEigen (StrainIncrement,&EigenVectors);
CalcElastic (StrainIncrement, InitialPoint, InitialPoint);    //Initial Point copied to the final Point, later Final point updated with the elastic stress
//cout<<"Stress after elastic step:"<<FinalPoint.stress[0]<<" "<<FinalPoint.stress[1]<<" "<<FinalPoint.stress[2]<<" "<<FinalPoint.stress[3]<<" "<<FinalPoint.stress[4]<<" "<<FinalPoint.stress[5]<<endl;


PurelyElastic=!(CheckIfPlastic (InitialPoint)); //CheckIfPlastic returns true if in the plastic part of the YL. If on the YL (within Tolerance INT_TOL returns false)
//cout<<"p="<<FinalPoint.GetMeanStress()<<"  q="<<FinalPoint.GetShearStress()<<"  Yield? "<<PurelyElastic<<endl; getchar();
//PurelyElastic=TRUE;

if (PurelyElastic)
{
	//FinalPoint.Copy(InitialPoint);
	//update with what we have - all updated after the procedure
	//cout<<"Purely Elastic Step"<<endl;
}
else  // we have elasto-plastic step
{
	ReturnImplicit (InitialPoint, Region);
}
//for (int i=0; i<3; i++) cout<<StrainIncrement[i]<<"  "<<endl;
//for (int i=0; i<3; i++) cout<<FinalPoint.stress[i]<<"  "<<endl;

//FinalPoint.SetSpecificVolume(ComputeNu(FinalPoint.stress,FinalPoint.state, FinalPoint.GetSuction()));

//EigenVectors.RotateToOrigin(InitialPoint->stress,&EigenVectors);
//EigenVectors.RotateToOrigin(InitialPoint->strain,&EigenVectors);
//EigenVectors.RotateToOrigin(InitialPoint->plastic_strain,&EigenVectors);
//getchar();
//*/
}

void ClassicMohrCoulomb::ReturnImplicit (BBMPoint* Point, int* Region)
{

BBMMatrix EigenVectors(3,3);
double EigenValues [3];
for (int i=0; i<6; i++) Point->stress[i]=-Point->stress[i];

//cout<<"Point stress initial:"<<Point->stress[0]<<"  "<<Point->stress[1]<<"  "<<Point->stress[2]<<"  "<<Point->stress[3]<<"  "<<Point->stress[4]<<"  "<<Point->stress[5]<<"  "<<endl;
Point->GetEigen(EigenValues,&EigenVectors);
Point->SetStressEigen(EigenValues);
EigenVectors.RotateToEigen(Point->strain,&EigenVectors);
EigenVectors.RotateToEigen(Point->plastic_strain,&EigenVectors);

double k=(1+SinPhi)/(1-SinPhi);
double m=(1+SinPsi)/(1-SinPsi);

BBMMatrix DEL (3,3), NUMERATOR (1,1), DENOMINATOR (1,1);

double K43G=K+4.0*G/3.0;
double K23G=K-2.0*G/3.0;

DEL.PutElement (1,1,K43G);
DEL.PutElement (1,2,K23G);
DEL.PutElement (1,3,K23G); //rest of the line are zeros and rightly so
DEL.PutElement (2,1,K23G);
DEL.PutElement (2,2,K43G);
DEL.PutElement (2,3,K23G); //yes, the matrix is symmetrical, but it is faster to put this 3 additional elements
DEL.PutElement (3,1,K23G); //than just mirror all, including zeros, which are there already
DEL.PutElement (3,2,K23G);
DEL.PutElement (3,3,K43G);

BBMMatrix APEX (3,1),SIGMA (3,1), SIGMAAPEX (3,1);

for (int i=0; i<3; i++)
	{
	APEX.PutElementZ(i,0,2*Cohesion*sqrt(k)/(k-1));
	SIGMA.PutElementZ(i,0,Point->stress[i]);
	SIGMA.Substract(&APEX,&SIGMAAPEX);
	}

BBMMatrix RP1 (3,1), A1 (3,1), B1(3,1);

A1.PutElement(1,1,k);A1.PutElement(3,1,-1.0);
B1.PutElement(1,1,m);B1.PutElement(3,1,-1.0);

DEL.Multiply(&B1,&RP1);
A1.Transpose(&DENOMINATOR);
DENOMINATOR.Multiply(&DEL,&DENOMINATOR);
DENOMINATOR.Multiply(&B1,&DENOMINATOR);
RP1.Multiply(1.0/DENOMINATOR.GetElement(1,1),&RP1);
//RP1.Print();
//first update direction done

BBMMatrix NI_II (3,1), R1 (3,1);
//vector product
R1.PutElement(1,1,1.0);R1.PutElement(2,1,1.0);R1.PutElement(3,1,k);

RP1.VectorProduct(&R1,&NI_II);
//NI_II.Print();
NI_II.Transpose(&NUMERATOR);
NUMERATOR.Multiply(&SIGMAAPEX,&NUMERATOR);

double pI_II=NUMERATOR.GetElement(1,1);


BBMMatrix NI_III (3,1), R2 (3,1);
//vector product
R2.PutElement(1,1,1.0);R2.PutElement(2,1,k);R2.PutElement(3,1,k);

RP1.VectorProduct(&R2,&NI_III);
NI_III.Transpose(&NUMERATOR);
NUMERATOR.Multiply(&SIGMAAPEX,&NUMERATOR);

double pI_III=NUMERATOR.GetElement(1,1);


BBMMatrix RP2 (3,1), A2 (3,1), B2(3,1);

A2.PutElement(2,1,k);A2.PutElement(3,1,-1.0);
B2.PutElement(2,1,m);B2.PutElement(3,1,-1.0);

DEL.Multiply(&B2,&RP2);
A2.Transpose(&DENOMINATOR);
DENOMINATOR.Multiply(&DEL,&DENOMINATOR);
DENOMINATOR.Multiply(&B2,&DENOMINATOR);
RP2.Multiply(1.0/DENOMINATOR.GetElement(1,1),&RP2);
//RP2.Print();
BBMMatrix N2 (3,1);

RP1.VectorProduct(&RP2,&N2);

R1.Transpose(&DENOMINATOR);
DENOMINATOR.Multiply(&N2, &DENOMINATOR);
N2.Transpose(&NUMERATOR);
NUMERATOR.Multiply(&SIGMAAPEX,&NUMERATOR);
double t1=NUMERATOR.GetElement(1,1)/DENOMINATOR.GetElement(1,1);

//first line finished

BBMMatrix RP3 (3,1), A3 (3,1), B3(3,1);

A3.PutElement(1,1,k);A3.PutElement(2,1,-1.0);
B3.PutElement(1,1,m);B3.PutElement(2,1,-1.0);

DEL.Multiply(&B3,&RP3);
A3.Transpose(&DENOMINATOR);
DENOMINATOR.Multiply(&DEL,&DENOMINATOR);
DENOMINATOR.Multiply(&B3,&DENOMINATOR);
RP3.Multiply(1.0/DENOMINATOR.GetElement(1,1),&RP3);
//RP3.Print();
BBMMatrix N3 (3,1);

RP1.VectorProduct(&RP3,&N3);
//RP3.Print();

R2.Transpose(&DENOMINATOR);
//R2.Print();
DENOMINATOR.Multiply(&N3, &DENOMINATOR);
N3.Transpose(&NUMERATOR);

//SIGMAAPEX.Print();
NUMERATOR.Multiply(&SIGMAAPEX,&NUMERATOR);
double t2=NUMERATOR.GetElement(1,1)/DENOMINATOR.GetElement(1,1);
//cout<<"t2="<<t2<<" Denominator="<<DENOMINATOR.GetElement(1,1)<<" N3="<<endl;
//N3.Print();

//Alternative calculations for SinPhi=0 [i.e. apex--> inf] is needed



if ((t1>0.0)&& (t2>0.0)) {
	*Region=4;
	for (int i=0; i<3; i++) Point->stress[i]=APEX.GetElement(1,1);
	}
else if (pI_II<0){
	*Region=2;
	R1.Multiply(t1,&SIGMA);
	SIGMA.Add(&APEX,&SIGMA);
	for (int i=0; i<3; i++) Point->stress[i]=SIGMA.GetElementZ(i,0);
	}
else if (pI_III<=0) {
	*Region=1;
	double f=k*Point->stress[0]-Point->stress[2]-2*Cohesion*sqrt(k);
	RP1.Multiply(f,&NUMERATOR);
	SIGMA.Substract(&NUMERATOR,&SIGMA);
	for (int i=0; i<3; i++) Point->stress[i]=SIGMA.GetElementZ(i,0);
	}
else
	{
		*Region=3;
		R2.Multiply(t2,&SIGMA);
		//SIGMA.Print();
		//APEX.Print();
		SIGMA.Add(&APEX,&SIGMA);
		//SIGMA.Print();
		for (int i=0; i<3; i++) Point->stress[i]=SIGMA.GetElementZ(i,0);
	}

EigenVectors.RotateToOrigin(Point->stress,&EigenVectors);
//cout<<"Region:"<<*Region<<endl;
//cout<<"Point stress:"<<Point->stress[0]<<"  "<<Point->stress[1]<<"  "<<Point->stress[2]<<"  "<<Point->stress[3]<<"  "<<Point->stress[4]<<"  "<<Point->stress[5]<<"  "<<endl;
for (int i=0; i<6; i++) Point->stress[i]=-Point->stress[i];
EigenVectors.RotateToOrigin(Point->strain,&EigenVectors);
EigenVectors.RotateToOrigin(Point->plastic_strain,&EigenVectors);

}

