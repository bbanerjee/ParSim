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

#include <CCA/Components/MPM/ConstitutiveModel/Models/MohrCoulombClassic.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/TensorUtils.h>

#include <Core/Exceptions/InvalidValue.h>
#include <Core/Util/DebugStream.h>

#include <chrono>
#include <cmath>
#include <iomanip>

using namespace Uintah;

static DebugStream dbg_doing("ClassicMC_doing", false);
static DebugStream dbg("ClassicMC", false);
static DebugStream dbg_unloading("ClassicMC_unloading", false);
static DebugStream dbg_yield("ClassicMC_yield", false);
static DebugStream dbg_dfdsigma("ClassicMC_dfdsigma", false);
static DebugStream dbg_RK("ClassicMC_RK", false);

constexpr long64 testParticleID = 6619136;
/**
 * Constructors
 */
MohrCoulombClassic::MohrCoulombClassic()
  : MohrCoulombBase()
{
}

MohrCoulombClassic::MohrCoulombClassic(double G,
                                       double K,
                                       double cohesion,
                                       double phi,
                                       double psi,
                                       double pMin)
  : MohrCoulombBase(G, K, cohesion, phi, psi, pMin)
{
}

MohrCoulombClassic::MohrCoulombClassic(const MohrCoulombClassic* cm)
  : MohrCoulombBase(cm)
{
}

/**
 Integrate
 */
MohrCoulombState
MohrCoulombClassic::integrate(const Vector7& strainIncrement,
                              const MohrCoulombState& initialState)
{
  dbg_doing << "Doing MohrCoulombClassic::integrate\n";

  Vector7 purelyElasticStrain, purelyplasticStrain;

  dbg << "Strain increment:" << strainIncrement.transpose() << "\n";

  bool onTheYieldLocus = checkYieldNormalized(initialState);

  if (dbg_yield.active()) {
    if (initialState.particleID == testParticleID) {
      dbg_yield << __LINE__ << ":MCClassic::Before return: p = " << initialState.meanStress()
          << " q = " << initialState.shearStress() 
          << " Initial state yield? " << std::boolalpha << onTheYieldLocus << "\n";
    }
  }

  MohrCoulombState finalState;
  calcElastic(strainIncrement, initialState, finalState);

  dbg << "Stress elastic" << finalState.stress.transpose() << "\n";

  bool elasticPlastic = checkYieldNormalized(finalState);

  if (dbg_yield.active()) {
    if (initialState.particleID == testParticleID) {
      dbg_yield << "p = " << finalState.meanStress()
          << " q = " << finalState.shearStress() 
          << " Trial state yield? " << std::boolalpha << elasticPlastic << "\n";
    }
  }

  bool unloading = false;
  if (elasticPlastic) {
    if (onTheYieldLocus) {
      unloading = checkGradient(initialState, finalState);

      Vector7 plasticStrainInc;
      if (unloading) {
        dbg << "\n\n Elasto-Plastic unloading=" << unloading << "\n";

        Vector7 elasticStrainInc;
        std::tie(elasticStrainInc, plasticStrainInc) =
          findIntersectionUnloading(strainIncrement, initialState);

        calcElastic(elasticStrainInc, initialState, finalState);

        double spVol =
          computeNu(finalState.stress, finalState.state, finalState.suction());
        finalState.specificVolume(spVol);

      } else {

        plasticStrainInc = strainIncrement;
      }


      calculatePlastic(plasticStrainInc, finalState);

    } else {

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

  if (dbg_yield.active()) {

    onTheYieldLocus = checkYieldNormalized(finalState);

    if (onTheYieldLocus) {
      std::ostringstream err;
      err << "**WARNING** The stress state of particle " << initialState.particleID
          << " after return is on or outside the yield surface\n"
          << "The value of the yield function is: " 
          << computeYieldNormalized(finalState.stress) << ". "
          << "The yield tolerance is:" << d_int.d_yieldTol << "\n";
      dbg_yield << err.str();
    }

    if (initialState.particleID == testParticleID) {
      dbg_yield << __LINE__ << "MCClassic::After return: p = " << finalState.meanStress()
          << " q = " << finalState.shearStress() 
          << " Final state yield? " << std::boolalpha << onTheYieldLocus << "\n\n";
    }
  }

  dbg << "strain = " << strainIncrement.transpose() << "\n"
      << "stress = " << finalState.stress.transpose() << "\n";

  double spVol =
    computeNu(finalState.stress, finalState.state, finalState.suction());
  finalState.specificVolume(spVol);

  return finalState;
}

/**
 * Modified implicit, should be accurate and equal to the explicit
 * take into account that different surfaces may be active during integration
 */
MohrCoulombState
MohrCoulombClassic::integrate(const Vector7& strainIncrement,
                              const MohrCoulombState& initialState,
                              RegionType& region)
{
  dbg_doing << "Doing MohrCoulombClassic::integrate region\n";

  MohrCoulombState finalState;
  calcElastic(strainIncrement, initialState, finalState);

  dbg << "Stress after elastic step: " << finalState.stress.transpose() << "\n";

  bool elasticPlastic = checkYieldNormalized(finalState);
  if (elasticPlastic) {
    finalState = doReturnImplicit(finalState, region);
  }

  dbg << " Strain inc = " << strainIncrement.transpose() << "\n";
  dbg << " Final stress = " << finalState.stress.transpose() << "\n";

  return finalState;
}

/**
 * Check yield
 *  check of the standard yield surface
 *  Note: value of function is normalised by the sum of eigenvalues
 *  Yielded -> true, Elastic -> false
 */
bool
MohrCoulombClassic::checkYieldNormalized(const MohrCoulombState& state) const
{
  dbg_doing << "Doing MohrCoulombClassic::checkYieldNormalized\n";

  double yieldFnValue = computeYieldNormalized(state.stress);

  dbg << "Check Yield: Mean Stress = " << state.meanStress()
      << " Shear Stress=" << state.shearStress()
      << " cohesion = " << d_yield.d_cohesion
      << " Yield Function = " << yieldFnValue << "\n";

  if (yieldFnValue > d_int.d_yieldTol) {
    return true;
  }
  return false;
}

double
MohrCoulombClassic::computeYieldNormalized(const Vector6& stress) const
{
  dbg_doing << "Doing MohrCoulombClassic::computeYieldNormalized\n";

  double yieldFnValue = -1.0;
  double meanStress = firstInvariant(stress) / 3.0;
  if (meanStress > -d_yield.d_pMin) {

    Vector3 eigenval;
    Matrix33 dummy;

    std::tie(eigenval, dummy) = getEigen(stress);
    double sigma_1 = eigenval(0);
    double sigma_3 = eigenval(2);

    yieldFnValue = (sigma_1 - sigma_3) -
                   (sigma_1 + sigma_3) * d_yield.d_sin_phi -
                   2.0 * d_yield.d_cohesion * d_yield.d_cos_phi;

    // Normalisation
    yieldFnValue /=
      (std::abs(sigma_1) + std::abs(sigma_3) + 2.0 * d_yield.d_cohesion);

    if (!std::isfinite(yieldFnValue)) {
      std::cout << "**ERROR** Yield function is not finite for stress state\n"
                << stress.transpose() << "\n"
                << " with eigenvals = " << eigenval.transpose() << "\n";
    }

  } else {

    yieldFnValue = -d_yield.d_pMin - meanStress;

  }


  return yieldFnValue;
}

/**
 * Eigenvalues are sorted from largest to smallest
 */
std::tuple<Vector3, Matrix33>
MohrCoulombClassic::getEigen(const Vector6& stress) const
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
  Vector3 eigenval  = solver.eigenvalues();
  Matrix33 eigenvec = solver.eigenvectors();

  double temp_val = eigenval(0);
  eigenval(0)     = eigenval(2);
  eigenval(2)     = temp_val;

  auto temp_vec   = eigenvec.col(0).eval();
  eigenvec.col(0) = eigenvec.col(2);
  eigenvec.col(2) = temp_vec;

  return std::make_tuple(eigenval, eigenvec);
}

Vector6
MohrCoulombClassic::rotateToOrigin(const Vector6& vec,
                                   const Matrix33& eigenVecs) const
{
  Matrix33 mat    = toMatrix33(vec);
  Matrix33 rotMat = eigenVecs * mat * eigenVecs.transpose();
  Vector6 rotVec  = toVector6(rotMat);

  return rotVec;
}

Vector6
MohrCoulombClassic::rotateToEigen(const Vector6& vec,
                                  const Matrix33& eigenVecs) const
{
  Matrix33 mat    = toMatrix33(vec);
  Matrix33 rotMat = eigenVecs.transpose() * mat * eigenVecs;
  Vector6 rotVec  = toVector6(rotMat);

  return rotVec;
}

/**
 * Actually compute the gradient of the yield function wrt stress
 *
 * The Mohr-Coulomb yield function in invariant space has the form
 *   f = R(theta) * q - p * sin(phi) - c * cos(phi)
 * where
 *   R(theta) = 1/sqrt(3) * sin(theta) - 1/3 * cos(theta) * sin(phi)
 * with
 *   theta = theta_c + pi/3
 * Also
 *   p = I1/3, q = sqrt(3 * J2), cos(3 * theta_c) = (r/q)^3,
 *   r^3 = 27/2 J3
 * where
 *   I1 = tr(stress), J2 = 1/2 * tr(stress_dev:stress_dev), J3 = det(stress_dev)
 *   stress_dev = stress - I1/3 * I
 */
std::tuple<Vector6, Vector6>
MohrCoulombClassic::computeDfDsigma(const Vector6& stress_vec) const
{
  dbg_doing << "Doing MohrCoulombClassic::computeDfDsigma\n";

  // Constants from TensorUtils.h
  using namespace Vaango::Tensor;

  Vector6 df_dsigma_vec = Vector6::Zero();
  Vector6 dg_dsigma_vec = Vector6::Zero();
  double meanStress = third * firstInvariant(stress_vec);

  if (meanStress > -d_yield.d_pMin) {

    // Convert the 6-vector into a 3x3 symmetric matrix for easier to
    // check calculations
    Matrix3 stress = toMatrix3(stress_vec);

    dbg_dfdsigma << "stress = " << stress << "\n";

    // Compute invariants
    double I1 = stress.Trace();
    double p = I1 * third;
    Matrix3 stress_dev = stress - Identity * p;
    Matrix3 ss = stress_dev * stress_dev;
    double J2 = ss.Trace() * half;
    double sqrt_J2 = std::max(std::sqrt(J2), TINY);
    double q = sqrt_three * sqrt_J2;
    double J3 = stress_dev.Determinant();
    double r_cubed = 27.0 / 2.0 * J3;
    double cos3theta_c = r_cubed / (q * q * q);
    if (cos3theta_c < -1) {
      cos3theta_c = -1.0;
    } else if (cos3theta_c > 1) {
      cos3theta_c = 1.0;
    }
    double theta_c = std::acos(cos3theta_c) * third;
    double theta = theta_c + pi * third;

    dbg_dfdsigma << std::setprecision(std::numeric_limits<double>::digits10)
                 << "I1 = " << I1 << " p = " << p << " J2 = " << J2 << " q = " << q 
                 << " J3 = " << J3 << " r^3 = " << r_cubed << "\n"
                 << "cos(3theta_c) = " << cos3theta_c << " theta_c = " << theta_c
                 << " theta = " << theta << "\n";

    // Compute R and dR_dtheta
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    double R_phi = one_sqrt_three * sin_theta - third * cos_theta * d_yield.d_sin_phi;
    double R_psi = one_sqrt_three * sin_theta - third * cos_theta * d_potential.d_sin_psi;
    double dR_phi_dtheta = one_sqrt_three * cos_theta + third * sin_theta * d_yield.d_sin_phi;
    double dR_psi_dtheta = one_sqrt_three * cos_theta + third * sin_theta * d_potential.d_sin_psi;

    dbg_dfdsigma << std::setprecision(std::numeric_limits<double>::digits10)
                 << "sin(theta) = " << sin_theta << " cos(theta) = " << cos_theta
                 << " R_phi = " << R_phi << " R_psi = " << R_psi << "\n"
                 << "dR_phi_dtheta = " << dR_phi_dtheta << " dR_psi_dtheta = " << dR_psi_dtheta << "\n";

    // Compute dq_dsigma and dtheta_dsigma
    double sin3theta_c = std::max(std::sin(3.0 * theta_c), 1.0e-6);
    Matrix3 dq_dsigma = stress_dev * (sqrt_three * half / sqrt_J2);
    Matrix3 dJ3_dsigma = ss - Identity * (two_third * J2);
    double q_cubed = std::max(q * q * q, TINY);

    /*  // Tricky when sin(3theta_c) = 0
    double fac_1 = -9.0 / 2.0 / sin3theta_c / q_cubed;
    dbg_dfdsigma << std::setprecision(std::numeric_limits<double>::digits10)
                 << "fac_1 = " << fac_1 << "\n";
    */

    Matrix3 dtheta_dsigma = 
      (dJ3_dsigma - dq_dsigma * (3.0 * J3 / q)) * (- 9.0 * half / (sin3theta_c * q_cubed));

    dbg_dfdsigma << std::setprecision(std::numeric_limits<double>::digits10)
                 << "sin(3theta_c) = " << sin3theta_c << " q^3 = " << q_cubed << "\n"
                 << "dq_dsigma = " << dq_dsigma << "\n"
                 << "dJ3_dsigma = " << dJ3_dsigma << "\n"
                 << "dtheta_dsigma = " << dtheta_dsigma << "\n";
                 

    // Compute df_dsigma and dg_dsigma    
    Matrix3 dR_phi_dsigma = dR_phi_dtheta * dtheta_dsigma;
    Matrix3 dR_psi_dsigma = dR_psi_dtheta * dtheta_dsigma;

    dbg_dfdsigma << std::setprecision(std::numeric_limits<double>::digits10)
                 << "dR_phi_dsigma = " << dR_phi_dsigma << "\n"
                 << "dR_phi_dsigma * q = " <<  dR_phi_dsigma * q << "\n"
                 << "R_phi * dq_dsigma = " <<  R_phi * dq_dsigma << "\n"
                 << "dp_dsigma * sin_phi = " <<  Identity * third * d_yield.d_sin_phi << "\n"; 

    Matrix3 df_dsigma = dR_phi_dsigma * q + dq_dsigma * R_phi - Identity * (third * d_yield.d_sin_phi);
    Matrix3 dg_dsigma = dR_psi_dsigma * q + dq_dsigma * R_psi - Identity * (third * d_potential.d_sin_psi);

    // Convert to vectors
    df_dsigma_vec = toVector6(df_dsigma);
    dg_dsigma_vec = toVector6(dg_dsigma);

  } else {
    double one_third = 1.0 / 3.0;

    df_dsigma_vec(0)     = -one_third;
    df_dsigma_vec(1)     = -one_third;
    df_dsigma_vec(2)     = -one_third;

    dg_dsigma_vec(0) = -one_third;
    dg_dsigma_vec(1) = -one_third;
    dg_dsigma_vec(2) = -one_third;
  }

  df_dsigma_vec /= df_dsigma_vec.norm();
  dg_dsigma_vec /= dg_dsigma_vec.norm();

  dbg_dfdsigma << std::setprecision(std::numeric_limits<double>::digits10)
               << "df_dsigma = " << df_dsigma_vec.transpose() << "\n"
               << "dg_dsigma = " << dg_dsigma_vec.transpose() << "\n";
              
  return std::make_tuple(df_dsigma_vec, dg_dsigma_vec);
}

/*
std::tuple<Vector6, Vector6>
MohrCoulombClassic::computeDfDsigmaInEigenSpace(const Vector6& stress) const
{
  dbg_doing << "Doing MohrCoulombClassic::computeDfDsigma\n";

  Vector6 df_dsigma = Vector6::Zero();
  Vector6 dg_dsigma = Vector6::Zero();

  double meanStress = firstInvariant(stress) / 3.0;
  if (meanStress > -d_yield.d_pMin) {

    double df_dsigma1 = 1.0 - d_yield.d_sin_phi;
    double df_dsigma3 = -1.0 - d_yield.d_sin_phi;
    double dg_dsigma1 = 1.0 - d_potential.d_sin_psi;
    double dg_dsigma3 = -1.0 - d_potential.d_sin_psi;

    df_dsigma(0) = df_dsigma1;
    df_dsigma(2) = df_dsigma3;

    dg_dsigma(0) = dg_dsigma1;
    dg_dsigma(2) = dg_dsigma3;

  } else {
    double one_third = 1.0 / 3.0;

    df_dsigma(0)     = -one_third;
    df_dsigma(1)     = -one_third;
    df_dsigma(2)     = -one_third;

    dg_dsigma(0) = -one_third;
    dg_dsigma(1) = -one_third;
    dg_dsigma(2) = -one_third;
  }

  df_dsigma /= df_dsigma.norm();
  dg_dsigma /= dg_dsigma.norm();

  return std::make_tuple(df_dsigma, dg_dsigma);
}
*/

/**
  This procedure calculate stress increment using Modified Euler method
  A - matrix with coefficients,
  C - x coefficients.
  BRes - result coefficients,
  B - error estimate,
 */
std::tuple<double, int>
MohrCoulombClassic::plasticRKME221(MohrCoulombState& state,
                                   const Vector7& epStrain) const

{
  constexpr int Order = 2;
  constexpr int Steps = 2;

  using RKMatrix = Eigen::Matrix<double, Steps, Steps>;
  using RKVector = Eigen::Matrix<double, Steps, 1>;

  RKMatrix A    = RKMatrix::Zero();
  RKVector C    = RKVector::Zero();
  RKVector BRes = RKVector::Zero();
  RKVector B    = RKVector::Zero();

  getParamRKME221(A, B, BRes, C);

  bool errorEstimate = false;
  auto result        = doRungeKuttaEig<Order, Steps>(
    A, B, BRes, C, state, epStrain, errorEstimate);
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
MohrCoulombClassic::plasticRK332(MohrCoulombState& state,
                                 const Vector7& epStrain) const

{
  constexpr int Steps = 3;
  constexpr int Order = 3;

  using RKMatrix = Eigen::Matrix<double, Steps, Steps>;
  using RKVector = Eigen::Matrix<double, Steps, 1>;

  RKMatrix A    = RKMatrix::Zero();
  RKVector C    = RKVector::Zero();
  RKVector BRes = RKVector::Zero();
  RKVector B    = RKVector::Zero();

  getParamRK332(A, B, BRes, C);

  bool errorEstimate = false;
  auto result        = doRungeKuttaEig<Order, Steps>(
    A, B, BRes, C, state, epStrain, errorEstimate);
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
MohrCoulombClassic::plasticRKBog432(MohrCoulombState& state,
                                    const Vector7& epStrain) const
{
  constexpr int Steps = 4;
  constexpr int Order = 3;

  using RKMatrix = Eigen::Matrix<double, Steps, Steps>;
  using RKVector = Eigen::Matrix<double, Steps, 1>;

  RKMatrix A    = RKMatrix::Zero();
  RKVector C    = RKVector::Zero();
  RKVector BRes = RKVector::Zero();
  RKVector B    = RKVector::Zero();

  getParamRKBog432(A, B, BRes, C);

  bool errorEstimate = false;
  auto result        = doRungeKuttaEig<Order, Steps>(
    A, B, BRes, C, state, epStrain, errorEstimate);
  return result;
}

/**
 This procedure calculate stress increment using 4-3 Runge - Kutta pair. The
 coeficients are given in  Ordinary and Partial Differential Equations routines
 in... by Lee & Schiesser, Chapman & Hall 2003
 The procedure consists of 5 stages and is 4th order accurate.
*/
std::tuple<double, int>
MohrCoulombClassic::plasticRK543(MohrCoulombState& state,
                                 const Vector7& epStrain) const
{
  constexpr int Steps = 5;
  constexpr int Order = 4;

  using RKMatrix = Eigen::Matrix<double, Steps, Steps>;
  using RKVector = Eigen::Matrix<double, Steps, 1>;

  RKMatrix A    = RKMatrix::Zero();
  RKVector C    = RKVector::Zero();
  RKVector BRes = RKVector::Zero();
  RKVector B    = RKVector::Zero();

  getParamRK543(A, B, BRes, C);

  bool errorEstimate = false;
  auto result        = doRungeKuttaEig<Order, Steps>(
    A, B, BRes, C, state, epStrain, errorEstimate);
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
MohrCoulombClassic::plasticRKEng654(MohrCoulombState& state,
                                    const Vector7& epStrain) const
{
  constexpr int Steps = 6;
  constexpr int Order = 5;

  using RKMatrix = Eigen::Matrix<double, Steps, Steps>;
  using RKVector = Eigen::Matrix<double, Steps, 1>;

  RKMatrix A    = RKMatrix::Zero();
  RKVector C    = RKVector::Zero();
  RKVector BRes = RKVector::Zero();
  RKVector B    = RKVector::Zero();

  getParamRKEng654(A, B, BRes, C);

  bool errorEstimate = false;
  auto result        = doRungeKuttaEig<Order, Steps>(
    A, B, BRes, C, state, epStrain, errorEstimate);
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
MohrCoulombClassic::plasticRKCK654(MohrCoulombState& state,
                                   const Vector7& epStrain) const
{
  constexpr int Steps = 6;
  constexpr int Order = 5;

  using RKMatrix = Eigen::Matrix<double, Steps, Steps>;
  using RKVector = Eigen::Matrix<double, Steps, 1>;

  RKMatrix A    = RKMatrix::Zero();
  RKVector C    = RKVector::Zero();
  RKVector BRes = RKVector::Zero();
  RKVector B    = RKVector::Zero();

  getParamRKCK654(A, B, BRes, C);

  bool errorEstimate = false;
  auto result        = doRungeKuttaEig<Order, Steps>(
    A, B, BRes, C, state, epStrain, errorEstimate);
  return result;
}

/**
  This procedure calculate stress increment using Runge - Kutta pair as given by
  DORMAND - PRINCE. The used pair is known also as DOPRI5 or RK5(4)7FM. The
  procedure consists of 7 stages and should be less
  efficient than the RKErr8544 that uses Bogacki - Shimpine pair.
 */
std::tuple<double, int>
MohrCoulombClassic::plasticRKDP754(MohrCoulombState& state,
                                   const Vector7& epStrain) const
{
  constexpr int Steps = 7;
  constexpr int Order = 5;

  using RKMatrix = Eigen::Matrix<double, Steps, Steps>;
  using RKVector = Eigen::Matrix<double, Steps, 1>;

  RKMatrix A    = RKMatrix::Zero();
  RKVector C    = RKVector::Zero();
  RKVector BRes = RKVector::Zero();
  RKVector B    = RKVector::Zero();

  getParamRKDP754(A, B, BRes, C);

  bool errorEstimate = false;
  auto result        = doRungeKuttaEig<Order, Steps>(
    A, B, BRes, C, state, epStrain, errorEstimate);
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
MohrCoulombClassic::plasticRKErr8544(MohrCoulombState& state,
                                     const Vector7& epStrain) const
{
  constexpr int Steps = 8;
  constexpr int Order = 5;

  using RKMatrix    = Eigen::Matrix<double, Steps, Steps>;
  using RKVector    = Eigen::Matrix<double, Steps, 1>;
  using RKVectorErr = Eigen::Matrix<double, Steps - 1, 1>;

  RKMatrix A            = RKMatrix::Zero();
  RKVector C            = RKVector::Zero();
  RKVector BRes         = RKVector::Zero();
  RKVector B            = RKVector::Zero();
  RKVectorErr ErrorCoef = RKVectorErr::Zero();

  getParamRKErr8544(A, B, BRes, C, ErrorCoef);

  bool errorEstimate = false;
  auto result        = doRungeKuttaEigErr<Order, Steps>(
    A, B, BRes, C, ErrorCoef, state, epStrain, errorEstimate);
  return result;
}

/**
 * Use table to integrate
 */
std::tuple<double, int>
MohrCoulombClassic::plasticExtrapol(MohrCoulombState& state,
                                    const Vector7& epStrain) const
{
  // six stresses, p0*, six absolute stresses, absolute p0*,
  // 6 plastic strains
  constexpr int TableRows = 20;
  constexpr int StepMax   = 15;
  using ApproxTableMatrix = Eigen::Matrix<double, TableRows, StepMax + 1>;
  using HTableMatrix      = Eigen::Matrix<double, StepMax, StepMax + 1>;

  Vector6 Zero = Vector6::Zero();

  std::chrono::time_point<std::chrono::system_clock> startTime, endTime;
  startTime = std::chrono::system_clock::now();

  HTableMatrix hSquareTable = HTableMatrix::Zero();
  for (int i = 1; i < StepMax + 1; i++) {
    for (int j = 0; j < i; j++) {
      hSquareTable(i, j) = (DivInt[i] / DivInt[j]) * (DivInt[i] / DivInt[j]);
      dbg << "DivInt[" << i << "] = " << DivInt[i] << " DivInt[" << j
          << "] = " << DivInt[j] << "\n";
      dbg << "hSquareTable[" << i << ", " << j << "]=" << hSquareTable(i, j)
          << "\n";
    }
  }

  MohrCoulombState state_new  = state;
  MohrCoulombState state_copy = state;
  MohrCoulombState state_old  = state;

  Vector6 initPlasticStrain = state.plasticStrain;
  Vector7 initStress        = Vector7::Zero();
  initStress.block<6, 1>(0, 0) = state.stress;
  initStress(6) = state.p0Star();

  ApproxTableMatrix approxTable    = ApproxTableMatrix::Zero();
  ApproxTableMatrix approxTableOld = ApproxTableMatrix::Zero();

  double dP0Star        = 0.0;
  Vector6 plasticStrain = Vector6::Zero();
  Vector7 dSigma        = Vector7::Zero();
  Vector7 absStress     = Vector7::Zero();
  Vector7 dError        = Vector7::Zero();

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
    curCol.block<7, 1>(7, 0)  = absStress;
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
        dbg << approxTable.col(i).transpose() << "\n";
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
      dSigma        = currData.block<7, 1>(0, 0);
      dP0Star       = currData(6);
      dError        = currData.block<7, 1>(0, 0) - prevData.block<7, 1>(0, 0);

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

    dbg << "Relative error after iteration " << loop
        << " is equal to: " << rError << "\n";

    approxTableOld = approxTable;

    // error less then requested...check for drift and add the error
    if (rError < d_int.d_integrationTol) {
      auto currData = approxTable.col(loop);
      dSigma        = currData.block<7, 1>(0, 0);
      dP0Star       = currData(6);

      state_copy = state_old;
      state_old.update(Zero, epStrain, dSigma.block<6, 1>(0, 0), dP0Star);
      state_old = state_new;

      switch (d_int.d_driftCorrection) {
        case DriftCorrection::CORRECTION_AT_END: {

          /*
          Vector3 eigenVal;
          Matrix33 eigenVec;
          std::tie(eigenVal, eigenVec) = getEigen(state_new.stress);
          state_new.setStressEigen(eigenVal);
          state_new.strain.block<6, 1>(0, 0) =
            rotateToEigen(state_new.strain.block<6, 1>(0, 0), eigenVec);
          state_new.plasticStrain =
            rotateToEigen(state_new.plasticStrain, eigenVec);
          */

          correctDriftEnd(state_new);

          /*
          state_new.stress = rotateToOrigin(state_new.stress, eigenVec);
          state_new.strain.block<6, 1>(0, 0) =
            rotateToOrigin(state_new.strain.block<6, 1>(0, 0), eigenVec);
          state_new.plasticStrain =
            rotateToOrigin(state_new.plasticStrain, eigenVec);
          */

        } break;
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
    absStress     = currData.block<7, 1>(7, 0);
    plasticStrain = currData.block<6, 1>(14, 0);
    dP0Star       = currData(6);
    state.update(plasticStrain, epStrain, dSigma.block<6, 1>(0, 0), dP0Star);

    endTime = std::chrono::system_clock::now();

    for (int i = 0; i < loop + 1; i++) {
      numIter += DivInt[i];
    }

    dbg << "Procedure has coverged after " << numIter << " iterations.\n";
  } else {
    loop--;

    auto currData = approxTable.col(loop);
    dSigma        = currData.block<7, 1>(0, 0);
    absStress     = currData.block<7, 1>(7, 0);
    plasticStrain = currData.block<6, 1>(14, 0);
    dP0Star       = currData(6);

    state.update(plasticStrain, epStrain, dSigma.block<6, 1>(0, 0), dP0Star);

    endTime = std::chrono::system_clock::now();

    for (int i = 0; i < loop + 1; i++) {
      numIter += DivInt[i];
      dbg << "Procedure has NOT CONVERGED after " << numIter << "iterations.\n";
    }
  }

  double time = std::chrono::duration<double>(endTime - startTime).count();
  dbg << "Calculation took: " << time << "s."
      << "\n";

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
MohrCoulombClassic::doRungeKuttaEig(
  const Eigen::Matrix<double, Steps, Steps>& AA,
  const Eigen::Matrix<double, Steps, 1>& BB,
  const Eigen::Matrix<double, Steps, 1>& BRes,
  const Eigen::Matrix<double, Steps, 1>& CC,
  MohrCoulombState& state,
  const Vector7& epStrain,
  bool errorEstimate) const
{
  dbg_doing << "Doing MohrCoulombClassic::doRungeKuttaEig\n";

  std::vector<MohrCoulombState> midStates(Steps);
  midStates[0] = state;

  Eigen::Matrix<double, 6, Steps> dSigma =
    Eigen::Matrix<double, 6, Steps>::Zero();
  Eigen::Matrix<double, 6, Steps> plasticStrain =
    Eigen::Matrix<double, 6, Steps>::Zero();
  Eigen::Matrix<double, 1, Steps> dP0Star =
    Eigen::Matrix<double, 1, Steps>::Zero();

  Vector7 absStress = Vector7::Zero();
  Vector7 reuseRes  = Vector7::Zero();

  double newStepSize = 1.0;
  double stepLength  = 1.0;
  double microStep   = 0.0;
  double totalSize   = 0.0;
  double frequency   = 100000.0 / Order; // how often display info about steps
  double methodPower = std::pow(2.0, Order) * d_int.d_integrationTol;
  double stepAccuracyCheck = 0.0;
  bool reuseStep           = false;

  Vector6 stressInc        = Vector6::Zero();
  Vector6 plasticStrainInc = Vector6::Zero();
  double p0StarInc         = 0.0;

  std::chrono::time_point<std::chrono::system_clock> startTime, endTime;
  startTime = std::chrono::system_clock::now();

  int stepNo    = 0;
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

    dbg_RK << "Step no. = " << stepNo << " Step Length = " << stepLength
           << " Current strain (0) = " << currentStrain(0) << "\n";

    // Make copies of the initial state
    for (auto& midState : midStates) {
      midState = state;
    }

    if (reuseStep) {
      dSigma.col(0)        = reuseRes.block<6, 1>(0, 0);
      plasticStrain.col(0) = plasticStrainInc.block<6, 1>(0, 0);
      dP0Star(0, 0) = reuseRes(6);

    } else {

      /*
      Vector3 eigenVal;
      Matrix33 eigenVec;
      std::tie(eigenVal, eigenVec) = getEigen(midStates[0].stress);

      midStates[0].setStressEigen(eigenVal);
      midStates[0].strain.block<6, 1>(0, 0) =
        rotateToEigen(midStates[0].strain.block<6, 1>(0, 0), eigenVec);
      midStates[0].plasticStrain =
        rotateToEigen(midStates[0].plasticStrain, eigenVec);
      substepStrain.block<6, 1>(0, 0) =
        rotateToEigen(substepStrain.block<6, 1>(0, 0), eigenVec);
      */

      std::tie(stressInc, p0StarInc) = calcPlastic(midStates[0], substepStrain);

      /*
      midStates[0].stress = rotateToOrigin(midStates[0].stress, eigenVec);
      midStates[0].strain.block<6, 1>(0, 0) =
        rotateToOrigin(midStates[0].strain.block<6, 1>(0, 0), eigenVec);
      midStates[0].plasticStrain =
        rotateToOrigin(midStates[0].plasticStrain, eigenVec);
      substepStrain.block<6, 1>(0, 0) =
        rotateToOrigin(substepStrain.block<6, 1>(0, 0), eigenVec);
      plasticStrainInc = rotateToOrigin(plasticStrainInc, eigenVec);
      stressInc        = rotateToOrigin(stressInc, eigenVec);
      */

      dbg_RK << "Stress midState[" << 0 << "] = " << midStates[0].stress.transpose() << "\n";
      dbg_RK << "Stress increment = " << stressInc.transpose() << "\n";

      dSigma.col(0)        = stressInc;
      plasticStrain.col(0) = substepStrain.block<6,1>(0,0);
      dP0Star(0, 0) = p0StarInc;
    }

    for (int rkloop = 1; rkloop < Steps; rkloop++) {
      stressInc        = Vector6::Zero();
      plasticStrainInc = Vector6::Zero();
      p0StarInc        = 0;

      currentStrain = substepStrain * CC[rkloop];

      for (int i = 0; i < rkloop; i++) {
        stressInc += (dSigma.col(i) * AA(rkloop, i));
        plasticStrainInc += (plasticStrain.col(i) * AA(rkloop, i));
        p0StarInc += (dP0Star(0, i) * AA(rkloop, i));
      }

      dbg_RK << "RK loop no. = " << rkloop 
             << " Current strain = " << currentStrain.transpose() << "\n";
      std::cout << "RK loop no. = " << rkloop 
                << " Current strain = " << currentStrain.transpose() << "\n";

      midStates[rkloop].update(
        plasticStrainInc, currentStrain, stressInc, p0StarInc);

      /*
      Vector3 eigenVal;
      Matrix33 eigenVec;
      std::tie(eigenVal, eigenVec) = getEigen(midStates[rkloop].stress);

      midStates[rkloop].setStressEigen(eigenVal);
      midStates[rkloop].strain.block<6, 1>(0, 0) =
        rotateToEigen(midStates[0].strain.block<6, 1>(0, 0), eigenVec);
      midStates[rkloop].plasticStrain =
        rotateToEigen(midStates[0].plasticStrain, eigenVec);
      substepStrain.block<6, 1>(0, 0) =
        rotateToEigen(substepStrain.block<6, 1>(0, 0), eigenVec);
      */

      dbg_RK << "Stress midState[" << rkloop << "] = " << midStates[rkloop].stress.transpose()
          << "\n";
      dbg_RK << "Strain increment = " << substepStrain.transpose() << "\n";

      Vector6 stressInc        = Vector6::Zero();
      std::tie(stressInc, p0StarInc) = calcPlastic(midStates[rkloop], substepStrain);

      /*
      midStates[rkloop].stress =
        rotateToOrigin(midStates[rkloop].stress, eigenVec);
      midStates[rkloop].strain.block<6, 1>(0, 0) =
        rotateToOrigin(midStates[rkloop].strain.block<6, 1>(0, 0), eigenVec);
      midStates[rkloop].plasticStrain =
        rotateToOrigin(midStates[rkloop].plasticStrain, eigenVec);
      substepStrain.block<6, 1>(0, 0) =
        rotateToOrigin(substepStrain.block<6, 1>(0, 0), eigenVec);
      plasticStrainInc = rotateToOrigin(plasticStrainInc, eigenVec);
      stressInc        = rotateToOrigin(stressInc, eigenVec);
      */

      dbg_RK << "Stress increment = " << stressInc.transpose() << "\n";

      dSigma.col(rkloop)        = stressInc;
      plasticStrain.col(rkloop) = substepStrain.block<6,1>(0,0);
      dP0Star(0, rkloop) = p0StarInc;
    }

    Vector6 sigBRes = dSigma * BRes;
    // Vector6 epsPBRes = plasticStrain * BRes;
    double p0BRes = dP0Star * BRes;
    Vector6 sigB  = dSigma * BB;
    double p0B    = dP0Star * BB;

    Vector7 stressAtYield, error;
    stressAtYield.block<6, 1>(0, 0) = sigBRes;
    error.block<6, 1>(0, 0)  = sigB;
    stressAtYield(6) = p0BRes;
    error(6)  = p0B;

    if (!errorEstimate) {
      error -= stressAtYield;
    }

    switch (d_int.d_tolMethod) {
      case ToleranceMethod::EPUS_RELATIVE_ERROR:
        rError = checkNorm(stressAtYield, stressAtYield(6), state, error);
        break;
      case ToleranceMethod::SLOAN:
        rError = checkNormSloan(stressAtYield, stressAtYield(6), state, error);
        break;
      default:
        std::ostringstream err;
        err << "**ERROR** Improper d_tolMethod in increment.dta"
            << "\n";
        throw InvalidValue(err.str(), __FILE__, __LINE__);
    }

    for (int i = 0; i < 7; i++) {
      if (!std::isfinite(stressAtYield(i))) {
        std::ostringstream err;
        err << "**ERROR** Results not a number. Try correcting issue, results "
               "may be "
               "incorrect."
            << " Results = \n"
            << stressAtYield << "\n";
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
      state.update(
        plasticStrainInc, substepStrain, stressAtYield.block<6, 1>(0, 0), stressAtYield(6));
      //******
      //******

      // Keeping a copy to check values
      MohrCoulombState trialState = state;

      switch (d_int.d_driftCorrection) {
        case DriftCorrection::CORRECTION_AT_END: {

          /*
          Vector3 eigenVal;
          Matrix33 eigenVec;
          std::tie(eigenVal, eigenVec) = getEigen(state.stress);
          state.setStressEigen(eigenVal);
          state.strain.block<6, 1>(0, 0) =
            rotateToEigen(state.strain.block<6, 1>(0, 0), eigenVec);
          state.plasticStrain =
            rotateToEigen(midStates[0].plasticStrain, eigenVec);
          */

          correctDriftEnd(state);

          /*
          state.stress = rotateToOrigin(state.stress, eigenVec);
          state.strain.block<6, 1>(0, 0) =
            rotateToOrigin(state.strain.block<6, 1>(0, 0), eigenVec);
          state.plasticStrain = rotateToOrigin(state.plasticStrain, eigenVec);
          */
        } break;
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
          rError = checkNorm(stressAtYield, stressAtYield(6), state, error);
          break;
        case ToleranceMethod::SLOAN:
          rError = checkNormSloan(stressAtYield, stressAtYield(6), state, error);
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

        state_old   = state;
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
        absStress.block<6, 1>(0, 0) =
          (state.stress - state_old.stress).cwiseAbs();
        absStress(6) = std::abs(state.p0Star() - state_old.p0Star());

        reuseStep         = false;
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
MohrCoulombClassic::doRungeKuttaEigErr(
  const Eigen::Matrix<double, Steps, Steps>& AA,
  const Eigen::Matrix<double, Steps, 1>& BB,
  const Eigen::Matrix<double, Steps, 1>& BRes,
  const Eigen::Matrix<double, Steps, 1>& CC,
  const Eigen::Matrix<double, Steps - 1, 1>& ErrCoef,
  MohrCoulombState& state,
  const Vector7& epStrain,
  bool errorEstimate) const
{
  constexpr double criticalStepSize = 1.0e-4;
  Vector7 substepStrain             = epStrain;
  double newStepSize                = 0;
  for (int i = 0; i < 6; i++) {
    newStepSize += std::abs(substepStrain[i]);
  }
  newStepSize = 0.01 / (newStepSize * Order);
  if (newStepSize > 1) {
    newStepSize = 1;
  }
  dbg << "In doRungeKuttaErr::newStepSize = " << newStepSize << "\n";

  double stepLength        = 1;
  double totalSize         = 0;
  bool reuseStep           = false;
  double microStep         = 0;
  double stepAccuracyCheck = 0;
  double frequency   = 100000.0 / Order; // how often display info about steps
  double methodPower = std::pow(2.0, Order) * d_int.d_integrationTol;

  std::vector<MohrCoulombState> midStates(Steps);

  Eigen::Matrix<double, 6, Steps> dSigma =
    Eigen::Matrix<double, 6, Steps>::Zero();
  Eigen::Matrix<double, 6, Steps> plasticStrain =
    Eigen::Matrix<double, 6, Steps>::Zero();
  Eigen::Matrix<double, 1, Steps> dP0Star =
    Eigen::Matrix<double, 1, Steps>::Zero();

  Vector7 absStress        = Vector7::Zero();
  Vector7 reuseRes         = Vector7::Zero();
  Vector6 stressInc        = Vector6::Zero();
  Vector6 plasticStrainInc = Vector6::Zero();
  double p0StarInc         = 0.0;

  std::chrono::time_point<std::chrono::system_clock> startTime, endTime;
  startTime = std::chrono::system_clock::now();

  int stepNo    = 0;
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

    double rError      = 0;
    double rErrorOther = 0;

    dbg << "In doRungeKuttaErr:: Step Length = " << stepLength
        << " Current strain (0) = " << currentStrain(0) << "\n";

    // Make copies of the initial state
    for (auto& midState : midStates) {
      midState = state;
    }

    // Below the main R-K loop to calculate the value of intermediate stresses;
    // values stored in dSigma_temp
    if (reuseStep) {
      dSigma.col(0)        = reuseRes.block<6, 1>(0, 0);
      plasticStrain.col(0) = plasticStrainInc.block<6, 1>(0, 0);
      dP0Star(0, 0) = reuseRes(6);

    } else {

      /*
      Vector3 eigenVal;
      Matrix33 eigenVec;
      std::tie(eigenVal, eigenVec) = getEigen(midStates[0].stress);

      midStates[0].setStressEigen(eigenVal);
      midStates[0].strain.block<6, 1>(0, 0) =
        rotateToEigen(midStates[0].strain.block<6, 1>(0, 0), eigenVec);
      midStates[0].plasticStrain =
        rotateToEigen(midStates[0].plasticStrain, eigenVec);
      substepStrain.block<6, 1>(0, 0) =
        rotateToEigen(substepStrain.block<6, 1>(0, 0), eigenVec);
      */

      stressInc        = Vector6::Zero();
      std::tie(stressInc, p0StarInc) = calcPlastic(midStates[0], substepStrain);

      /*
      midStates[0].stress = rotateToOrigin(midStates[0].stress, eigenVec);
      midStates[0].strain.block<6, 1>(0, 0) =
        rotateToOrigin(midStates[0].strain.block<6, 1>(0, 0), eigenVec);
      midStates[0].plasticStrain =
        rotateToOrigin(midStates[0].plasticStrain, eigenVec);
      substepStrain.block<6, 1>(0, 0) =
        rotateToOrigin(substepStrain.block<6, 1>(0, 0), eigenVec);
      plasticStrainInc = rotateToOrigin(plasticStrainInc, eigenVec);
      stressInc        = rotateToOrigin(stressInc, eigenVec);
      */

      dSigma.col(0)        = stressInc;
      plasticStrain.col(0) = substepStrain.block<6,1>(0,0);
      dP0Star(0, 0) = p0StarInc;
    }

    for (int rkloop = 1; rkloop < Steps; rkloop++) {
      stressInc        = Vector6::Zero();
      plasticStrainInc = Vector6::Zero();
      p0StarInc        = 0;

      currentStrain = substepStrain * CC[rkloop];

      for (int i = 0; i < rkloop; i++) {
        stressInc += (dSigma.col(i) * AA(rkloop, i));
        plasticStrainInc += (plasticStrain.col(i) * AA(rkloop, i));
        p0StarInc += (dP0Star(0, i) * AA(rkloop, i));
      }

      midStates[rkloop].update(
        plasticStrainInc, currentStrain, stressInc, p0StarInc);

      /*
      Vector3 eigenVal;
      Matrix33 eigenVec;
      std::tie(eigenVal, eigenVec) = getEigen(midStates[rkloop].stress);

      midStates[rkloop].setStressEigen(eigenVal);
      midStates[rkloop].strain.block<6, 1>(0, 0) =
        rotateToEigen(midStates[0].strain.block<6, 1>(0, 0), eigenVec);
      midStates[rkloop].plasticStrain =
        rotateToEigen(midStates[0].plasticStrain, eigenVec);
      substepStrain.block<6, 1>(0, 0) =
        rotateToEigen(substepStrain.block<6, 1>(0, 0), eigenVec);
      */

      std::tie(stressInc, p0StarInc) = calcPlastic(midStates[rkloop], substepStrain);

      /*
      midStates[rkloop].stress =
        rotateToOrigin(midStates[rkloop].stress, eigenVec);
      midStates[rkloop].strain.block<6, 1>(0, 0) =
        rotateToOrigin(midStates[rkloop].strain.block<6, 1>(0, 0), eigenVec);
      midStates[rkloop].plasticStrain =
        rotateToOrigin(midStates[rkloop].plasticStrain, eigenVec);
      substepStrain.block<6, 1>(0, 0) =
        rotateToOrigin(substepStrain.block<6, 1>(0, 0), eigenVec);
      plasticStrainInc = rotateToOrigin(plasticStrainInc, eigenVec);
      stressInc        = rotateToOrigin(stressInc, eigenVec);
      */

      dSigma.col(rkloop)        = stressInc;
      plasticStrain.col(rkloop) = substepStrain.block<6,1>(0,0);
      dP0Star(0, rkloop) = p0StarInc;
    }

    Vector6 sigBRes = dSigma * BRes;
    // Vector6 epsPBRes = plasticStrain * BRes;
    double p0BRes = dP0Star * BRes;
    Vector6 sigB  = dSigma * BB;
    double p0B    = dP0Star * BB;

    Vector7 stressAtYield, error;
    stressAtYield.block<6, 1>(0, 0) = sigBRes;
    error.block<6, 1>(0, 0)  = sigB;
    stressAtYield(6) = p0BRes;
    error(6)  = p0B;

    // error estimate calculated in case we
    // have lower order solution instead of
    // error estimate
    if (!errorEstimate) {
      error -= stressAtYield;
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
        rError      = checkNorm(stressAtYield, stressAtYield(6), state, error);
        rErrorOther = checkNorm(stressAtYield, stressAtYield(6), state, errorOther);
        break;
      case ToleranceMethod::SLOAN:
        rError      = checkNormSloan(stressAtYield, stressAtYield(6), state, error);
        rErrorOther = checkNormSloan(stressAtYield, stressAtYield(6), state, errorOther);
        break;
      default:
        std::ostringstream err;
        err << "**ERROR** Improper d_tolMethod in increment.dta"
            << "\n";
        throw InvalidValue(err.str(), __FILE__, __LINE__);
    }

    for (int i = 0; i < 7; i++) {
      if (!std::isfinite(stressAtYield(i))) {
        std::cout << "**ERROR** Results not a number. Try correcting issue, "
                     "results may be "
                     "incorrect. Results = \n"
                  << stressAtYield << "\n";
        stressAtYield(i) = 0;
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
      state.update(
        plasticStrainInc, substepStrain, stressAtYield.block<6, 1>(0, 0), stressAtYield(6));
      //******
      //******

      // Keeping a copy to check values
      MohrCoulombState trialState = state;

      switch (d_int.d_driftCorrection) {
        case DriftCorrection::CORRECTION_AT_END: {

          /*
          Vector3 eigenVal;
          Matrix33 eigenVec;
          std::tie(eigenVal, eigenVec) = getEigen(state.stress);
          state.setStressEigen(eigenVal);
          state.strain.block<6, 1>(0, 0) =
            rotateToEigen(state.strain.block<6, 1>(0, 0), eigenVec);
          state.plasticStrain =
            rotateToEigen(midStates[0].plasticStrain, eigenVec);
          */

          correctDriftEnd(state);

          /*
          state.stress = rotateToOrigin(state.stress, eigenVec);
          state.strain.block<6, 1>(0, 0) =
            rotateToOrigin(state.strain.block<6, 1>(0, 0), eigenVec);
          state.plasticStrain = rotateToOrigin(state.plasticStrain, eigenVec);
          */

        } break;
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
          temp        = checkNorm(stressAtYield, stressAtYield(6), state, error);
          rErrorOther = checkNorm(stressAtYield, stressAtYield(6), state, errorOther);
          if (temp > rError) {
            rError = temp;
          }
          if (rErrorOther > rError) {
            rError = rErrorOther;
          }
          break;
        case ToleranceMethod::SLOAN:
          temp        = checkNormSloan(stressAtYield, stressAtYield(6), state, error);
          rErrorOther = checkNormSloan(stressAtYield, stressAtYield(6), state, errorOther);
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

        state_old   = state;
        double temp = double(stepNo) / frequency;
        if (std::modf(temp, &temp) == 0) {
          std::cout << "Step number : " << stepNo << "\n";
          std::cout << "Total size done is : " << totalSize
                    << " of whole step. Current stepLength = " << stepLength
                    << "\n";
          std::cout << "Stress state :\n" << state.stress << "\n";
        }

      } else {
        absStress.block<6, 1>(0, 0) =
          (state.stress - state_old.stress).cwiseAbs();
        absStress(6) = std::abs(state.p0Star() - state_old.p0Star());

        // this may be done only after successful re-evaluation of the substep
        reuseStep         = false;
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
MohrCoulombClassic::plasticMidpoint(MohrCoulombState& state,
                                    const Vector7& epStrain,
                                    Vector7& absStress,
                                    int numIter) const
{
  double h = 1.0 / (double)numIter;

  absStress                 = Vector7::Zero();
  Vector7 currentStrain     = epStrain * h;
  Vector7 halfCurrentStrain = currentStrain * 0.5;

  dbg << "Step Length = " << h << "\n";
  dbg << "Current strain [0] = " << currentStrain[0] << "\n";

  MohrCoulombState state_new, state_mid;
  Vector6 dSigma        = Vector6::Zero();
  Vector6 plasticStrain = Vector6::Zero();
  double dP0Star        = 0.0;
  for (int loop = 0; loop < numIter; loop++) {
    state_new = state;
    state_mid = state;

    /*
    Vector3 eigenVal;
    Matrix33 eigenVec;
    std::tie(eigenVal, eigenVec) = getEigen(state_new.stress);
    state_new.setStressEigen(eigenVal);
    state_mid.setStressEigen(eigenVal);

    state_new.strain.block<6, 1>(0, 0) =
      rotateToEigen(state_new.strain.block<6, 1>(0, 0), eigenVec);
    state_mid.strain.block<6, 1>(0, 0) =
      rotateToEigen(state_new.strain.block<6, 1>(0, 0), eigenVec);
    state_new.plasticStrain = rotateToEigen(state_new.plasticStrain, eigenVec);
    state_mid.plasticStrain = rotateToEigen(state_new.plasticStrain, eigenVec);
    halfCurrentStrain.block<6, 1>(0, 0) =
      rotateToEigen(halfCurrentStrain.block<6, 1>(0, 0), eigenVec);
    */

    std::tie(dSigma, dP0Star) = calcPlastic(state_new, halfCurrentStrain);

    state_mid.update(plasticStrain, halfCurrentStrain, dSigma, dP0Star);

    std::tie(dSigma, dP0Star) = calcPlastic(state_mid, halfCurrentStrain);

    /*
    plasticStrain = rotateToOrigin(plasticStrain, eigenVec);
    dSigma        = rotateToOrigin(dSigma, eigenVec);
    */

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
 * Implicit return to yield surface
 */
MohrCoulombState
MohrCoulombClassic::doReturnImplicit(const MohrCoulombState& state,
                                     RegionType& region) const
{
  MohrCoulombState state_new = state;
  state_new.stress *= -1.0; // tension negative

  dbg << "stress initial: " << state_new.stress.transpose() << "\n";

  /*
  Vector3 eigenVal;
  Matrix33 eigenVec;
  std::tie(eigenVal, eigenVec) = getEigen(state_new.stress);
  state_new.setStressEigen(eigenVal);

  state_new.strain.block<6, 1>(0, 0) =
    rotateToEigen(state_new.strain.block<6, 1>(0, 0), eigenVec);
  state_new.plasticStrain = rotateToEigen(state_new.plasticStrain, eigenVec);
  */

  double K                = d_elastic.d_K;
  double G                = d_elastic.d_G;
  Matrix66 elasticTangent = calculateElasticTangentMatrix(K, G);
  Matrix33 Delastic       = elasticTangent.block<3, 3>(0, 0);

  double k = (1.0 + d_yield.d_sin_phi) / (1.0 - d_yield.d_sin_phi);
  double m = (1.0 + d_potential.d_sin_psi) / (1.0 - d_potential.d_sin_psi);

  Vector3 apex  = Vector3::Ones();
  double sqrt_k = std::sqrt(k);
  double kRatio = sqrt_k / (k - 1);
  apex *= (2.0 * d_yield.d_cohesion * kRatio);

  Vector3 sigma     = state_new.stress.block<3, 1>(0, 0);
  Vector3 sigmaApex = sigma - apex;

  //  First set
  Vector3 df_dsigma = Vector3::Zero();
  df_dsigma(0)      = k;
  df_dsigma(2)      = -1.0;

  Vector3 dg_dsigma = Vector3::Zero();
  dg_dsigma(0)      = m;
  dg_dsigma(2)      = -1.0;

  double denominator = df_dsigma.transpose() * Delastic * dg_dsigma;
  Vector3 RP1        = (Delastic * dg_dsigma) / denominator;

  dbg << "RP1 = " << RP1.transpose() << "\n";

  Vector3 R1 = Vector3::Ones();
  R1(2)      = k;

  Vector3 R2 = Vector3::Ones();
  R2(1)      = k;
  R2(2)      = k;

  auto NI_I1  = RP1.cross(R1);
  auto NI_III = RP1.cross(R2);
  dbg << "NI_I1 = " << NI_I1.transpose() << "\n";
  dbg << "NI_III = " << NI_III.transpose() << "\n";

  double pI_II  = NI_I1.transpose() * sigmaApex;
  double pI_III = NI_III.transpose() * sigmaApex;

  //  Second set
  df_dsigma(0) = 0.0;
  df_dsigma(1) = k;
  df_dsigma(2) = -1.0;

  dg_dsigma(0) = 0.0;
  dg_dsigma(1) = m;
  dg_dsigma(2) = -1.0;

  denominator = df_dsigma.transpose() * Delastic * dg_dsigma;
  Vector3 RP2 = (Delastic * dg_dsigma) / denominator;

  dbg << "RP2 = " << RP2.transpose() << "\n";

  auto N2 = RP1.cross(RP2);

  denominator      = R1.transpose() * N2;
  double numerator = N2.transpose() * sigmaApex;

  double t1 = numerator / denominator;

  //  Third set
  df_dsigma(0) = k;
  df_dsigma(1) = 0;
  df_dsigma(2) = -1.0;

  dg_dsigma(0) = m;
  dg_dsigma(1) = 0;
  dg_dsigma(2) = -1.0;

  denominator = df_dsigma.transpose() * Delastic * dg_dsigma;
  Vector3 RP3 = (Delastic * dg_dsigma) / denominator;

  dbg << "RP3 = " << RP3.transpose() << "\n";

  auto N3 = RP1.cross(RP3);

  denominator = R2.transpose() * N3;
  numerator   = N3.transpose() * sigmaApex;

  double t2 = numerator / denominator;

  // Alternative calculations for d_sin_phi=0 [i.e. apex--> inf] is needed
  if ((t1 > 0.0) && (t2 > 0.0)) {

    region = RegionType::APEX_REGION_4;
    state_new.stress.block<3, 1>(0, 0) = apex;

  } else if (pI_II < 0) {

    region = RegionType::APEX_REGION_2;
    state_new.stress.block<3, 1>(0, 0) = R1 * t1 + apex;

  } else if (pI_III <= 0) {

    region   = RegionType::APEX_REGION_1;
    double f = k * state_new.stress(0) - state_new.stress(2) -
               2.0 * d_yield.d_cohesion * sqrt_k;

    state_new.stress.block<3, 1>(0, 0) = sigma - RP1 * f;

  } else {

    region = RegionType::APEX_REGION_3;

    state_new.stress.block<3, 1>(0, 0) = R2 * t2 + apex;
  }

  /*
  state_new.stress = rotateToOrigin(state_new.stress, eigenVec);
  state_new.strain.block<6, 1>(0, 0) =
    rotateToOrigin(state_new.strain.block<6, 1>(0, 0), eigenVec);
  state_new.plasticStrain = rotateToOrigin(state_new.plasticStrain, eigenVec);
  */

  state_new.stress *= -1.0;

  dbg << "region: " << region << "\n";
  dbg << "stress: " << state_new.stress.transpose() << "\n";

  return state_new;
}
