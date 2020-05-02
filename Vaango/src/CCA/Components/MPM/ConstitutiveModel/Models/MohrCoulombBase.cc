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

#include <CCA/Components/MPM/ConstitutiveModel/Models/MohrCoulombBase.h>

#include <Core/Exceptions/InvalidValue.h>
#include <Core/Util/DebugStream.h>

#include <chrono>
#include <cmath>
#include <iomanip>

using namespace Uintah;

static DebugStream dbg_doing("BaseMC_doing", false);
static DebugStream dbg("BaseMC", false);
static DebugStream dbg_params("BaseMC_params", false);
static DebugStream dbg_unloading("BaseMC_unloading", false);
static DebugStream dbg_alpha("BaseMC_alpha", false);
static DebugStream dbg_mod("BaseMC_mod", false);
static DebugStream dbg_calcplastic("BaseMC_calcplastic", false);
static DebugStream dbg_estimate("BaseMC_estimate", false);

constexpr long64 testParticleID = 6619136;

/**
 * Constructors
 */
MohrCoulombBase::MohrCoulombBase()
{
  double G = 10000.0;
  double K = 20000.0;
  double cohesion = 0;
  double phi = 30;
  setModelParameters(G, K, cohesion, phi, phi, -1);
  d_int.setDefaults(d_yield);
}

MohrCoulombBase::MohrCoulombBase(double G, double K, double cohesion,
                                 double phi, double psi, double pMin)
{
  setModelParameters(G, K, cohesion, phi, psi, pMin);
  d_int.setDefaults(d_yield);
}

MohrCoulombBase::MohrCoulombBase(const MohrCoulombBase* cm)
{
  d_elastic = cm->d_elastic;
  d_yield = cm->d_yield;
  d_potential = cm->d_potential;
  d_nonAssociated = cm->d_nonAssociated;
}

/**
 * Initialization methods
 */
void
MohrCoulombBase::setModelParameters(double G, double K, double cohesion,
                                     double phi, double psi, double pMin)
{
  d_elastic.set(G, K);
  d_yield.set(cohesion, phi, pMin);
  d_potential.set(psi);
  if (phi != psi) {
    d_nonAssociated = true;
  } else {
    d_nonAssociated = false;
  }
}

void
MohrCoulombBase::setIntegrationParameters(int maxIterPegasus, double alfaCheck,
                                double alfaChange, double alfaRatio,
                                double yieldTol, double integrationTolerance,
                                double betaFactor, double minMeanStress,
                                double suctionTol,
                                SolutionAlgorithm solutionAlgorithm,
                                ToleranceMethod toleranceMethod,
                                DriftCorrection driftCorrection)
{
  d_int.d_maxIter = maxIterPegasus;
  d_int.d_alfaCheck = alfaCheck;
  d_int.d_alfaChange = alfaChange;
  d_int.d_alfaRatio = alfaRatio;
  d_int.d_yieldTol = yieldTol;
  d_int.d_integrationTol = integrationTolerance;
  d_int.d_betaFactor = betaFactor;
  d_int.d_minMeanStress = minMeanStress;
  d_int.d_suctionTol = suctionTol;

  d_int.d_solutionAlgorithm = solutionAlgorithm;
  d_int.d_tolMethod = toleranceMethod;
  d_int.d_driftCorrection = driftCorrection;
}

/**
 * Calc elastic state
 */
void
MohrCoulombBase::calcElastic(const Vector7& strainInc,
                              const MohrCoulombState& initialState,
                              MohrCoulombState& finalState) const
{
  finalState = initialState;

  Vector6 stressInc =
    calcStressIncElast(initialState.specificVolume(), initialState.stress,
                       initialState.strain, strainInc);

  dbg << __LINE__ << ":BaseMC::Stress increment is: " << stressInc.transpose() << "\n";

  Vector6 plasticStrainInc = Vector6::Zero();
  double p0StarInc = 0.0;
  finalState.update(plasticStrainInc, strainInc, stressInc, p0StarInc);
}

Vector6
MohrCoulombBase::calcStressIncElast(double nu0, const Vector6& stress_old,
                                     const Vector7& strain_old,
                                     const Vector7& strain_inc) const
{
  double K = d_elastic.d_K;
  double G = d_elastic.d_G;

  double K43G = K + 4.0 * G / 3.0;
  double K23G = K - 2.0 * G / 3.0;

  Vector6 stress_inc;
  stress_inc(0) = K43G * strain_inc(0) + K23G * (strain_inc(1) + strain_inc(2));
  stress_inc(1) = K43G * strain_inc(1) + K23G * (strain_inc(0) + strain_inc(2));
  stress_inc(2) = K43G * strain_inc(2) + K23G * (strain_inc(0) + strain_inc(1));
  stress_inc(3) = 2 * G * strain_inc(3);
  stress_inc(4) = 2 * G * strain_inc(4);
  stress_inc(5) = 2 * G * strain_inc(5);

  return stress_inc;
}

/*
 *  Compute elastic tangent stiffness
 */
Matrix67
MohrCoulombBase::calculateElasticTangentMatrix(
  const MohrCoulombState& state) const
{
  double K = d_elastic.d_K;
  double G = d_elastic.d_G;
  Matrix66 DEP_66 = calculateElasticTangentMatrix(K, G);

  Matrix67 DEP = Matrix67::Zero();
  DEP.block<6, 6>(0, 0) = DEP_66;
  return DEP;
}

Matrix66
MohrCoulombBase::calculateElasticTangentMatrix(double K, double G) const
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

/*
 *  Check the gradient between initial and final state
 *
 *  The normalisation is important to catch the unloading due to shear stress
 *  12 13 23. As the q is always positive, large unloading - loading values
 *  of 12 13 or 23 component lead to larger q, the shear stress change is
 *  taken as positive and there is no unloading. This fix, though not the
 *  most elegant, should work.
 */
bool
MohrCoulombBase::checkGradient(const MohrCoulombState& initialState,
                                const MohrCoulombState& finalState) const
{
  dbg_doing << "Doing MohrCoulombBase::checkGradient\n";

  Vector7 strainInc = finalState.strain - initialState.strain;

  dbg << __LINE__ << ":BaseMC::strain_final = " << finalState.strain.transpose()
      << "\n      strain_initial = " << initialState.strain.transpose()
      << "\n      Strain inc = " << strainInc.transpose() << "\n";

  // Scale up the strain increment
  /*
  double max_strain_inc = strainInc.block<6,1>(0,0).cwiseAbs().maxCoeff();
  strainInc.block<6,1>(0,0) /= (max_strain_inc * 1.0e-10);
  strainInc(6) = finalState.suction() - initialState.suction();
  */

  Matrix67 tangentElastic = calculateElasticTangentMatrix(initialState);

  Vector6 stressInc = tangentElastic * strainInc;

  dbg << __LINE__ << ":BaseMC::Strain inc = " << strainInc.transpose() << "\n";
  dbg << "BaseMC::TangentElastic: \n" << tangentElastic << "\n";
  dbg << "BaseMC::Stress inc = " << stressInc.transpose() << "\n";
  dbg << "BaseMC::Suction Increment = " << strainInc(6) << "\n";

  Vector6 df; // not initialised; will contain values of derivatives of the
              // yield locus function
  double cosinus = computeDsigmaDotDf(initialState.stress, stressInc, df, 0, 0);

  if (cosinus > -d_int.d_yieldTol) {
    return false;
  }

  dbg << __LINE__ << ":BaseMC::Unloading occuring... cosinus = " << cosinus << "\n";

  return true; //(negative cosinus means unloading occurs)
}

/**
  This procedure finds a gradient to the yield locus at given point.
  It is used later to determine whether we have partly elastic unloading or only
  plastic step.
  The gradient may be found only when we are on "main" plastic yield locus;
  Otherwise results may be erroneous.

  Parameters: state[], s[] (stress), ds[] (stress increment), suction,
  Gradient = gradient[1], gradient[2] etc - there the result will be stored
  First the C(s) part must be calculated, then the gradient.

  returns cosinus of the angle...
 */
double
MohrCoulombBase::computeDsigmaDotDf(const Vector6& s, const Vector6& ds,
                               Vector6& df_dSigma, double /*suction*/,
                               double /*dsuction*/) const
{
  dbg_doing << "Doing MohrCoulombBase::findGradient\n";

  // compute gradient of yield function
  Vector6 dg_dSigma;
  std::tie(df_dSigma, dg_dSigma) = computeDfDsigma(s);

  dbg << __LINE__ << ":BaseMC::df_dsigma = " << df_dSigma.transpose() << "\n";

  // Compute norm of ds and df_dsigma
  double dfLength = df_dSigma.norm();
  double dsLength = ds.norm();

  dbg << __LINE__ << ":BaseMC::||df|| = " << dfLength 
      << " ||dsig|| = " << dsLength << "\n";
  
  // Calculate dot product of the two normalized vectors
  double cosin = df_dSigma.dot(ds)/ (dfLength * dsLength);
  
  if (dbg.active()) {
    if (cosin < -d_int.d_yieldTol) {
      dbg << __LINE__ << ":BaseMC::Check if no problem in the findGradient\n";
      dbg << "cosin = " << cosin << "\n";
    }
  }

  return cosin;
}

/**
 Find intersection with yield surface during unloading
 */
std::tuple<Vector7, Vector7>
MohrCoulombBase::findIntersectionUnloading(
  const Vector7& strainIncrement, const MohrCoulombState& initialState)
{
  dbg_doing << "Doing MohrCoulombBase::findIntersectionUnloading\n";

  dbg_alpha << "ParticleID = " << initialState.particleID << "\n";

  /*
  if (initialState.particleID == testParticleID) {
    dbg_alpha.setActive(true);
  } else {
    dbg_alpha.setActive(false);
  }
  */
  double alpha = findYieldAlpha(initialState.state, initialState.stress,
                                initialState.strain, strainIncrement);

  Vector7 elasticStrainInc = strainIncrement * alpha;
  Vector7 plasticStrainInc = strainIncrement - elasticStrainInc;
  
  return std::make_tuple(elasticStrainInc, plasticStrainInc);
}

std::tuple<Vector7, Vector7>
MohrCoulombBase::findIntersection(const Vector7& strainIncrement,
                                   const MohrCoulombState& initialState) const
{
  dbg_doing << "Doing MohrCoulombBase::findIntersection\n";

  dbg_mod << "ParticleID = " << initialState.particleID << "\n";

  /*
  if (initialState.particleID == testParticleID) {
    dbg_mod.setActive(true);
  } else {
    dbg_mod.setActive(false);
  }
  */

  double alpha = findYieldModified(initialState.state, initialState.stress,
                                   initialState.strain, strainIncrement);
  Vector7 elasticStrainInc = strainIncrement * alpha;
  Vector7 plasticStrainInc = strainIncrement - elasticStrainInc;

  return std::make_tuple(elasticStrainInc, plasticStrainInc);
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
MohrCoulombBase::findYieldAlpha(const Vector3& state, const Vector6& s0,
                                 const Vector7& eps0, const Vector7& dep)
{
  dbg_doing << "Doing MohrCoulombBase::findYieldAlpha\n";
  dbg_doing << "BaseMC::Call of this procedure is generally rare "
               "and most likely improper\n";

  dbg_params << "BaseMC::Parameters are set to:\n";
  dbg_params << "BaseMC::Maximum number of iteration d_maxIter: " << d_int.d_maxIter << "\n";
  dbg_params << "BaseMC::Value of d_alfaCheck - when the additional iter. is to "
                "enter: " << d_int.d_alfaCheck << "\n";
  dbg_params << "BaseMC::Value of change of alfa in the step: " << d_int.d_alfaChange << "\n";
  dbg_params << "BaseMC::alfa old/alfa ratio: " << d_int.d_alfaRatio << "\n";

  Vector6 sIni = s0;
  dbg_alpha << __LINE__ << ":BaseMC::FindYieldAlpha:. sIni = " << sIni.transpose() << "\n";

  double f0 = computeYieldNormalized(sIni);
  dbg_alpha << __LINE__ << ":BaseMC:: f0 (first) = " << f0 << "\n";

  double yieldTol_old = d_int.d_yieldTol;
  double delta = 0.0;
  if (f0 < 0) {
    delta = 0;
    d_int.d_yieldTol = -0.9 * f0;
  } else {
    dbg_alpha << __LINE__ << ":BaseMC::findYieldAlpha : case not coded for yet\n";

    delta = (f0 + d_int.d_yieldTol) / 2;
    d_int.d_yieldTol = 0.9 * (d_int.d_yieldTol - f0) / 2;
  }

  f0 -= delta;
  dbg_alpha << __LINE__ << ":BaseMC:: f0 (second) = " << f0 << " delta = " << delta << "\n";


  double alfa0 = 0;
  double alfa1 = 1;
  Vector7 epsFin = dep * alfa1;
  dbg_alpha << __LINE__ << ":BaseMC::dEps = " << dep.transpose() << "\n";
  dbg_alpha << __LINE__ << ":BaseMC::epsFin = " << epsFin.transpose() << "\n";

  Vector6 dsigma = calcStressIncElast(state(2), sIni, eps0, epsFin);
  dbg_alpha << __LINE__ << ":BaseMC::dsigma = " << dsigma.transpose() << "\n";
  Vector6 sFin = sIni + dsigma;
  dbg_alpha << __LINE__ << ":BaseMC::sFin = " << sFin.transpose() << "\n";

  double f1 = computeYieldNormalized(sFin);
  dbg_alpha << __LINE__ << ":BaseMC:: f1 (first) = " << f1 << "\n";

  f1 -= delta;
  dbg_alpha << __LINE__ << ":BaseMC:: f1 (second) = " << f1 << " delta = " << delta << "\n";

  dbg_alpha << __LINE__ << ":BaseMC::Procedure findYieldAlpha. Initial values of f0 = " << f0
      << " f1 = " << f1 << "\n";
  dbg_alpha << "BaseMC::f0 should lie on the yield locus within tolerance: " << yieldTol_old
      << "\n";

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

      dbg_alpha << __LINE__ << ":BaseMC::In iteration " << iter << " alfa = " << alfa
          << " and f = " << fAlfa << "\n";

    // if fAlfa within tolerance, we have the solution
    if (std::abs(fAlfa) < d_int.d_yieldTol) {

      alphaYield = alfa0 + alfa * (alfa1 - alfa0);

      dbg_alpha << __LINE__ << ":BaseMC:Solution in findYieldAlpha procedure was found after " << iter
          << " iterations."
          << "\n";
      dbg_alpha << "BaseMC:Value of the yield function is equal to: " << fAlfa << "\n";
      dbg_alpha << "BaseMC:Modified value of tolerance is: " << d_int.d_yieldTol << "\n";
      dbg_alpha << "BaseMC:Value of delta is: " << delta << "\n";
      dbg_alpha << "BaseMC:Alfa is equal to: " << alphaYield << "\n";

      if (iter > 50) {
        std::cout << "**WARNING** Large number of iterations!!! Solution is "
                     "however correct... "
                  << "\n";
      }

      d_int.d_yieldTol = yieldTol_old;

      dbg_alpha << __LINE__ << ":BaseMC::Yield tolerance is set back to: " << d_int.d_yieldTol << "\n";

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
        dbg_alpha << __LINE__ << ":BaseMC::Problematic iteration entered !!!" << "\n";

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
        dbg_alpha << __LINE__ << ":BaseMC::Problematic iteration entered !!!" << "\n";

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
MohrCoulombBase::findYieldModified(const Vector3& state, const Vector6& stress_old,
                                    const Vector7& strain_old,
                                    const Vector7& strain_inc) const
{
  dbg_doing << "Doing MohrCoulombBase::findYieldModified\n";

  dbg_params << __LINE__ << ":BaseMC::Parameters are set to:\n";
  dbg_params << "BaseMC::Maximum number of iteration d_maxIter: " << d_int.d_maxIter << "\n";
  dbg_params << "BaseMC::Value of d_alfaCheck - when the additional iter. is to "
         "enter: " << d_int.d_alfaCheck << "\n";
  dbg_params << "BaseMC::Value of change of alfa in the step: " << d_int.d_alfaChange << "\n";
  dbg_params << "BaseMC::alfa old/alfa ratio: " << d_int.d_alfaRatio << "\n";

  double alfa0 = 0.0;
  double alfa1 = 1.0;

  double alfa = alfa1;
  Vector7 deps = strain_inc * alfa;

  Vector6 dsigma = calcStressIncElast(state(2), stress_old, strain_old, deps);
  Vector6 stress_alfa = stress_old + dsigma;

  double f0 = computeYieldNormalized(stress_old);
  double f1 = computeYieldNormalized(stress_alfa);

  dbg_mod << __LINE__ << ":BaseMC::Procedure findYieldModified. Initial values of f0 = " << f0
      << " f1 = " << f1 << "\n";
  dbg_mod << "BaseMC::Value of f0 should be negative, and value of f1 should be "
         "positive.\n";
  dbg_mod << "BaseMC::Values should be larger than tolerance for yield: "
      << d_int.d_yieldTol << "\n";

  double alfa_old = 0;
  for (int iter = 0; iter < d_int.d_maxIter; iter++) {

    bool problems = false;
    alfa_old = alfa;

    // Interpolate linearly between the two yield surfaces to find the f = 0 location
    double tt = f0 / (f0 - f1);
    alfa = (1.0 - tt) * alfa0 + tt * alfa1;

    Vector7 deps_alfa = strain_inc * alfa;
    Vector6 dsigma_alfa = calcStressIncElast(state(2), stress_old, strain_old, deps_alfa);
    stress_alfa = stress_old + dsigma_alfa;

    double f_alfa = computeYieldNormalized(stress_alfa);

    dbg_mod << __LINE__ << ":BaseMC::In iteration " << iter << " alfa = " << alfa
            << " and f = " << f_alfa << " alfa0 = " << alfa0 << "\n";

    // if the difference is below numerical the accuracy, we have the solution
    if ((alfa1 - alfa0) < TINY) {
      f_alfa = 0.0;
    }

    if (std::abs(f_alfa) < d_int.d_yieldTol) {

      // if f_alfa within tolerance, we have the solution

      if (iter > 50) {
        std::cout << "**WARNING** Large number of iterations in "
                     "findYieldModified procedure!!! "
                     "Solution is however correct..."
                  << "\n";
      }
      dbg_mod << __LINE__ 
          << ":BaseMC::Solution alfa = " << alfa 
          << " in findYieldModified procedure was found after "
          << iter << " iterations." << "\n";

      return alfa;
    }

    if (f_alfa > 0) {

      // if f_alfa > 0 - we are yielding - max alfa set to current alfa
      f_alfa = computeYieldNormalized(stress_alfa);
      if (((1 - alfa) < d_int.d_alfaCheck) &&
          ((1 - alfa_old) / (1 - alfa)) < d_int.d_alfaRatio) {
        problems = true;
      }

      alfa1 = alfa;
      f1 = f_alfa;

      if (problems) {

        dbg_mod << __LINE__ << ":BaseMC::Problematic iteration entered !!!\n";

        alfa = (1.0 - d_int.d_alfaChange) * alfa0 + d_int.d_alfaChange * alfa1;
        deps_alfa = strain_inc * alfa;

        dsigma_alfa = calcStressIncElast(state(2), stress_old, strain_old, deps_alfa);
        stress_alfa = stress_old + dsigma_alfa;

        f_alfa = computeYieldNormalized(stress_alfa);

        if (f_alfa > 0) {
          // if f_alfa > 0 - we are yielding - max alfa set to current alfa
          alfa1 = alfa;
          f1 = f_alfa;
        } else {
          // if f_alfa < 0 - we are elastic - minimum alfa is set to current alfa
          alfa0 = alfa;
          f0 = f_alfa;
        }
      }

    } else {

      // if f_alfa < 0 - we are elastic - minimum alfa is set to current alfa
      f_alfa = computeYieldNormalized(stress_alfa);
      if ((alfa < d_int.d_alfaCheck) && (alfa_old / alfa) < d_int.d_alfaRatio) {
        problems = true;
      }
      alfa0 = alfa;
      f0 = f_alfa;

      if (problems) {

        dbg_mod << __LINE__ << ":BaseMC::Problematic iteration entered !!!\n";

        alfa = (1.0 - d_int.d_alfaChange) * alfa0 + d_int.d_alfaChange * alfa1;
        deps_alfa = strain_inc * alfa;

        dsigma_alfa = calcStressIncElast(state(2), stress_old, strain_old, deps_alfa);
        stress_alfa = stress_old + dsigma_alfa;

        f_alfa = computeYieldNormalized(stress_alfa);
        if (f_alfa > 0) {
          // if f_alfa > 0 - we are yielding - max alfa set to current alfa
          alfa1 = alfa;
          f1 = f_alfa;
        } else {
          // if f_alfa < 0 - we are elastic - minimum alfa is set to current alfa
          alfa0 = alfa;
          f0 = f_alfa;
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
  err << "Stress:" << stress_old << "\n";
  err << "Strain:" << strain_inc << "\n";
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
MohrCoulombBase::computeNu(const Vector6& s, const Vector3& state,
                            double suction) const
{
  // does nothing for SMC
  return 1.0;
}

/**
  Calculate stress increment during plastic deformation
 */
double
MohrCoulombBase::calculatePlastic(const Vector7& plasticStrainInc,
                                   MohrCoulombState& state) const
{
  dbg_doing << "Doing MohrCoulombBase::calculatePlastic\n";

  double time;
  int numIter;

  switch (d_int.d_solutionAlgorithm) {
    case SolutionAlgorithm::RUNGE_KUTTA_SECOND_ORDER_MODIFIED_EULER: {
      // calculate using Modified Euler scheme
      std::tie(time, numIter) = plasticRKME221(state, plasticStrainInc);
      break;
    }
    case SolutionAlgorithm::RUNGE_KUTTA_THIRD_ORDER_NYSTROM: {
      // calculate using 3rd order RK scheme (Nystrom)
      std::tie(time, numIter) = plasticRK332(state, plasticStrainInc);
      break;
    }
    case SolutionAlgorithm::RUNGE_KUTTA_THIRD_ORDER_BOGACKI: {
      // calculate using 3rd order RK scheme (Bogacki - Shampine)
      std::tie(time, numIter) = plasticRKBog432(state, plasticStrainInc);
      break;
    }
    case SolutionAlgorithm::RUNGE_KUTTA_FOURTH_ORDER: {
      // calculate using 4th order RK scheme
      std::tie(time, numIter) = plasticRK543(state, plasticStrainInc);
      break;
    }
    case SolutionAlgorithm::RUNGE_KUTTA_FIFTH_ORDER_ENGLAND: {
      // calculate using 5th order RK scheme (England)
      std::tie(time, numIter) = plasticRKEng654(state, plasticStrainInc);
      break;
    }
    case SolutionAlgorithm::RUNGE_KUTTA_FIFTH_ORDER_CASH: {
      // calculate using 5th order RK scheme (Cash - Karp)
      std::tie(time, numIter) = plasticRKCK654(state, plasticStrainInc);
      break;
    }
    case SolutionAlgorithm::RUNGE_KUTTA_FIFTH_ORDER_DORMAND: {
      // calculate using 5th order RK scheme (Dormand - Prince)
      std::tie(time, numIter) = plasticRKDP754(state, plasticStrainInc);
      break;
    }
    case SolutionAlgorithm::RUNGE_KUTTA_FIFTH_ORDER_BOGACKI: {
      // calculate using 5th order RK scheme (Bogacki - Shampine)
      std::tie(time, numIter) = plasticRKErr8544(state, plasticStrainInc);
      break;
    }
    case SolutionAlgorithm::EXTRAPOLATION_BULIRSCH: {
      // calculate using extrapolation method (Bulirsch - Stoer)
      std::tie(time, numIter) = plasticExtrapol(state, plasticStrainInc);
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
void
MohrCoulombBase::getParamRKME221(Eigen::Matrix<double, 2, 2>& A,
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
void
MohrCoulombBase::getParamRK332(Eigen::Matrix<double, 3, 3>& A,
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
void
MohrCoulombBase::getParamRKBog432(Eigen::Matrix<double, 4, 4>& A,
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
void
MohrCoulombBase::getParamRK543(Eigen::Matrix<double, 5, 5>& A,
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
void
MohrCoulombBase::getParamRKEng654(Eigen::Matrix<double, 6, 6>& A,
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
void
MohrCoulombBase::getParamRKCK654(Eigen::Matrix<double, 6, 6>& A,
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
void
MohrCoulombBase::getParamRKDP754(Eigen::Matrix<double, 7, 7>& A,
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
void
MohrCoulombBase::getParamRKErr8544(Eigen::Matrix<double, 8, 8>& A,
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
 * Compute the stress and p0star increment,
 */
std::tuple<Vector6, double>
MohrCoulombBase::calcPlastic(const MohrCoulombState& state,
                             const Vector7& strainInc) const 
{
  std::cout << "\n\n** Particle ID = " << state.particleID << "\n\n";
  // No p0Star calculation yet
  double dP0Star = 1.0;

  if (!state.checkIfFinite()) {
    std::ostringstream err;
    err << "**Error** in the calcPlastic Procedure. "
        << "Point internal values are not finite.\n";
    throw InvalidValue(err.str(), __FILE__, __LINE__);
  }

  // Save stress
  auto stress = state.stress;

  // Save elastic tangent modulus
  double K = d_elastic.d_K;
  double G = d_elastic.d_G;
  Matrix66 elasticTangent = calculateElasticTangentMatrix(K, G);

  // Compute trial stress
  Vector6 stress_trial = stress + elasticTangent * strainInc.block<6,1>(0,0);

  dbg_calcplastic << "stress = " << toMatrix3(stress) << "\n"
                  << "stress_trial = " << toMatrix3(stress_trial) << "\n";

  // Project trial stress on to yield surface
  Vector6 stress_closest_Pn = 
    projectTrialStressToYieldSurface(strainInc.block<6,1>(0,0), 
                                     stress, elasticTangent, 
                                     stress_trial);

  if (dbg_calcplastic.active()) {
    double f_new = computeYieldNormalized(stress_closest_Pn);
    dbg_calcplastic << "f_final = " << f_new << "\n"
                    << "stress_closest = " << toMatrix3(stress_closest_Pn) << "\n";
  }

  Vector6 dSigma = stress_closest_Pn - stress;

  return std::make_tuple(dSigma, dP0Star);
}

/**
 * Find intersection with yield surface along projection direction
 */
Vector6 
MohrCoulombBase::projectTrialStressToYieldSurface(const Vector6& strainInc,
                                                  const Vector6& stress_old, 
                                                  const Matrix66& elasticTangent, 
                                                  const Vector6& stress_trial) const 
{
  // Update the normals to the yield surface
  Vector6 df_dsigma = Vector6::Zero();
  Vector6 dg_dsigma = Vector6::Zero();
  std::tie(df_dsigma, dg_dsigma) = computeDfDsigma(stress_old);

  // Update projection direction
  Vector6 proj_direction = elasticTangent * dg_dsigma;
  proj_direction /= proj_direction.norm();

  dbg_calcplastic << "df_dsigma = " << df_dsigma.transpose() << "\n"
                  << "dg_dsigma = " << dg_dsigma.transpose() << "\n"
                  << "P_n = " << proj_direction.transpose() << "\n";

  bool foundInitialAlpha = false;
  double alpha, f_trial, f_alpha;
  Vector6 sigma_alpha;

  std::tie(foundInitialAlpha, alpha, f_trial, f_alpha, sigma_alpha) = 
    estimateInitialBisectionParameter(stress_old, stress_trial, proj_direction);

  bool successfulBisection = true;

  if (foundInitialAlpha) {

    // If elastic for some reason
    if (f_trial < 0.0) {
      return stress_trial;
    }

    dbg_calcplastic << "foundInitialAlpha::alpha = " << alpha 
                    << " f_trial = " << f_trial << " f_alpha = " << f_alpha << "\n";

    std::tie(successfulBisection, alpha, f_alpha, sigma_alpha) = 
      findIntersectionWithBisection(alpha, f_alpha, stress_trial, proj_direction);

  } else {

    std::cout << "**WARNING** f_trial = " << f_trial << " and f_alpha = " << f_alpha 
              << " have the same sign.  Cannot use bisection.\n";

    if ((stress_trial - stress_old).norm() > 1.0 && 
        (std::abs(f_trial - f_alpha) > 1.0e-5)) {

      Vector6 strainInc_step1 = strainInc * 0.5;
      Vector6 stress_trial_step1 = stress_old + elasticTangent * strainInc_step1;

      std::cout << "\n\nalpha = "  << alpha << "\n";
      std::cout << "df_dsigma = " << df_dsigma.transpose() << "\n"
                << "P_n = " << proj_direction.transpose() << "\n";
      std::cout << "s0     = " << toMatrix3(stress_old) << "\n"
                << "strial = " << toMatrix3(stress_trial_step1) << "\n"
                << "salpha = " << toMatrix3(stress_trial_step1 - alpha * proj_direction) << "\n\n";

      Vector6 stress_step1 = 
        projectTrialStressToYieldSurface(strainInc_step1, stress_old, elasticTangent,
                                         stress_trial_step1);

      dbg_calcplastic << "After step1: s0     = " << toMatrix3(stress_old) << "\n"
                      << "             strial = " << toMatrix3(stress_trial_step1) << "\n"
                      << "             snew   = " << toMatrix3(stress_step1) << "\n";

      Vector6 stress_trial_step2 = stress_step1 + elasticTangent * strainInc_step1;
      Vector6 stress_step2 = 
        projectTrialStressToYieldSurface(strainInc_step1, stress_step1, elasticTangent,
                                         stress_trial_step2);

      dbg_calcplastic << "After step2: s0     = " << toMatrix3(stress_step1) << "\n"
                      << "             strial = " << toMatrix3(stress_trial_step2) << "\n"
                      << "             snew   = " << toMatrix3(stress_step2) << "\n";


      sigma_alpha = stress_step2;

    } else {

      sigma_alpha = stress_trial;

    }
  }

  if (!successfulBisection) {
    std::cout << "**WARNING** Bisection unsuccessful.\n"
              << "            Using first-order stress projection on to yield surface instead\n";
    std::cout << "            Final stress state is not on the yield surface.\n"
                 "            Some error may arise and results may be incorrect.\n";
    std::cout << " sigma_alpha = " << sigma_alpha.transpose() << "\n";
    std::cout << " f_alpha = " << f_alpha << " alpha = " << alpha << "\n"
              << " df_dsigma = " << df_dsigma.transpose() << "\n"
              << " proj_dir = " << proj_direction.transpose() << "\n"
              << " stress_old = " << toMatrix3(stress_old) << "\n"
              << " stress_trial = " << toMatrix3(stress_trial) << "\n"
              << " stress_new = " << toMatrix3(sigma_alpha) << "\n";

    auto dSigma = firstOrderStressUpdate(strainInc, elasticTangent, 
                                         df_dsigma, dg_dsigma);
    sigma_alpha = stress_old + dSigma;
  }

  return sigma_alpha;
}

/**
 * Returns status: success = true, failure = false
 *         alpha, f_trial, f_alpha, sigma_alpha
 */
std::tuple<bool, double, double, double, Vector6>
MohrCoulombBase::estimateInitialBisectionParameter(const Vector6& stress_old,
                                                   const Vector6& stress_trial,
                                                   const Vector6& proj_direction) const
{
  double f_trial = computeYieldNormalized(stress_trial);

  // Check if elastic for some reason
  if (f_trial < 0.0) {
    return std::make_tuple(true, 0.0, f_trial, 0.0, stress_trial);
  }

  double alpha = (stress_trial - stress_old).norm();
  Vector6 sigma_alpha = stress_trial - alpha * proj_direction;
  double f_alpha = computeYieldNormalized(sigma_alpha);

  dbg_estimate << "\n\nestimate alpha: stress_old   = " << toMatrix3(stress_old) << "\n";
  dbg_estimate << "estimate alpha: stress_trial = " << toMatrix3(stress_trial) << "\n";
  dbg_estimate << "estimate alpha: stress_alpha = " << toMatrix3(sigma_alpha) << "\n";
  dbg_estimate << "estimate alpha: Pn = " << proj_direction.transpose() << "\n";
  dbg_estimate << "estimate alpha: alpha = " << alpha 
                  << " f_trial = " << f_trial << " f_alpha = " << f_alpha << "\n\n";

  // Make sure the yield functions have opposite signs for bisection to work
  int numIter = 0;
  while (std::signbit(f_trial) == std::signbit(f_alpha)) {
    alpha *= 1.1;
    sigma_alpha = stress_trial - alpha * proj_direction;
    f_alpha = computeYieldNormalized(sigma_alpha);

    if (numIter > 10 || !std::isfinite(alpha) || !std::isfinite(f_alpha)) {

      dbg_estimate << "alpha = " << alpha 
                      << " f_trial = " << f_trial << " f_alpha = " << f_alpha << "\n";

      return std::make_tuple(false, alpha, f_trial, f_alpha, sigma_alpha);
    }

    ++numIter;

  } // end while signbit

  dbg_estimate << "alpha = " << alpha 
                  << " f_trial = " << f_trial << " f_alpha = " << f_alpha << "\n";

  return std::make_tuple(true, alpha, f_trial, f_alpha, sigma_alpha);
}

/**
 *  Do bisection algorithm to find intersection with yield surface
 *  Returns: solved, alpha, f_alpha, sigma_alpha
 */
std::tuple<bool, double, double, Vector6>
MohrCoulombBase::findIntersectionWithBisection(double alpha_in, double f_alpha_in,
                                               const Vector6& stress_trial,
                                               const Vector6& proj_direction) const
{
  double alpha_max = 0;
  double alpha_min = alpha_in;
  double f_min = f_alpha_in;

  bool solved = false;
  double alpha, f_alpha;
  Vector6 sigma_alpha;
   
  int numIter = 0;
  while (numIter < d_int.d_maxIter) {

    alpha = (alpha_min + alpha_max) / 2.0;

    sigma_alpha = stress_trial - alpha * proj_direction;

    f_alpha = computeYieldNormalized(sigma_alpha);

    dbg_calcplastic << " f_alpha = " << f_alpha << " alpha = " << alpha << "\n";

    if (std::abs(f_alpha) < d_int.d_yieldTol ||
        std::abs(alpha_min - alpha_max) < d_int.d_yieldTol) {
      solved = true;
      break;
    }

    if (std::signbit(f_alpha) == std::signbit(f_min)) {
      alpha_min = alpha;
      f_min = f_alpha;
    } else {
      alpha_max = alpha;
    }

    ++numIter;
  }

  return std::make_tuple(solved, alpha, f_alpha, sigma_alpha);
}
 
Vector6 
MohrCoulombBase::firstOrderStressUpdate(const Vector6& strainInc,
                                        const Matrix66& elasticTangent, 
                                        const Vector6& df_dsigma,
                                        const Vector6& dg_dsigma) const
{
  auto numerator = df_dsigma.transpose() * elasticTangent;
  double denominator = numerator * dg_dsigma;
  if (denominator < TINY) {
    std::cout << "**WARNING** denominator of plastic multiplier is very small."
                 " Some error may arise and results may be incorrect\n";
  }

  auto plasticStrainInc = (dg_dsigma * numerator * strainInc) / denominator;
  Vector6 dSigma = elasticTangent * (strainInc - plasticStrainInc);

  for (int i = 0; i < 6; i++) {
    if (!std::isfinite(dSigma(i))) {
      std::ostringstream err;
      err << "calcPlastic: dSigma not finite. dSigma = \n" << dSigma << "\n";
      throw InvalidValue(err.str(), __FILE__, __LINE__);
    }
  }

  return dSigma;
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
MohrCoulombBase::checkNorm(const Vector7& dSigma, double dP0Star,
                            const MohrCoulombState& initialState,
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
MohrCoulombBase::checkNormSloan(const Vector7& dSigma, double dP0Star,
                                 const MohrCoulombState& initialState,
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
MohrCoulombBase::correctDriftBeg(MohrCoulombState& state,
                                  const MohrCoulombState& stateOld) const
{
  dbg_doing << "Doing MohrCoulombBase::correctDriftBeg\n";

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

      dbg << __LINE__ << ":BaseMC::Delta Epsilon Plastic:\n" << dEps_p.transpose() << "\n";

      auto dSigma = (elasticTangent * dg_dsigma) * (-lambda);

      Vector7 zeros = Vector7::Zero();
      state.update(dEps_p, zeros, dSigma, 0);
    }

    if (numIter > 10) {
      dbg << __LINE__ << ":**WARNING** BaseMC::Drift Correction Procedure failed.\n";
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
MohrCoulombBase::correctDriftEnd(MohrCoulombState& state) const
{
  dbg_doing << "Doing MohrCoulombBase::correctDriftEnd\n";

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

      dbg << __LINE__ << ":BaseMC::Drift Correction, Iteration = " << numIter
          << " Function Value = " << fValue << "\n";

      std::tie(df_dsigma, dg_dsigma) = computeDfDsigma(state.stress);

      Matrix66 elasticTangent = calculateElasticTangentMatrix(K, G);

      auto numerator = df_dsigma.transpose() * elasticTangent;
      double denominator = numerator * dg_dsigma;

      double lambda = fValue / denominator;

      Vector6 dEps_p = df_dsigma * lambda;

      dbg << __LINE__ << ":BaseMC::Delta Epsilon Plastic:\n" << dEps_p << "\n";

      auto dSigma = (elasticTangent * dg_dsigma) * (-lambda);

      Vector7 zeros = Vector7::Zero();
      state.update(dEps_p, zeros, dSigma, 0);
    }
    if (numIter > 10) {
      dbg << __LINE__ << ":**WARNING** BaseMC::Drift Correction Procedure failed\n";
      correctDrift = false;
    }
  } while (correctDrift);
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
MohrCoulombBase::getTangentMatrix(const MohrCoulombState& state) const
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
MohrCoulombBase::calculateElastoPlasticTangentMatrix(
  const MohrCoulombState& state) const
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
MohrCoulombBase::computeG1(
  const MohrCoulombState& initialState, RetentionModel retentionModel,
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
MohrCoulombBase::computeG2(
  const MohrCoulombState& initialState, RetentionModel retentionModel,
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
