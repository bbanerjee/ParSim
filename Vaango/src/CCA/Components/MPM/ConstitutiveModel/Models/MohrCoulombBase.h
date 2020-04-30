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

#ifndef __MPM_CONSTITUTIVEMODEL_MODELS_BASE_MOHRCOULOMB__
#define __MPM_CONSTITUTIVEMODEL_MODELS_BASE_MOHRCOULOMB__

#include <CCA/Components/MPM/ConstitutiveModel/Models/MohrCoulombState.h>

#include <cmath>
#include <vector>

namespace Uintah {

// used to check for accuracies below machine accuracy, and to accept anything
// less as ok
constexpr double TINY = 1.0e-14;

// *WARNING** prevents algorithm to have more than 1e4 steps; cost: accuracy
// reduction in some unusual cases, but the whole thing will keep on going
constexpr double criticalStepSize = 1.0e-4;

/**
 * Integration intervals
 */

// constexpr int DivInt00[15] =
// {2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768};
// constexpr int DivInt01[15] = {2,4,6,8,10,12,14,16,18,20,22,24,26,28,30};
// constexpr int DivInt02[15] = {2,4,6,8,12,16,24,32,48,64,96,128,192,256,384};
// constexpr int DivInt03[15] = {2,4,6,8,12,16,24,32,48,64,96,128,192,256,384};
// constexpr int DivInt04[15] =
// {12,16,24,32,48,64,96,128,160,192,256,320,384,448,512};

// fastest so far
constexpr int DivInt[15] = { 32,  48,  64,  96,  128, 160, 192, 256,
                             320, 384, 448, 512, 608, 736, 992 };

// N+N-4
// constexpr int DivInt06[19] =
// {2,4,6,8,10,14,20,28,38,52,72,100,138,190,262,362,500,690,952};
// N+N-5
// constexpr int DivInt07[21] =
// {2,4,6,8,10,12,16,22,30,40,52,68,90,120,170,222,290,380,500,670,892};
// constexpr int DivInt08[15] =
// {28,32,40,52,64,78,94,120,154,200,240,290,330,380,440};
// constexpr int DivInt09[15] =
// {32,40,52,64,78,94,120,154,200,240,290,330,380,440,520};
// constexpr int DivInt10[15] =
// {32,36,40,46,52,60,70,82,96,112,130,150,176,220,380};
// constexpr int DivInt11[15] =
// {32,36,40,52,68,92,114,154,200,240,290,330,380,440,520};
// n=n-1*1.1, doesn't converge too often
// constexpr int DivInt12[15] =
// {32,36,40,44,50,56,62,68,76,84,92,102,112,124,136};
// n=n-1*1.2
// constexpr int DivInt13[15] =
// {32,38,46,56,66,80,96,114,138,166,198,238,286,344,412};


class MohrCoulombBase
{

public:
  MohrCoulombBase();
  MohrCoulombBase(double G, double K, double cohesion, double phi, double psi, double pMin = -1.0);
  MohrCoulombBase(const MohrCoulombBase* cm);
  MohrCoulombBase& operator=(const MohrCoulombBase&) = default;
  MohrCoulombBase(const MohrCoulombBase&) = default;
  virtual ~MohrCoulombBase() = default;

  void setModelParameters(double G, double K, double cohesion, double phi,
                          double psi, double pMin = -1.0);

  void setIntegrationParameters(int maxIterPegasus, double alfaCheck,
                                double alfaChange, double alfaRatio,
                                double yieldTol, double integrationTolerance,
                                double betaFactor, double minMeanStress,
                                double yieldLocTolerance,
                                SolutionAlgorithm solutionAlgorithm,
                                ToleranceMethod toleranceMethod,
                                DriftCorrection driftCorrection);

  virtual MohrCoulombState integrate(const Vector7& strainIncrement,
                                     const MohrCoulombState& initialState) = 0;

  virtual MohrCoulombState integrate(const Vector7& strainIncrement,
                             const MohrCoulombState& initialState, 
                             RegionType& region) = 0;

  Matrix67 getTangentMatrix(const MohrCoulombState& state) const;

  Vector7 computeG1(const MohrCoulombState& initialState,
                    RetentionModel retentionModel,
                    const std::vector<double>& retentionParameters) const;

  double computeG2(const MohrCoulombState& initialState,
                   RetentionModel retentionModel,
                   const std::vector<double>& retentionParameters) const;

protected:

  virtual bool checkYieldNormalized(const MohrCoulombState& state) const = 0;

  virtual double computeYieldNormalized(const Vector6& stress) const = 0;

  virtual std::tuple<Vector6, Vector6> computeDfDsigma(const Vector6& stress) const = 0;

  virtual 
  std::tuple<double, int> plasticRKME221(MohrCoulombState& state,
                                         const Vector7& epStrain) const = 0;

  virtual 
  std::tuple<double, int> plasticRK332(MohrCoulombState& state,
                                       const Vector7& epStrain) const = 0;

  virtual 
  std::tuple<double, int> plasticRKBog432(MohrCoulombState& point,
                                          const Vector7& epStrain) const = 0;

  virtual 
  std::tuple<double, int> plasticRK543(MohrCoulombState& state,
                                       const Vector7& epStrain) const = 0;

  virtual 
  std::tuple<double, int> plasticRKEng654(MohrCoulombState& state,
                                          const Vector7& epStrain) const = 0;

  virtual 
  std::tuple<double, int> plasticRKCK654(MohrCoulombState& state,
                                         const Vector7& epStrain) const = 0;

  virtual 
  std::tuple<double, int> plasticRKDP754(MohrCoulombState& state,
                                         const Vector7& epStrain) const = 0;

  virtual 
  std::tuple<double, int> plasticRKErr8544(MohrCoulombState& state,
                                           const Vector7& epStrain) const = 0;

  virtual 
  std::tuple<double, int> plasticExtrapol(MohrCoulombState& state,
                                          const Vector7& epStrain) const = 0;

  virtual
  double plasticMidpoint(MohrCoulombState& state, const Vector7& epStrain,
                         Vector7& absStress, int numIter) const = 0;


  // Elastic parameters
  struct Elastic
  {
    double d_G;  // Shear modulus
    double d_K;  // Bulk modulus
    double d_E;  // Young's modulus
    double d_nu; // Poisson's ratio

    void set(double G, double K)
    {
      d_G = G;
      d_K = K;
      d_E = 9.0 * K * G / (G + 3.0 * K);
      d_nu = (3.0 * K - 2.0 * G) / (2.0 * G + 6.0 * K);
    }
  };

  // Mohr - Coulomb parameters
  struct Yield
  {
    double d_pMin;
    double d_cohesion;
    double d_phi;
    double d_sin_phi;
    double d_cos_phi;

    // Rounded Mohr-Coulomb parameters
    double d_alpha;  // For the M parameter
    double d_alpha4; // alpha^4;

    void set(double c, double phi, double pMin = -1)
    {
      d_cohesion = c;
      d_phi = phi * M_PI / 180.0;
      d_sin_phi = std::sin(d_phi);
      d_cos_phi = std::cos(d_phi);
      d_alpha = (3.0 - d_sin_phi) / (3.0 + d_sin_phi);
      d_alpha4 = d_alpha * d_alpha * d_alpha * d_alpha;
      if (pMin != -1) {
        d_pMin = pMin;
      } else {
        d_pMin = c * d_cos_phi / d_sin_phi;
      }
    }
  };

  // If the flow rule is not associated;
  struct PotentialFn
  {
    double d_psi;
    double d_sin_psi;
    double d_cos_psi;

    void set(double psi)
    {
      d_psi = psi * M_PI / 180.0;
      d_sin_psi = std::sin(d_psi);
      d_cos_psi = std::cos(d_psi);
    }
  };

  // Integration parameters
  struct IntegrationParameters
  {

    // Elastic parameters, used in the Pegasus algorithm
    int d_maxIter;
    double d_alfaCheck;
    double d_alfaChange;
    double d_alfaRatio;

    // Yield tolerance
    double
      d_yieldTol; // tolerance for the Yield locus (relatively defined, used
                  // to check if the stress state is on the YL)
    double d_integrationTol;
    double d_betaFactor; // safety factor
    double d_minMeanStress;
    double d_suctionTol; // used in checking suction yield locus. advisible not
                         // to modify...

    // Algorithms
    DriftCorrection d_driftCorrection;
    ToleranceMethod d_tolMethod;
    SolutionAlgorithm d_solutionAlgorithm;

    void setDefaults()
    {
      d_maxIter = 200;
      d_alfaCheck = 1;
      d_alfaChange = 0.05;
      d_alfaRatio = 10;

      d_yieldTol = 1e-6;
      d_integrationTol = 0.01;
      d_betaFactor = 0.9;
      d_minMeanStress = -1.0e8;
      d_suctionTol = 1e-8;

      d_driftCorrection = DriftCorrection::CORRECTION_AT_END;
      d_tolMethod = ToleranceMethod::SLOAN;
      d_solutionAlgorithm =
        SolutionAlgorithm::RUNGE_KUTTA_SECOND_ORDER_MODIFIED_EULER;
    }

    void setDefaults(const Yield& yield)
    {
      setDefaults();
      if (yield.d_sin_phi > 0) {
        d_minMeanStress = yield.d_cohesion * yield.d_cos_phi / yield.d_sin_phi;
      }
    }
  };

  // Model parameters
  bool d_nonAssociated;
  Elastic d_elastic;
  Yield d_yield;
  PotentialFn d_potential;

  // Integration parameters
  IntegrationParameters d_int;

  void calcElastic(const Vector7& strain, const MohrCoulombState& initialPoint,
                   MohrCoulombState& finalPoint) const;

  Vector6 calcStressIncElast(double nu0, const Vector6& s0, const Vector7& eps0,
                             const Vector7& deps) const;

  Matrix67 calculateElasticTangentMatrix(const MohrCoulombState& state) const;
  Matrix66 calculateElasticTangentMatrix(double K, double G) const;

  Matrix67 calculateElastoPlasticTangentMatrix(
    const MohrCoulombState& state) const;

  bool checkGradient(const MohrCoulombState& initialState,
                     const MohrCoulombState& finalState) const;

  double computeDsigmaDotDf(const Vector6& s, const Vector6& ds, Vector6& dF,
                      double suction, double dsuction) const;


  inline double firstInvariant(const Vector6& s) const
  {
    double I1 = s(0) + s(1) + s(2);
    return I1;
  }

  inline double secondInvariant(const Vector6& s) const
  {
    double I2 = s(0) * s(1) + s(1) * s(2) + s(2) * s(0) - s(3) * s(3) -
                s(4) * s(4) - s(5) * s(5);
    return I2;
  }

  inline double thirdInvariant(const Vector6& s) const
  {
    double I3 = s(0) * s(1) * s(2) + 2 * s(3) * s(4) * s(5) -
                s(0) * s(5) * s(5) - s(1) * s(4) * s(4) - s(2) * s(3) * s(3);
    return I3;
  }

  inline double firstDevInvariant(const Vector6& s) const { return 0.0; }

  inline double secondDevInvariant(const Vector6& s) const
  {
    double J2 = ((s(0) - s(1)) * (s(0) - s(1)) + (s(0) - s(2)) * (s(0) - s(2)) +
                 (s(1) - s(2)) * (s(1) - s(2))) /
                  6.0 +
                (s(3) * s(3) + s(4) * s(4) + s(5) * s(5));
    if (std::abs(J2) < TINY) {
      J2 = TINY;
    }
    return J2;
  }

  inline double thirdDevInvariant(const Vector6& s) const
  {
    double I1 = firstInvariant(s);
    double I2 = secondInvariant(s);
    double I3 = thirdInvariant(s);
    double J3 = I1 * I1 * I1 * 2.0 / 27.0 - I1 * I2 / 3.0 + I3;

    return J3;
  }

  inline std::tuple<double, double> vonMisesStress(const Vector6& s) const
  {
    double J2 = secondDevInvariant(s);
    double vmStress = std::sqrt(3.0 * J2);
    if (std::abs(vmStress) < TINY) {
      vmStress = TINY;
    }
    return std::make_tuple(J2, vmStress);
  }

  std::tuple<Vector7, Vector7>
  findIntersectionUnloading(const Vector7& strainIncrement,
                            const MohrCoulombState& initialState);

  std::tuple<Vector7, Vector7>
  findIntersection(const Vector7& strainIncrement,
                        const MohrCoulombState& initialState) const;

  double findYieldAlpha(const Vector3& state, const Vector6& s0,
                        const Vector7& eps0, const Vector7& deps);

  double findYieldModified(const Vector3& state, const Vector6& s0,
                           const Vector7& eps0, const Vector7& deps) const;

  double computeNu(const Vector6& s, const Vector3& state,
                   double suction) const;

  double calculatePlastic(const Vector7& purelyPlasticStrain,
                          MohrCoulombState& state) const;

  std::tuple<Vector6, double> calcPlastic(const MohrCoulombState& state,
                                          const Vector7& strainInc) const;

  Vector6 projectTrialStressToYieldSurface(const Vector6& strainInc,
                                           const Vector6& stress_old, 
                                           const Matrix66& elasticTangent, 
                                           const Vector6& df_dsigma,
                                           const Vector6& dg_dsigma,
                                           const Vector6& stress_trial, 
                                           const Vector6& proj_direction) const;

  std::tuple<bool, double, double, double, Vector6>
    estimateInitialBisectionParameter(const Vector6& stress_old,
                                      const Vector6& stress_trial,
                                      const Vector6& proj_direction) const;

  std::tuple<bool, double, double, Vector6>
    findIntersectionWithBisection(double alpha_in, double f_alpha_in,
                                  const Vector6& stress_trial,
                                  const Vector6& proj_direction) const;

  Vector6 firstOrderStressUpdate(const Vector6& strainInc,
                                 const Matrix66& elasticTangent, 
                                 const Vector6& df_dsigma,
                                 const Vector6& dg_dsigma) const;

  void getParamRKME221(Eigen::Matrix<double, 2, 2>& A,
                       Eigen::Matrix<double, 2, 1>& B,
                       Eigen::Matrix<double, 2, 1>& BRes,
                       Eigen::Matrix<double, 2, 1>& C) const;

  void getParamRK332(Eigen::Matrix<double, 3, 3>& A,
                     Eigen::Matrix<double, 3, 1>& B,
                     Eigen::Matrix<double, 3, 1>& BRes,
                     Eigen::Matrix<double, 3, 1>& C) const;

  void getParamRKBog432(Eigen::Matrix<double, 4, 4>& A,
                        Eigen::Matrix<double, 4, 1>& B,
                        Eigen::Matrix<double, 4, 1>& BRes,
                        Eigen::Matrix<double, 4, 1>& C) const;

  void getParamRK543(Eigen::Matrix<double, 5, 5>& A,
                     Eigen::Matrix<double, 5, 1>& B,
                     Eigen::Matrix<double, 5, 1>& BRes,
                     Eigen::Matrix<double, 5, 1>& C) const;

  void getParamRKEng654(Eigen::Matrix<double, 6, 6>& A,
                        Eigen::Matrix<double, 6, 1>& B,
                        Eigen::Matrix<double, 6, 1>& BRes,
                        Eigen::Matrix<double, 6, 1>& C) const;

  void getParamRKCK654(Eigen::Matrix<double, 6, 6>& A,
                       Eigen::Matrix<double, 6, 1>& B,
                       Eigen::Matrix<double, 6, 1>& BRes,
                       Eigen::Matrix<double, 6, 1>& C) const;

  void getParamRKDP754(Eigen::Matrix<double, 7, 7>& A,
                       Eigen::Matrix<double, 7, 1>& B,
                       Eigen::Matrix<double, 7, 1>& BRes,
                       Eigen::Matrix<double, 7, 1>& C) const;

  void getParamRKErr8544(Eigen::Matrix<double, 8, 8>& A,
                         Eigen::Matrix<double, 8, 1>& B,
                         Eigen::Matrix<double, 8, 1>& BRes,
                         Eigen::Matrix<double, 8, 1>& C,
                         Eigen::Matrix<double, 7, 1>& ErrorCoef) const;

  double checkNorm(const Vector7& dSigma, double dP0Star,
                   const MohrCoulombState& initialState,
                   const Vector7& dError) const;

  double checkNormSloan(const Vector7& dSigma, double dP0Star,
                        const MohrCoulombState& initialState,
                        const Vector7& dError) const;

  void correctDriftBeg(MohrCoulombState& state,
                       const MohrCoulombState& stateOld) const;

  void correctDriftEnd(MohrCoulombState& state) const;

};

} // end namespace Uintah

#endif //__MPM_CONSTITUTIVEMODEL_MODELS_BASE_MOHRCOULOMB__
