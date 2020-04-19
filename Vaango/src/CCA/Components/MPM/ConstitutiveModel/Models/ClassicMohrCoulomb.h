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
#ifndef __MPM_CONSTITUTIVEMODEL_MODELS_SHENG_MOHRCOULOMB__
#define __MPM_CONSTITUTIVEMODEL_MODELS_SHENG_MOHRCOULOMB__

#include "StateMohrCoulomb.h"

#include <cmath>
#include <vector>

namespace Uintah {

constexpr double TINY = 1.0e-14;

enum class DriftCorrection
{
  NO_CORRECTION = 1,
  CORRECTION_AT_BEGIN = 2,
  CORRECTION_AT_END = 3
};

std::ostream&
operator<<(std::ostream& out, const DriftCorrection& dc)
{
  static std::array<const char*, 3> names = {
    { "NO_CORRECTION", "CORRECTION_AT_BEGIN", "CORRECTION_AT_END" }
  };
  out << names[static_cast<int>(dc)];
  return out;
}

enum class ToleranceMethod
{
  EPUS_RELATIVE_ERROR = 0,
  SLOAN = 1
};

enum class SolutionAlgorithm
{
  RUNGE_KUTTA_SECOND_ORDER_MODIFIED_EULER = 1,
  RUNGE_KUTTA_THIRD_ORDER_NYSTROM = 2,
  RUNGE_KUTTA_THIRD_ORDER_BOGACKI = 3,
  RUNGE_KUTTA_FOURTH_ORDER = 4,
  RUNGE_KUTTA_FIFTH_ORDER_ENGLAND = 5,
  RUNGE_KUTTA_FIFTH_ORDER_CASH = 6,
  RUNGE_KUTTA_FIFTH_ORDER_DORMAND = 7,
  RUNGE_KUTTA_FIFTH_ORDER_BOGACKI = 8,
  EXTRAPOLATION_BULIRSCH = 9
};

std::ostream&
operator<<(std::ostream& out, const SolutionAlgorithm& sa)
{
  static std::array<const char*, 9> names = {
    { "RUNGE_KUTTA_SECOND_ORDER_MODIFIED_EULER",
      "RUNGE_KUTTA_THIRD_ORDER_NYSTROM", "RUNGE_KUTTA_THIRD_ORDER_BOGACKI",
      "RUNGE_KUTTA_FOURTH_ORDER", "RUNGE_KUTTA_FIFTH_ORDER_ENGLAND",
      "RUNGE_KUTTA_FIFTH_ORDER_CASH", "RUNGE_KUTTA_FIFTH_ORDER_DORMAND",
      "RUNGE_KUTTA_FIFTH_ORDER_BOGACKI", "EXTRAPOLATION_BULIRSCH" }
  };
  out << names[static_cast<int>(sa)];
  return out;
}

enum class RetentionModel
{
  STATE_SURFACE = 1,
  VAN_GENUCHTEN = 2,
  GALLIPOLI = 3
};


/** 
 * This ia a classic Mohr-Coulomb model
 * integration works based on the principal stress space
 * note that the strain must be rotated into the same coordinate system as
 * stress
*/

class ClassicMohrCoulomb
{

public:
  ClassicMohrCoulomb();
  ClassicMohrCoulomb(double G, double K, double cohesion, double phi, double psi);

  ClassicMohrCoulomb(const ClassicMohrCoulomb&) = delete;
  ClassicMohrCoulomb& operator=(const ClassicMohrCoulomb&) = delete;

  ~ClassicMohrCoulomb() = default;

  void setModelParameters(double G, double K, double cohesion, double phi,
                          double psi);

  void setIntegrationParameters(int maxIterPegasus, double integrationTolerance,
                                double betaFactor, double yieldLocTolerance,
                                SolutionAlgorithm solutionAlgorithm,
                                ToleranceMethod toleranceMethod,
                                DriftCorrection driftCorrection);

  StateMohrCoulomb integrate(const Vector7& strainIncrement,
                             const StateMohrCoulomb& initialState);

  void Integrate(double* StrainIncrement, const StateMohrCoulomb& state_init);
  void Integrate(double* StrainIncrement, double SuctionIncrement,
                 const StateMohrCoulomb& state_init, double* StressIncrement,
                 double P0StarIncrement, Vector6& plasticStrainIncrement);
  // void IntegrateMC (double* StrainIncrement,const StateMohrCoulomb& state_init);
  // void IntegrateMCI (double* StrainIncrement,const StateMohrCoulomb& state_init);
  // int IntegrateMCIA (double* StrainIncrement,const StateMohrCoulomb& state_init);
  void IntegrateMCIClassic(double* StrainIncrement, const StateMohrCoulomb& state_init,
                           int* Region);
  void IntegrateConst(double* StrainIncrement, const StateMohrCoulomb& state_init,
                      int StepNo, int Method);
  double CalcStressElast(double nu0, double* s0, double* eps0, double* deps,
                         double* ds);
  double CalcElastic(double* Strain, const StateMohrCoulomb& state_init,
                     BBMPoint* FinalPoint);
  void CalcStressElastM(double* deps, double* ds);
  void FindElStrGradPQ(double nu0, double* s0, double* eps0, double* deps,
                       double* ds);
  bool CheckGradient(const StateMohrCoulomb& state_init, BBMPoint* FinalPoint);
  void FindElStrGrad(double nu0, double* s0, double* eps0, double* deps,
                     double* ds);
  double CalculatePZero(const StateMohrCoulomb& state);
  bool CheckYield(const StateMohrCoulomb& state);
  bool CheckYieldNormalised(const StateMohrCoulomb& state);
  void CheckYield(const Vector3& state, double* s, double suction, double* FValue);
  void CheckYieldNormalised(const Vector3& state, double* s, double suction,
                            double* FValue);
  bool CheckIfPlastic(const StateMohrCoulomb& state);
  double ComputeYieldFunction(const StateMohrCoulomb& state);
  double ComputeYieldFunctionNN(const StateMohrCoulomb& state);
  // bool CheckYield (double *state, double* s, double suction);
  void FindYieldOriginal(const Vector3& state, double* s0, double* eps0, double* deps,
                         double* a);
  // double FindYieldOriginal (double *state, double*s0, double* eps0, double*
  // deps);
  void FindYieldModified(const Vector3& state, double* s0, double* eps0, double* deps,
                         double* a);
  double ComputeNu(double* s, const Vector3& state, double suction);
  double FindGradient(const Vector3& state, double* s, double* ds, double* dF,
                      double suction, double dsuction);
  double FindGradientPQ(const StateMohrCoulomb& state, double* ds, double* dF,
                        double dsuction);
  void MoveYieldaBit(const Vector3& state, double* s, double* ds, double* eps0,
                     double* deps, double* gradient, double F0);
  void FindYieldYield(const Vector3& state, double* s0, double* eps0, double* deps,
                      double* a);
  void FindIntersectionUnloading(double* StrainIncrement,
                                 const StateMohrCoulomb& state_init, double* PurelyElastic,
                                 double* PurelyPlastic);
  void FindIntersection(double* StrainIncrement, const StateMohrCoulomb& state_init,
                        double* PurelyElasticStrain,
                        double* PurelyPlasticStrain);
  void PaintLocus(const Vector3& state, double suction, int Max);
  bool FindIntersectionPlastic(double* PurelyPlasticStrain,
                               BBMPoint* FinalPoint,
                               double* IntersectionStrain);
  bool CheckStayInCorner(BBMPoint* FinalPoint, double* IntersectionStrain);
  void ComputeG1(const StateMohrCoulomb& state_init, int RetentionModel,
                 double* RetentionParameters, double* G1);
  void ComputeG2(const StateMohrCoulomb& state_init, int RetentionModel,
                 double* RetentionParameters, double* G2);

  double Getk();
  double GetLambdaZero();
  double Getr();
  double GetBeta();
  double GetKappaP();
  double Getpc();
  void read();
  void write();

  // *********************************** Plastic Procedures below
  // ******************************************************

  void GetTangentMatrixPQ(const StateMohrCoulomb& state, Matrix67& DEP);
  void CalculateElastoPlasticTangentMatrixPQ(const StateMohrCoulomb& state, Matrix67& DEP);
  void CalculateElasticTangentMatrixPQ(const StateMohrCoulomb& state, Matrix67& DEP);
  void GetTangentMatrix(const StateMohrCoulomb& state, Matrix67& DEP);
  void CalculateElastoPlasticTangentMatrix(const StateMohrCoulomb& state, Matrix67& DEP);
  void CalculateElasticTangentMatrix(const StateMohrCoulomb& state, Matrix67& DEP);
  void GetDerivative(double MeanStress, double ShearStress, double suction,
                     double PZero, const Vector3& state, double* deriv);
  double GetLambda(double* deriv, double stresspq[3], double strainpq[3]);
  double PlasticEuler(
    const StateMohrCoulomb& state, const Vector7& epStrain, Vector7& absStress,
    int NumberIterations); // returns elapsed time of computations
  double RungeKutta(double A[][8], double* B, double* BRes, double* C,
                    const StateMohrCoulomb& state, const Vector7& epStrain, Vector7& absStress,
                    int* NumberIter, double MethodOrder, int MethodSteps,
                    bool ErrorEstimate);
  double RungeKuttaEqualStep(double A[][8], double* B, double* BRes, double* C,
                             const StateMohrCoulomb& state, const Vector7& epStrain,
                             Vector7& absStress, double* RelError,
                             int NumberIter, double MethodOrder,
                             int MethodSteps, bool ErrorEstimate);
  double RungeKuttaExtrapol(double A[][8], double* B, double* BRes, double* C,
                            const StateMohrCoulomb& state, const Vector7& epStrain,
                            Vector7& absStress, int* NumberIter,
                            double MethodOrder, int MethodSteps,
                            bool ErrorEstimate);
  // Runge Kutta schemes
  double CalculatePlastic(double* PurelyPlasticStrain, const StateMohrCoulomb& state);
  double CalculatePlasticMC(double* PurelyPlasticStrain, const StateMohrCoulomb& state);
  int CalculatePlasticMC_NEW(double* PurelyPlasticStrain, const StateMohrCoulomb& state);
  double CalculatePlasticConst(double* PurelyPlasticStrain, const StateMohrCoulomb& state,
                               int StepNo);
  double PlasticRKErr8544(const StateMohrCoulomb& state, const Vector7& epStrain, Vector7& absStress,
                          int* NumberIter); // Bogacki - Shimpine
  double PlasticRKDP754(const StateMohrCoulomb& state, const Vector7& epStrain, Vector7& absStress,
                        int* NumberIter); // Dormand Prince
  double PlasticRKCK654(const StateMohrCoulomb& state, const Vector7& epStrain, Vector7& absStress,
                        int* NumberIter); // Cash - Karp
  double PlasticRKEng654(const StateMohrCoulomb& state, const Vector7& epStrain, Vector7& absStress,
                         int* NumberIter); // England as given by Sloan
  double PlasticRK543(const StateMohrCoulomb& state, const Vector7& epStrain, Vector7& absStress,
                      int* NumberIter); // 4th order with 3rd ord estimate
  double PlasticRK332(const StateMohrCoulomb& state, const Vector7& epStrain, Vector7& absStress,
                      int* NumberIter); // 3rd order R-K scheme
  double PlasticRKBog432(
    const StateMohrCoulomb& state, const Vector7& epStrain, Vector7& absStress,
    int* NumberIter); // Bogacki - Shimpine 3rd order Runge Kutta scheme
  double PlasticRKME221(const StateMohrCoulomb& state, const Vector7& epStrain, Vector7& absStress,
                        int* NumberIter); // Modified Euler
  double PlasticRKNoExTry(const StateMohrCoulomb& state, const Vector7& epStrain, Vector7& absStress,
                          int* NumberIter); // using not in an extrapolation way
  // Extrapolation Schemes
  double PlasticExtrapol(const StateMohrCoulomb& state, const Vector7& epStrain, Vector7& absStress,
                         int* NumberIter);
  double RKExtrapolation(double A[][8], double* B, double* BRes, double* C,
                         const StateMohrCoulomb& state, const Vector7& epStrain, Vector7& absStress,
                         const StateMohrCoulomb& state_old, double* RelError, int* NumberIter,
                         double MethodOrder, int MethodSteps,
                         bool ErrorEstimate);
  double PlasticMidpoint(const StateMohrCoulomb& state, const Vector7& epStrain, Vector7& absStress,
                         int* NumberIter);
  double PlasticMidpointGallipoli(const StateMohrCoulomb& state, const Vector7& epStrain,
                                  Vector7& absStress, int* NumberIter);
  double CheckNorm(double* DSigma, double DPZeroStar, const StateMohrCoulomb& state_init,
                   double* DError); // returns RError
  double CheckNormSloan(double* DSigma, double DPZeroStar,
                        const StateMohrCoulomb& state_init,
                        double* DError); // returns RError
  void CorrectDrift(const StateMohrCoulomb& state);
  void CorrectDriftOnce(const StateMohrCoulomb& state);
  void CorrectDriftMC(const StateMohrCoulomb& state);
  void CorrectDriftBeg(BBMPoint* EndPoint, const StateMohrCoulomb& stateOld);
  void ReturnImplicit(const StateMohrCoulomb& state, int* Region);

  double ComputeYieldFunctionEigen(double Eigen[3]);
  double ComputeYieldFunctionEigenNorm(double Eigen[3]);
  double FindPercentage(double* Strain, const StateMohrCoulomb& state_init);
  double FindPercentageBis(double* Strain, const StateMohrCoulomb& state_init);
  void CorrectDriftEdge(const StateMohrCoulomb& state_init, int Region);
  void CorrectDriftEdgeJC(const StateMohrCoulomb& state_init, int Region);
  void CorrectDriftTip(const StateMohrCoulomb& state_init);

  // Used Procedures before, not updated anymore, though, mostly, working.
  // void DriftCorrect (const StateMohrCoulomb& state, const Vector7& epStrain, Vector7& dSigma,
  // double* Lambda, double& dP0Star, double* FValue);
  // void CorrectDriftBeg (const StateMohrCoulomb& state, const Vector7& epStrain, Vector7& dSigma,
  // double* Lambda, double& dP0Star, double* FValue);
  // double Plastic (const StateMohrCoulomb& state, const Vector7& epStrain, Vector7& absStress, int*
  // NumberIter);		//returns elapsed time of computations
  // double PlasticNewSlow (const StateMohrCoulomb& state, const Vector7& epStrain, double*
  // AbsStress, int* NumberIterations);	//returns elapsed time of computations
  // double PlasticRKErr6 (const StateMohrCoulomb& state, const Vector7& epStrain, Vector7& absStress,
  // int* NumberIterations);	//returns elapsed time of computations
  // double PlasticRKErr75 (const StateMohrCoulomb& state, const Vector7& epStrain, double*
  // AbsStress, int* NumberIter);

  // double PlasticRK5Err4_2 (const StateMohrCoulomb& state, const Vector7& epStrain, double*
  // AbsStress, int* NumberIter);
  // double PlasticRK4Err3  (const StateMohrCoulomb& state, const Vector7& epStrain, double*
  // AbsStress, int* NumberIter);
  // double PlasticRK4Err3v2  (const StateMohrCoulomb& state, const Vector7& epStrain, double*
  // AbsStress, int* NumberIter);
  // double PlasticRK3Err2  (const StateMohrCoulomb& state, const Vector7& epStrain, double*
  // AbsStress, int* NumberIter);
  // double PlasticRK3Err2v2  (const StateMohrCoulomb& state, const Vector7& epStrain, double*
  // AbsStress, int* NumberIter);
  // double PlasticRK3Err2v3  (const StateMohrCoulomb& state, const Vector7& epStrain, double*
  // AbsStress, int* NumberIter);
  // double PlasticRKSloan (const StateMohrCoulomb& state, const Vector7& epStrain, double*
  // AbsStress, int* NumberIterations);	//returns elapsed time of computations
  // double PlasticMidpointC (const StateMohrCoulomb& state, const Vector7& epStrain, double*
  // AbsStress, int* NumberIter);
  // double PlasticMidpointCN (const StateMohrCoulomb& state, const Vector7& epStrain, double*
  // AbsStress, int* NumberIter);
  // double PlasticMidpointC4 (const StateMohrCoulomb& state, const Vector7& epStrain, double*
  // AbsStress, int* NumberIter);
  // double PlasticMidpointC6 (const StateMohrCoulomb& state, const Vector7& epStrain, double*
  // AbsStress, int* NumberIter);

protected:

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
    double d_cohesion;
    double d_phi;
    double d_sin_phi;
    double d_cos_phi;

    void set(double c, double phi)
    {
      d_cohesion = c;
      d_phi = phi * M_PI / 180.0;
      d_sin_phi = std::sin(phi);
      d_cos_phi = std::cos(phi);
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
      d_sin_psi = std::sin(psi);
      d_cos_psi = std::cos(psi);
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
    double d_yieldTol; 
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

  bool checkYieldNormalized(const StateMohrCoulomb& state) const;

private:

  double computeYieldNormalized(const Vector6& stress) const;

  void calcElastic(const Vector7& strain, const StateMohrCoulomb& initialPoint,
                   StateMohrCoulomb& finalPoint) const;

  Vector6 calcStressIncElast(double nu0, const Vector6& s0, const Vector7& eps0,
                             const Vector7& deps) const;

  bool checkGradient(const StateMohrCoulomb& initialState,
                     const StateMohrCoulomb& finalState) const;

  Matrix67 calculateElasticTangentMatrix(const StateMohrCoulomb& state) const;
  Matrix66 calculateElasticTangentMatrix(double K, double G) const;

  Matrix67 calculateElastoPlasticTangentMatrix(
    const StateMohrCoulomb& state) const;

  double findGradient(const Vector6& s, const Vector6& ds, Vector6& dF,
                      double suction, double dsuction) const;

  std::tuple<Vector6, Vector6> computeDfDsigma(const Vector6& stress) const;

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

  void findIntersectionUnloading(const Vector7& strainIncrement,
                                 const StateMohrCoulomb& initialState,
                                 Vector7& purelyElastic,
                                 Vector7& purelyPlastic);

  void findIntersection(const Vector7& strainIncrement,
                        const StateMohrCoulomb& initialState,
                        Vector7& elasticStrainInc,
                        Vector7& plasticStrainInc) const;

  double findYieldAlpha(const Vector3& state, const Vector6& s0,
                        const Vector7& eps0, const Vector7& deps);

  double findYieldModified(const Vector3& state, const Vector6& s0,
                           const Vector7& eps0, const Vector7& deps) const;

  double computeNu(const Vector6& s, const Vector3& state,
                   double suction) const;

  double calculatePlastic(const Vector7& purelyPlasticStrain,
                          StateMohrCoulomb& state) const;

  int calcPlastic(const StateMohrCoulomb& state, const Vector7& epStrainInc,
                  Vector6& dSigma, Vector6& dEps_p, double& dP0Star) const;

  std::tuple<double, int> plasticRKME221(StateMohrCoulomb& state,
                                         const Vector7& epStrain) const;

  std::tuple<double, int> plasticRK332(StateMohrCoulomb& state,
                                       const Vector7& epStrain) const;

  std::tuple<double, int> plasticRKBog432(StateMohrCoulomb& point,
                                          const Vector7& epStrain) const;

  std::tuple<double, int> plasticRK543(StateMohrCoulomb& state,
                                       const Vector7& epStrain) const;

  std::tuple<double, int> plasticRKEng654(StateMohrCoulomb& state,
                                          const Vector7& epStrain) const;

  std::tuple<double, int> plasticRKCK654(StateMohrCoulomb& state,
                                         const Vector7& epStrain) const;

  std::tuple<double, int> plasticRKDP754(StateMohrCoulomb& state,
                                         const Vector7& epStrain) const;

  std::tuple<double, int> plasticRKErr8544(StateMohrCoulomb& state,
                                           const Vector7& epStrain) const;

  template <int Order, int Steps>
  std::tuple<double, int> doRungeKutta(
    const Eigen::Matrix<double, Steps, Steps>& AA,
    const Eigen::Matrix<double, Steps, 1>& BB,
    const Eigen::Matrix<double, Steps, 1>& BRes,
    const Eigen::Matrix<double, Steps, 1>& CC, StateMohrCoulomb& state,
    const Vector7& epStrain, bool errorEstimate) const;

  template <int Order, int Steps>
  std::tuple<double, int> doRungeKuttaErr(
    const Eigen::Matrix<double, Steps, Steps>& AA,
    const Eigen::Matrix<double, Steps, 1>& BB,
    const Eigen::Matrix<double, Steps, 1>& BRes,
    const Eigen::Matrix<double, Steps, 1>& CC,
    const Eigen::Matrix<double, Steps - 1, 1>& ErrCoef, StateMohrCoulomb& state,
    const Vector7& epStrain, bool errorEstimate) const;

  double checkNorm(const Vector7& dSigma, double dP0Star,
                   const StateMohrCoulomb& initialState,
                   const Vector7& dError) const;

  double checkNormSloan(const Vector7& dSigma, double dP0Star,
                        const StateMohrCoulomb& initialState,
                        const Vector7& dError) const;

  void correctDriftBeg(StateMohrCoulomb& state,
                       const StateMohrCoulomb& stateOld) const;

  void correctDriftEnd(StateMohrCoulomb& state) const;

  double plasticMidpoint(StateMohrCoulomb& state, const Vector7& epStrain,
                         Vector7& absStress, int numIter) const;

  std::tuple<double, int> plasticExtrapol(StateMohrCoulomb& state,
                                          const Vector7& epStrain) const;
};

} // end namespace Uintah

#endif //__MPM_CONSTITUTIVEMODEL_MODELS_SHENG_MOHRCOULOMB__
