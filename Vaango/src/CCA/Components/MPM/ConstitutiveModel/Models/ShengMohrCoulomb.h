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

#include "BBMMatrix.h"
#include "BBMPoint.h"


#include <cmath>

/* this Mohr-Coulomb like model uses a rounded Mohr-Coulomb surface, see
   eq 13 in Sheng D, Sloan SW & Yu HS Computational Mechanics 26:185-196 (2000)
   Springer
*/

namespace Uintah {


class ShengMohrCoulomb {

public:

  // used to check for accuracies below machine accuracy, and to accept anything
  // less as ok
  constexpr double TINY = 1.0e-14;

  enum class DriftCorrection
  {
    NO_CORRECTION = 1,
    CORRECTION_AT_END = 3;
  };

  enum class ToleranceMethod
  {
    EPUS_RELATIVE_ERROR = 0,
    SLOAN = 1
  };

  enum class SolutionAlgorithm
  {
    MODIFIED_EULER = 1,
    RUNGE_KUTTA_332 = 2,
    NYSTROM = 3
  };

  ShengMohrCoulomb();
  ShengMohrCoulomb(double G, double K, double cohesion, double phi, double psi);

  ShengMohrCoulomb(const ShengMohrCoulomb& ) = delete;
  ShengMohrCoulomb& operator=(const ShengMohrCoulomb& ) = delete;
  ~ShengMohrCoulomb() = default;

  void setModelParameters(double G, double K, double cohesion, double phi, 
                          double psi);

  void setIntegrationParameters(int maxIterPegasus,
                                double integrationTolerance, 
                                double betaFactor, 
                                double yieldLocTolerance,
                                SolutionAlgorithm solutionAlgorithm,
                                ToleranceMethod toleranceMethod,
                                DriftCorrection driftCorrection) ;

  void integrate(double *strainIncrement, BBMPoint *initialPoint);
  void integrate(double *strainIncrement, double suctionIncrement,
                 BBMPoint *initialPoint, double *stressIncrement,
                 double P0StarIncrement, double *plasticStrainIncrement);
  void integrateConst(double *strainIncrement, BBMPoint *initialPoint,
                      int stepNo, int method);
  double calcStressElast(double nu0, double *s0, double *eps0, double *deps,
                         double *ds);
  double calcElastic(double *strain, BBMPoint *initialPoint,
                     BBMPoint *finalPoint);
  void calcStressElastM(double *deps, double *ds);
  void findElStrGradPQ(double nu0, double *s0, double *eps0, double *deps,
                       double *ds);
  bool checkGradient(BBMPoint *initialPoint, BBMPoint *finalPoint);
  void findElStrGrad(double nu0, double *s0, double *eps0, double *deps,
                     double *ds);
  double calculatepZero(BBMPoint *point);
  bool checkYield(BBMPoint *point);
  bool checkYieldNormalised(BBMPoint *point);
  void checkYield(double *state, double *s, double suction, double *fValue);
  void checkYieldNormalised(double *state, double *s, double suction,
                            double *fValue);
  bool checkIfPlastic(BBMPoint *point);
  double computeYieldFunction(BBMPoint *point);
  double computeYieldFunctionNN(BBMPoint *point);
  // bool checkYield (double *state, double* s, double suction);
  void findYieldOriginal(double *state, double *s0, double *eps0, double *deps,
                         double *a);
  // double findYieldOriginal (double *state, double*s0, double* eps0, double*
  // deps);
  void findYieldModified(double *state, double *s0, double *eps0, double *deps,
                         double *a);
  double ComputeNu(double *s, double *state, double suction);
  double findGradient(double *state, double *s, double *ds, double *dF,
                      double suction, double dsuction);
  double findGradientPQ(BBMPoint *point, double *ds, double *dF,
                        double dsuction);
  void moveYieldaBit(double *state, double *s, double *ds, double *eps0,
                     double *deps, double *gradient, double F0);
  void findYieldYield(double *state, double *s0, double *eps0, double *deps,
                      double *a);
  void findIntersectionUnloading(double *strainIncrement,
                                 BBMPoint *initialPoint, double *purelyElastic,
                                 double *purelyPlastic);
  void findIntersection(double *strainIncrement, BBMPoint *initialPoint,
                        double *purelyElasticStrain,
                        double *purelyPlasticStrain);
  void paintLocus(double *state, double suction, int Max);
  void computeG1(BBMPoint *initialPoint, int retentionModel,
                 double *retentionParameters, double *G1);
  void computeG2(BBMPoint *initialPoint, int retentionModel,
                 double *retentionParameters, double *G2);

  double getk();
  double getLambdaZero();
  double getr();
  double getBeta();
  double getKappaP();
  double getpc();
  void read();
  void write();

  // *********************************** Plastic Procedures below
  // ******************************************************

  void getTangentMatrixPQ(BBMPoint *point, BBMMatrix *DEP);
  void calculateElastoPlasticTangentMatrixPQ(BBMPoint *point, BBMMatrix *DEP);
  void calculateElasticTangentMatrixPQ(BBMPoint *point, BBMMatrix *DEP);
  void getTangentMatrix(BBMPoint *point, BBMMatrix *DEP);
  void calculateElastoPlasticTangentMatrix(BBMPoint *point, BBMMatrix *DEP);
  void calculateElasticTangentMatrix(BBMPoint *point, BBMMatrix *DEP);
  void getDerivative(double meanStress, double shearStress, double suction,
                     double pZero, double *state, double *deriv);
  double getLambda(double *deriv, double stresspq[3], double strainpq[3]);
  double
  plasticEuler(BBMPoint *point, double *epStrain, double *absStress,
               int numberIterations); // returns elapsed time of computations
  double doRungeutta(double A[][8], double *B, double *BRes, double *C,
                    BBMPoint *point, double *epStrain, double *absStress,
                    int *numberIter, double methodOrder, int methodSteps,
                    bool errorEstimate);
  double doRungeuttaEqualStep(double A[][8], double *B, double *BRes, double *C,
                             BBMPoint *point, double *epStrain,
                             double *absStress, double *RelError,
                             int numberIter, double methodOrder,
                             int methodSteps, bool errorEstimate);
  double doRungeuttaExtrapol(double A[][8], double *B, double *BRes, double *C,
                            BBMPoint *point, double *epStrain,
                            double *absStress, int *numberIter,
                            double methodOrder, int methodSteps,
                            bool errorEstimate);
  // Runge Kutta schemes
  double calculatePlastic(double *purelyPlasticStrain, BBMPoint *point);
  double calculatePlasticConst(double *purelyPlasticStrain, BBMPoint *point,
                               int stepNo);
  double plasticRKErr8544(BBMPoint *point, double *epStrain, double *absStress,
                          int *numberIter); // Bogacki - Shimpine
  double plasticRKDP754(BBMPoint *point, double *epStrain, double *absStress,
                        int *numberIter); // Dormand Prince
  double plasticRKCK654(BBMPoint *point, double *epStrain, double *absStress,
                        int *numberIter); // Cash - Karp
  double plasticRKEng654(BBMPoint *point, double *epStrain, double *absStress,
                         int *numberIter); // England as given by Sloan
  double plasticRK543(BBMPoint *point, double *epStrain, double *absStress,
                      int *numberIter); // 4th order with 3rd ord estimate
  double plasticRK332(BBMPoint *point, double *epStrain, double *absStress,
                      int *numberIter); // 3rd order R-K scheme
  double plasticRKBog432(
      BBMPoint *point, double *epStrain, double *absStress,
      int *numberIter); // Bogacki - Shimpine 3rd order Runge Kutta scheme
  double plasticRKME221(BBMPoint *point, double *epStrain, double *absStress,
                        int *numberIter); // Modified Euler
  double plasticRKNoExTry(BBMPoint *point, double *epStrain, double *absStress,
                          int *numberIter); // using not in an extrapolation way
  // Extrapolation Schemes
  double plasticExtrapol(BBMPoint *point, double *epStrain, double *absStress,
                         int *numberIter);
  double doRKExtrapolation(double A[][8], double *B, double *BRes, double *C,
                         BBMPoint *point, double *epStrain, double *absStress,
                         BBMPoint *OldPoint, double *RelError, int *numberIter,
                         double methodOrder, int methodSteps,
                         bool errorEstimate);
  double plasticMidpoint(BBMPoint *point, double *epStrain, double *absStress,
                         int *numberIter);
  double plasticMidpointGallipoli(BBMPoint *point, double *epStrain,
                                  double *absStress, int *numberIter);
  double checkNorm(double *DSigma, double dpZeroStar, BBMPoint *initialPoint,
                   double *DError); // returns RError
  double checkNormSloan(double *DSigma, double dpZeroStar,
                        BBMPoint *initialPoint,
                        double *DError); // returns RError
  void correctDrift(BBMPoint *point);
  void correctDriftBeg(BBMPoint *endPoint, BBMPoint *pointOld);

  // Used Procedures before, not updated anymore, though, mostly, working.
  // void driftCorrect (BBMPoint point, double* epStrain, BBMMatrix* dSigma,
  // double* Lambda, double* dpZeroStar, double* fValue);
  // void correctDriftBeg (BBMPoint point, double* epStrain, BBMMatrix* dSigma,
  // double* Lambda, double* dpZeroStar, double* fValue);
  // double Plastic (BBMPoint* point, double* epStrain, double* absStress, int*
  // numberIter);		//returns elapsed time of computations
  // double PlasticNewSlow (BBMPoint* point, double* epStrain, double*
  // absStress, int* numberIterations);	//returns elapsed time of computations
  // double plasticRKErr6 (BBMPoint* point, double* epStrain, double* absStress,
  // int* numberIterations);	//returns elapsed time of computations
  // double plasticRKErr75 (BBMPoint* point, double* epStrain, double*
  // absStress, int* numberIter);

  // double plasticRK5Err4_2 (BBMPoint* point, double* epStrain, double*
  // absStress, int* numberIter);
  // double plasticRK4Err3  (BBMPoint* point, double* epStrain, double*
  // absStress, int* numberIter);
  // double plasticRK4Err3v2  (BBMPoint* point, double* epStrain, double*
  // absStress, int* numberIter);
  // double plasticRK3Err2  (BBMPoint* point, double* epStrain, double*
  // absStress, int* numberIter);
  // double plasticRK3Err2v2  (BBMPoint* point, double* epStrain, double*
  // absStress, int* numberIter);
  // double plasticRK3Err2v3  (BBMPoint* point, double* epStrain, double*
  // absStress, int* numberIter);
  // double plasticRKSloan (BBMPoint* point, double* epStrain, double*
  // absStress, int* numberIterations);	//returns elapsed time of computations
  // double plasticMidpointC (BBMPoint* point, double* epStrain, double*
  // absStress, int* numberIter);
  // double plasticMidpointCN (BBMPoint* point, double* epStrain, double*
  // absStress, int* numberIter);
  // double plasticMidpointC4 (BBMPoint* point, double* epStrain, double*
  // absStress, int* numberIter);
  // double plasticMidpointC6 (BBMPoint* point, double* epStrain, double*
  // absStress, int* numberIter);

private:

  // Elastic parameters
  struct Elastic {
    double d_G;  // Shear modulus
    double d_K;  // Bulk modulus
    double d_E;  // Young's modulus
    double d_nu; // Poisson's ratio

    void set(double G, double K) {
      d_G = G; d_K = K; 
      d_E = 9.0 * K * G / (G + 3.0 * K);
      d_nu = (3.0 * K - 2.0 * G) / (2.0 * G + 6.0 * K);
    }
  };

  // Mohr - Coulomb parameters
  struct Yield {
    double d_cohesion;
    double d_phi;
    double d_sin_phi;
    double d_cos_phi;

    // Rounded Mohr-Coulomb parameters
    double d_alpha;  // For the M parameter
    double d_alpha4; // alpha^4;

    void set(double c, double phi) {
      d_cohesion = c; 
      d_phi = phi * M_PI / 180.0; 
      d_sin_phi = std::sin(phi); d_cos_phi = std::cos(phi);
      d_alpha = (3.0 - d_sin_phi) / (3.0 + d_sin_phi);
      d_alpha4 = alpha * alpha * alpha * alpha;
    }
  };

  // If the flow rule is not associated;
  struct PotentialFn {
    double d_psi;
    double d_sin_psi;
    double d_cos_psi;

    void set(double psi) {
      d_psi = psi * M_PI / 180.0; 
      d_sin_psi = std::sin(psi); d_cos_psi = std::cos(psi);
    }
  };

  // Integration parameters
  struct IntegrationParameters {

    // Elastic parameters, used in the Pegasus algorithm
    int    d_maxIter;
    double d_alfaCheck;
    double d_alfaChange;
    double d_alfaRatio;

    // Yield tolerance
    double d_yieldTol;      // tolerance for the Yield locus (relatively defined, used
                            // to check if the stress state is on the YL)
    double d_integrationTol;
    double d_betaFactor;    // safety factor
    double d_minMeanStress; 
    double d_suctionTol;    // used in checking suction yield locus. advisible not to modify...

    // Algorithms
    DriftCorrection d_driftCorrection;
    ToleranceMethod d_tolMethod;
    SolutionAlgorithm d_solutionAlgorithm;

    void setDefaults() {
      d_maxIter    = 200;
      d_alfaCheck  = 1;
      d_alfaChange = 0.05;
      d_alfaRatio  = 10;

      d_yieldTol       = 1e-6; 
      d_integrationTol = 0.01; 
      d_betaFactor     = 0.9;   
      d_minMeanStress  = -1.0e8; 
      d_suctionTol     = 1e-8; 

      d_driftCorrection   = DriftCorrection::CORRECTION_AT_END;
      d_tolMethod         = ToleranceMethod::SLOAN; 
      d_solutionAlgorithm = SolutionAlgorithm::MODIFIED_EULER;
    }

    void setDefaults(const Yield& yield) {
      setDefaults();
      if (yield.sin_phi > 0) {
        d_minMeanStress = yield.cohesion * yield.cos_phi / yield.sin_phi;
      } 
    }
  }

  // Model parameters
  bool        d_nonAssociated;
  Elastic     d_elastic;
  Yield       d_yield;
  PotentialFn d_potential;

  // Integration parameters
  IntegrationParameters d_integration;

  int calcPlastic(BBMPoint point, double *epStrain, BBMMatrix *dSigma,
                  double *plasticStrain, double *dpZeroStar, double fValue,
                  double *dS, double *dLambda);
  int calcPlasticPQ(BBMPoint point, double *epStrain, BBMMatrix *dSigma,
                    double *plasticStrain, double *dpZeroStar, double fValue,
                    double *dS, double *dLambda);
  int calcPlasticFaster(BBMPoint point, double *epStrain, BBMMatrix *dSigma,
                        double *plasticStrain, double *dpZeroStar,
                        double fValue, double *dS, double *dLambda);
};

} // end namespace Uintah

#endif //__MPM_CONSTITUTIVEMODEL_MODELS_SHENG_MOHRCOULOMB__
