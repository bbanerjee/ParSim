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

/** 
 * This ia a classic Mohr-Coulomb model
 * integration works based on the principal stress space
 * note that the strain must be rotated into the same coordinate system as
 * stress
*/

class ClassicMohrCoulomb : public ShengMohrCoulomb
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

protected:

  bool checkYieldNormalized(const StateMohrCoulomb& state) const override;

  double computeYieldNormalized(const Vector6& stress) const override;

private:


  Matrix67 calculateElastoPlasticTangentMatrix(
    const StateMohrCoulomb& state) const;

  std::tuple<Vector6, Vector6> computeDfDsigma(const Vector6& stress) const override;

  std::tuple<double, int> plasticRKME221(StateMohrCoulomb& state,
                                         const Vector7& epStrain) const override;

  std::tuple<double, int> plasticRK332(StateMohrCoulomb& state,
                                       const Vector7& epStrain) const override;

  std::tuple<double, int> plasticRKBog432(StateMohrCoulomb& point,
                                          const Vector7& epStrain) const override;

  std::tuple<double, int> plasticRK543(StateMohrCoulomb& state,
                                       const Vector7& epStrain) const override;

  std::tuple<double, int> plasticRKEng654(StateMohrCoulomb& state,
                                          const Vector7& epStrain) const override;

  std::tuple<double, int> plasticRKCK654(StateMohrCoulomb& state,
                                         const Vector7& epStrain) const override;

  std::tuple<double, int> plasticRKDP754(StateMohrCoulomb& state,
                                         const Vector7& epStrain) const override;

  std::tuple<double, int> plasticRKErr8544(StateMohrCoulomb& state,
                                           const Vector7& epStrain) const override;

  std::tuple<double, int> plasticExtrapol(StateMohrCoulomb& state,
                                          const Vector7& epStrain) const override;

  template <int Order, int Steps>
  std::tuple<double, int> doRungeKuttaEig(const Eigen::Matrix<double, Steps, Steps>& AA,
                                          const Eigen::Matrix<double, Steps, 1>& BB,
                                          const Eigen::Matrix<double, Steps, 1>& BRes,
                                          const Eigen::Matrix<double, Steps, 1>& CC,
                                          StateMohrCoulomb& state, const Vector7& epStrain,
                                          bool errorEstimate) const;

  template <int Order, int Steps>
  std::tuple<double, int> doRungeKuttaEigErr(
    const Eigen::Matrix<double, Steps, Steps>& AA,
    const Eigen::Matrix<double, Steps, 1>& BB,
    const Eigen::Matrix<double, Steps, 1>& BRes,
    const Eigen::Matrix<double, Steps, 1>& CC,
    const Eigen::Matrix<double, Steps - 1, 1>& ErrCoef, StateMohrCoulomb& state,
    const Vector7& epStrain, bool errorEstimate) const;

  double plasticMidpoint(StateMohrCoulomb& state, const Vector7& epStrain,
                         Vector7& absStress, int numIter) const override;

  std::tuple<Vector3, Matrix33> getEigen(const Vector6& stress) const;
  Vector6 rotateToOrigin(const Vector6& vec, const Matrix33& eigenVecs);
  Vector6 rotateToEigen(const Vector6& vec, const Matrix33& eigenVecs);
  Matrix33 toMatrix33(const Vector6& vec);
  Vector6 toVector6(const Matrix33& mat);
};

} // end namespace Uintah

#endif //__MPM_CONSTITUTIVEMODEL_MODELS_SHENG_MOHRCOULOMB__
