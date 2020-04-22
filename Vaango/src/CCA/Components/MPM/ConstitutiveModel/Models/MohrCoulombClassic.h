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
#ifndef __MPM_CONSTITUTIVEMODEL_MODELS_CLASSIC_MOHRCOULOMB__
#define __MPM_CONSTITUTIVEMODEL_MODELS_CLASSIC_MOHRCOULOMB__

#include <CCA/Components/MPM/ConstitutiveModel/Models/MohrCoulombBase.h>

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

namespace Uintah {

/** 
 * This ia a classic Mohr-Coulomb model
 * integration works based on the principal stress space
 * note that the strain must be rotated into the same coordinate system as
 * stress
*/

class MohrCoulombClassic : public MohrCoulombBase
{

public:
  MohrCoulombClassic();
  MohrCoulombClassic(double G, double K, double cohesion, double phi, double psi);
  MohrCoulombClassic(const MohrCoulombClassic* cm);

  MohrCoulombClassic(const MohrCoulombClassic&) = delete;
  MohrCoulombClassic& operator=(const MohrCoulombClassic&) = delete;

  ~MohrCoulombClassic() = default;

  MohrCoulombState integrate(const Vector7& strainIncrement,
                             const MohrCoulombState& initialState) override;

  MohrCoulombState integrate(const Vector7& strainIncrement,
                             const MohrCoulombState& initialState, 
                             RegionType& region) override;
protected:

  bool checkYieldNormalized(const MohrCoulombState& state) const override;

  double computeYieldNormalized(const Vector6& stress) const override;

private:


  Matrix67 calculateElastoPlasticTangentMatrix(
    const MohrCoulombState& state) const;

  std::tuple<Vector6, Vector6> computeDfDsigma(const Vector6& stress) const override;

  std::tuple<double, int> plasticRKME221(MohrCoulombState& state,
                                         const Vector7& epStrain) const override;

  std::tuple<double, int> plasticRK332(MohrCoulombState& state,
                                       const Vector7& epStrain) const override;

  std::tuple<double, int> plasticRKBog432(MohrCoulombState& point,
                                          const Vector7& epStrain) const override;

  std::tuple<double, int> plasticRK543(MohrCoulombState& state,
                                       const Vector7& epStrain) const override;

  std::tuple<double, int> plasticRKEng654(MohrCoulombState& state,
                                          const Vector7& epStrain) const override;

  std::tuple<double, int> plasticRKCK654(MohrCoulombState& state,
                                         const Vector7& epStrain) const override;

  std::tuple<double, int> plasticRKDP754(MohrCoulombState& state,
                                         const Vector7& epStrain) const override;

  std::tuple<double, int> plasticRKErr8544(MohrCoulombState& state,
                                           const Vector7& epStrain) const override;

  std::tuple<double, int> plasticExtrapol(MohrCoulombState& state,
                                          const Vector7& epStrain) const override;

  template <int Order, int Steps>
  std::tuple<double, int> doRungeKuttaEig(const Eigen::Matrix<double, Steps, Steps>& AA,
                                          const Eigen::Matrix<double, Steps, 1>& BB,
                                          const Eigen::Matrix<double, Steps, 1>& BRes,
                                          const Eigen::Matrix<double, Steps, 1>& CC,
                                          MohrCoulombState& state, const Vector7& epStrain,
                                          bool errorEstimate) const;

  template <int Order, int Steps>
  std::tuple<double, int> doRungeKuttaEigErr(
    const Eigen::Matrix<double, Steps, Steps>& AA,
    const Eigen::Matrix<double, Steps, 1>& BB,
    const Eigen::Matrix<double, Steps, 1>& BRes,
    const Eigen::Matrix<double, Steps, 1>& CC,
    const Eigen::Matrix<double, Steps - 1, 1>& ErrCoef, MohrCoulombState& state,
    const Vector7& epStrain, bool errorEstimate) const;

  double plasticMidpoint(MohrCoulombState& state, const Vector7& epStrain,
                         Vector7& absStress, int numIter) const override;

  MohrCoulombState doReturnImplicit(const MohrCoulombState& state, 
                                    RegionType& region) const;

  std::tuple<Vector3, Matrix33> getEigen(const Vector6& stress) const;
  Vector6 rotateToOrigin(const Vector6& vec, const Matrix33& eigenVecs) const;
  Vector6 rotateToEigen(const Vector6& vec, const Matrix33& eigenVecs) const;
  Matrix33 toMatrix33(const Vector6& vec) const;
  Vector6 toVector6(const Matrix33& mat) const;
};

} // end namespace Uintah

#endif //__MPM_CONSTITUTIVEMODEL_MODELS_CLASSIC_MOHRCOULOMB__
