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

#ifndef __MPM_CONSTITUTIVEMODEL_MODELS_SHENG_MOHRCOULOMB__
#define __MPM_CONSTITUTIVEMODEL_MODELS_SHENG_MOHRCOULOMB__

#include <CCA/Components/MPM/ConstitutiveModel/RockSoilModels/MohrCoulombBase.h>

#include <cmath>
#include <vector>

/* this Mohr-Coulomb like model uses a rounded Mohr-Coulomb surface, see
   eq 13 in Sheng D, Sloan SW & Yu HS Computational Mechanics 26:185-196 (2000)
   Springer
*/

namespace Uintah {

class MohrCoulombSheng : public MohrCoulombBase
{

public:
  MohrCoulombSheng();
  MohrCoulombSheng(double G, double K, double cohesion, double phi, double psi, double pMin = -1);
  MohrCoulombSheng(const MohrCoulombSheng* cm);
  MohrCoulombSheng(const MohrCoulombSheng&) = default;
  MohrCoulombSheng& operator=(const MohrCoulombSheng&) = default;
  virtual ~MohrCoulombSheng() = default;

  MohrCoulombState integrate(const Vector7& strainIncrement,
                             const MohrCoulombState& initialState) override;

  MohrCoulombState integrate(const Vector7& strainIncrement,
                             const MohrCoulombState& initialState, 
                             RegionType& region) override;

protected:

  bool checkYieldNormalized(const MohrCoulombState& state) const override;

  double computeYieldNormalized(const Vector6& stress) const override;

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

private:

  template <int Order, int Steps>
  std::tuple<double, int> doRungeKutta(
    const Eigen::Matrix<double, Steps, Steps>& AA,
    const Eigen::Matrix<double, Steps, 1>& BB,
    const Eigen::Matrix<double, Steps, 1>& BRes,
    const Eigen::Matrix<double, Steps, 1>& CC, MohrCoulombState& state,
    const Vector7& epStrain, bool errorEstimate) const;

  template <int Order, int Steps>
  std::tuple<double, int> doRungeKuttaErr(
    const Eigen::Matrix<double, Steps, Steps>& AA,
    const Eigen::Matrix<double, Steps, 1>& BB,
    const Eigen::Matrix<double, Steps, 1>& BRes,
    const Eigen::Matrix<double, Steps, 1>& CC,
    const Eigen::Matrix<double, Steps - 1, 1>& ErrCoef, MohrCoulombState& state,
    const Vector7& epStrain, bool errorEstimate) const;

  double plasticMidpoint(MohrCoulombState& state, const Vector7& epStrain,
                         Vector7& absStress, int numIter) const override;

};

} // end namespace Uintah

#endif //__MPM_CONSTITUTIVEMODEL_MODELS_SHENG_MOHRCOULOMB__
