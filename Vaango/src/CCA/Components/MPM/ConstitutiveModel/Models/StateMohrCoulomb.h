/*
 * The MIT License
 *
 * copyright (c) 1997-2019 Center for the Simulation of Accidental Fires and
 * Explosions (CSAFE), and  Scientific Computing and Imaging Institute (SCI),
 * University of Utah.
 * copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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

#ifndef __MPM_CONSTITUTIVEMODEL_MODELS_STATE_MOHRCOULOMB__
#define __MPM_CONSTITUTIVEMODEL_MODELS_STATE_MOHRCOULOMB__

#include <Eigen/Dense>

namespace Uintah {

using Vector3 = Eigen::Matrix<double, 3, 1, Eigen::DontAlign>;
using Vector6 = Eigen::Matrix<double, 6, 1, Eigen::DontAlign>;
using Vector7 = Eigen::Matrix<double, 7, 1, Eigen::DontAlign>;

using Matrix6 = Eigen::Matrix<double, 6, 6, Eigen::DontAlign>;
using Matrix7 = Eigen::Matrix<double, 7, 7, Eigen::DontAlign>;

class StateMohrCoulomb
{

public:

  Vector6 stress;         // stress state
  Vector7 strain;         // strain state + suction
  Vector6 plasticStrain;  // plastic strain, no suction included

  Vector3 state;          // state [0] - p0* - preconsolidation stress at suction==0
                          // state [1] - s0  - max suction experienced
                          // state [2] -       specific volume  

  Vector6 microStress;
  Vector7 microStrain;
  Vector6 microPlasticStrain;
  Vector3 microState;

  StateMohrCoulomb();
  StateMohrCoulomb(const StateMohrCoulomb& state) = default;
  StateMohrCoulomb& operator=(const StateMohrCoulomb& state) = default;
  ~StateMohrCoulomb() = default;

  void read ();
  void write ();
  void update (double* PlasticStrainInc, double* Strain_Incr, double *Stress_Incr, double dPZeroStar);

  double getMeanStress ();
  double getShearStress ();
  double getPStar ();
  double getSpecVol ();
  double getSuction ();
  double getYieldSuction ();
  double getTheta ();
  double getThetaDeg ();
  double getThetaDeg_0 ();
  void updatePStar (double value); //sets P0Star parameter
  void setSuction (double value); //sets suction
  void setSpecificVolume (double value); //sets specific volume
  void setPStar (double value);
  void setYieldSuction (double value);


  //Invariants functions
  double getFirstInvariant ();
  double getSecondInvariant ();
  double getThirdInvariant ();
  double getFirstDevInvariant();
  double getSecondDevInvariant();
  double getThirdDevInvariant();

  //Eigenvalues
  void getEigen(double Eigen[3]);
  void getEigen(double Eigen[3],BBMMatrix* EigenVectors);
  void getEigen(BBMMatrix* EigenValues, BBMMatrix* EigenVectors);
  void setStressEigen (double Eigen[3]);

  bool checkIfFinite ();
};

} // end namespace Uintah

#endif //__MPM_CONSTITUTIVEMODEL_MODELS_STATE_MOHRCOULOMB__
