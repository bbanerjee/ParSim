/*
 * The MIT License
 *
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef __BORJA_PRESSURE_TEMPLATED_MODEL_H__
#define __BORJA_PRESSURE_TEMPLATED_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/EquationOfStateT.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Borja.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

////////////////////////////////////////////////////////////////////////////
/*!
  \class EOS_BorjaT

  \brief The Borja model for calculating pressure

  Reference:Borja, R.I. and Tamagnini, C.(1998) Cam-Clay plasticity Part III:
  Extension of the infinitesimal model to include finite strains,
  Computer Methods in Applied Mechanics and Engineering, 155 (1-2),
  pp. 73-95.

  The pressure is given by

  p = p0 beta exp[(epse_v - epse_v0)/kappatilde]

  where

  p0 = constant
  beta = 1 + 3/2 alpha/kappatilde epse_s^2
  alpha = constant
  kappatilde = constant
  epse_s = sqrt(2/3) ||epse||
  epse_v = tr(epse)
  epse_v0 = constant
  epse = elastic strain tensor

*/
////////////////////////////////////////////////////////////////////////////

class EOS_BorjaT : public EquationOfStateT<EOS_BorjaT, ModelState_Borja>
{

public:
  // constructors
  EOS_BorjaT(Uintah::ProblemSpecP& ps);
  EOS_BorjaT(const EOS_BorjaT* cm);
  EOS_BorjaT&
  operator=(const EOS_BorjaT& cm) = delete;

  // Special operator for computing internal energy
  double
  operator()(double eta) const;

  // destructor
  ~EOS_BorjaT();

  void
  outputProblemSpec(Uintah::ProblemSpecP& ps);

  /*! Get parameters */
  std::map<std::string, double>
  getParameters() const
  {
    std::map<std::string, double> params;
    params["p0"]         = d_p0;
    params["alpha"]      = d_alpha;
    params["kappatilde"] = d_kappatilde;
    params["epse_v0"]    = d_epse_v0;
    return params;
  }

  /////////////////////////////////////////////////////////////////////////
  /*! Calculate the pressure using a equation of state */
  /////////////////////////////////////////////////////////////////////////
  double
  computePressure(const Uintah::MPMMaterial* matl,
                  const ModelState_Borja* state,
                  const Uintah::Matrix3& deformGrad,
                  const Uintah::Matrix3& rateOfDeformation,
                  const double& delT);

  double
  eval_dp_dJ(const Uintah::MPMMaterial* matl,
             const double& delF,
             const ModelState_Borja* state);

  ////////////////////////////////////////////////////////////////////////
  /*! Calculate the derivative of p with respect to epse_v
      where epse_v = tr(epse)
            epse = total elastic strain */
  ////////////////////////////////////////////////////////////////////////
  double
  computeDpDepse_v(const ModelState_Borja* state) const;

  ////////////////////////////////////////////////////////////////////////
  /*! Calculate the derivative of p with respect to epse_s
      where epse_s = sqrt{2}{3} ||ee||
            ee = epse - 1/3 tr(epse) I
            epse = total elastic strain */
  ////////////////////////////////////////////////////////////////////////
  double
  computeDpDepse_s(const ModelState_Borja* state) const;

  // Calculate rate of temperature change due to compression/expansion
  double
  computeIsentropicTemperatureRate(const double T,
                                   const double rho_0,
                                   const double rho_cur,
                                   const double Dtrace);

  // Compute pressure (option 1)
  double
  computePressure(const double& rho_orig, const double& rho_cur);

  // Compute pressure (option 2)
  void
  computePressure(const double& rho_orig,
                  const double& rho_cur,
                  double& pressure,
                  double& dp_drho,
                  double& csquared);

  // Compute bulk modulus
  void
  setInitialBulkModulus();
  double
  computeInitialBulkModulus();
  double
  computeBulkModulus(const ModelState_Borja* state);
  double
  computeBulkModulus(const double& rho_orig, const double& rho_cur);

  // Compute strain energy
  double
  computeStrainEnergy(const ModelState_Borja* state);
  double
  computeStrainEnergy(const double& rho_orig, const double& rho_cur);

  // Compute density given pressure
  double
  computeDensity(const double& rho_orig, const double& pressure);

  double
  computeElasticVolumetricStrain(const double& pp, const double& p0);
  double
  computeExpElasticVolumetricStrain(const double& pp,
                                    const double& p0);
  double
  computeDerivExpElasticVolumetricStrain(const double& pp,
                                         const double& p0,
                                         double& exp_eps_e_v);

private:
  double d_p0;         // Reference pressure
  double d_alpha;      // Pressure-shear coupling constant
  double d_kappatilde; // Reference compressibility
  double d_epse_v0;    // Volumetric strain at reference pressure

  //  Pressure computation
  double
  evalPressure(const double& epse_v, const double& epse_s) const;

  //  Pressure derivative computation
  double
  evalDpDepse_v(const double& epse_v, const double& epse_s) const;

  //  Shear derivative computation
  double
  evalDpDepse_s(const double& epse_v, const double& epse_s) const;
};

} // End namespace Uintah

#endif // __BORJA_PRESSURE_TEMPLATED_MODEL_H__
