/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#ifndef __BORJA_PRESSURE_MODEL_H__
#define __BORJA_PRESSURE_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/MPMEquationOfState.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_CamClay.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

////////////////////////////////////////////////////////////////////////////
/*!
  \class BorjaEOS

  \brief The Borja model for calculating pressure

  Reference:Borja, R.I. and Tamagnini, C.(1998) Cam-Clay plasticity Part III:
  Extension of the infinitesimal model to include finite strains,
  Computer Methods in Applied Mechanics and Engineering, 155 (1-2),
  pp. 73-95.

  The pressure is given by

  p = p0 beta exp[(epse_v - epse_v0)/kappahat]

  where

  p0 = constant
  beta = 1 + 3/2 alpha/kappahat epse_s^2
  alpha = constant
  kappahat = constant
  epse_s = sqrt(2/3) ||epse||
  epse_v = tr(epse)
  epse_v0 = constant
  epse = elastic strain tensor

*/
////////////////////////////////////////////////////////////////////////////

class BorjaEOS : public MPMEquationOfState
{

private:
  double d_p0;         // Reference pressure
  double d_alpha;      // Pressure-shear coupling constant
  double d_kappatilde; // Reference compressibility
  double d_kappahat;   // Large deformation compressibility
  double d_epse_v0;    // Volumetric strain at reference pressure

public:
  // constructors
  BorjaEOS(Uintah::ProblemSpecP& ps);
  BorjaEOS(const BorjaEOS* cm);
  BorjaEOS&
  operator=(const BorjaEOS& cm) = delete;

  // Special operator for computing internal energy
  double
  operator()(double eta) const;

  // destructor
  ~BorjaEOS() override;

  void
  outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  EOSMaterialType
  materialType() const override
  {
    return EOSMaterialType::ROCK_SOIL;
  }

  /*! Get parameters */
  std::map<std::string, double>
  getParameters() const override
  {
    std::map<std::string, double> params;
    params["p0"]         = d_p0;
    params["alpha"]      = d_alpha;
    params["kappatilde"] = d_kappatilde;
    params["kappahat"]   = d_kappahat;
    params["epse_v0"]    = d_epse_v0;
    return params;
  }

  /////////////////////////////////////////////////////////////////////////
  /*! Calculate the pressure using a equation of state */
  /////////////////////////////////////////////////////////////////////////
  double
  computePressure(const Uintah::MPMMaterial* matl,
                  const ModelStateBase* state,
                  const Uintah::Matrix3& deformGrad,
                  const Uintah::Matrix3& rateOfDeformation,
                  const double& delT) override;

  // Compute the bulk modulus
  double
  computeBulkModulus(const ModelStateBase* state) const override;

  // Compute the volumetric strain energy
  double
  computeStrainEnergy(const ModelStateBase* state) override;

  double
  eval_dp_dJ(const Uintah::MPMMaterial* matl,
             const double& delF,
             const ModelStateBase* state) override;

  ////////////////////////////////////////////////////////////////////////
  /*! Calculate the derivative of p with respect to epse_v
      where epse_v = tr(epse)
            epse = total elastic strain */
  ////////////////////////////////////////////////////////////////////////
  double
  computeDpDepse_v(const ModelStateBase* state) const override;

  ////////////////////////////////////////////////////////////////////////
  /*! Calculate the derivative of p with respect to epse_s
      where epse_s = sqrt{2}{3} ||ee||
            ee = epse - 1/3 tr(epse) I
            epse = total elastic strain */
  ////////////////////////////////////////////////////////////////////////
  double
  computeDpDepse_s(const ModelStateBase* state) const override;

  // Calculate rate of temperature change due to compression/expansion
  double
  computeIsentropicTemperatureRate(const double T,
                                   const double rho_0,
                                   const double rho_cur,
                                   const double Dtrace) override;

  // Compute pressure (option 1)
  double
  computePressure(const double& rho_orig, const double& rho_cur) const override;

  // Compute pressure (option 2)
  void
  computePressure(const double& rho_orig,
                  const double& rho_cur,
                  double& pressure,
                  double& dp_drho,
                  double& csquared) override;

  // Compute bulk modulus
  void
  setInitialBulkModulus();
  double
  computeInitialBulkModulus() const override;
  double
  computeBulkModulus(const double& rho_orig,
                     const double& rho_cur) const override;

  // Compute strain energy
  double
  computeStrainEnergy(const double& rho_orig, const double& rho_cur) override;

  // Compute density given pressure
  double
  computeDensity(const double& rho_orig, const double& pressure) override;

  double
  computeElasticVolumetricStrain(const double& pp, const double& p0) override;
  double
  computeExpElasticVolumetricStrain(const double& pp,
                                    const double& p0) override;
  double
  computeDerivExpElasticVolumetricStrain(const double& pp,
                                         const double& p0,
                                         double& exp_eps_e_v) override;

private:
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

} // End namespace Vaango

#endif // __BORJA_PRESSURE_MODEL_H__
