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

#ifndef __MODELS_WATER_EOS_MODEL_H__
#define __MODELS_WATER_EOS_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/MPMEquationOfState.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

////////////////////////////////////////////////////////////////////////////
/*!
  \class WaterEOS
  \brief Murnaghan equation of state for pure water.

  The equation of state is given by
  \f[
    p = p_0 + K_0/n*(J^{-n} - 1)
  \f]
  where \n
  \f$p\f$ = pressure\n
  \f$p_0\f$ = reference pressure\n
  \f$K_0\f$ = reference bulk modulus\n
  \f$n\f$ = pressure derivative of reference bulk modulus
  \f$J = \rho_0/\rho\f$ = ratio of mass densities

  \warning For use only with Arena
*/
////////////////////////////////////////////////////////////////////////////

class WaterEOS : public MPMEquationOfState
{

private:
  double d_p0;
  double d_K0;
  double d_n;

  // Prevent copying of this class
  // copy constructor
  // WaterEOS(const WaterEOS &cm);
  WaterEOS&
  operator=(const WaterEOS& cm);

public:
  // constructors
  WaterEOS();
  WaterEOS(Uintah::ProblemSpecP& ps);
  WaterEOS(const WaterEOS* cm);

  // destructor
  ~WaterEOS() override;

  void
  outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  EOSMaterialType
  materialType() const override
  {
    return EOSMaterialType::FLUID;
  }

  /*! Get parameters */
  std::map<std::string, double>
  getParameters() const override
  {
    std::map<std::string, double> params;
    params["p0"] = d_p0;
    params["K0"] = d_K0;
    params["n"]  = d_n;
    params["Kw"] = d_bulkModulus;
    return params;
  }

  //////////
  // Calculate the pressure using a equation of state
  double
  computePressure(const Uintah::MPMMaterial* matl,
                  const ModelStateBase* state,
                  const Uintah::Matrix3& deformGrad,
                  const Uintah::Matrix3& rateOfDeformation,
                  const double& delT) override;
  double
  computePressure(const double& rho_orig, const double& rho_cur) const override;
  void
  computePressure(const double& rho_orig,
                  const double& rho_cur,
                  double& pressure,
                  double& dp_drho,
                  double& csquared) override;

  //////////
  // Calculate the derivative of the pressure
  double
  eval_dp_dJ(const Uintah::MPMMaterial* matl,
             const double& detF,
             const ModelStateBase* state) override;

  // Compute bulk modulus
  double
  computeInitialBulkModulus() const override;
  double
  computeBulkModulus(const double& pressure) const;
  double
  computeBulkModulus(const double& rho_orig,
                     const double& rho_cur) const override;
  double
  computeBulkModulus(const ModelStateBase* state) const override;

  // Compute strain energy
  double
  computeStrainEnergy(const double& rho_orig, const double& rho_cur) override;
  double
  computeStrainEnergy(const ModelStateBase* state) override;

  // Compute density given pressure
  double
  computeDensity(const double& rho_orig, const double& pressure) override;

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
  computeDpDepse_s(const ModelStateBase* state) const override
  {
    return 0.0;
  };

  ////////////////////////////////////////////////////////////////////////
  /**
   * Function: computeElasticVolumetricStrain
   *
   * Purpose:
   *   Compute the volumetric strain given a pressure (p)
   *
   * Inputs:
   *   pp  = current pressure
   *   p0 = initial pressure
   *
   * Returns:
   *   eps_e_v = current elastic volume strain
   */
  ////////////////////////////////////////////////////////////////////////
  double
  computeElasticVolumetricStrain(const double& pp, const double& p0) override;

  ////////////////////////////////////////////////////////////////////////
  /**
   * Function: computeExpElasticVolumetricStrain
   *
   * Purpose:
   *   Compute the exponential of volumetric strain given a pressure (p)
   *
   * Inputs:
   *   pp  = current pressure
   *   p0 = initial pressure
   *
   * Returns:
   *   exp(eps_e_v) = exponential of the current elastic volume strain
   */
  ////////////////////////////////////////////////////////////////////////
  double
  computeExpElasticVolumetricStrain(const double& pp,
                                    const double& p0) override;

  ////////////////////////////////////////////////////////////////////////
  /**
   * Function: computeDerivExpElasticVolumetricStrain
   *
   * Purpose:
   *   Compute the pressure drivative of the exponential of
   *   the volumetric strain at a given pressure (p)
   *
   * Inputs:
   *   pp  = current pressure
   *   p0 = initial pressure
   *
   * Outputs:
   *   exp_eps_e_v = exp(eps_e_v) = exponential of elastic volumeric strain
   *
   * Returns:
   *   deriv = d/dp[exp(eps_e_v)] = derivative of the exponential of
   *                                current elastic volume strain
   */
  ////////////////////////////////////////////////////////////////////////
  double
  computeDerivExpElasticVolumetricStrain(const double& pp,
                                         const double& p0,
                                         double& exp_eps_e_v) override;
};

} // End namespace Vaango

#endif // __MODELS_WATER_EOS_MODEL_H__
