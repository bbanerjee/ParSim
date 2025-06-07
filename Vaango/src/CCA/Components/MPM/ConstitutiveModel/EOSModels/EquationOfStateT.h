/*
 * The MIT License
 *
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#ifndef __EQUATION_OF_STATE_TEMPLATED_H__
#define __EQUATION_OF_STATE_TEMPLATED_H__

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateT.h>

#include <Core/Math/Matrix3.h>

namespace Vaango {

////////////////////////////////////////////////////////////////////////////
/*!
  \class EquationOfStateT
  \brief Base class for solid equations of state
  \author Biswajit Banerjee, \n
*/
////////////////////////////////////////////////////////////////////////////

template <typename DerivedT, typename StateT>
class EquationOfStateT
{

public:
  ~EquationOfStateT() = default;

  void
  outputProblemSpec(Uintah::ProblemSpecP& ps)
  {
    derived()->l_outputProblemSpec(ps);
  }

  std::map<std::string, double>
  getParameters() const
  {
    return derived()->l_getParameters();
  }

  void
  setBulkModulus(const double& bulk)
  {
    d_bulkModulus = bulk;
  }

  ////////////////////////////////////////////////////////////////////////
  /*! Calculate the hydrostatic component of stress (pressure)
      using an equation of state */
  ////////////////////////////////////////////////////////////////////////
  double
  computePressure(const Uintah::MPMMaterial* matl,
                  const StateT* state,
                  const Uintah::Matrix3& F, /* def grad */
                  const Uintah::Matrix3& D, /* rate of def */
                  const double& delT)
  {
    return derived()->l_computePressure(matl, state, F, D, delT);
  }

  ////////////////////////////////////////////////////////////////////////
  /*! Calculate the pressure without considering internal energy
      (option 1)*/
  ////////////////////////////////////////////////////////////////////////
  double
  computePressure(const double& rho_orig, const double& rho_cur) const
  {
    return derived()->l_computePressure(rho_orig, rho_cur);
  }

  ////////////////////////////////////////////////////////////////////////
  /*! Calculate the pressure without considering internal energy
      (option 2).  Also compute dp/drho and c^2. */
  ////////////////////////////////////////////////////////////////////////
  void
  computePressure(const double& rho_orig,
                  const double& rho_cur,
                  double& pressure,
                  double& dp_drho,
                  double& csquared)
  {
    derived()->l_computePressure(
      rho_orig, rho_cur, pressure, dp_drho, csquared);
  }

  /*! Calculate the derivative of \f$p(J)\f$ wrt \f$J\f$
      where \f$J = det(F) = rho_0/rho\f$ */
  double
  eval_dp_dJ(const Uintah::MPMMaterial* matl,
             const double& delF,
             const StateT* state)
  {
    return derived()->l_eval_dp_dJ(matl, delF, state);
  }

  ////////////////////////////////////////////////////////////////////////
  /*! Calculate the derivative of p with respect to epse_v
      where epse_v = tr(epse)
            epse = total elastic strain */
  ////////////////////////////////////////////////////////////////////////
  double
  computeDpDepse_v(const StateT* state) const
  {
    return derived()->l_computeDpDepse_v(state);
  }

  ////////////////////////////////////////////////////////////////////////
  /*! Calculate the derivative of p with respect to epse_s
      where epse_s = sqrt{2}{3} ||ee||
            ee = epse - 1/3 tr(epse) I
            epse = total elastic strain */
  ////////////////////////////////////////////////////////////////////////
  double
  computeDpDepse_s(const StateT* state) const
  {
    return derived()->l_computeDpDepse_s(state);
  }

  // Calculate rate of temperature change due to compression/expansion
  double
  computeIsentropicTemperatureRate(const double T,
                                   const double rho_0,
                                   const double rho_cur,
                                   const double Dtrace)
  {
    return derived()->l_computeIsentropicTemperatureRate(
      T, rho_0, rho_cur, Dtrace);
  }

  ////////////////////////////////////////////////////////////////////////
  /*! Calculate the tangent bulk modulus */
  ////////////////////////////////////////////////////////////////////////
  double
  computeInitialBulkModulus()
  {
    return derived()->l_computeInitialBulkModulus();
  }

  double
  computeBulkModulus(const StateT* state)
  {
    return derived()->l_computeBulkModulus(state);
  }

  double
  computeBulkModulus(const double& rho_orig, const double& rho_cur)
  {
    return derived()->l_computeBulkModulus(rho_orig, rho_cur);
  }

  ////////////////////////////////////////////////////////////////////////
  /*! Calculate the accumulated strain energy */
  ////////////////////////////////////////////////////////////////////////
  double
  computeStrainEnergy(const StateT* state)
  {
    return derived()->l_computeStrainEnergy(state);
  }

  double
  computeStrainEnergy(const double& rho_orig, const double& rho_cur)
  {
    return derived()->l_computeStrainEnergy(rho_orig, rho_cur);
  }

  ////////////////////////////////////////////////////////////////////////
  /*! Calculate the mass density given a pressure */
  ////////////////////////////////////////////////////////////////////////
  double
  computeDensity(const double& rho_orig, const double& pressure)
  {
    return derived()->l_computeDensity(rho_orig, pressure);
  }

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
  computeElasticVolumetricStrain(const double& pp, const double& p0)
  {
    return derived()->l_computeElasticVolumetricStrain(pp, p0);
  }

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
  computeExpElasticVolumetricStrain(const double& pp, const double& p0)
  {
    return derived()->l_computeExpElasticVolumetricStrain(pp, p0);
  }

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
                                         double& exp_eps_e_v)
  {
    return derived()->l_computeDerivExpElasticVolumetricStrain(
      pp, p0, exp_eps_e_v);
  }

protected:
  double d_bulkModulus;

private:
  EquationOfStateT() { d_bulkModulus = 0.0; }

  DerivedT*
  derived()
  {
    return static_cast<DerivedT*>(this);
  }

  const DerivedT*
  derived() const
  {
    return static_cast<const DerivedT*>(this);
  }

  friend DerivedT;
};
} // End namespace Vaango

#endif // __EQUATION_OF_STATE_TEMPLATED_H__
