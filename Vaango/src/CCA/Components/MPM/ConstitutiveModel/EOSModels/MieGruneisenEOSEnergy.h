/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef __MIE_GRUNEISEN_EOS_ENERGY_MODEL_H__
#define __MIE_GRUNEISEN_EOS_ENERGY_MODEL_H__

#include "MPMEquationOfState.h"
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

////////////////////////////////////////////////////////////////////////////
/*!
  \class MieGruneisenEOSEnergy

  \brief A Mie-Gruneisen equation of state model with three-parameter Us-up fit

   Modified by Jim Guilkey to be energy based.

*/
////////////////////////////////////////////////////////////////////////////

class MieGruneisenEOSEnergy : public MPMEquationOfState
{
public:
  struct CMData
  {
    double C_0;
    double Gamma_0;
    double S_1;
    double S_2;
    double S_3;
    double rho_0;
  };

  // constructors
  MieGruneisenEOSEnergy(Uintah::ProblemSpecP& ps);
  MieGruneisenEOSEnergy(const MieGruneisenEOSEnergy* cm);
  MieGruneisenEOSEnergy&
  operator=(const MieGruneisenEOSEnergy& cm) = delete;

  // Special operator for computing internal energy
  double
  operator()(double eta) const;

  // destructor
  ~MieGruneisenEOSEnergy() override;

  void
  outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  EOSMaterialType
  materialType() const override
  {
    return EOSMaterialType::ALL;
  }

  std::map<std::string, double>
  getParameters() const override;

  /////////////////////////////////////////////////////////////////////////
  /*! Calculate the pressure using a equation of state */
  /////////////////////////////////////////////////////////////////////////
  double
  computePressure(const Uintah::MPMMaterial* matl,
                  const ModelStateBase* state,
                  const Uintah::Matrix3& deformGrad,
                  const Uintah::Matrix3& rateOfDeformation,
                  const double& delT) override;

  // Calculate rate of temperature change due to compression/expansion
  double
  computeIsentropicTemperatureRate(const double T,
                                   const double rho_0,
                                   const double rho_cur,
                                   const double Dtrace) override;

  double
  eval_dp_dJ(const Uintah::MPMMaterial* matl,
             const double& delF,
             const ModelStateBase* state) override;

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
  double
  computeInitialBulkModulus() const override;
  double
  computeBulkModulus(const double& rho_orig,
                     const double& rho_cur) const override;
  double
  computeBulkModulus(const ModelStateBase* state) const override;

  // Compute strain energy
  double
  computeStrainEnergy(const ModelStateBase* state) override;
  double
  computeStrainEnergy(const double& rho_orig, const double& rho_cur) override;

  // Compute density given pressure
  double
  computeDensity(const double& rho_orig, const double& pressure) override;

  double
  computeDpDepse_v(const Vaango::ModelStateBase*) const override;
  double
  computeDpDepse_s(const Vaango::ModelStateBase*) const override;
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
  double
  eval_dp_dJ(double rho_0, double rho) const;

  typedef double (MieGruneisenEOSEnergy::*pFuncPtr)(const double&,
                                                    const double&) const;
  typedef double (MieGruneisenEOSEnergy::*dpdJFuncPtr)(const double&,
                                                       const double&) const;

  // Find root of p(eta) - p0 = 0 using Ridder's method
  double
  findEtaRidder(pFuncPtr pFunc,
                const double& rho_orig,
                const double& p0,
                double& etamin,
                double& etamax,
                const double& tolerance,
                const int& maxIter) const;

  // Find root of p(eta) - p0 = 0 using Newton's method
  double
  findEtaNewton(pFuncPtr pFunc,
                const dpdJFuncPtr dpdJFunc,
                const double& rho_orig,
                const double& p0,
                const double& J0,
                const double& tolerance,
                const int& maxIter) const;

  // Compute p for compressive volumetric deformations
  double
  pCompression(const double& rho_orig, const double& eta) const;

  // Compute dp/dJ for compressive volumetric deformations
  double
  dpdJCompression(const double& rho_orig, const double& eta) const;

  // Compute p for tensile volumetric deformations
  double
  pTension(const double& rho_orig, const double& eta) const;

  // Compute dp/dJ for tensile volumetric deformations
  double
  dpdJTension(const double& rho_orig, const double& eta) const;

private:
  CMData d_const;
  double d_J_min;
};

} // End namespace Vaango

#endif // __MIE_GRUNEISEN_EOS_ENERGY_MODEL_H__
