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

#include "DefaultHypoElasticEOS.h"
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <cmath>

using namespace Uintah;
using Vaango::ModelStateBase;

DefaultHypoElasticEOS::DefaultHypoElasticEOS()
{
  d_bulkModulus = -1.0;
}

DefaultHypoElasticEOS::DefaultHypoElasticEOS(ProblemSpecP& ps)
{
  ps->require("bulk_modulus", d_bulkModulus);
}

DefaultHypoElasticEOS::DefaultHypoElasticEOS(const DefaultHypoElasticEOS* cm)
{
  d_bulkModulus = cm->d_bulkModulus;
}

DefaultHypoElasticEOS::~DefaultHypoElasticEOS() = default;

void
DefaultHypoElasticEOS::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP eos_ps = ps->appendChild("equation_of_state");
  eos_ps->setAttribute("type", "default_hypo");
  eos_ps->appendElement("bulk_modulus", d_bulkModulus);
}

std::map<std::string, double> 
DefaultHypoElasticEOS::getParameters() const 
{
  std::map<std::string, double> params;
  params["bulk_modulus"] = d_bulkModulus;
  return params;
}

//////////
// Calculate the pressure using the elastic constitutive equation
double
DefaultHypoElasticEOS::computePressure(const MPMMaterial*,
                                       const ModelStateBase* state,
                                       const Matrix3&,
                                       const Matrix3& rateOfDeformation,
                                       const double& delT)
{
  // Get the state data
  double kappa = state->bulkModulus;
  double p_n = state->pressure;

  // Calculate pressure increment
  double delp = rateOfDeformation.Trace() * (kappa * delT);

  // Calculate pressure
  double p = p_n + delp;
  return p;
}

/* **WARNING** Not hypoelastic.  Hypoelastic is history-dependent. */
double
DefaultHypoElasticEOS::eval_dp_dJ(const MPMMaterial* matl, const double& detF,
                                  const ModelStateBase* state)
{
  return (state->bulkModulus / detF);
}

/* **WARNING** Not hypoelastic.  Hypoelastic is history-dependent. */
// Compute pressure (option 1)
//  (assume linear relation holds and bulk modulus does not change
//   with deformation)
double
DefaultHypoElasticEOS::computePressure(const double& rho_orig,
                                       const double& rho_cur)
{
  if (d_bulkModulus < 0.0) {
    throw ParameterNotFound(
      "Please initialize bulk modulus in EOS before computing pressure",
      __FILE__, __LINE__);
  }

  double J = rho_orig / rho_cur;
  double p = d_bulkModulus * (1.0 - 1.0 / J);
  return p;
}

/* **WARNING** Not hypoelastic.  Hypoelastic is history-dependent. */
// Compute pressure (option 2)  (assume small strain relation holds)
//  (assume linear relation holds and bulk modulus does not change
//   with deformation)
void
DefaultHypoElasticEOS::computePressure(const double& rho_orig,
                                       const double& rho_cur, double& pressure,
                                       double& dp_drho, double& csquared)
{
  if (d_bulkModulus < 0.0) {
    throw ParameterNotFound(
      "Please initialize bulk modulus in EOS before computing pressure",
      __FILE__, __LINE__);
  }

  double J = rho_orig / rho_cur;
  pressure = d_bulkModulus * (1.0 - 1.0 / J);
  dp_drho = -d_bulkModulus / rho_orig;
  csquared = d_bulkModulus / rho_cur;
}

// Compute bulk modulus
double 
DefaultHypoElasticEOS::computeInitialBulkModulus()
{
  return d_bulkModulus;
}

double
DefaultHypoElasticEOS::computeBulkModulus(const double& rho_orig,
                                          const double& rho_cur)
{
  return d_bulkModulus;
}

double 
DefaultHypoElasticEOS::computeBulkModulus(const ModelStateBase* state)
{
  return d_bulkModulus;
}

// Compute strain energy
double 
DefaultHypoElasticEOS::computeStrainEnergy(const ModelStateBase* state)
{
  return state->energy;
}

/* **WARNING** Hypoelastic is history-dependent. */
double
DefaultHypoElasticEOS::computeStrainEnergy(const double& rho_orig,
                                           const double& rho_cur)
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute strain energy of hypoelastic material"
         " using only the initial and current densities.  Please change"
         " the equation_of_state if you need this functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}

/* **WARNING** Hypoelastic is history-dependent. */
// Compute density given pressure
double
DefaultHypoElasticEOS::computeDensity(const double& rho_orig,
                                      const double& pressure)
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute density of hypoelastic material"
         " using only the initial density and current pressure."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}

double 
DefaultHypoElasticEOS::computeDpDepse_v(const Vaango::ModelStateBase*) const
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute dp/deps_v of hypoelastic material"
         " unless the incremental pressure and volume strains are provided." 
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}
double 
DefaultHypoElasticEOS::computeDpDepse_s(const Vaango::ModelStateBase*) const
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute dp/deps_s of hypoelastic material"
         " unless the incremental pressure and volume strains are provided."  
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}
double 
DefaultHypoElasticEOS::computeElasticVolumetricStrain(const double& pp,
                                      const double& p0)
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute volume strain of hypoelastic material."
         " It should be provided as an input."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}

double 
DefaultHypoElasticEOS::computeExpElasticVolumetricStrain(const double& pp,
                                         const double& p0)
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute exp(volume strain) of hypoelastic material."
         " It should be provided as an input."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}
double 
DefaultHypoElasticEOS::computeDerivExpElasticVolumetricStrain(const double& pp,
                                              const double& p0,
                                              double& exp_eps_e_v)
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute derivative of exp(volume strain) of "
         " hypoelastic material. It should be provided as an input."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}
