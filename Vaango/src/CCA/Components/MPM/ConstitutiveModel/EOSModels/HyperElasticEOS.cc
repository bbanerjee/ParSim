/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2022 Parresia Research Limited, New Zealand
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

#include "HyperElasticEOS.h"
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <cmath>

using namespace Uintah;
using namespace Vaango;
using std::ostringstream;

HyperElasticEOS::HyperElasticEOS()
{
  d_bulkModulus = -1.0;
}

HyperElasticEOS::HyperElasticEOS(ProblemSpecP& ps)
{
  ps->require("bulk_modulus", d_bulkModulus);
}

HyperElasticEOS::HyperElasticEOS(const HyperElasticEOS* cm)
{
  d_bulkModulus = cm->d_bulkModulus;
}

HyperElasticEOS::~HyperElasticEOS() = default;

void
HyperElasticEOS::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP eos_ps = ps->appendChild("equation_of_state");
  eos_ps->setAttribute("type", "default_hyper");
  eos_ps->appendElement("bulk_modulus", d_bulkModulus);
}

std::map<std::string, double>
HyperElasticEOS::getParameters() const
{
  std::map<std::string, double> params;
  params["bulk_modulus"] = d_bulkModulus;
  return params;
}

//////////
// Calculate the pressure using the elastic constitutive equation
double
HyperElasticEOS::computePressure(const MPMMaterial* matl,
                                 const ModelStateBase* state,
                                 const Matrix3&,
                                 const Matrix3& rateOfDeformation,
                                 const double& delT)
{
  double rho_0 = matl->getInitialDensity();
  double rho   = state->density;
  double J     = rho_0 / rho;
  double kappa = state->bulkModulus;

  double p = 0.5 * kappa * (J - 1.0 / J);
  return p;
}

double
HyperElasticEOS::eval_dp_dJ(const MPMMaterial* matl,
                            const double& detF,
                            const ModelStateBase* state)
{
  double J     = detF;
  double kappa = state->bulkModulus;

  double dpdJ = 0.5 * kappa * (1.0 + 1.0 / (J * J));
  return dpdJ;
}

// Compute pressure (option 1)
double
HyperElasticEOS::computePressure(const double& rho_orig,
                                 const double& rho_cur) const
{
  double J = rho_orig / rho_cur;
  double p = 0.5 * d_bulkModulus * (J - 1.0 / J);
  return p;
}

// Compute pressure (option 2)
void
HyperElasticEOS::computePressure(const double& rho_orig,
                                 const double& rho_cur,
                                 double& pressure,
                                 double& dp_drho,
                                 double& csquared)
{
  double J     = rho_orig / rho_cur;
  pressure     = 0.5 * d_bulkModulus * (J - 1.0 / J);
  double dp_dJ = 0.5 * d_bulkModulus * (1.0 + 1.0 / (J * J));
  dp_drho      = -0.5 * d_bulkModulus * (1.0 + J * J) / rho_orig;
  csquared     = dp_dJ / rho_cur;
}

// Compute bulk modulus
double
HyperElasticEOS::computeInitialBulkModulus() const
{
  return computeBulkModulus(1.0, 1.0);
}

double
HyperElasticEOS::computeBulkModulus(const double& rho_orig,
                                    const double& rho_cur) const
{
  double J    = rho_orig / rho_cur;
  double bulk = 0.5 * d_bulkModulus * (1.0 + 1.0 / (J * J));
  return bulk;
}

double
HyperElasticEOS::computeBulkModulus(const ModelStateBase* state) const
{
  return computeBulkModulus(state->initialDensity, state->density);
}

double
HyperElasticEOS::computeStrainEnergy(const ModelStateBase* state)
{
  return computeStrainEnergy(state->initialDensity, state->density);
}

// Compute strain energy
double
HyperElasticEOS::computeStrainEnergy(const double& rho_orig,
                                     const double& rho_cur)
{
  double J = rho_orig / rho_cur;
  double U = 0.5 * d_bulkModulus * (0.5 * (J * J - 1.0) - log(J));
  return U;
}

// Compute density given pressure (tension +ve)
double
HyperElasticEOS::computeDensity(const double& rho_orig, const double& pressure)
{
  double numer1    = d_bulkModulus * d_bulkModulus + pressure * pressure;
  double sqrtNumer = sqrt(numer1);
  double rho       = rho_orig / d_bulkModulus * (-pressure + sqrtNumer);
  if (rho < 0) {
    ostringstream desc;
    desc << "Value of pressure (" << pressure
         << ") is beyond the range of validity of model"
         << "\n"
         << "  density = " << rho << "\n";
    throw InvalidValue(desc.str(), __FILE__, __LINE__);
  }
  return rho;
}

double
HyperElasticEOS::computeDpDepse_v(const Vaango::ModelStateBase*) const
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute dp/deps_e_v of hyperelastic material"
         " unless the elastic part of J is provided."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}

double
HyperElasticEOS::computeDpDepse_s(const Vaango::ModelStateBase*) const
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute dp/deps_e_s of hyperelastic material"
         " unless the elastic part of J is provided."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}

double
HyperElasticEOS::computeElasticVolumetricStrain(const double& pp,
                                                const double& p0)
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute volume strain of hyperelastic material."
         " It should be provided as an input."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}

double
HyperElasticEOS::computeExpElasticVolumetricStrain(const double& pp,
                                                   const double& p0)
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute exp(volume strain) of hyperelastic material."
         " It should be provided as an input."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}

double
HyperElasticEOS::computeDerivExpElasticVolumetricStrain(const double& pp,
                                                        const double& p0,
                                                        double& exp_eps_e_v)
{
  std::ostringstream err;
  err << "**ERROR** Cannot compute derivative of exp(volume strain) of "
         " hyperelastic material. It should be provided as an input."
         " Please change the equation_of_state if you need this "
         " functionality.\n";
  throw InternalError(err.str(), __FILE__, __LINE__);

  return -1;
}
