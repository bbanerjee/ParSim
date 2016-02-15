/*
 * The MIT License
 *
 * Copyright (c) 2015-2016 Parresia Research Limited, New Zealand
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


#include <CCA/Components/MPM/ConstitutiveModel/Models/Pressure_Air.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_MasonSand.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <cmath>

using namespace Uintah;
using namespace Vaango;

Pressure_Air::Pressure_Air()
{
  d_p0 = 101325.0;  // Hardcoded (SI units).  *TODO* Get as input with ProblemSpec later.
  d_gamma = 1.4;
  d_bulkModulus = d_gamma*d_p0;  
} 

Pressure_Air::Pressure_Air(Uintah::ProblemSpecP&)
{
  d_p0 = 101325.0;  // Hardcoded (SI units).  *TODO* Get as input with ProblemSpec later.
  d_gamma = 1.4;
  d_bulkModulus = d_gamma*d_p0;  
} 
         
Pressure_Air::Pressure_Air(const Pressure_Air* cm)
{
  d_p0 = cm->d_p0;
  d_gamma = cm->d_gamma;
  d_bulkModulus = cm->d_bulkModulus;
} 
         
Pressure_Air::~Pressure_Air()
{
}

void Pressure_Air::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP eos_ps = ps->appendChild("pressure_model");
  eos_ps->setAttribute("type","air");
}
         
//////////
// Calculate the pressure using the elastic constitutive equation
double 
Pressure_Air::computePressure(const Uintah::MPMMaterial* matl,
                              const ModelStateBase* state_input,
                              const Uintah::Matrix3& ,
                              const Uintah::Matrix3& rateOfDeformation,
                              const double& delT)
{
  const ModelState_MasonSand* state = dynamic_cast<const ModelState_MasonSand*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_MasonSand.";
    throw SCIRun::InternalError(out.str(), __FILE__, __LINE__);
  }

  double rho_0 = matl->getInitialDensity();
  double rho = state->density;
  double p = computePressure(rho_0, rho);
  return p;
}

// Compute pressure (option 1)
double 
Pressure_Air::computePressure(const double& rho_orig,
                              const double& rho_cur)
{
  double J = rho_orig/rho_cur;
  double eps_v = -std::log(J);
  double p = d_p0*(std::exp(d_gamma*eps_v) - 1.0);
  return p;
}

// Compute pressure (option 2)
void 
Pressure_Air::computePressure(const double& rho_orig,
                              const double& rho_cur,
                              double& pressure,
                              double& dp_drho,
                              double& csquared)
{
  double J = rho_orig/rho_cur;
  double eps_v = -std::log(J);
  pressure = d_p0*(std::exp(d_gamma*eps_v) - 1.0);
  double dp_dJ = -d_gamma*pressure/J;
  dp_drho = -dp_dJ*rho_orig/(rho_cur*rho_cur);
  csquared = dp_dJ/rho_cur;
}

// Compute derivative of pressure 
double 
Pressure_Air::eval_dp_dJ(const Uintah::MPMMaterial* matl,
                         const double& detF, 
                         const ModelStateBase* state_input)
{
  const ModelState_MasonSand* state = dynamic_cast<const ModelState_MasonSand*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_MasonSand.";
    throw SCIRun::InternalError(out.str(), __FILE__, __LINE__);
  }

  double J = detF;
  double eps_v = -std::log(J);
  double p = d_p0*(std::exp(d_gamma*eps_v) - 1.0);
  double dpdJ = -d_gamma*p/J;
  return dpdJ;
}

// Compute bulk modulus
double 
Pressure_Air::computeInitialBulkModulus()
{
  d_bulkModulus = d_gamma*d_p0;  
  return d_bulkModulus;
}

double 
Pressure_Air::computeBulkModulus(const double& pressure)
{
  d_bulkModulus = d_gamma*(pressure + d_p0);
  return d_bulkModulus;
}

double 
Pressure_Air::computeBulkModulus(const double& rho_orig,
                                 const double& rho_cur)
{
  double p = computePressure(rho_orig, rho_cur);
  d_bulkModulus = computeBulkModulus(p);
  return d_bulkModulus;
}

double 
Pressure_Air::computeBulkModulus(const ModelStateBase* state_input)
{
  const ModelState_MasonSand* state = dynamic_cast<const ModelState_MasonSand*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_MasonSand.";
    throw SCIRun::InternalError(out.str(), __FILE__, __LINE__);
  }

  double p = -state->I1/3.0;
  d_bulkModulus = computeBulkModulus(p);
  return d_bulkModulus;
}

// Compute strain energy
double 
Pressure_Air::computeStrainEnergy(const double& rho_orig,
                                  const double& rho_cur)
{
  throw InternalError("ComputeStrainEnergy has not been implemented yet for Air.",
                      __FILE__, __LINE__);
  return 0.0;
}

double 
Pressure_Air::computeStrainEnergy(const ModelStateBase* state)
{
  throw InternalError("ComputeStrainEnergy has not been implemented yet for Air.",
                      __FILE__, __LINE__);
  return 0.0;
}


// Compute density given pressure (tension +ve)
double 
Pressure_Air::computeDensity(const double& rho_orig,
                             const double& pressure)
{
  throw InternalError("ComputeDensity has not been implemented yet for Air.",
                      __FILE__, __LINE__);
  return 0.0;
}

//  Calculate the derivative of p with respect to epse_v
double 
Pressure_Air::computeDpDepse_v(const ModelStateBase* state_input) const
{
  const ModelState_MasonSand* state = dynamic_cast<const ModelState_MasonSand*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_MasonSand.";
    throw SCIRun::InternalError(out.str(), __FILE__, __LINE__);
  }

  double p = -state->I1/3.0;
  double dp_depse_v = d_gamma*(p + d_p0);
  return dp_depse_v;
}
