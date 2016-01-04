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


#include <CCA/Components/MPM/ConstitutiveModel/Models/Pressure_Granite.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_Arenisca3.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <cmath>

using namespace Uintah;
using namespace Vaango;

Pressure_Granite::Pressure_Granite()
{
  d_p0 = 101325;  // Hardcoded (SI units).  *TODO* Get as input with ProblemSpec later.
  d_K0 = 40e9;
  d_n = 4.0;
  d_bulk = d_K0;  
} 

Pressure_Granite::Pressure_Granite(Uintah::ProblemSpecP&)
{
  d_p0 = 101325;  // Hardcoded (SI units).  *TODO* Get as input with ProblemSpec later.
  d_K0 = 40e9;
  d_n = 4.0;
  d_bulk = d_K0;  
} 
         
Pressure_Granite::Pressure_Granite(const Pressure_Granite* cm)
{
  d_p0 = cm->d_p0;
  d_K0 = cm->d_K0;
  d_n = cm->d_n;
  d_bulk = cm->d_bulk;
} 
         
Pressure_Granite::~Pressure_Granite()
{
}

void Pressure_Granite::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP eos_ps = ps->appendChild("pressure_model");
  eos_ps->setAttribute("type","granite");
}
         
//////////
// Calculate the pressure using the elastic constitutive equation
double 
Pressure_Granite::computePressure(const Uintah::MPMMaterial* matl,
                                const ModelStateBase* state_input,
                                const Uintah::Matrix3& ,
                                const Uintah::Matrix3& rateOfDeformation,
                                const double& delT)
{
  const ModelState_Arenisca3* state = dynamic_cast<const ModelState_Arenisca3*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Arenisca3.";
    throw SCIRun::InternalError(out.str(), __FILE__, __LINE__);
  }

  double rho_0 = matl->getInitialDensity();
  double rho = state->density;
  double p = computePressure(rho_0, rho);
  return p;
}

// Compute pressure (option 1)
double 
Pressure_Granite::computePressure(const double& rho_orig,
                                const double& rho_cur)
{
  double J = rho_orig/rho_cur;
  double p = d_p0 + d_K0/d_n*(std::pow(J, -d_n) - 1);
  return p;
}

// Compute pressure (option 2)
void 
Pressure_Granite::computePressure(const double& rho_orig,
                                const double& rho_cur,
                                double& pressure,
                                double& dp_drho,
                                double& csquared)
{
  double J = rho_orig/rho_cur;
  pressure = d_p0 + d_K0/d_n*(std::pow(J, -d_n) - 1);
  double dp_dJ = -d_K0*std::pow(J, -(d_n+1));
  dp_drho = d_K0/rho_cur*std::pow(J, -d_n);
  csquared = dp_dJ/rho_cur;
}

// Compute derivative of pressure 
double 
Pressure_Granite::eval_dp_dJ(const Uintah::MPMMaterial* matl,
                           const double& detF, 
                           const ModelStateBase* state_input)
{
  const ModelState_Arenisca3* state = dynamic_cast<const ModelState_Arenisca3*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Arenisca3.";
    throw SCIRun::InternalError(out.str(), __FILE__, __LINE__);
  }

  double J = detF;
  double dpdJ = -d_K0*std::pow(J, -(d_n+1));
  return dpdJ;
}

// Compute bulk modulus
double 
Pressure_Granite::computeInitialBulkModulus()
{
  return d_K0;  
}

double 
Pressure_Granite::computeBulkModulus(const double& pressure)
{
  double bulk = d_K0;
  if (pressure > 0.0) {
    bulk += d_n*(pressure - d_p0);
  }
  return bulk;
}

double 
Pressure_Granite::computeBulkModulus(const double& rho_orig,
                                   const double& rho_cur)
{
  double p = computePressure(rho_orig, rho_cur);
  double bulk = computeBulkModulus(p);
  return bulk;
}

double 
Pressure_Granite::computeBulkModulus(const ModelStateBase* state_input)
{
  const ModelState_Arenisca3* state = dynamic_cast<const ModelState_Arenisca3*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Arenisca3.";
    throw SCIRun::InternalError(out.str(), __FILE__, __LINE__);
  }

  double p = -state->I1/3.0;
  double bulk = computeBulkModulus(p);
  return bulk;
}

// Compute strain energy
double 
Pressure_Granite::computeStrainEnergy(const double& rho_orig,
                                    const double& rho_cur)
{
  throw InternalError("ComputeStrainEnergy has not been implemented yet for Granite.",
    __FILE__, __LINE__);
  return 0.0;
}

double 
Pressure_Granite::computeStrainEnergy(const ModelStateBase* state)
{
  throw InternalError("ComputeStrainEnergy has not been implemented yet for Granite.",
    __FILE__, __LINE__);
  return 0.0;
}


// Compute density given pressure (tension +ve)
double 
Pressure_Granite::computeDensity(const double& rho_orig,
                               const double& pressure)
{
  throw InternalError("ComputeDensity has not been implemented yet for Granite.",
    __FILE__, __LINE__);
  return 0.0;
}

//  Calculate the derivative of p with respect to epse_v
double 
Pressure_Granite::computeDpDepse_v(const ModelStateBase* state_input) const
{
  const ModelState_Arenisca3* state = dynamic_cast<const ModelState_Arenisca3*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Arenisca3.";
    throw SCIRun::InternalError(out.str(), __FILE__, __LINE__);
  }

  double p = -state->I1/3.0;
  double dp_depse_v = d_K0 + d_n*(p - d_p0);
  return dp_depse_v;
}
