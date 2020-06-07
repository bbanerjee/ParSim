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

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Default.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/Pressure_Hyperelastic.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <cmath>

using namespace Uintah;
using namespace Vaango;

using std::ostringstream;
using std::endl;

Pressure_Hyperelastic::Pressure_Hyperelastic()
{
  d_bulkModulus = -1.0;
}

Pressure_Hyperelastic::Pressure_Hyperelastic(Uintah::ProblemSpecP&)
{
  d_bulkModulus = -1.0;
}

Pressure_Hyperelastic::Pressure_Hyperelastic(const Pressure_Hyperelastic* cm)
{
  d_bulkModulus = cm->d_bulkModulus;
}

Pressure_Hyperelastic::~Pressure_Hyperelastic() = default;

void
Pressure_Hyperelastic::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP eos_ps = ps->appendChild("pressure_model");
  eos_ps->setAttribute("type", "default_Hyper");
}

//////////
// Calculate the pressure using the elastic constitutive equation
double
Pressure_Hyperelastic::computePressure(const Uintah::MPMMaterial* matl,
                                       const ModelStateBase* state_input,
                                       const Uintah::Matrix3&,
                                       const Uintah::Matrix3& rateOfDeformation,
                                       const double& delT)
{
  const ModelState_Default* state =
    static_cast<const ModelState_Default*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Default.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  double rho_0 = matl->getInitialDensity();
  double rho = state->density;
  double J = rho_0 / rho;
  double kappa = state->bulkModulus;

  double p = 0.5 * kappa * (J - 1.0 / J);
  return p;
}

double
Pressure_Hyperelastic::eval_dp_dJ(const Uintah::MPMMaterial* matl,
                                  const double& detF,
                                  const ModelStateBase* state_input)
{
  const ModelState_Default* state =
    static_cast<const ModelState_Default*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Default.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  double J = detF;
  double kappa = state->bulkModulus;

  double dpdJ = 0.5 * kappa * (1.0 + 1.0 / (J * J));
  return dpdJ;
}

// Compute pressure (option 1)
double
Pressure_Hyperelastic::computePressure(const double& rho_orig,
                                       const double& rho_cur)
{
  /*
  if (d_bulkModulus < 0.0) {
    throw InternalError("Please initialize bulk modulus in EOS before computing
  pressure",
                            __FILE__, __LINE__);
  }
  */

  double J = rho_orig / rho_cur;
  double p = 0.5 * d_bulkModulus * (J - 1.0 / J);
  return p;
}

// Compute pressure (option 2)
void
Pressure_Hyperelastic::computePressure(const double& rho_orig,
                                       const double& rho_cur, double& pressure,
                                       double& dp_drho, double& csquared)
{
  /*
  if (d_bulkModulus < 0.0) {
    throw InternalError("Please initialize bulk modulus in EOS before computing
  pressure",
                            __FILE__, __LINE__);
  }
  */

  double J = rho_orig / rho_cur;
  pressure = 0.5 * d_bulkModulus * (J - 1.0 / J);
  double dp_dJ = 0.5 * d_bulkModulus * (1.0 + 1.0 / J * J);
  dp_drho = -0.5 * d_bulkModulus * (1.0 + J * J) / rho_orig;
  csquared = dp_dJ / rho_cur;
}

// Compute bulk modulus
double
Pressure_Hyperelastic::computeInitialBulkModulus()
{
  return d_bulkModulus; // return -1 if uninitialized
}

double
Pressure_Hyperelastic::computeBulkModulus(const double& rho_orig,
                                          const double& rho_cur)
{
  /*
  if (d_bulkModulus < 0.0) {
    throw InternalError("Please initialize bulk modulus in EOS before computing
  modulus",
                            __FILE__, __LINE__);
  }
  */

  double J = rho_orig / rho_cur;
  double bulk = 0.5 * d_bulkModulus * (1.0 + 1.0 / J * J);
  return bulk;
}

// Compute strain energy
double
Pressure_Hyperelastic::computeStrainEnergy(const double& rho_orig,
                                           const double& rho_cur)
{
  /*
  if (d_bulkModulus < 0.0) {
    throw InternalError("Please initialize bulk modulus in EOS before computing
  energy",
                            __FILE__, __LINE__);
  }
  */

  double J = rho_orig / rho_cur;
  double U = 0.5 * d_bulkModulus * (0.5 * (J * J - 1.0) - log(J));
  return U;
}

// Compute density given pressure (tension +ve)
double
Pressure_Hyperelastic::computeDensity(const double& rho_orig,
                                      const double& pressure)
{
  /*
  if (d_bulkModulus < 0.0) {
    throw InternalError("Please initialize bulk modulus in EOS before computing
  density",
                            __FILE__, __LINE__);
  }
  */
  double numer1 = d_bulkModulus * d_bulkModulus + pressure * pressure;
  double sqrtNumer = sqrt(numer1);
  double rho = rho_orig / d_bulkModulus * (-pressure + sqrtNumer);
  if (rho < 0) {
    ostringstream desc;
    desc << "Value of pressure (" << pressure
         << ") is beyond the range of validity of model" << endl
         << "  density = " << rho << endl;
    throw InvalidValue(desc.str(), __FILE__, __LINE__);
  }
  return rho;
}
