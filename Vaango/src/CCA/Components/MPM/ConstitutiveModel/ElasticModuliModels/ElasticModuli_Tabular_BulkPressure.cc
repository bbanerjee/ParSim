/*
 * The MIT License
 *
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

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_Tabular_BulkPressure.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Tabular.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ProblemSetupException.h>

//#define DEBUG_INTERPOLATION

using namespace Vaango;

// Construct a default elasticity model.
ElasticModuli_Tabular_BulkPressure::ElasticModuli_Tabular_BulkPressure(Uintah::ProblemSpecP& ps)
  : d_bulk(ps)
{
  ps->require("G0", d_shear.G0);
  ps->require("nu", d_shear.nu);

  checkInputParameters();
}

//--------------------------------------------------------------
// Check that the input parameters are reasonable
//--------------------------------------------------------------
void
ElasticModuli_Tabular_BulkPressure::checkInputParameters()
{
  std::ostringstream warn;
  if (d_shear.G0 <= 0.0) {
    warn << "G0 must be positive. G0 = " << d_shear.G0 << std::endl;
    throw Uintah::ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
}

// Construct a copy of a elasticity model.
ElasticModuli_Tabular_BulkPressure::ElasticModuli_Tabular_BulkPressure(const ElasticModuli_Tabular_BulkPressure* model)
{
  d_bulk = model->d_bulk;
  d_shear = model->d_shear;
}

void
ElasticModuli_Tabular_BulkPressure::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP elasticModuli_ps = ps->appendChild("elastic_moduli_model");
  elasticModuli_ps->setAttribute("type", "tabular_bulk_pressure");

  d_bulk.table.outputProblemSpec(elasticModuli_ps);

  elasticModuli_ps->appendElement("G0", d_shear.G0);
  elasticModuli_ps->appendElement("nu", d_shear.nu);
}

// Compute the elastic moduli
ElasticModuli
ElasticModuli_Tabular_BulkPressure::getInitialElasticModuli() const
{
  double K = computeBulkModulus(0, 0);
  double G = computeShearModulus(K);
  return ElasticModuli(K, G);
}

ElasticModuli
ElasticModuli_Tabular_BulkPressure::getCurrentElasticModuli(const ModelStateBase* state_input) const
{
  const ModelState_Tabular* state =
    static_cast<const ModelState_Tabular*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  // Make sure the quantities are positive in compression
  double p_bar = -(state->stressTensor).Trace() / 3.0;
  double ev_p_bar = -(state->plasticStrainTensor).Trace();

  // Compute the elastic moduli
  double K = computeBulkModulus(p_bar, ev_p_bar);
  double G = computeShearModulus(K);

  #ifdef DEBUG_INTERPOLATION
    if (K < 0 || !std::isfinite(K)) {
      std::cout << "p = " << p_bar << " ev_p = " << ev_p_bar << "\n";
      std::cout << " K = " << K << " G = " << G << std::endl;
    }
  #endif

  return ElasticModuli(K, G);
}

double 
ElasticModuli_Tabular_BulkPressure::computeBulkModulus(const double& pressure,
                                                        const double& plasticVolStrain) const
{
  double p_bar = (pressure < 0) ? 0 : pressure;
  double ev_p_bar = (plasticVolStrain < 0) ? 0 : plasticVolStrain;

  #ifdef DEBUG_INTERPOLATION
  std::cout << "p = (" << pressure << ", " << p_bar 
            << ") ev_p_bar = (" << plasticVolStrain << "," << ev_p_bar << ")\n";
  #endif

  DoubleVec1D K;
  try {
    K = 
      d_bulk.table.interpolate<2>({{ev_p_bar, p_bar}});
  } catch (Uintah::InvalidValue& e) {
    std::ostringstream out;
    out << "**WARNING** In computeBulkModulus (Low):"
        << " pressure = " << pressure
        << " plasticVolStrain = " << plasticVolStrain << "\n"
        << e.message();
    throw Uintah::InvalidValue(out.str(), __FILE__, __LINE__);
  }

  #ifdef DEBUG_INTERPOLATION
    if (K[0] < 0 || !std::isfinite(K[0])) {
      std::cout <<std::setprecision(16) 
                << "p = " << pressure
                << " ev_p = " << plasticVolStrain
                << " K = " << K[0] << "\n";
    }
  #endif

  return K[0];
}

double 
ElasticModuli_Tabular_BulkPressure::computeShearModulus(const double& K) const
{
  double nu = d_shear.nu;
  double G = (nu > -1.0 && nu < 0.5) 
             ? 1.5*K*(1.0 - 2.0*nu)/(1.0 + nu) 
             : d_shear.G0;
  return G;
}

/* Get the elastic moduli and their derivatives with respect to a single
   plastic internal variable */
std::pair<ElasticModuli, ElasticModuli>
ElasticModuli_Tabular_BulkPressure::getElasticModuliAndDerivatives(const ModelStateBase* state_input) const
{
  const ModelState_Tabular* state =
    dynamic_cast<const ModelState_Tabular*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  // Make sure the quantities are positive in compression
  double p_bar = -(state->stressTensor).Trace() / 3.0;
  double ev_p_bar = -(state->plasticStrainTensor).Trace();

  // Compute the elastic moduli
  double epsilon = 1.0e-8;
  double K = computeBulkModulus(p_bar, ev_p_bar);
  double K_min = computeBulkModulus(p_bar, ev_p_bar-epsilon);
  double K_max = computeBulkModulus(p_bar, ev_p_bar+epsilon);
  double G = computeShearModulus(K);
  double G_min = computeShearModulus(K_min);
  double G_max = computeShearModulus(K_max);

  // Compute derivatives (negative sign to convert to compression -ve convention)
  double dK_deps_p = -(K_max - K_min)/(2*epsilon);
  double dG_deps_p = -(G_max - G_min)/(2*epsilon);

  return std::make_pair(ElasticModuli(K, G),
                        ElasticModuli(dK_deps_p, dG_deps_p));
}

/*! Compute derivatives of moduli with respect to internal variables */
std::vector<ElasticModuli> 
ElasticModuli_Tabular_BulkPressure::computeDModuliDIntVar(const ModelStateBase* state) const
{
  auto K_dK = getElasticModuliAndDerivatives(state);
  std::vector<ElasticModuli> derivs;
  derivs.push_back(K_dK.second);
  return derivs;
}

/*! Compute moduli and derivatives of moduli with respect to internal variables */
std::pair<ElasticModuli, std::vector<ElasticModuli>>
ElasticModuli_Tabular_BulkPressure::computeModuliAndDModuliDIntVar(const ModelStateBase* state) const
{
  auto K_dK = getElasticModuliAndDerivatives(state);
  auto moduli = K_dK.first;
  std::vector<ElasticModuli> derivs;
  derivs.push_back(K_dK.second);
  return std::make_pair(moduli, derivs);
}
