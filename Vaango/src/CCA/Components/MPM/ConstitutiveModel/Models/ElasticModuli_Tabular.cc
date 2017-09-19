/*
 * The MIT License
 *
 * Copyright (c) 2015 Parresia Research Limited, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/Models/ElasticModuli_Tabular.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_Tabular.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ProblemSetupException.h>

using namespace Vaango;

// Construct a default elasticity model.
ElasticModuli_Tabular::ElasticModuli_Tabular(Uintah::ProblemSpecP& ps)
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
ElasticModuli_Tabular::checkInputParameters()
{
  std::ostringstream warn;

  /* TODO : Add checks for table parameters */
  if (d_shear.G0 <= 0.0) {
    warn << "G0 must be positive. G0 = " << d_shear.G0 << std::endl;
    throw Uintah::ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
}

// Construct a copy of a elasticity model.
ElasticModuli_Tabular::ElasticModuli_Tabular(const ElasticModuli_Tabular* model)
{
  d_bulk = model->d_bulk;
  d_shear = model->d_shear;
}

void
ElasticModuli_Tabular::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP elasticModuli_ps = ps->appendChild("elastic_moduli_model");
  elasticModuli_ps->setAttribute("type", "tabular");

  d_bulk.table.outputProblemSpec(elasticModuli_ps);

  elasticModuli_ps->appendElement("G0", d_shear.G0);
  elasticModuli_ps->appendElement("nu", d_shear.nu);
}

// Compute the elastic moduli
ElasticModuli
ElasticModuli_Tabular::getInitialElasticModuli() const
{
  double K = computeBulkModulus(1.0e-6, 0);
  double G = computeShearModulus(K);
  return ElasticModuli(K, G);
}

ElasticModuli
ElasticModuli_Tabular::getCurrentElasticModuli(const ModelStateBase* state_input)
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
  double ev_e_bar = -(state->elasticStrainTensor).Trace();
  double ev_p_bar = -(state->plasticStrainTensor).Trace();

  // Compute the elastic moduli
  //std::cout << "ev_e = " << ev_e_bar << " ev_p = " << ev_p_bar;
  double K = computeBulkModulus(ev_e_bar, ev_p_bar);
  double G = computeShearModulus(K);
  //std::cout << " K = " << K << " G = " << G << std::endl;

  return ElasticModuli(K, G);
}

double 
ElasticModuli_Tabular::computeBulkModulus(const double& elasticVolStrain,
                                          const double& plasticVolStrain) const
{
  double epsilon = 1.0e-6;
  DoubleVec1D pressure_lo;
  DoubleVec1D pressure_hi;
  try {
    pressure_lo = 
      d_bulk.table.interpolate<2>({{plasticVolStrain, elasticVolStrain-epsilon}});
    pressure_hi = 
      d_bulk.table.interpolate<2>({{plasticVolStrain, elasticVolStrain+epsilon}});
  } catch (Uintah::InvalidValue& e) {
    std::ostringstream out;
    out << "**ERROR** In ElasticModuli_Tabular::computeBulkModulus:"
        << " elasticVolStrain = " << elasticVolStrain
        << " plasticVolStrain = " << plasticVolStrain
        << e.message();
    throw Uintah::InvalidValue(out.str(), __FILE__, __LINE__);
  }
  double K = (pressure_hi[0] - pressure_lo[0])/(2*epsilon);
  return K;
}

double 
ElasticModuli_Tabular::computeShearModulus(const double& K) const
{
  double nu = d_shear.nu;
  double G = (nu > -1.0 && nu < 0.5) 
             ? 1.5*K*(1.0 - 2.0*nu)/(1.0 + nu) 
             : d_shear.G0;
  return G;
}
