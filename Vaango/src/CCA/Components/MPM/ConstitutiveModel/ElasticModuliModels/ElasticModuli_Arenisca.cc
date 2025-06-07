/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_Arenisca.h>
#include <Core/Exceptions/InternalError.h>

using namespace Uintah;
using namespace Vaango;

// Construct a default elasticity model.
/* From Arenisca3:
 * If the user has specified a nonzero G1 and G2, these are used to define a
 * pressure
 * dependent poisson ratio, which is used to adjust the shear modulus along with
 * the
 * bulk modulus.  The high pressure limit has nu=G1+G2;
 * if ((d_cm.G1!=0.0)&&(d_cm.G2!=0.0)){
 *     // High pressure bulk modulus:
 *     double nu = d_cm.G1+d_cm.G2;
 *     shear = 1.5*bulk*(1.0-2.0*nu)/(1.0+nu);
 * }
 */

ElasticModuli_Arenisca::ElasticModuli_Arenisca(Uintah::ProblemSpecP& ps)
{
  ps->require("B0", d_bulk.B0); // Tangent Elastic Bulk Modulus Parameter
  ps->require("B1", d_bulk.B1); // Tangent Elastic Bulk Modulus Parameter
  ps->require("B2", d_bulk.B2); // Tangent Elastic Bulk Modulus Parameter
  ps->require("B3", d_bulk.B3); // Tangent Elastic Bulk Modulus Parameter
  ps->require("B4", d_bulk.B4); // Tangent Elastic Bulk Modulus Parameter
  ps->getWithDefault("B01", d_bulk.B01,
                     0.0); // Tangent Elastic Bulk Modulus Parameter

  ps->require("G0", d_shear.G0); // Tangent Elastic Shear Modulus Parameter
  ps->require("G1", d_shear.G1); // Tangent Elastic Shear Modulus Parameter
  ps->require("G2", d_shear.G2); // Tangent Elastic Shear Modulus Parameter
  ps->require("G3", d_shear.G3); // Tangent Elastic Shear Modulus Parameter
  ps->require("G4", d_shear.G4); // Tangent Elastic Shear Modulus Parameter
}

// Construct a copy of a elasticity model.
ElasticModuli_Arenisca::ElasticModuli_Arenisca(
  const ElasticModuli_Arenisca* model)
{
  d_bulk = model->d_bulk;
  d_shear = model->d_shear;
}

// Destructor of elasticity model.
ElasticModuli_Arenisca::~ElasticModuli_Arenisca() = default;

void
ElasticModuli_Arenisca::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP elasticModuli_ps = ps->appendChild("elastic_moduli_model");
  elasticModuli_ps->setAttribute("type", "arenisca");

  elasticModuli_ps->appendElement("B0", d_bulk.B0);
  elasticModuli_ps->appendElement("B01", d_bulk.B01);
  elasticModuli_ps->appendElement("B1", d_bulk.B1);
  elasticModuli_ps->appendElement("B2", d_bulk.B2);
  elasticModuli_ps->appendElement("B3", d_bulk.B3);
  elasticModuli_ps->appendElement("B4", d_bulk.B4);
  elasticModuli_ps->appendElement("G0", d_shear.G0);
  elasticModuli_ps->appendElement("G1",
                                  d_shear.G1); // Low pressure Poisson ratio
  elasticModuli_ps->appendElement(
    "G2", d_shear.G2); // Pressure-dependent Poisson ratio term
  elasticModuli_ps->appendElement("G3", d_shear.G3); // Not used
  elasticModuli_ps->appendElement("G4", d_shear.G4); // Not used
}

// Compute the elasticity
ElasticModuli
ElasticModuli_Arenisca::getInitialElasticModuli() const
{
  return ElasticModuli(d_bulk.B0, d_shear.G0);
}

ElasticModuli
ElasticModuli_Arenisca::getElasticModuliUpperBound() const
{
  return ElasticModuli(d_bulk.B0 + d_bulk.B1, d_shear.G0);
}

ElasticModuli
ElasticModuli_Arenisca::getElasticModuliLowerBound() const
{
  return ElasticModuli(d_bulk.B0, d_shear.G0);
}

ElasticModuli
ElasticModuli_Arenisca::getCurrentElasticModuli(
  const ModelStateBase* state_input) const
{
  const ModelState_Arenisca3* state =
    static_cast<const ModelState_Arenisca3*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Arenisca3.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  // Make sure the quantities are positive in compression
  double I1_bar = -state->I1;
  double ev_p_bar = -(state->plasticStrainTensor).Trace();

  double KK = d_bulk.B0;
  double GG = d_shear.G0;
  if (I1_bar > 0.0) {
    double exp_B2_by_I1 = std::exp(-d_bulk.B2 / I1_bar);
    KK += d_bulk.B01 * I1_bar + d_bulk.B1 * exp_B2_by_I1;
    double nu = d_shear.G1 + d_shear.G2 * exp_B2_by_I1;
    GG = (nu > 0.0) ? 1.5 * KK * (1.0 - 2.0 * nu) / (1.0 + nu) : GG;
  }

  if (ev_p_bar > 0) {
    KK -= d_bulk.B3 * std::exp(-d_bulk.B4 / ev_p_bar);
  }

  return ElasticModuli(KK, GG);
}

/*! Compute derivatives of moduli with respect to internal variables */
std::vector<ElasticModuli> 
ElasticModuli_Arenisca::computeDModuliDIntVar(const ModelStateBase* state) const
{
  std::vector<ElasticModuli> derivs;
  return derivs;
}

/*! Compute moduli and derivatives of moduli with respect to internal variables */
std::pair<ElasticModuli, std::vector<ElasticModuli>>
ElasticModuli_Arenisca::computeModuliAndDModuliDIntVar(const ModelStateBase* state) const
{
  ElasticModuli moduli = getCurrentElasticModuli(state);
  std::vector<ElasticModuli> derivs;
  return std::make_pair(moduli, derivs);
}
