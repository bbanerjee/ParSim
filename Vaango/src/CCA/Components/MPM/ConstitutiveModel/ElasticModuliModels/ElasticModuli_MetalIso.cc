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

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_MetalIso.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/MPMEquationOfStateFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/MPMEquationOfState.h>
#include <CCA/Components/MPM/ConstitutiveModel/ShearModulusModels/ShearModulusModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/ShearModulusModels/ShearModulusModel.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/ProblemSetupException.h>

#include <iostream>

using namespace Vaango;

ElasticModuli_MetalIso::ElasticModuli_MetalIso(Uintah::ProblemSpecP& ps)
{
  d_eos = MPMEquationOfStateFactory::create(ps);

  // Don't allow unphysical EOS models for metals
  if (d_eos->materialType() != EOSMaterialType::ALL) {  
    std::ostringstream err;
    err << "**ERROR** Fluid/soil/rock equations of state cannot be used for metal elasticity."
           " Please correct the input file.\n";
    throw Uintah::ProblemSetupException(err.str(), __FILE__, __LINE__);
  }

  d_shear = Vaango::ShearModulusModelFactory::create(ps, d_eos.get());

  d_Km = d_eos->computeInitialBulkModulus();
  d_Gm = d_shear->computeInitialShearModulus();

  if (!(d_Km > 0.0 && d_Gm > 0.0)) {
    std::ostringstream err;
    err << "**ERROR** The model state in the elastic modulus model"
           " has not been initialized correctly. Please contact the developers.\n";
    throw Uintah::InternalError(err.str(), __FILE__, __LINE__);
  }

  double nu = (3.0 * d_Km - 2.0 * d_Gm)/(6.0 * d_Km + 2.0 * d_Gm);
  if (nu < -1.0 || nu > 0.5) {
    std::ostringstream err;
    err << "**ERROR** The Poisson's ratio (" << nu << ") of the material is unphysical."
           " K = " << d_Km << " and G = " << d_Gm << " Please correct input data.\n";
    throw Uintah::ProblemSetupException(err.str(), __FILE__, __LINE__);
  }

}

ElasticModuli_MetalIso::ElasticModuli_MetalIso(MPMEquationOfState* eos, 
                                               ShearModulusModel* shear)
{
  d_eos = MPMEquationOfStateFactory::createCopy(eos);
  d_shear = ShearModulusModelFactory::createCopy(shear);
  if (d_eos == nullptr || d_shear == nullptr) {
    std::ostringstream err;
    err << "**ERROR** The model state in the elastic modulus model"
           " for ElasticModuli_MetalIso has not been initialized correctly. "
           " Please contact the developers.\n";
    throw Uintah::InternalError(err.str(), __FILE__, __LINE__);
  }

  d_Km = d_eos->computeInitialBulkModulus();
  d_Gm = d_shear->computeInitialShearModulus();

  if (!(d_Km > 0.0 && d_Gm > 0.0)) {
    std::ostringstream err;
    err << "**ERROR** The model state in the elastic modulus model"
           " has not been initialized correctly. Please contact the developers.\n";
    throw Uintah::InternalError(err.str(), __FILE__, __LINE__);
  }

  double nu = (3.0 * d_Km - 2.0 * d_Gm)/(6.0 * d_Km + 2.0 * d_Gm);
  if (nu < -1.0 || nu > 0.5) {
    std::ostringstream err;
    err << "**ERROR** The Poisson's ratio (" << nu << ") of the material is unphysical."
           " K = " << d_Km << " and G = " << d_Gm << " Please correct input data.\n";
    throw Uintah::ProblemSetupException(err.str(), __FILE__, __LINE__);
  }
}

// Construct a copy of a elasticity model.
ElasticModuli_MetalIso::ElasticModuli_MetalIso(
  const ElasticModuli_MetalIso* model)
{
  d_eos = MPMEquationOfStateFactory::createCopy(model->d_eos.get());
  d_shear = ShearModulusModelFactory::createCopy(model->d_shear.get());
  
  d_Km = model->d_Km;
  d_Gm = model->d_Gm;

  d_rho0 = model->d_rho0;
}

void
ElasticModuli_MetalIso::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP elasticModuli_ps = ps->appendChild("elastic_moduli_model");
  elasticModuli_ps->setAttribute("type", "metal_iso");

  d_eos->outputProblemSpec(ps);
  d_shear->outputProblemSpec(ps);
}

std::map<std::string, double> 
ElasticModuli_MetalIso::getParameters() const
{
  std::map<std::string, double> params = d_eos->getParameters();
  std::map<std::string, double> shear_params = d_shear->getParameters();
  params.insert(shear_params.begin(), shear_params.end());
  return params;
}

// Compute the elasticity
ElasticModuli
ElasticModuli_MetalIso::getInitialElasticModuli() const
{
  return ElasticModuli(d_Km, d_Gm);
}

ElasticModuli
ElasticModuli_MetalIso::getElasticModuliUpperBound() const
{
  ModelStateBase state;
  state.initialDensity = 1000;
  state.density = 1.0e2 * 1000 ;
  state.initialBulkModulus = d_Km;
  state.initialShearModulus = d_Gm;
  state.temperature = 0.01;
  state.meltingTemp = 1.0e9;
  state.pressure = 1.0e6;
  state.porosity = 0.0;
  double Km = d_eos->computeBulkModulus(&state);
  double Gm = d_shear->computeShearModulus(&state);
  return ElasticModuli(Km, Gm);
}

ElasticModuli
ElasticModuli_MetalIso::getElasticModuliLowerBound() const
{
  ModelStateBase state;
  state.initialDensity = 1000;
  state.density = 0.01 * 1000;
  state.initialShearModulus = d_Gm;
  state.meltingTemp = 1.0e4;
  state.temperature = state.meltingTemp - 1.0;
  state.pressure = 0.01;
  state.porosity = 0.0;
  double Km = d_eos->computeBulkModulus(&state);
  double Gm = d_shear->computeShearModulus(&state);
  return ElasticModuli(Km, Gm);
}

ElasticModuli
ElasticModuli_MetalIso::getCurrentElasticModuli(const ModelStateBase* state) const
{
  double Km = d_eos->computeBulkModulus(state);
  double Gm = d_shear->computeShearModulus(state);

  double kappa_m = (4.0 * Gm + 3.0 * Km)/(4.0 * Gm);
  double gamma_m = 5.0 * (4.0 * Gm + 3.0 * Km)/(8.0 * Gm + 9.0 * Km);

  double phi = state->porosity;
  double Kfac = (1.0 - phi)/(1.0 + (kappa_m - 1.0) * phi);
  double Gfac = (1.0 - phi)/(1.0 + (gamma_m - 1.0) * phi);

  double K = Km * Kfac;
  double G = Gm * Gfac;
  return ElasticModuli(K, G);
}

/*! Compute derivatives of moduli with respect to internal variables */
std::vector<ElasticModuli> 
ElasticModuli_MetalIso::computeDModuliDIntVar(const ModelStateBase* state) const
{
  std::vector<ElasticModuli> derivs;

  double dK_dep = 0.0;
  double dG_dep = 0.0;
  derivs.push_back(ElasticModuli(dK_dep, dG_dep));

  double Km = d_eos->computeBulkModulus(state);
  double Gm = d_shear->computeShearModulus(state);

  double kappa_m = (4.0 * Gm + 3.0 * Km)/(4.0 * Gm);
  double gamma_m = 5.0 * (4.0 * Gm + 3.0 * Km)/(8.0 * Gm + 9.0 * Km);

  double phi = state->porosity;
  double Kfac = 1.0 + (kappa_m - 1.0) * phi;
  double Gfac = 1.0 + (gamma_m - 1.0) * phi;

  double dK_dphi = - (Km * kappa_m)/(Kfac * Kfac);
  double dG_dphi = - (Gm * gamma_m)/(Gfac * Gfac);

  derivs.push_back(ElasticModuli(dK_dphi, dG_dphi));
  return derivs;
}

/*! Compute moduli and derivatives of moduli with respect to internal variables */
std::pair<ElasticModuli, std::vector<ElasticModuli>>
ElasticModuli_MetalIso::computeModuliAndDModuliDIntVar(const ModelStateBase* state) const
{
  std::vector<ElasticModuli> derivs;

  double dK_dep = 0.0;
  double dG_dep = 0.0;
  derivs.push_back(ElasticModuli(dK_dep, dG_dep));

  double Km = d_eos->computeBulkModulus(state);
  double Gm = d_shear->computeShearModulus(state);

  double kappa_m = (4.0 * Gm + 3.0 * Km)/(4.0 * Gm);
  double gamma_m = 5.0 * (4.0 * Gm + 3.0 * Km)/(8.0 * Gm + 9.0 * Km);

  double phi = state->porosity;
  double Kfac = 1.0 + (kappa_m - 1.0) * phi;
  double Gfac = 1.0 + (gamma_m - 1.0) * phi;

  double K = Km * (1.0 - phi) / Kfac;
  double G = Gm * (1.0 - phi) / Gfac;

  double dK_dphi = - (Km * kappa_m)/(Kfac * Kfac);
  double dG_dphi = - (Gm * gamma_m)/(Gfac * Gfac);

  derivs.push_back(ElasticModuli(dK_dphi, dG_dphi));
  return std::make_pair(ElasticModuli(K, G), derivs);
}
