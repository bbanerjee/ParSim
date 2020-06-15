/*
 * The MIT License
 *
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

#ifndef ___ELASTIC_MODULUS_ISOTROPIC_METAL_TEMPLATED_H__
#define ___ELASTIC_MODULUS_ISOTROPIC_METAL_TEMPLATED_H__

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuliT.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <iostream>
#include <memory>

namespace Vaango {

/*! \class ElasticModuli_MetalIsoT
 *  \brief An isotropic elastic metal model with porosity
 *  \author Biswajit Banerjee,
 *
*/
template <typename PressureT, typename ShearT>
class ElasticModuli_MetalIsoT : public ElasticModuliT<ElasticModuli_MetalIsoT, ModelState_MetalT>
{
public:

  ElasticModuli_MetalIsoT(Uintah::ProblemSpecP& ps)
  {
    d_eos = std::make_unique<EquationOfStateT<PressureT, ModelState_MetalT>>(ps);
    d_shear = std::make_unique<ShearModulusT<ShearT, ModelState_MetalT, PressureT>(ps);

    d_Km = d_eos->computeInitialBulkModulus();
    d_Gm = d_shear->computeInitialShearModulus();

    if (!(d_Km > 0.0 && d_Gm > 0.0)) {
      std::ostringstream err;
      err << "**ERROR** The model state in the elastic modulus model"
             " has not been initialized correctly. Please contact the developers.\n";
      throw Uintah::InternalError(err.str(), __FILE__, __LINE__);
    }

  }

  ElasticModuli_MetalIsoT(const ElasticModuli_MetalIsoT* emm)
  {
    d_eos = model->d_eos;
    d_shear = model->d_shear;
    d_Km = model->d_Km;
    d_Gm = model->d_Gm;
    d_rho0 = model->d_rho0; 
  }

  ~ElasticModuli_MetalIsoT() = default;

  ElasticModuli_MetalIsoT& operator=(const ElasticModuli_MetalIsoT& emm) = delete;

  void outputProblemSpec(Uintah::ProblemSpecP& ps)
  {
    Uintah::ProblemSpecP elasticModuli_ps = ps->appendChild("elastic_moduli_model");
    elasticModuli_ps->setAttribute("type", "metal_iso");

    d_eos->outputProblemSpec(ps);
    d_shear->outputProblemSpec(ps);
  }

  /*! Get parameters */
  std::map<std::string, double> getParameters() const
  {
    std::map<std::string, double> params = d_eos->getParameters();
    std::map<std::string, double> shear_params = d_shear->getParameters();
    params.insert(shear_params.begin(), shear_params.end());
    return params;
  }

  /*! Compute the elasticity */
  ElasticModuli getInitialElasticModuli() const
  {
    return ElasticModuli(d_Km, d_Gm);
  }

  ElasticModuli getCurrentElasticModuli(const ModelState_MetalT* state)
  {
    double Km = d_eos->computeBulkModulus(state);
    double Gm = d_shear->computeShearModulus(state);
    return ElasticModuli(Km, Gm);
  }

  ElasticModuli getElasticModuliLowerBound() const
  {
    ModelState_MetalT state;
    state.initialDensity = 1000;
    state.density = 0.01 * 1000 ;
    state.initialBulkModulus = d_Km;
    state.initialShearModulus = d_Gm;
    state.meltingTemp = 1.0e4;
    state.temperature = state.meltingTemp - 1.0;
    state.pressure = 0.01;
    double Km = d_eos->computeBulkModulus(&state);
    double Gm = d_shear->computeShearModulus(&state);
    return ElasticModuli(Km, Gm);
  }

  ElasticModuli getElasticModuliUpperBound() const
  {
    ModelState_MetalT state;
    state.initialDensity = 1000;
    state.density = 1.0e2 * 1000 ;
    state.initialBulkModulus = d_Km;
    state.initialShearModulus = d_Gm;
    state.temperature = 0.01;
    state.meltingTemp = 1.0e9;
    state.pressure = 1.0e6;
    double Km = d_eos->computeBulkModulus(&state);
    double Gm = d_shear->computeShearModulus(&state);
    return ElasticModuli(Km, Gm);
  }

private:

  /* For tangent bulk modulus parameters */
  std::unique_ptr<EquationOfStateT<PressureT, ModelState_MetalT>> d_eos;

  /* For tangent shear modulus parameters */
  std::unique_ptr<ShearModulusT<ShearT, ModelState_MetalT, PressureT>> d_shear;

  double d_Km, d_Gm, d_rho0;

};
} // End namespace Vaango

#endif // __ELASTIC_MODULUS_ISOTROPIC_METAL_TEMPLATED_H__
