/*
 * The MIT License
 *
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#ifndef ___ELASTIC_MODULUS_ISOTROPIC_METAL_H__
#define ___ELASTIC_MODULUS_ISOTROPIC_METAL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuliModel.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <memory>

namespace Vaango {

/*! \class ElasticModuli_MetalIso
 *  \brief An isotropic elastic metal model with porosity
 *  \author Biswajit Banerjee,
 *
*/
class ShearModulusModel;
class MPMEquationOfState;

class ElasticModuli_MetalIso : public ElasticModuliModel
{
public:

  ElasticModuli_MetalIso(Uintah::ProblemSpecP& ps);
  ElasticModuli_MetalIso(MPMEquationOfState* eos, ShearModulusModel* shear);
  ElasticModuli_MetalIso(const ElasticModuli_MetalIso* emm);
  ~ElasticModuli_MetalIso() override = default;
  ElasticModuli_MetalIso& operator=(const ElasticModuli_MetalIso& emm) = delete;

  void outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  /*! Get parameters */
  std::map<std::string, double> getParameters() const override;

  /*! Compute the elasticity */
  ElasticModuli getInitialElasticModuli() const override;
  ElasticModuli getCurrentElasticModuli(const ModelStateBase* state) const override;
  ElasticModuli getElasticModuliLowerBound() const override;
  ElasticModuli getElasticModuliUpperBound() const override;

  /*! Compute derivatives of moduli with respect to internal variables */
  std::vector<ElasticModuli> computeDModuliDIntVar(const ModelStateBase* state) const override;

  /*! Compute moduli and derivatives of moduli with respect to internal variables */
  std::pair<ElasticModuli, std::vector<ElasticModuli>>
  computeModuliAndDModuliDIntVar(const ModelStateBase* state) const override;

private:

  /* For tangent bulk modulus parameters */
  std::unique_ptr<MPMEquationOfState> d_eos;

  /* For tangent shear modulus parameters */
  ShearModulusModel* d_shear;

  double d_Km, d_Gm, d_rho0;

};
} // End namespace Vaango

#endif // __ELASTIC_MODULUS_ISOTROPIC_METAL_H__
