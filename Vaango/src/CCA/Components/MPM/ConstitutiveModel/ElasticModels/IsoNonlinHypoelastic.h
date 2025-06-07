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

#ifndef __ISOTROPIC_NONLINEAR_HYPOELASTICITY_MODEL_H__
#define __ISOTROPIC_NONLINEAR_HYPOELASTICITY_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuliModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/InternalVariableModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/DeformationState.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/TensorUtils.h>
#include <Core/ProblemSpec/ProblemSpec.h>

namespace Vaango {

/*
 * Isotropic nonlinear hypoelasticity model
 *
 * Computes Cauchy stress given deformation rate
 * Computes derivative of stress wrt internal variables
 */
class IsoNonlinHypoelastic
{

public:

  IsoNonlinHypoelastic();
  IsoNonlinHypoelastic(const ElasticModuliModel* elastic,
                       const InternalVariableModel* intvar);
  IsoNonlinHypoelastic(const IsoNonlinHypoelastic& model) = delete;
  IsoNonlinHypoelastic&
  operator=(const IsoNonlinHypoelastic& model) = delete;

  virtual ~IsoNonlinHypoelastic() = default;

  Uintah::Matrix3
  computeStress(double delT,
                const Uintah::Matrix3& stress_old,
                const Uintah::DeformationState* deformState,
                const ModelStateBase* modelState) const; 

  std::vector<Uintah::Matrix3>
  computeDStressDIntVar(double delT,
                        const std::vector<Uintah::Matrix3>& derivStress_old,
                        const Uintah::DeformationState* deformState,
                        const ModelStateBase* modelState) const;

  std::pair<Uintah::Matrix3, std::vector<Uintah::Matrix3>>
  computeStressAndDStressDIntVar(double delT,
                                 const Uintah::Matrix3& stress_old,
                                 const std::vector<Uintah::Matrix3>& derivStress_old,
                                 const Uintah::DeformationState* deformState,
                                 const ModelStateBase* modelState) const;

  Tensor::Matrix6Mandel
  computeElasticTangentModulus(const ModelStateBase* state) const;

private:

  const ElasticModuliModel* d_elastic;
  const InternalVariableModel* d_intvar;

};
} // End namespace Vaango

#endif // __ISOTROPIC_NONLINEAR_HYPOELASTICITY_MODEL_H__
