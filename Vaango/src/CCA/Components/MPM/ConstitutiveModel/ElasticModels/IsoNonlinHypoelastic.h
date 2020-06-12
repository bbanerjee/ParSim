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

#ifndef __ISOTROPIC_NONLINEAR_HYPOELASTICITY_MODEL_H__
#define __ISOTROPIC_NONLINEAR_HYPOELASTICITY_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
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

  IsoNonlinHypoElastic();
  IsoNonlinHypoelastic(const ElasticModuliModel* elastic);
  IsoNonlinHypoelastic(const PressureModel* pressure,
                       const ShearModulusModel* shear);
  IsoNonlinHypoelastic(const MPMEquationOfState* eos,
                       const ShearModulusModel* shear);
  IsoNonlinHypoElastic(const IsoNonlinHypoelastic& model) = delete;
  IsoNonlinHypoelastic&
  operator=(const IsoNonlinHypoelastic& model) = delete;

  virtual ~IsoNonlinHypoelastic();

  Uintah::Matrix3
  computeStress(const Matrix3& deformRate,
                const ModelStateBase* state); 

  Uintah::Matrix3
  computeDStressDIntVar(const InternalVariableModel* intvar,
                        const ModelStateBase* state);

private:

  ElasticModuliModel* d_elastic;
  PressureModel* d_pressure;
  MPMEquationOfState* d_eos;
  ShearModulusModel* d_shear;

};
} // End namespace Uintah

#endif // __ISOTROPIC_NONLINEAR_HYPOELASTICITY_MODEL_H__
