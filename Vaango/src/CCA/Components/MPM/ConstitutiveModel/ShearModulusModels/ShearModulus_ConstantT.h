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

#ifndef __TEMPLATED_CONSTANT_SHEAR_MODEL_H__
#define __TEMPLATED_CONSTANT_SHEAR_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/EOS_NullT.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_MetalT.h>
#include <CCA/Components/MPM/ConstitutiveModel/ShearModulusModels/ShearModulusT.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

/*!
  \class ShearModulus_ConstantT

  \brief The Constant model for calculating shear stress
 *
*/
class ShearModulus_ConstantT
  : public ShearModulusT<ShearModulus_ConstantT, ModelState_MetalT, EOS_NullT>
{
public:
  /*! Construct a constant shear modulus model. */
  // ShearModulus_ConstantT(Uintah::ProblemSpecP& ps, EOS_NullT* eos);
  ShearModulus_ConstantT(Uintah::ProblemSpecP& ps);

  /*! Construct a copy of constant shear modulus model. */
  ShearModulus_ConstantT(const ShearModulus_ConstantT* smm);

  ShearModulus_ConstantT&
  operator=(const ShearModulus_ConstantT& smm) = delete;

  /*! Destructor of constant shear modulus model.   */
  ~ShearModulus_ConstantT() = default;

  void
  l_outputProblemSpec(Uintah::ProblemSpecP& ps);

  /*! Get parameters */
  ParameterDict
  l_getParameters() const;

  /*! Compute the shear modulus */
  double
  l_computeInitialShearModulus();

  double
  l_computeShearModulus(const ModelState_MetalT* state);

  double
  l_computeShearModulus(const ModelState_MetalT* state) const;

  /*! Compute the shear strain energy */
  double
  l_computeStrainEnergy(const ModelState_MetalT* state);

  /////////////////////////////////////////////////////////////////////////
  /*
    Compute q = 3 mu epse_s
       where mu = shear modulus
             epse_s = sqrt{2/3} ||ee||
             ee = deviatoric part of elastic strain = epse - 1/3 epse_v I
             epse = total elastic strain
             epse_v = tr(epse)
  */
  /////////////////////////////////////////////////////////////////////////
  double
  l_computeQ(const ModelState_MetalT* state) const;

  /////////////////////////////////////////////////////////////////////////
  /*
    Compute dq/depse_s
  */
  /////////////////////////////////////////////////////////////////////////
  double
  l_computeDqDepse_s(const ModelState_MetalT* state) const;

  /////////////////////////////////////////////////////////////////////////
  /*
    Compute dq/depse_v
  */
  /////////////////////////////////////////////////////////////////////////
  double
  l_computeDqDepse_v(const ModelState_MetalT* state) const;
};
} // End namespace Vaango

#endif // __TEMPLATED_CONSTANT_SHEAR_MODEL_H__
