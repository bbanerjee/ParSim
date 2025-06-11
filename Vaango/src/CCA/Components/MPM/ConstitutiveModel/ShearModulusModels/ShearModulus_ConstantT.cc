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

#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/EOS_NullT.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_MetalT.h>
#include <CCA/Components/MPM/ConstitutiveModel/ShearModulusModels/ShearModulus_ConstantT.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/ProblemSpec/ProblemSpec.h>

using namespace Vaango;

/*! Construct a constant shear modulus model. */
ShearModulus_ConstantT::ShearModulus_ConstantT(Uintah::ProblemSpecP& ps)
{
  d_eos = nullptr;
  ps->require("shear_modulus", d_shearModulus);
}

/*! Construct a copy of constant shear modulus model. */
ShearModulus_ConstantT::ShearModulus_ConstantT(
  const ShearModulus_ConstantT* smm)
{
  d_eos          = smm->d_eos;
  d_shearModulus = smm->d_shearModulus;
}

void
ShearModulus_ConstantT::l_outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP shear_ps = ps->appendChild("shear_modulus_model");
  shear_ps->setAttribute("type", "constant_shear");
  shear_ps->appendElement("shear_modulus", d_shearModulus);
}

/*! Get parameters */
ParameterDict
ShearModulus_ConstantT::l_getParameters() const
{
  ParameterDict params;
  params["shear_modulus"] = d_shearModulus;
  return params;
}

/*! Compute the shear modulus */
double
ShearModulus_ConstantT::l_computeInitialShearModulus()
{
  return d_shearModulus;
}

double
ShearModulus_ConstantT::l_computeShearModulus([[maybe_unused]] const ModelState_MetalT* state)
{
  return d_shearModulus;
}

double
ShearModulus_ConstantT::l_computeShearModulus(
  [[maybe_unused]] const ModelState_MetalT* state) const
{
  return d_shearModulus;
}

/*! Compute the shear strain energy */
double
ShearModulus_ConstantT::l_computeStrainEnergy(const ModelState_MetalT* state)
{
  return state->energy;
}

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
ShearModulus_ConstantT::l_computeQ(const ModelState_MetalT* state) const
{
  return state->q;
}

/////////////////////////////////////////////////////////////////////////
/*
  Compute dq/depse_s
*/
/////////////////////////////////////////////////////////////////////////
double
ShearModulus_ConstantT::l_computeDqDepse_s([[maybe_unused]] const ModelState_MetalT* state) const
{
  return 3.0 * d_shearModulus;
}

/////////////////////////////////////////////////////////////////////////
/*
  Compute dq/depse_v
*/
/////////////////////////////////////////////////////////////////////////
double
ShearModulus_ConstantT::l_computeDqDepse_v([[maybe_unused]] const ModelState_MetalT* state) const
{
  return 0.0;
}
