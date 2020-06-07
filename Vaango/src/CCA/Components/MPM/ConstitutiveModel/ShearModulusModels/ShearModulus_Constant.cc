/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Default.h>
#include <CCA/Components/MPM/ConstitutiveModel/ShearModulusModels/ShearModulus_Constant.h>
#include <Core/Exceptions/InternalError.h>

using namespace Uintah;
using namespace Vaango;

// Construct a shear modulus model.
ShearModulus_Constant::ShearModulus_Constant(Uintah::ProblemSpecP& ps,
                                             PressureModel* eos)
{
  d_eos = eos;

  ps->require("shear_modulus", d_shear);
}

// Construct a copy of a shear modulus model.
ShearModulus_Constant::ShearModulus_Constant(const ShearModulus_Constant* smm)
{
  d_eos = smm->d_eos;

  d_shear = smm->d_shear;
}

// Destructor of shear modulus model.
ShearModulus_Constant::~ShearModulus_Constant() = default;

void
ShearModulus_Constant::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP shear_ps = ps->appendChild("elastic_shear_modulus_model");
  shear_ps->setAttribute("type", "constant_shear");

  shear_ps->appendElement("shear_modulus", d_shear);
}

// Compute the shear modulus
double
ShearModulus_Constant::computeInitialShearModulus()
{
  return d_shear;
}

double
ShearModulus_Constant::computeShearModulus(const ModelStateBase* state_input)
{
  const ModelState_Default* state =
    static_cast<const ModelState_Default*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Default.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  d_shear = state->initialShearModulus;
  return d_shear;
}

double
ShearModulus_Constant::computeShearModulus(
  const ModelStateBase* state_input) const
{
  const ModelState_Default* state =
    static_cast<const ModelState_Default*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Default.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  return state->initialShearModulus;
}
