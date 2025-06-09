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

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Arenisca3.h>
#include "ModelState_Arenisca3.h"
using namespace Vaango;

ModelState_Arenisca3::ModelState_Arenisca3()
  : ModelStateBase()
  , I1(0.0)
  , sqrt_J2(0.0)
  , plasticStrainTensor(Uintah::Matrix3(0.0)) // The tensor form of plastic strain
  , kappa(0.0)            // The cap kappa parameter
  , capX(0.0)             // The cap hydrostatic compressive strength X
  , zeta(0.0)             // The back stress parameter
{
}


ModelState_Arenisca3&
ModelState_Arenisca3::operator=(const ModelState_Arenisca3& state)
{
  if (this == &state)
    return *this;

  // Call base class assignment operator to handle base part
  ModelStateBase::operator=(state);

  // Copy derived class specific members
  this->I1 = state.I1;
  this->sqrt_J2 = state.sqrt_J2;
  this->plasticStrainTensor = state.plasticStrainTensor;
  this->kappa = state.kappa;
  this->capX = state.capX;
  this->zeta = state.zeta;

  return *this;
}

ModelState_Arenisca3&
Vaango::ModelState_Arenisca3::operator=(
  const ModelState_Arenisca3&& state) noexcept
{
  if (this == &state)
    return *this;

  // Call base class assignment operator to handle base part
  ModelStateBase::operator=(std::move(state));

  // Move derived class specific members
  this->I1 = std::move(state.I1);
  this->sqrt_J2 = std::move(state.sqrt_J2);
  this->plasticStrainTensor = std::move(state.plasticStrainTensor);
  this->kappa = std::move(state.kappa);
  this->capX = std::move(state.capX);
  this->zeta = std::move(state.zeta);

  return *this;
}