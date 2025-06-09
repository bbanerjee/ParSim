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

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_TabularCap.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/Constants.h>
#include <Core/Exceptions/InternalError.h>
#include <iostream>
#include "ModelState_TabularCap.h"

using namespace Vaango;

ModelState_TabularCap::ModelState_TabularCap() 
 : ModelState_Tabular()
 , capX(0)
 , I1_min(0)
 , I1_max(0)
 , sqrtJ2_max(0)
{
}

ModelState_TabularCap&
Vaango::ModelState_TabularCap::operator=(const ModelState_TabularCap& state)
{
  if (this == &state) {
    return *this;
  }

  // Call base class assignment operator to handle base part
  ModelState_Tabular::operator=(state);

  // Copy derived class specific members
  this->capX        = state.capX;
  this->I1_min      = state.I1_min;
  this->I1_max      = state.I1_max;
  this->sqrtJ2_max  = state.sqrtJ2_max;
  this->closest     = state.closest;
  this->tangent     = state.tangent;
  this->yield_f_pts = state.yield_f_pts;

  return *this;
}

ModelState_TabularCap&
Vaango::ModelState_TabularCap::operator=(
  const ModelState_TabularCap&& state) noexcept
{
  if (this == &state) {
    return *this;
  }

  // Call base class assignment operator to handle base part
  ModelState_Tabular::operator=(std::move(state));

  // Move derived class specific members
  this->capX        = std::move(state.capX);
  this->I1_min      = std::move(state.I1_min);
  this->I1_max      = std::move(state.I1_max);
  this->sqrtJ2_max  = std::move(state.sqrtJ2_max);
  this->closest     = std::move(state.closest);
  this->tangent     = std::move(state.tangent);
  this->yield_f_pts = std::move(state.yield_f_pts);

  return *this;
}
