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

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Tabular.h>
#include <Core/Exceptions/InternalError.h>
#include <iostream>
#include "ModelState_Tabular.h"

using namespace Vaango;

const Uintah::Matrix3 ModelState_Tabular::Identity(1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                                                 0.0, 0.0, 1.0);
const double ModelState_Tabular::sqrtTwo = std::sqrt(2.0);
const double ModelState_Tabular::sqrtThree = std::sqrt(3.0);

ModelState_Tabular::ModelState_Tabular()
  : ModelStateBase()
  , particleID(0)
  , I1(0)
  , J2(0)
  , sqrt_J2(0)
  , zz(0)
  , rr(0)
  , ep_v(0)
  , ep_eq(0)
  , ep_cum_eq(0)
  , stressTensor(0)
  , deviatoricStressTensor(0)
  , elasticStrainTensor(0)
  , plasticStrainTensor(0)
{
}

ModelState_Tabular&
Vaango::ModelState_Tabular::operator=(const ModelState_Tabular& state)
{
  if (this == &state) {
    return *this;
  }

  // Call base class assignment operator to handle base part
  ModelStateBase::operator=(state);

  // Copy derived class specific members
  this->particleID = state.particleID;
  this->I1 = state.I1;
  this->J2 = state.J2;
  this->sqrt_J2 = state.sqrt_J2;
  this->zz = state.zz;
  this->rr = state.rr;
  this->ep_v = state.ep_v;
  this->ep_eq = state.ep_eq;
  this->ep_cum_eq = state.ep_cum_eq;
  this->stressTensor = state.stressTensor;
  this->deviatoricStressTensor = state.deviatoricStressTensor;
  this->elasticStrainTensor = state.elasticStrainTensor;
  this->plasticStrainTensor = state.plasticStrainTensor;

  return *this;
}

ModelState_Tabular&
Vaango::ModelState_Tabular::operator=(const ModelState_Tabular&& state) noexcept
{
  if (this == &state) {
    return *this;
  }

  // Call base class assignment operator to handle base part
  ModelStateBase::operator=(std::move(state));

  // Move derived class specific members
  this->particleID = std::move(state.particleID);
  this->I1 = std::move(state.I1);
  this->J2 = std::move(state.J2);
  this->sqrt_J2 = std::move(state.sqrt_J2);
  this->zz = std::move(state.zz);
  this->rr = std::move(state.rr);
  this->ep_v = std::move(state.ep_v);
  this->ep_eq = std::move(state.ep_eq);
  this->ep_cum_eq = std::move(state.ep_cum_eq);
  this->stressTensor = std::move(state.stressTensor);
  this->deviatoricStressTensor = std::move(state.deviatoricStressTensor);
  this->elasticStrainTensor = std::move(state.elasticStrainTensor);
  this->plasticStrainTensor = std::move(state.plasticStrainTensor);

  return *this;
}

void
ModelState_Tabular::updateStressInvariants()
{
  // Compute the first invariant of the total stress
  I1 = stressTensor.Trace();

  // Compute the deviatoric part of the total stress tensor
  deviatoricStressTensor = stressTensor - Identity * (I1 / 3.0);

  // Compute the second invariant of the deviatoric total stress
  J2 = 0.5 * deviatoricStressTensor.Contract(deviatoricStressTensor);
  J2 = (J2 < 1e-16 * (I1 * I1 + J2)) ? 0.0 : J2;
  sqrt_J2 = std::sqrt(J2);

  // Compute the Lode coordinates (r, z) of the effective stress
  rr = sqrtTwo * sqrt_J2;
  zz = I1 / sqrtThree;
}

void
ModelState_Tabular::updatePlasticStrainInvariants()
{
  // Compute volumetric strain
  ep_v = plasticStrainTensor.Trace();

  // Compute equivalent plastic strain
  Uintah::Matrix3 devPlasticStrain =
    plasticStrainTensor - Identity * (ep_v / 3.0);
  ep_eq = std::sqrt(2.0 / 3.0 * devPlasticStrain.Contract(devPlasticStrain));
  ep_cum_eq = ep_eq;
}
