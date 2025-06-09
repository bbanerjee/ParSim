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

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Arena.h>
#include <Core/Exceptions/InternalError.h>
#include <iostream>
#include "ModelState_Arena.h"

using namespace Vaango;

const Uintah::Matrix3 ModelState_Arena::Identity(1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                                                 0.0, 0.0, 1.0);
const double ModelState_Arena::sqrtTwo = std::sqrt(2.0);
const double ModelState_Arena::sqrtThree = std::sqrt(3.0);

ModelState_Arena::ModelState_Arena()
  : ModelStateBase()
{
  particleID = 0;

  bulkModulus = 0.0;
  shearModulus = 0.0;

  capX = 0.0;
  kappa = 0.0;
  pbar_w = 0.0;

  stressTensor = Uintah::Matrix3(0.0);
  deviatoricStressTensor = Uintah::Matrix3(0.0);
  I1_eff = 0.0;
  J2 = 0.0;
  sqrt_J2 = 0.0;
  rr = 0.0;
  zz_eff = 0.0;

  elasticStrainTensor = Uintah::Matrix3(0.0);
  plasticStrainTensor = Uintah::Matrix3(0.0);
  ep_v = 0.0;
  dep_v = 0.0;
  ep_cum_eq = 0.0;
  ep_eq = 0.0;

  phi0 = 0.0;
  Sw0 = 0.0;

  porosity = 0.0;
  saturation = 0.0;

  p3 = 0.0;
  t_grow = 1.0e10;
  coherence = 1.0;
}


ModelState_Arena&
ModelState_Arena::operator=(const ModelState_Arena& state)
{
  if (this == &state)
    return *this;

  // Call base class assignment operator to handle base part
  ModelStateBase::operator=(state);

  // Copy derived class specific members
  this->particleID = state.particleID;

  this->bulkModulus = state.bulkModulus;
  this->shearModulus = state.shearModulus;

  this->capX = state.capX;
  this->kappa = state.kappa;
  this->pbar_w = state.pbar_w;

  this->stressTensor = state.stressTensor;
  this->deviatoricStressTensor = state.deviatoricStressTensor;
  this->I1_eff = state.I1_eff;
  this->J2 = state.J2;
  this->sqrt_J2 = state.sqrt_J2;
  this->rr = state.rr;
  this->zz_eff = state.zz_eff;

  this->elasticStrainTensor = state.elasticStrainTensor;
  this->plasticStrainTensor = state.plasticStrainTensor;
  this->ep_v = state.ep_v;
  this->dep_v = state.dep_v;
  this->ep_cum_eq = state.ep_cum_eq;
  this->ep_eq = state.ep_eq;

  this->phi0 = state.phi0;
  this->Sw0 = state.Sw0;
  this->porosity = state.porosity;
  this->saturation = state.saturation;

  this->p3 = state.p3;
  this->t_grow = state.t_grow;
  this->coherence = state.coherence;

  this->yieldParams = state.yieldParams;

  return *this;
}

ModelState_Arena&
Vaango::ModelState_Arena::operator=(const ModelState_Arena&& state) noexcept
{
  if (this == &state)
    return *this;

  // Call base class assignment operator to handle base part
  ModelStateBase::operator=(std::move(state));

  // Move derived class specific members
  this->particleID = std::move(state.particleID);

  this->bulkModulus = std::move(state.bulkModulus);
  this->shearModulus = std::move(state.shearModulus);

  this->capX = std::move(state.capX);
  this->kappa = std::move(state.kappa);
  this->pbar_w = std::move(state.pbar_w);

  this->stressTensor = std::move(state.stressTensor);
  this->deviatoricStressTensor = std::move(state.deviatoricStressTensor);
  this->I1_eff = std::move(state.I1_eff);
  this->J2 = std::move(state.J2);
  this->sqrt_J2 = std::move(state.sqrt_J2);
  this->rr = std::move(state.rr);
  this->zz_eff = std::move(state.zz_eff);

  this->elasticStrainTensor = std::move(state.elasticStrainTensor);
  this->plasticStrainTensor = std::move(state.plasticStrainTensor);
  this->ep_v = std::move(state.ep_v);
  this->dep_v = std::move(state.dep_v);
  this->ep_cum_eq = std::move(state.ep_cum_eq);
  this->ep_eq = std::move(state.ep_eq);

  this->phi0 = std::move(state.phi0);
  this->Sw0 = std::move(state.Sw0);
  this->porosity = std::move(state.porosity);
  this->saturation = std::move(state.saturation);

  this->p3 = std::move(state.p3);
  this->t_grow = std::move(state.t_grow);
  this->coherence = std::move(state.coherence);

  this->yieldParams = std::move(state.yieldParams);

  return *this;
}

void
ModelState_Arena::updateStressInvariants()
{
  // Compute the first invariant of the total stress
  double I1 = stressTensor.Trace(); // Pa

  // Compute the deviatoric part of the total stress tensor
  deviatoricStressTensor = stressTensor - Identity * (I1 / 3.0); // Pa

  // Compute the second invariant of the deviatoric total stress
  J2 = 0.5 * deviatoricStressTensor.Contract(deviatoricStressTensor); // Pa^2
  J2 = (J2 < 1e-16 * (I1 * I1 + J2)) ? 0.0 : J2;
  sqrt_J2 = std::sqrt(J2);

  // Compute I1_eff for partially saturated Arena model
  I1_eff = I1 + (pbar_w * 3.0);

  // Compute the Lode coordinates (r, z) of the effective stress
  rr = sqrtTwo * sqrt_J2;
  zz_eff = I1_eff / sqrtThree;

#ifdef TEST_EFFECT_OF_J2_SIGN
  // Compute the third invariant of the deviatoric total stress
  double J3 = deviatoricStressTensor(0, 0) * deviatoricStressTensor(1, 1) *
                deviatoricStressTensor(2, 2) +
              2.0 * deviatoricStressTensor(0, 1) *
                deviatoricStressTensor(1, 2) * deviatoricStressTensor(2, 0) -
              (deviatoricStressTensor(0, 0) * deviatoricStressTensor(1, 2) *
                 deviatoricStressTensor(1, 2) +
               deviatoricStressTensor(1, 1) * deviatoricStressTensor(2, 0) *
                 deviatoricStressTensor(2, 0) +
               deviatoricStressTensor(2, 2) * deviatoricStressTensor(0, 1) *
                 deviatoricStressTensor(0, 1));

  // Change the sign of sqrtJ2 and rr based on sign of J3
  sqrt_J2 = std::copysign(sqrt_J2, J3);
  rr = std::copysign(rr, J3);
#endif
}

void
ModelState_Arena::updatePlasticStrainInvariants()
{
  // Compute volumetric strain
  ep_v = plasticStrainTensor.Trace();

  // Compute equivalent plastic strain
  Uintah::Matrix3 devPlasticStrain =
    plasticStrainTensor - Identity * (ep_v / 3.0);
  ep_eq = std::sqrt(2.0 / 3.0 * devPlasticStrain.Contract(devPlasticStrain));
}
