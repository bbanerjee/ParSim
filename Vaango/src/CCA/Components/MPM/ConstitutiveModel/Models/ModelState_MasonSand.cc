/*
 * The MIT License
 *
 * Copyright (c) 2015-2016 Parresia Research Limited, New Zealand
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


#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_MasonSand.h>
#include <Core/Exceptions/InternalError.h>
#include <iostream>

using namespace Vaango;

const Uintah::Matrix3 ModelState_MasonSand::Identity(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
const double ModelState_MasonSand::sqrtTwo = std::sqrt(2.0);
const double ModelState_MasonSand::sqrtThree = std::sqrt(3.0);

ModelState_MasonSand::ModelState_MasonSand()
  : ModelState_Default()
{
  bulkModulus = 0.0;
  shearModulus = 0.0;

  capX = 0.0;
  kappa = 0.0;
  zeta = 0.0;

  stressTensor = Uintah::Matrix3(0.0);
  deviatoricStressTensor = Uintah::Matrix3(0.0);
  I1 = 0.0;
  J2 = 0.0;
  sqrt_J2 = 0.0;
  rr = 0.0;
  zz = 0.0;

  plasticStrainTensor = Uintah::Matrix3(0.0);
  ep_v = 0.0;
  dep_v = 0.0;
  ev_0 = 0.0;
  p3 = 0.0;

  phi0 = 0.0;
  Sw0 = 0.0;

  porosity = 0.0;
  saturation = 0.0;
}

ModelState_MasonSand::ModelState_MasonSand(const ModelState_MasonSand& state)
{
  bulkModulus = state.bulkModulus;
  shearModulus = state.shearModulus;

  capX = state.capX;
  kappa = state.kappa;
  zeta = state.zeta;

  stressTensor = state.stressTensor;
  deviatoricStressTensor = state.deviatoricStressTensor;
  I1 = state.I1;
  J2 = state.J2;
  sqrt_J2 = state.sqrt_J2;
  rr = state.rr;
  zz = state.zz;

  plasticStrainTensor = state.plasticStrainTensor;
  ep_v = state.ep_v;
  dep_v = state.dep_v;
  ev_0 = state.ev_0;
  p3 = state.p3;

  phi0 = state.phi0;
  Sw0 = state.Sw0;
  porosity = state.porosity;
  saturation = state.saturation;

  yieldParams = state.yieldParams;
}

ModelState_MasonSand::ModelState_MasonSand(const ModelState_MasonSand* state)
{
  bulkModulus = state->bulkModulus;
  shearModulus = state->shearModulus;

  capX = state->capX;
  kappa = state->kappa;
  zeta = state->zeta;

  stressTensor = state->stressTensor;
  deviatoricStressTensor = state->deviatoricStressTensor;
  I1 = state->I1;
  J2= state->J2;
  sqrt_J2 = state->sqrt_J2;
  rr = state->rr;
  zz = state->zz;

  plasticStrainTensor = state->plasticStrainTensor;
  ep_v = state->ep_v;
  dep_v = state->dep_v;
  ev_0 = state->ev_0;
  p3 = state->p3;

  phi0 = state->phi0;
  Sw0 = state->Sw0;
  porosity = state->porosity;
  saturation = state->saturation;

  yieldParams = state->yieldParams;
}

ModelState_MasonSand::~ModelState_MasonSand()
{
}

ModelState_MasonSand&
ModelState_MasonSand::operator=(const ModelState_MasonSand& state)
{
  if (this == &state) return *this;

  bulkModulus = state.bulkModulus;
  shearModulus = state.shearModulus;

  capX = state.capX;
  kappa = state.kappa;
  zeta = state.zeta;

  stressTensor = state.stressTensor;
  deviatoricStressTensor = state.deviatoricStressTensor;
  I1 = state.I1;
  J2 = state.J2;
  sqrt_J2 = state.sqrt_J2;
  rr = state.rr;
  zz = state.zz;

  plasticStrainTensor = state.plasticStrainTensor;
  ep_v = state.ep_v;
  dep_v = state.dep_v;
  ev_0 = state.ev_0;
  p3 = state.p3;

  phi0 = state.phi0;
  Sw0 = state.Sw0;
  porosity = state.porosity;
  saturation = state.saturation;

  yieldParams = state.yieldParams;

  return *this;
}

ModelState_MasonSand*
ModelState_MasonSand::operator=(const ModelState_MasonSand* state)
{
  if (this == state) return this;

  bulkModulus = state->bulkModulus;
  shearModulus = state->shearModulus;

  capX = state->capX;
  kappa = state->kappa;
  zeta = state->zeta;

  stressTensor = state->stressTensor;
  deviatoricStressTensor = state->deviatoricStressTensor;
  I1 = state->I1;
  J2 = state->J2;
  sqrt_J2 = state->sqrt_J2;
  rr = state->rr;
  zz = state->zz;

  plasticStrainTensor = state->plasticStrainTensor;
  ep_v = state->ep_v;
  dep_v = state->dep_v;
  ev_0 = state->ev_0;
  p3 = state->p3;

  phi0 = state->phi0;
  Sw0 = state->Sw0;
  porosity = state->porosity;
  saturation = state->saturation;

  yieldParams = state->yieldParams;

  return this;
}

void 
ModelState_MasonSand::updateStressInvariants()
{
  // Compute the first invariant
  I1 = stressTensor.Trace();  //Pa

  // Compute the deviatoric part of the tensor
  deviatoricStressTensor = stressTensor - Identity*(I1/3.0);  //Pa

  // Compute the second invariant
  J2 = 0.5*deviatoricStressTensor.Contract(deviatoricStressTensor);  //Pa^2
  J2 = (J2 < 1e-16*(I1*I1+J2)) ? 0.0 : J2;
  sqrt_J2 = std::sqrt(J2);

  // Compute the Lode coordinates (r, z)
  rr = sqrtTwo*sqrt_J2;
  zz = I1/sqrtThree;
}

void 
ModelState_MasonSand::updateVolumetricPlasticStrain()
{
  // Compute volumetric strain
  ep_v = plasticStrainTensor.Trace();
}
