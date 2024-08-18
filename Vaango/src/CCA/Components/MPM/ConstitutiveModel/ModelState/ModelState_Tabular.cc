/*
 * The MIT License
 *
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

using namespace Vaango;

const Uintah::Matrix3 ModelState_Tabular::Identity(1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                                                 0.0, 0.0, 1.0);
const double ModelState_Tabular::sqrtTwo = std::sqrt(2.0);
const double ModelState_Tabular::sqrtThree = std::sqrt(3.0);


ModelState_Tabular::ModelState_Tabular() 
 : ModelStateBase()
 , particleID(0)
 , I1(0) , J2(0) , sqrt_J2(0) , zz(0) , rr(0)
 , ep_v(0) , ep_eq(0), ep_cum_eq(0)
 , stressTensor(0) , deviatoricStressTensor(0)
 , elasticStrainTensor(0) , plasticStrainTensor(0)
{
}

ModelState_Tabular::ModelState_Tabular(const ModelState_Tabular* state)
{
  *this = *state;
}

ModelState_Tabular*
ModelState_Tabular::operator=(const ModelState_Tabular* state)
{
  if (this == state)
    return this;
  
  *this = *state;

  return this;
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
