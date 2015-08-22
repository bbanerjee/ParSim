/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015 Parresia Research Limited, New Zealand
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


#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_CamClay.h>
using namespace Vaango;

ModelState_CamClay::ModelState_CamClay()
  : ModelState_Default()
{
  p_c = 0.0;
  p_c0 = 0.0;
  p = 0.0;
  q = 0.0;
  epse_v = 0.0;
  epse_s = 0.0;
  epse_v_tr = 0.0;
  epse_s_tr = 0.0;
  elasticStrainTensor = Uintah::Matrix3(0.0);
  elasticStrainTensorTrial = Uintah::Matrix3(0.0);
}

ModelState_CamClay::ModelState_CamClay(const ModelState_CamClay& state)
{
  p_c = state.p_c;
  p_c0 = state.p_c0;
  p = state.p;
  q = state.q;
  epse_v = state.epse_v;
  epse_s = state.epse_s;
  epse_v_tr = state.epse_v_tr;
  epse_s_tr = state.epse_s_tr;
  elasticStrainTensor = state.elasticStrainTensor;
  elasticStrainTensorTrial = state.elasticStrainTensorTrial;
}

ModelState_CamClay::ModelState_CamClay(const ModelState_CamClay* state)
{
  p_c = state->p_c;
  p_c0 = state->p_c0;
  p = state->p;
  q = state->q;
  epse_v = state->epse_v;
  epse_s = state->epse_s;
  epse_v_tr = state->epse_v_tr;
  epse_s_tr = state->epse_s_tr;
  elasticStrainTensor = state->elasticStrainTensor;
  elasticStrainTensorTrial = state->elasticStrainTensorTrial;
}

ModelState_CamClay::~ModelState_CamClay()
{
}

ModelState_CamClay&
ModelState_CamClay::operator=(const ModelState_CamClay& state)
{
  if (this == &state) return *this;
  p_c = state.p_c;
  p_c0 = state.p_c0;
  p = state.p;
  q = state.q;
  epse_v = state.epse_v;
  epse_s = state.epse_s;
  epse_v_tr = state.epse_v_tr;
  epse_s_tr = state.epse_s_tr;
  elasticStrainTensor = state.elasticStrainTensor;
  elasticStrainTensorTrial = state.elasticStrainTensorTrial;
  return *this;
}

ModelState_CamClay*
ModelState_CamClay::operator=(const ModelState_CamClay* state)
{
  if (this == state) return this;
  p_c = state->p_c;
  p_c0 = state->p_c0;
  p = state->p;
  q = state->q;
  epse_v = state->epse_v;
  epse_s = state->epse_s;
  epse_v_tr = state->epse_v_tr;
  epse_s_tr = state->epse_s_tr;
  elasticStrainTensor = state->elasticStrainTensor;
  elasticStrainTensorTrial = state->elasticStrainTensorTrial;
  return this;
}
