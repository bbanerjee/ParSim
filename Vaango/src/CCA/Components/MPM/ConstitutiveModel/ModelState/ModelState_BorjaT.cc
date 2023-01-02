/*
 * The MIT License
 *
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_BorjaT.h>

using namespace Vaango;

ModelState_BorjaT::ModelState_BorjaT()
  : ModelStateT<ModelState_BorjaT>()
{
  p_c = 0.0;
  p_c0 = 0.0;
  epse_v = 0.0;
  epse_s = 0.0;
  epse_v_tr = 0.0;
  epse_s_tr = 0.0;
  elasticStrainTensor = Uintah::Matrix3(0.0);
  elasticStrainTensorTrial = Uintah::Matrix3(0.0);
}

ModelState_BorjaT::ModelState_BorjaT(const ModelState_BorjaT& state)
{
  *this = &state;
}

ModelState_BorjaT::ModelState_BorjaT(const ModelState_BorjaT* state)
{
  *this = state;
}

ModelState_BorjaT::~ModelState_BorjaT() = default;

ModelState_BorjaT&
ModelState_BorjaT::operator=(const ModelState_BorjaT& state)
{
  if (this == &state)
    return *this;

  *this = &state;
  return *this;

}

ModelState_BorjaT*
ModelState_BorjaT::operator=(const ModelState_BorjaT* state)
{
  if (this == state) {
    return this;
  }
  this->copyState(state);
  return this;
}

void 
ModelState_BorjaT::copyLocalState(const ModelState_BorjaT* state)
{
  p_c = state->p_c;
  p_c0 = state->p_c0;
  epse_v = state->epse_v;
  epse_s = state->epse_s;
  epse_v_tr = state->epse_v_tr;
  epse_s_tr = state->epse_s_tr;
  elasticStrainTensor = state->elasticStrainTensor;
  elasticStrainTensorTrial = state->elasticStrainTensorTrial;
}
