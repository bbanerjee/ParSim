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

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_SoilModelBrannon.h>

using namespace Vaango;

ModelState_SoilModelBrannon::ModelState_SoilModelBrannon()
  : ModelState_Default()
{
  kappa = 0.0;
  CR = 0.0;
  maxX = 0.0;
  eps_v = 0.0;
  delta_eps_v = 0.0;
  scale_eps_v = 0.0;
}

ModelState_SoilModelBrannon::ModelState_SoilModelBrannon(
  const ModelState_SoilModelBrannon& state)
{
  kappa = state.kappa;
  CR = state.CR;
  maxX = state.maxX;
  eps_v = state.eps_v;
  delta_eps_v = state.delta_eps_v;
  scale_eps_v = state.scale_eps_v;
}

ModelState_SoilModelBrannon::ModelState_SoilModelBrannon(
  const ModelState_SoilModelBrannon* state)
{
  kappa = state->kappa;
  CR = state->CR;
  maxX = state->maxX;
  eps_v = state->eps_v;
  delta_eps_v = state->delta_eps_v;
  scale_eps_v = state->scale_eps_v;
}

ModelState_SoilModelBrannon::~ModelState_SoilModelBrannon() = default;

ModelState_SoilModelBrannon&
ModelState_SoilModelBrannon::operator=(const ModelState_SoilModelBrannon& state)
{
  if (this == &state)
    return *this;
  kappa = state.kappa;
  CR = state.CR;
  maxX = state.maxX;
  eps_v = state.eps_v;
  delta_eps_v = state.delta_eps_v;
  scale_eps_v = state.scale_eps_v;
  return *this;
}

ModelState_SoilModelBrannon*
ModelState_SoilModelBrannon::operator=(const ModelState_SoilModelBrannon* state)
{
  if (this == state)
    return this;
  kappa = state->kappa;
  CR = state->CR;
  maxX = state->maxX;
  eps_v = state->eps_v;
  delta_eps_v = state->delta_eps_v;
  scale_eps_v = state->scale_eps_v;
  return this;
}
