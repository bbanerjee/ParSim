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

#include "ConstantCp.h"
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>

using namespace Uintah;
using Vaango::ModelStateBase;

// Construct a specific heat model.
ConstantCp::ConstantCp()
{
}

ConstantCp::ConstantCp(ProblemSpecP&)
{
}

// Construct a copy of a specific heat model.
ConstantCp::ConstantCp(const ConstantCp*)
{
}

// Destructor of specific heat model.
ConstantCp::~ConstantCp() = default;

void
ConstantCp::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP cm_ps = ps->appendChild("specific_heat_model");
  cm_ps->setAttribute("type", "constant_Cp");
}

// Compute the specific heat
double
ConstantCp::computeSpecificHeat(const ModelStateBase* state)
{
  return state->specificHeat;
}
