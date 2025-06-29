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

#include "NoneCheck.h"
#include <Core/ProblemSpec/ProblemSpec.h>
#include <cmath>
#include <vector>

using namespace Uintah;


NoneCheck::NoneCheck()
{
}

NoneCheck::NoneCheck(ProblemSpecP&)
{
}

NoneCheck::NoneCheck(const NoneCheck*)
{
}

NoneCheck::~NoneCheck() = default;

void
NoneCheck::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP stability_ps = ps->appendChild("stability_check");
  stability_ps->setAttribute("type", "none");
}

bool
NoneCheck::checkStability(const Matrix3&, const Matrix3&,
                          const TangentModulusTensor&, Vector&)
{
  return false;
}

bool 
NoneCheck::checkStability([[maybe_unused]] const Matrix3& cauchyStress,
                          [[maybe_unused]] const Matrix3& deformRate,
                          [[maybe_unused]] const Vaango::Tensor::Matrix6Mandel& C_e,
                          [[maybe_unused]] const Vaango::Tensor::Vector6Mandel& P_vec,
                          [[maybe_unused]] const Vaango::Tensor::Vector6Mandel& N_vec,
                          [[maybe_unused]] double H,
                          [[maybe_unused]] Vector& direction)
{
  return false;
}
