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

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_TabularCap.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/Constants.h>
#include <Core/Exceptions/InternalError.h>
#include <iostream>

using namespace Vaango;

ModelState_TabularCap::ModelState_TabularCap() 
 : ModelState_Tabular()
 , capX(0)
 , I1_min(0)
 , I1_max(0)
 , sqrtJ2_max(0)
{
}

ModelState_TabularCap::ModelState_TabularCap(const ModelState_TabularCap* state)
{
  *this = *state;
}

ModelState_TabularCap*
ModelState_TabularCap::operator=(const ModelState_TabularCap* state)
{
  if (this == state)
    return this;
  
  *this = *state;

  return this;
}

