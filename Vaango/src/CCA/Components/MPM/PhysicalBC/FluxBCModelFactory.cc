/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2025 Biswajit Banerjee, Parresia Research Ltd, NZ
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
#include <CCA/Components/MPM/PhysicalBC/FluxBCModelFactory.h>

#include <CCA/Components/MPM/PhysicalBC/AutoCycleFluxBC.h>
#include <CCA/Components/MPM/PhysicalBC/FluxBCModel.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <memory>
#include <string>

using namespace Uintah;

std::unique_ptr<FluxBCModel>
FluxBCModelFactory::create(const MaterialManager* ss,
                           const MPMLabel* labels,
                           const MPMFlags* flags)
{
  if (flags->d_doAutoCycleBC) {
    return std::make_unique<AutoCycleFluxBC>(ss, labels, flags);
  } else {
    return std::make_unique<FluxBCModel>(ss, labels, flags);
  }

  return nullptr;
}
