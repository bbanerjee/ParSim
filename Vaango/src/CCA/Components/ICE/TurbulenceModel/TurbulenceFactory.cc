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

#include <CCA/Components/ICE/TurbulenceModel/TurbulenceFactory.h>

#include <CCA/Components/ICE/TurbulenceModel/DynamicModel.h>
#include <CCA/Components/ICE/TurbulenceModel/SmagorinskyModel.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>

using namespace Uintah;

TurbulenceFactory::TurbulenceFactory() {}

TurbulenceFactory::~TurbulenceFactory() {}

std::unique_ptr<Turbulence>
TurbulenceFactory::create(ProblemSpecP& ps, MaterialManagerP& materialManager)
{
  ProblemSpecP turb_ps = ps->findBlock("turbulence");

  if (turb_ps) {
    std::string turbulence_model;
    if (!turb_ps->getAttribute("model", turbulence_model)) {
      throw ProblemSetupException(
        "No model for turbulence", __FILE__, __LINE__);
    }
    if (turbulence_model == "Smagorinsky") {
      return std::make_unique<Smagorinsky_Model>(turb_ps, materialManager);
    } else if (turbulence_model == "Germano") {
      return std::make_unique<DynamicModel>(turb_ps, materialManager);
    } else {
      std::ostringstream warn;
      warn << "ERROR ICE: Unknown turbulence model (" << turbulence_model
           << " )\n"
           << "Valid models are:\n"
           << "Smagorinsky\n"
           << "Germano\n"
           << endl;
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }
  }
  return nullptr;
}
