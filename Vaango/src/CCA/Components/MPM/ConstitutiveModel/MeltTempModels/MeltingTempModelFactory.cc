/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include "MeltingTempModelFactory.h"
#include "BPSMeltTemp.h"
#include "ConstantMeltTemp.h"
#include "LinearMeltTemp.h"
#include "SCGMeltTemp.h"
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Parallel/Parallel.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <string>

using namespace Uintah;

/// Create an instance of a Melting Temperature Model
std::unique_ptr<MeltingTempModel>
MeltingTempModelFactory::create(ProblemSpecP& ps)
{
  ProblemSpecP child = ps->findBlock("melting_temp_model");
  if (!child) {
    proc0cout
      << "**WARNING** Creating default (constant melting temperature) model"
      << std::endl;
    return (std::make_unique<ConstantMeltTemp>());
    // throw ProblemSetupException("MPM::ConstitutiveModel:Cannot find melting
    // temp model.", __FILE__, __LINE__);
  }
  string mat_type;
  if (!child->getAttribute("type", mat_type)) {
    throw ProblemSetupException(
      "MPM::ConstitutiveModel:No type for melting temp model.",
      __FILE__,
      __LINE__);
  }

  if (mat_type == "constant_Tm") {
    return (std::make_unique<ConstantMeltTemp>(child));
  } else if (mat_type == "linear_Tm") {
    return (std::make_unique<LinearMeltTemp>(child));
  } else if (mat_type == "scg_Tm") {
    return (std::make_unique<SCGMeltTemp>(child));
  } else if (mat_type == "bps_Tm") {
    return (std::make_unique<BPSMeltTemp>(child));
  } else {
    proc0cout
      << "**WARNING** Creating default (constant melting temperature) model"
      << std::endl;
    return (std::make_unique<ConstantMeltTemp>(child));
    // throw ProblemSetupException("MPM::ConstitutiveModel:Unknown Melting Temp
    // Model ("+mat_type+")",
    //                            __FILE__, __LINE__);
  }
}

std::unique_ptr<MeltingTempModel>
MeltingTempModelFactory::createCopy(const MeltingTempModel* mtm)
{
  if (dynamic_cast<const ConstantMeltTemp*>(mtm)) {
    return (std::make_unique<ConstantMeltTemp>(
      dynamic_cast<const ConstantMeltTemp*>(mtm)));
  } else if (dynamic_cast<const LinearMeltTemp*>(mtm)) {
    return (std::make_unique<LinearMeltTemp>(
      dynamic_cast<const LinearMeltTemp*>(mtm)));
  } else if (dynamic_cast<const SCGMeltTemp*>(mtm)) {
    return (
      std::make_unique<SCGMeltTemp>(dynamic_cast<const SCGMeltTemp*>(mtm)));
  } else if (dynamic_cast<const BPSMeltTemp*>(mtm)) {
    return (
      std::make_unique<BPSMeltTemp>(dynamic_cast<const BPSMeltTemp*>(mtm)));
  } else {
    proc0cout << "**WARNING** Creating copy of default (constant melting "
                 "temperature) model"
              << std::endl;
    return (std::make_unique<ConstantMeltTemp>(
      dynamic_cast<const ConstantMeltTemp*>(mtm)));
    // throw ProblemSetupException("Cannot create copy of unknown melting temp
    // model", __FILE__, __LINE__);
  }
}
