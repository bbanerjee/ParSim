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

#include "StabilityCheckFactory.h"
#include "AcousticTensorCheck.h"
#include "BeckerCheck.h"
#include "DruckerBeckerCheck.h"
#include "DruckerCheck.h"
#include "NoneCheck.h"
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Parallel/Parallel.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <string>

using namespace Uintah;

/// Create an instance of a stabilty check method
/*! Available checks are : loss of ellipticity of the acoustic tensor */
std::unique_ptr<StabilityCheck>
StabilityCheckFactory::create(ProblemSpecP& ps)
{
  ProblemSpecP child = ps->findBlock("stability_check");
  if (!child) {
    proc0cout << "**WARNING** Creating default action (no stability check)"
              << "\n";
    return (std::make_unique<NoneCheck>());
    throw ProblemSetupException(
      "Cannot find stability check criterion.", __FILE__, __LINE__);
  }

  string mat_type;
  if (!child->getAttribute("type", mat_type)) {
    throw ProblemSetupException(
      "No type for stability check criterion.", __FILE__, __LINE__);
  }

  if (mat_type == "drucker") {
    return (std::make_unique<DruckerCheck>(child));
  } else if (mat_type == "acoustic") {
    return (std::make_unique<AcousticTensorCheck>(child));
  } else if (mat_type == "becker") {
    return (std::make_unique<BeckerCheck>(child));
  } else if (mat_type == "drucker_becker") {
    return (std::make_unique<DruckerBeckerCheck>(child));
  } else if (mat_type == "none") {
    return (std::make_unique<NoneCheck>(child));
  } else {
    proc0cout << "**WARNING** Creating default action (no stability check)"
              << "\n";
    return (std::make_unique<NoneCheck>(child));
    // throw ProblemSetupException("Unknown Stability Check ("+mat_type+")",
    // __FILE__, __LINE__);
  }
}

std::unique_ptr<StabilityCheck>
StabilityCheckFactory::createCopy(const StabilityCheck* sc)
{
  if (dynamic_cast<const DruckerCheck*>(sc)) {
    return (
      std::make_unique<DruckerCheck>(dynamic_cast<const DruckerCheck*>(sc)));
  }

  else if (dynamic_cast<const AcousticTensorCheck*>(sc)) {
    return (std::make_unique<AcousticTensorCheck>(
      dynamic_cast<const AcousticTensorCheck*>(sc)));
  }

  else if (dynamic_cast<const BeckerCheck*>(sc)) {
    return (
      std::make_unique<BeckerCheck>(dynamic_cast<const BeckerCheck*>(sc)));
  }

  else if (dynamic_cast<const DruckerBeckerCheck*>(sc)) {
    return (std::make_unique<DruckerBeckerCheck>(
      dynamic_cast<const DruckerBeckerCheck*>(sc)));
  } else if (dynamic_cast<const NoneCheck*>(sc)) {
    return (std::make_unique<NoneCheck>(dynamic_cast<const NoneCheck*>(sc)));
  }

  else {
    proc0cout
      << "**WARNING** Creating copy of default action (no stability check)"
      << "\n";
    return (std::make_unique<NoneCheck>(dynamic_cast<const NoneCheck*>(sc)));
    //  return 0;
  }
}
