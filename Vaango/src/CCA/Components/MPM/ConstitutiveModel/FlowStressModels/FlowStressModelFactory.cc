/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#include "FlowStressModelFactory.h"
#include "JohnsonCookFlow.h"
#include "LinearHardeningFlow.h"
#include "MTSFlow.h"
#include "PTWFlow.h"
#include "SCGFlow.h"
#include "ZAFlow.h"
#include "ZAPolymerFlow.h"
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <fstream>
#include <iostream>
#include <string>
using std::cerr;
using std::ifstream;
using std::ofstream;

using namespace Uintah;

std::unique_ptr<FlowStressModel>
FlowStressModelFactory::create(ProblemSpecP& ps)
{
  ProblemSpecP child = ps->findBlock("flow_model");
  if (!child) {
    throw ProblemSetupException(
      "Cannot find flow_model tag", __FILE__, __LINE__);
  }
  string mat_type;
  if (!child->getAttribute("type", mat_type)) {
    throw ProblemSetupException("No type for flow_model", __FILE__, __LINE__);
  }
  if (mat_type == "linear") {
    return (std::make_unique<LinearHardeningFlow>(child));
  } else if (mat_type == "johnson_cook") {
    return (std::make_unique<JohnsonCookFlow>(child));
  } else if (mat_type == "zerilli_armstrong") {
    return (std::make_unique<ZAFlow>(child));
  } else if (mat_type == "zerilli_armstrong_polymer") {
    return (std::make_unique<ZAPolymerFlow>(child));
  } else if (mat_type == "mechanical_threshold_stress") {
    return (std::make_unique<MTSFlow>(child));
  } else if (mat_type == "steinberg_cochran_guinan") {
    return (std::make_unique<SCGFlow>(child));
  } else if (mat_type == "preston_tonks_wallace") {
    return (std::make_unique<PTWFlow>(child));
  } else {
    throw ProblemSetupException(
      "Unknown flow Model (" + mat_type + ")", __FILE__, __LINE__);
  }
}

std::unique_ptr<FlowStressModel>
FlowStressModelFactory::createCopy(const FlowStressModel* pm)
{
  if (dynamic_cast<const LinearHardeningFlow*>(pm)) {
    return (std::make_unique<LinearHardeningFlow>(
      dynamic_cast<const LinearHardeningFlow*>(pm)));
  }

  else if (dynamic_cast<const JohnsonCookFlow*>(pm)) {
    return (std::make_unique<JohnsonCookFlow>(
      dynamic_cast<const JohnsonCookFlow*>(pm)));
  }

  else if (dynamic_cast<const ZAFlow*>(pm)) {
    return (std::make_unique<ZAFlow>(dynamic_cast<const ZAFlow*>(pm)));
  }

  else if (dynamic_cast<const ZAPolymerFlow*>(pm)) {
    return (
      std::make_unique<ZAPolymerFlow>(dynamic_cast<const ZAPolymerFlow*>(pm)));
  }

  else if (dynamic_cast<const MTSFlow*>(pm)) {
    return (std::make_unique<MTSFlow>(dynamic_cast<const MTSFlow*>(pm)));
  }

  else if (dynamic_cast<const SCGFlow*>(pm)) {
    return (std::make_unique<SCGFlow>(dynamic_cast<const SCGFlow*>(pm)));
  }

  else if (dynamic_cast<const PTWFlow*>(pm)) {
    return (std::make_unique<PTWFlow>(dynamic_cast<const PTWFlow*>(pm)));
  }

  else {
    throw ProblemSetupException(
      "Cannot create copy of unknown flow model", __FILE__, __LINE__);
  }
}
