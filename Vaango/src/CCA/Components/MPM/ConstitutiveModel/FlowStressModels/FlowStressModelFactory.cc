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

#include "FlowStressModelFactory.h"
#include "LinearHardeningFlow.h"
#include "JohnsonCookFlow.h"
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

FlowStressModel*
FlowStressModelFactory::create(ProblemSpecP& ps)
{
  ProblemSpecP child = ps->findBlock("flow_model");
  if (!child)
    throw ProblemSetupException("Cannot find flow_model tag", __FILE__,
                                __LINE__);
  string mat_type;
  if (!child->getAttribute("type", mat_type))
    throw ProblemSetupException("No type for flow_model", __FILE__, __LINE__);
  if (mat_type == "linear")
    return (scinew LinearHardeningFlow(child));
  else if (mat_type == "johnson_cook")
    return (scinew JohnsonCookFlow(child));
  else if (mat_type == "zerilli_armstrong")
    return (scinew ZAFlow(child));
  else if (mat_type == "zerilli_armstrong_polymer")
    return (scinew ZAPolymerFlow(child));
  else if (mat_type == "mechanical_threshold_stress")
    return (scinew MTSFlow(child));
  else if (mat_type == "steinberg_cochran_guinan")
    return (scinew SCGFlow(child));
  else if (mat_type == "preston_tonks_wallace")
    return (scinew PTWFlow(child));
  else {
    throw ProblemSetupException("Unknown flow Model (" + mat_type + ")",
                                __FILE__, __LINE__);
  }
}

FlowStressModel*
FlowStressModelFactory::createCopy(const FlowStressModel* pm)
{
  if (dynamic_cast<const LinearHardeningFlow*>(pm))
    return (scinew LinearHardeningFlow(dynamic_cast<const LinearHardeningFlow*>(pm)));

  else if (dynamic_cast<const JohnsonCookFlow*>(pm))
    return (scinew JohnsonCookFlow(dynamic_cast<const JohnsonCookFlow*>(pm)));

  else if (dynamic_cast<const ZAFlow*>(pm))
    return (scinew ZAFlow(dynamic_cast<const ZAFlow*>(pm)));

  else if (dynamic_cast<const ZAPolymerFlow*>(pm))
    return (scinew ZAPolymerFlow(dynamic_cast<const ZAPolymerFlow*>(pm)));

  else if (dynamic_cast<const MTSFlow*>(pm))
    return (scinew MTSFlow(dynamic_cast<const MTSFlow*>(pm)));

  else if (dynamic_cast<const SCGFlow*>(pm))
    return (scinew SCGFlow(dynamic_cast<const SCGFlow*>(pm)));

  else if (dynamic_cast<const PTWFlow*>(pm))
    return (scinew PTWFlow(dynamic_cast<const PTWFlow*>(pm)));

  else {
    throw ProblemSetupException("Cannot create copy of unknown flow model",
                                __FILE__, __LINE__);
  }
}
