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

#include <CCA/Components/SwitchingCriteria/SwitchingCriteriaFactory.h>

#include <CCA/Components/SwitchingCriteria/DDT1Criterion.h>
#include <CCA/Components/SwitchingCriteria/None.h>
#include <CCA/Components/SwitchingCriteria/SimpleBurn.h>
#include <CCA/Components/SwitchingCriteria/SteadyBurn.h>
#include <CCA/Components/SwitchingCriteria/SteadyState.h>
#include <CCA/Components/SwitchingCriteria/TimestepNumber.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>
#include <string>

namespace Uintah {

std::shared_ptr<SwitchingCriteria>
SwitchingCriteriaFactory::create(ProblemSpecP& ps,
                                 [[maybe_unused]] const ProcessorGroup* world)
{
  string criteria("");
  ProblemSpecP switch_ps = ps->findBlock("SwitchCriteria");
  if (switch_ps) {
    std::map<std::string, string> attributes;
    switch_ps->getAttributes(attributes);
    criteria = attributes["type"];
  } else {
    return nullptr;
  }

  if (criteria == "none" || criteria == "None" || criteria == "NONE") {
    return std::make_shared<None>();
  } else if (criteria == "timestep" || criteria == "Timestep" ||
             criteria == "TIMESTEP") {
    return std::make_shared<TimestepNumber>(switch_ps);
  } else if (criteria == "SimpleBurn" || criteria == "Simple_Burn" ||
             criteria == "simpleBurn" || criteria == "simple_Burn") {
    return std::make_shared<SimpleBurnCriteria>(switch_ps);
  } else if (criteria == "SteadyBurn" || criteria == "Steady_Burn" ||
             criteria == "steadyBurn" || criteria == "steady_Burn") {
    return std::make_shared<SteadyBurnCriteria>(switch_ps);
  } else if (criteria == "SteadyState" || criteria == "steadystate") {
    return std::make_shared<SteadyState>(switch_ps);
  } else if (criteria == "DDT1") {
    return std::make_shared<DDT1Criterion>(switch_ps);
  } else {
    std::ostringstream warn;
    warn << "\n ERROR:\n Unknown switching criteria (" << criteria << ")\n";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  return nullptr;
}

} // namespace Uintah