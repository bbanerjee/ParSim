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

#include <CCA/Components/MPM/ThermalContact/ThermalContactFactory.h>

#include <CCA/Components/MPM/ThermalContact/NullThermalContact.h>
#include <CCA/Components/MPM/ThermalContact/STThermalContact.h>
#include <Core/Malloc/Allocator.h>
#include <string>
using std::cerr;

using namespace Uintah;

std::unique_ptr<ThermalContact>
ThermalContactFactory::create(const ProblemSpecP& ps,
                              const MaterialManagerP& mat_manager,
                              const MPMLabel* labels,
                              const MPMFlags* flags)
{
  ProblemSpecP mpm_ps =
    ps->findBlockWithOutAttribute("MaterialProperties")->findBlock("MPM");

  for (ProblemSpecP child = mpm_ps->findBlock("thermal_contact"); child != 0;
       child              = child->findNextBlock("thermal_contact")) {
    return (
      std::make_unique<STThermalContact>(child, mat_manager, labels, flags));
  }

  ProblemSpecP child;
  return (
    std::make_unique<NullThermalContact>(child, mat_manager, labels, flags));
}
