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

#include <StandAlone/tools/puda/util.h>
#include <Core/Parallel/Parallel.h>


#include <cstdlib>
#include <iostream>

namespace Uintah {
void
findTimestep_loopLimits(const bool tslow_set,
                        const bool tsup_set,
                        const std::vector<double> times,
                        unsigned long& time_step_lower,
                        unsigned long& time_step_upper)
{
  if (!tslow_set) {
    time_step_lower = 0;
  } else if (time_step_lower >= times.size()) {
    std::cerr << "\n";
    std::cerr << "ERROR: 'timesteplow' must be between 0 and "
              << times.size() - 1 << ".  You had " << time_step_lower << ".\n";
    std::cerr << "\n";
    Uintah::Parallel::exitAll(2);
  }
  if (!tsup_set) {
    time_step_upper = times.size() - 1;
  } else if (time_step_upper >= times.size()) {
    std::cerr << "\n";
    std::cerr << "Error: 'timestephigh' must be between 0 and "
              << times.size() - 1 << ".  You had " << time_step_upper << ".\n";
    std::cerr << "\n";
    Uintah::Parallel::exitAll(2);
  }
}

} // namespace Uintah