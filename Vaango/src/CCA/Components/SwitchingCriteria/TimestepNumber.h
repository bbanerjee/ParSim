/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 1997-2012 The University of Utah
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

#ifndef Packages_Uintah_CCA_Components_Switching_TimestepNumber_h
#define Packages_Uintah_CCA_Components_Switching_TimestepNumber_h

#include <CCA/Ports/SwitchingCriteria.h>

#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {

class ProcessorGroup;
class DataWarehouse;

class TimestepNumber : public SwitchingCriteria
{
public:
  // this function has a switch for all known SwitchingCriteria

  TimestepNumber(ProblemSpecP& ps);
  virtual ~TimestepNumber();

  virtual void
  problemSetup(const ProblemSpecP& ps,
               const ProblemSpecP& restart_prob_spec,
               MaterialManagerP& mat_manager);

  virtual void
  scheduleSwitchTest(const LevelP& level, SchedulerP& sched);

  void
  switchTest(const ProcessorGroup*,
             const PatchSubset* patches,
             const MaterialSubset* matls,
             DataWarehouse*,
             DataWarehouse*);

private:
  int d_timestep{0};
  const VarLabel* d_timestep_label{nullptr};
  MaterialManagerP d_mat_manager{nullptr};
};
} // End namespace Uintah

#endif
