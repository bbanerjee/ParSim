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

#include <CCA/Components/SwitchingCriteria/None.h>

#include <CCA/Ports/Scheduler.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/VarTypes.h>

namespace Uintah {

DebugStream switching_dbg("SwitchingCriteria",
                          "SwitchingCriteria",
                          "Switching criteria debug stream",
                          false);

void
None::problemSetup([[maybe_unused]] const ProblemSpecP& ps,
                   [[maybe_unused]] const ProblemSpecP& restart_prob_spec,
                   MaterialManagerP& mat_manager)
{
  d_mat_manager = mat_manager;
}

void
None::scheduleSwitchTest(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level, switching_dbg, "Switching Criteria:None::scheduleSwitchTest");

  Task* t = scinew Task("switchTest", this, &None::switchTest);

  t->computes(d_switch_label);
  sched->addTask(t, level->eachPatch(), d_mat_manager->allMaterials());
}

void
None::switchTest([[maybe_unused]] const ProcessorGroup* group,
                 [[maybe_unused]] const PatchSubset* patches,
                 [[maybe_unused]] const MaterialSubset* matls,
                 [[maybe_unused]] DataWarehouse* old_dw,
                 DataWarehouse* new_dw)
{
  switching_dbg << "  Doing Switching Criteria:None::switchTest" << std::endl;
  double sw = 0;
  max_vartype switch_condition(sw);
  const Level* allLevels = 0;
  new_dw->put(switch_condition, d_switch_label, allLevels);
}

} // namespace Uintah