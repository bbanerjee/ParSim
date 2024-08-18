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

#include <CCA/Components/SwitchingCriteria/TimestepNumber.h>

#include <CCA/Ports/Scheduler.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Parallel/Parallel.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <iostream>
#include <string>

namespace Uintah {

extern DebugStream switching_dbg;

TimestepNumber::TimestepNumber(ProblemSpecP& ps)
{
  ps->require("timestep", d_timestep);
  proc0cout
    << "Switching criteria: \tTimestep Number: switch components on timestep "
    << d_timestep << std::endl;

  d_timestep_label =
    VarLabel::create(timeStep_name, timeStep_vartype::getTypeDescription());
}

TimestepNumber::~TimestepNumber()
{
  VarLabel::destroy(d_timestep_label);
}

void
TimestepNumber::problemSetup([[maybe_unused]] const ProblemSpecP& ps,
                             [[maybe_unused]] const ProblemSpecP& restart_prob_spec,
                             MaterialManagerP& mat_manager)
{
  d_mat_manager = mat_manager;
}

void
TimestepNumber::scheduleSwitchTest(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level,
                switching_dbg,
                "Switchinng Criteria:TimestepNumber::scheduleSwitchTest");

  Task* t = scinew Task("switchTest", this, &TimestepNumber::switchTest);

  t->requires(Task::OldDW, d_timestep_label);
  t->computes(d_switch_label);
  sched->addTask(t, level->eachPatch(), d_mat_manager->allMaterials());
}

void
TimestepNumber::switchTest([[maybe_unused]] const ProcessorGroup* group,
                           [[maybe_unused]] const PatchSubset* patches,
                           [[maybe_unused]] const MaterialSubset* matls,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw)
{
  switching_dbg << "Doing Switch Criteria:TimestepNumber";

  timeStep_vartype timeStep_var(0);
  if (old_dw->exists(d_timestep_label)) {
    old_dw->get(timeStep_var, d_timestep_label);
  } else if (new_dw->exists(d_timestep_label)) {
    new_dw->get(timeStep_var, d_timestep_label);
  }
  int timeStep = timeStep_var;

  double sw = (timeStep == d_timestep);

  switching_dbg << " is it time to switch components: " << sw << std::endl;

  max_vartype switch_condition(sw);

  const Level* allLevels = nullptr;
  new_dw->put(switch_condition, d_switch_label, allLevels);
}

} // namespace Uintah