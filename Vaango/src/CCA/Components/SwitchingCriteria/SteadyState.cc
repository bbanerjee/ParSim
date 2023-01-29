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

#include <CCA/Components/SwitchingCriteria/SteadyState.h>

#include <CCA/Ports/Scheduler.h>

#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Parallel/Parallel.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>
#include <string>

namespace Uintah {

extern DebugStream switching_dbg;

SteadyState::SteadyState(ProblemSpecP& ps)
{
  ps->require("material", d_material);
  ps->require("num_steps", d_numSteps);

  proc0cout << "material = " << d_material << std::endl;
  proc0cout << "num_steps  = " << d_numSteps << std::endl;

  d_heatRateCCLabel =
    VarLabel::create("heatRate_CC", CCVariable<double>::getTypeDescription());

  d_heatFluxSumLabel =
    VarLabel::create("heatFluxSum", sum_vartype::getTypeDescription());

  d_heatFluxSumTimeDerivativeLabel =
    VarLabel::create("heatFluxSumTimeDerivative",
                     sum_vartype::getTypeDescription());

  // delta t
  VarLabel* nonconstDelT =
    VarLabel::create(delT_name, delt_vartype::getTypeDescription());
  nonconstDelT->isReductionTask(false);
  d_delTLabel = nonconstDelT;
}

SteadyState::~SteadyState()
{
  VarLabel::destroy(d_heatRateCCLabel);
  VarLabel::destroy(d_heatFluxSumLabel);
  VarLabel::destroy(d_heatFluxSumTimeDerivativeLabel);
  VarLabel::destroy(d_delTLabel);
}

void
SteadyState::problemSetup([[maybe_unused]] const ProblemSpecP& ps,
                          [[maybe_unused]] const ProblemSpecP& restart_prob_spec,
                          MaterialManagerP& mat_manager)
{
  d_mat_manager = mat_manager;
}

void
SteadyState::scheduleInitialize(const LevelP& level, SchedulerP& sched)
{

  Task* t = scinew Task("SteadyState::actuallyInitialize",
                        this,
                        &SteadyState::initialize);

  t->computes(d_heatFluxSumLabel);
  t->computes(d_heatFluxSumTimeDerivativeLabel);
  t->computes(d_switch_label);

  sched->addTask(t, level->eachPatch(), d_mat_manager->allMaterials());
}

void
SteadyState::initialize([[maybe_unused]] const ProcessorGroup*,
                        [[maybe_unused]] const PatchSubset* patches,
                        [[maybe_unused]] const MaterialSubset* matls,
                        [[maybe_unused]] DataWarehouse*,
                        DataWarehouse* new_dw)
{
  proc0cout << "Initializing heatFluxSum and heatFluxSumTimeDerivative" << std::endl;

  new_dw->put(max_vartype(0.0), d_heatFluxSumLabel);
  new_dw->put(max_vartype(0.0), d_heatFluxSumTimeDerivativeLabel);
  new_dw->put(max_vartype(0.0), d_switch_label);
}

void
SteadyState::scheduleSwitchTest(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level,
                switching_dbg,
                "Switching Criteria:SteadyState::scheduleSwitchTest");

  Task* t = scinew Task("switchTest", this, &SteadyState::switchTest);

  std::unique_ptr<MaterialSubset> container = std::make_unique<MaterialSubset>();
  container->addReference();
  container->add(d_material);

  t->requires(Task::NewDW, d_heatRateCCLabel, container.get(), Ghost::None);
  t->requires(Task::OldDW, d_heatFluxSumLabel);
  t->requires(Task::OldDW, d_delTLabel);

  t->computes(d_heatFluxSumLabel);
  t->computes(d_heatFluxSumTimeDerivativeLabel);
  t->computes(d_switch_label);

  sched->addTask(t, level->eachPatch(), d_mat_manager->allMaterials());

  scheduleDummy(level, sched);

  container->removeReference();
}

void
SteadyState::switchTest([[maybe_unused]] const ProcessorGroup* group,
                        const PatchSubset* patches,
                        [[maybe_unused]] const MaterialSubset* matls,
                        DataWarehouse* old_dw,
                        DataWarehouse* new_dw)
{
  double sw          = 0;
  double heatFluxSum = 0;

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches,
              patch,
              switching_dbg,
              "Doing Switching Criteria:SimpleBurnCriteria::switchTest");

    constCCVariable<double> heatFlux;
    new_dw->get(heatFlux, d_heatRateCCLabel, 0, patch, Ghost::None, 0);

    for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
      heatFluxSum += heatFlux[*iter];
    }
  }

  new_dw->put(max_vartype(heatFluxSum), d_heatFluxSumLabel);
  proc0cout << "heatFluxSum = " << heatFluxSum << std::endl;

  max_vartype oldHeatFluxSum;
  old_dw->get(oldHeatFluxSum, d_heatFluxSumLabel);
  proc0cout << "oldHeatFluxSum = " << oldHeatFluxSum << std::endl;

  delt_vartype delT;
  old_dw->get(delT, d_delTLabel, getLevel(patches));

  double dH_dt = (heatFluxSum - oldHeatFluxSum) / delT;
  max_vartype heatFluxSumTimeDerivative(dH_dt);
  proc0cout << "heatFluxSumTimeDerivative = " << heatFluxSumTimeDerivative
            << std::endl;

  new_dw->put(heatFluxSumTimeDerivative, d_heatFluxSumTimeDerivativeLabel);

  max_vartype switch_condition(sw);

  const Level* allLevels = 0;
  new_dw->put(switch_condition, d_switch_label, allLevels);
}

void
SteadyState::scheduleDummy(const LevelP& level, SchedulerP& sched)
{
  Task* t = scinew Task("SteadyState::dummy", this, &SteadyState::dummy);
  t->requires(Task::OldDW, d_switch_label, level.get_rep());
  sched->addTask(t, level->eachPatch(), d_mat_manager->allMaterials());
}

void
SteadyState::dummy([[maybe_unused]] const ProcessorGroup* group,
                   [[maybe_unused]] const PatchSubset* patches,
                   [[maybe_unused]] const MaterialSubset* matls,
                   DataWarehouse* old_dw,
                   [[maybe_unused]] DataWarehouse* new_dw)
{
  max_vartype old_sw(1.23);
  old_dw->get(old_sw, d_switch_label);
  proc0cout << "old_sw = " << old_sw << std::endl;
}

} // namespace Uintah