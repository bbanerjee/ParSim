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

#include <CCA/Components/Examples/HeatEquation.h>

#include <CCA/Components/Examples/ExamplesLabel.h>

#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Scheduler.h>
#include <Core/Grid/EmptyMaterial.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>

using namespace Uintah;

HeatEquation::HeatEquation(const ProcessorGroup* myworld, const MaterialManagerP& mat_manager)
  : SimulationCommon(myworld, mat_manager)
{
  d_temperature_label =
    VarLabel::create("temperature", NCVariable<double>::getTypeDescription());
  d_residual_label =
    VarLabel::create("residual", sum_vartype::getTypeDescription());
}

HeatEquation::~HeatEquation()
{
  VarLabel::destroy(d_temperature_label);
  VarLabel::destroy(d_residual_label);
}

void
HeatEquation::problemSetup(const ProblemSpecP& params,
                           [[maybe_unused]] const ProblemSpecP& restart_prob_spec,
                           [[maybe_unused]] GridP& grid)
{
  ProblemSpecP heateqn = params->findBlock("HeatEquation");
  heateqn->require("delt", d_delT);
  heateqn->require("maxresidual", d_maxresidual);
  d_mymat = std::make_shared<EmptyMaterial>();
  d_materialManager->registerEmptyMaterial(d_mymat);
}

void
HeatEquation::scheduleInitialize(const LevelP& level, SchedulerP& sched)
{
  Task* task = scinew Task("initialize", this, &HeatEquation::initialize);
  task->computes(d_temperature_label);
  sched->addTask(task, level->eachPatch(), d_materialManager->allMaterials());
}

void
HeatEquation::scheduleComputeStableTimestep(const LevelP& level,
                                            SchedulerP& sched)
{
  Task* task = scinew Task("computeStableTimestep",
                           this,
                           &HeatEquation::computeStableTimestep);
  task->computes(getDelTLabel(), level.get_rep());
  sched->addTask(task, level->eachPatch(), d_materialManager->allMaterials());
}

void
HeatEquation::scheduleTimeAdvance(const LevelP& level, SchedulerP& sched)
{
  Task* task = scinew Task("timeAdvance", this, &HeatEquation::timeAdvance);

  task->needs(Task::OldDW, d_temperature_label, Ghost::AroundNodes, 1);
  task->computes(d_temperature_label);
  task->computes(d_residual_label);
  sched->addTask(task, level->eachPatch(), d_materialManager->allMaterials());
}

void
HeatEquation::computeStableTimestep(const ProcessorGroup*,
                                    const PatchSubset* patches,
                                    const MaterialSubset*,
                                    DataWarehouse*,
                                    DataWarehouse* new_dw)
{
  new_dw->put(delt_vartype(d_delT), getDelTLabel(), getLevel(patches));
}

void
HeatEquation::initialize(const ProcessorGroup*,
                         const PatchSubset* patches,
                         const MaterialSubset* matls,
                         DataWarehouse*,
                         DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    for (int m = 0; m < matls->size(); m++) {
      int matl = matls->get(m);

      NCVariable<double> temperature;
      new_dw->allocateAndPut(temperature, d_temperature_label, matl, patch);
      temperature.initialize(0);

      if (patch->getBCType(Patch::xminus) != Patch::Neighbor) {
        IntVector l, h;
        patch->getFaceNodes(Patch::xminus, 0, l, h);

        for (NodeIterator iter(l, h); !iter.done(); iter++) {
          temperature[*iter] = 1;
        }
      }
    }
  }
}

void
HeatEquation::timeAdvance([[maybe_unused]] const ProcessorGroup* pg,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw)
{

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    for (int m = 0; m < matls->size(); m++) {
      int matl = matls->get(m);
      constNCVariable<double> temperature;

      old_dw->get(temperature,
                  d_temperature_label,
                  matl,
                  patch,
                  Ghost::AroundNodes,
                  1);

      NCVariable<double> newtemperature;

      new_dw->allocateAndPut(newtemperature, d_temperature_label, matl, patch);
      newtemperature.copyPatch(temperature,
                               newtemperature.getLow(),
                               newtemperature.getHigh());

      double residual = 0;
      IntVector l     = patch->getNodeLowIndex();
      IntVector h     = patch->getNodeHighIndex();

      l +=
        IntVector(patch->getBCType(Patch::xminus) == Patch::Neighbor ? 0 : 1,
                  patch->getBCType(Patch::yminus) == Patch::Neighbor ? 0 : 1,
                  patch->getBCType(Patch::zminus) == Patch::Neighbor ? 0 : 1);
      h -= IntVector(patch->getBCType(Patch::xplus) == Patch::Neighbor ? 0 : 1,
                     patch->getBCType(Patch::yplus) == Patch::Neighbor ? 0 : 1,
                     patch->getBCType(Patch::zplus) == Patch::Neighbor ? 0 : 1);

      delt_vartype dt;
      old_dw->get(dt, getDelTLabel());
      Vector dx = patch->getLevel()->dCell();
      Vector diffusion_number(1. / (dx.x() * dx.x()),
                              1. / (dx.y() * dx.y()),
                              1. / (dx.z() * dx.z()));

      double k = .5;

      std::cout << "dx = " << dx << std::endl;
      diffusion_number = diffusion_number * k * dt;
      std::cout << "diffusion_number = " << diffusion_number << std::endl;

      for (NodeIterator iter(l, h); !iter.done(); iter++) {
        newtemperature[*iter] =
          (1. / 6) * (temperature[*iter + IntVector(1, 0, 0)] +
                      temperature[*iter + IntVector(-1, 0, 0)] +
                      temperature[*iter + IntVector(0, 1, 0)] +
                      temperature[*iter + IntVector(0, -1, 0)] +
                      temperature[*iter + IntVector(0, 0, 1)] +
                      temperature[*iter + IntVector(0, 0, -1)]);
        double diff = newtemperature[*iter] - temperature[*iter];
        residual += diff * diff;
      }

      new_dw->put(sum_vartype(residual), d_residual_label);
    }
  }
}
