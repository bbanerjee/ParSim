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

#include <CCA/Components/Examples/ParticleTest1.h>

#include <CCA/Components/Examples/ExamplesLabel.h>

#include <CCA/Ports/Scheduler.h>
#include <Core/Grid/EmptyMaterial.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>

using namespace Uintah;

ParticleTest1::ParticleTest1(const ProcessorGroup* myworld,
                             const MaterialManagerP& mat_manager)
  : SimulationCommon(myworld, mat_manager)
{
  d_labels = scinew ExamplesLabel();
}

ParticleTest1::~ParticleTest1()
{
  delete d_labels;
}

void
ParticleTest1::problemSetup(const ProblemSpecP& params,
                            const ProblemSpecP& restart_prob_spec,
                            GridP& grid,
                            const std::string& input_ups_dir)
{
  dynamic_cast<Scheduler*>(getPort("scheduler"))
    ->setPositionVar(d_labels->pXLabel);
  ProblemSpecP pt1 = params->findBlock("ParticleTest1");
  pt1->getWithDefault("doOutput", d_doOutput, 0);
  pt1->getWithDefault("doGhostCells", d_doGhostCells, 0);

  d_mymat = std::make_shared<EmptyMaterial>();
  d_materialManager->registerEmptyMaterial(d_mymat);
}

void
ParticleTest1::scheduleInitialize(const LevelP& level, SchedulerP& sched)
{
  Task* task = scinew Task("initialize", this, &ParticleTest1::initialize);
  task->computes(d_labels->pXLabel);
  task->computes(d_labels->pMassLabel);
  task->computes(d_labels->pParticleIDLabel);
  sched->addTask(task, level->eachPatch(), d_materialManager->allMaterials());
}

void
ParticleTest1::scheduleComputeStableTimestep(const LevelP& level,
                                             SchedulerP& sched)
{
  Task* task = scinew Task(
    "computeStableTimestep", this, &ParticleTest1::computeStableTimestep);
  task->computes(getDelTLabel(), level.get_rep());
  sched->addTask(task, level->eachPatch(), d_materialManager->allMaterials());
}

void
ParticleTest1::scheduleTimeAdvance(const LevelP& level, SchedulerP& sched)
{
  const MaterialSet* matls = d_materialManager->allMaterials();

  Task* task = scinew Task("timeAdvance", this, &ParticleTest1::timeAdvance);

  // set this in problemSetup.  0 is no ghost cells, 1 is all with 1 ghost
  // atound-node, and 2 mixes them
  if (d_doGhostCells == 0) {
    task->requires(Task::OldDW, d_labels->pParticleIDLabel, Ghost::None, 0);
    task->requires(Task::OldDW, d_labels->pXLabel, Ghost::None, 0);
    task->requires(Task::OldDW, d_labels->pMassLabel, Ghost::None, 0);
  }

  else if (d_doGhostCells == 1) {
    task->requires(Task::OldDW, d_labels->pXLabel, Ghost::AroundNodes, 1);
    task->requires(Task::OldDW, d_labels->pMassLabel, Ghost::AroundNodes, 1);
    task->requires(
      Task::OldDW, d_labels->pParticleIDLabel, Ghost::AroundNodes, 1);
  } else if (d_doGhostCells == 2) {
    task->requires(Task::OldDW, d_labels->pXLabel, Ghost::None, 0);
    task->requires(Task::OldDW, d_labels->pMassLabel, Ghost::AroundNodes, 1);
    task->requires(Task::OldDW, d_labels->pParticleIDLabel, Ghost::None, 0);
  }

  task->computes(d_labels->pXLabel_preReloc);
  task->computes(d_labels->pMassLabel_preReloc);
  task->computes(d_labels->pParticleIDLabel_preReloc);
  sched->addTask(task, level->eachPatch(), d_materialManager->allMaterials());

  d_labels->d_particleState.clear();
  d_labels->d_particleState_preReloc.clear();
  for (int m = 0; m < matls->size(); m++) {
    std::vector<const VarLabel*> vars;
    std::vector<const VarLabel*> vars_preReloc;

    vars.push_back(d_labels->pMassLabel);
    vars.push_back(d_labels->pParticleIDLabel);

    vars_preReloc.push_back(d_labels->pMassLabel_preReloc);
    vars_preReloc.push_back(d_labels->pParticleIDLabel_preReloc);
    d_labels->d_particleState.push_back(vars);
    d_labels->d_particleState_preReloc.push_back(vars_preReloc);
  }

  sched->scheduleParticleRelocation(level,
                                    d_labels->pXLabel_preReloc,
                                    d_labels->d_particleState_preReloc,
                                    d_labels->pXLabel,
                                    d_labels->d_particleState,
                                    d_labels->pParticleIDLabel,
                                    matls);
}

void
ParticleTest1::computeStableTimestep(const ProcessorGroup* /*pg*/,
                                     const PatchSubset* patches,
                                     const MaterialSubset* /*matls*/,
                                     DataWarehouse*,
                                     DataWarehouse* new_dw)
{
  new_dw->put(delt_vartype(1), getDelTLabel(), getLevel(patches));
}

void
ParticleTest1::initialize(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse* /*old_dw*/,
                          DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    Point low          = patch->cellPosition(patch->getCellLowIndex());
    Point high         = patch->cellPosition(patch->getCellHighIndex());
    for (int m = 0; m < matls->size(); m++) {
      srand(1);
      int numParticles = 10;
      int matl         = matls->get(m);

      ParticleVariable<Point> px;
      ParticleVariable<double> pMass;
      ParticleVariable<long64> pids;

      ParticleSubset* subset =
        new_dw->createParticleSubset(numParticles, matl, patch);
      new_dw->allocateAndPut(px, d_labels->pXLabel, subset);
      new_dw->allocateAndPut(pMass, d_labels->pMassLabel, subset);
      new_dw->allocateAndPut(pids, d_labels->pParticleIDLabel, subset);

      for (int i = 0; i < numParticles; i++) {
        Point pos(
          (((float)rand()) / RAND_MAX * (high.x() - low.x() - 1) + low.x()),
          (((float)rand()) / RAND_MAX * (high.y() - low.y() - 1) + low.y()),
          (((float)rand()) / RAND_MAX * (high.z() - low.z() - 1) + low.z()));
        px[i]    = pos;
        pids[i]  = patch->getID() * numParticles + i;
        pMass[i] = ((float)rand()) / RAND_MAX * 10;
      }
    }
  }
}

void
ParticleTest1::timeAdvance(const ProcessorGroup*,
                           const PatchSubset* patches,
                           const MaterialSubset* matls,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    for (int m = 0; m < matls->size(); m++) {
      int matl               = matls->get(m);
      ParticleSubset* pset   = old_dw->getParticleSubset(matl, patch);
      ParticleSubset* delset = scinew ParticleSubset(0, matl, patch);

      // Get the arrays of particle values to be changed
      constParticleVariable<Point> px;
      ParticleVariable<Point> pxnew;
      constParticleVariable<double> pMass;
      ParticleVariable<double> pMassnew;
      constParticleVariable<long64> pids;
      ParticleVariable<long64> pidsnew;

      old_dw->get(pMass, d_labels->pMassLabel, pset);
      old_dw->get(px, d_labels->pXLabel, pset);
      old_dw->get(pids, d_labels->pParticleIDLabel, pset);

      new_dw->allocateAndPut(pMassnew, d_labels->pMassLabel_preReloc, pset);
      new_dw->allocateAndPut(pxnew, d_labels->pXLabel_preReloc, pset);
      new_dw->allocateAndPut(
        pidsnew, d_labels->pParticleIDLabel_preReloc, pset);

      // every timestep, move down the +x axis, and decay the mass a little bit
      for (auto i = 0u; i < pset->numParticles(); i++) {
        Point pos(px[i].x() + .25, px[i].y(), px[i].z());
        pxnew[i]    = pos;
        pidsnew[i]  = pids[i];
        pMassnew[i] = pMass[i] * .9;
        if (d_doOutput) {
          std::cout << " Patch " << patch->getID() << ": ID " << pidsnew[i]
                    << ", pos " << pxnew[i] << ", mass " << pMassnew[i]
                    << std::endl;
        }
      }
      new_dw->deleteParticles(delset);
    }
  }
}
