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

#include <CCA/Components/Examples/Poisson2.h>

#include <CCA/Components/Examples/ExamplesLabel.h>

#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Scheduler.h>
#include <Core/Grid/BoundaryConditions/BCDataArray.h>
#include <Core/Grid/BoundaryConditions/BoundCond.h>
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

Poisson2::Poisson2(const ProcessorGroup* myworld,
                   const MaterialManagerP& mat_manager)
  : SimulationCommon(myworld, mat_manager)
{
  d_phi_label =
    VarLabel::create("phi", NCVariable<double>::getTypeDescription());
  d_residual_label =
    VarLabel::create("residual", sum_vartype::getTypeDescription());
}

Poisson2::~Poisson2()
{
  VarLabel::destroy(d_phi_label);
  VarLabel::destroy(d_residual_label);
}

void
Poisson2::problemSetup(const ProblemSpecP& params,
                       [[maybe_unused]] const ProblemSpecP& restart_prob_spec,
                       [[maybe_unused]] GridP& grid,
                       [[maybe_unused]] const std::string& input_ups_dir)
{
  ProblemSpecP poisson = params->findBlock("Poisson");
  poisson->require("delt", d_delt);
  poisson->require("maxresidual", d_maxresidual);
  d_mymat = std::make_shared<EmptyMaterial>();
  d_materialManager->registerEmptyMaterial(d_mymat);
}

void
Poisson2::scheduleInitialize(const LevelP& level, SchedulerP& sched)
{
  Task* task = scinew Task("initialize", this, &Poisson2::initialize);
  task->computes(d_phi_label);
  sched->addTask(task, level->eachPatch(), d_materialManager->allMaterials());
}

void
Poisson2::scheduleComputeStableTimestep(const LevelP& level, SchedulerP& sched)
{
  Task* task = scinew Task(
    "computeStableTimestep", this, &Poisson2::computeStableTimestep);
  task->computes(getDelTLabel(), level.get_rep());
  sched->addTask(task, level->eachPatch(), d_materialManager->allMaterials());
}

void
Poisson2::scheduleTimeAdvance(const LevelP& level, SchedulerP& sched)
{
  Task* task = scinew Task(
    "timeAdvance", this, &Poisson2::timeAdvance, level, sched.get_rep());
  task->hasSubScheduler();
  task->needs(Task::OldDW, d_phi_label, Ghost::AroundNodes, 1);
  task->computes(d_phi_label);
  LoadBalancer* lb                = sched->getLoadBalancer();
  const PatchSet* perproc_patches = lb->getPerProcessorPatchSet(level);
  sched->addTask(task, perproc_patches, d_materialManager->allMaterials());
}

void
Poisson2::computeStableTimestep(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset*,
                                DataWarehouse*,
                                DataWarehouse* new_dw)
{
  new_dw->put(delt_vartype(d_delt), getDelTLabel(), getLevel(patches));
}

void
Poisson2::initialize(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset* matls,
                     DataWarehouse*,
                     DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    for (int m = 0; m < matls->size(); m++) {
      int matl = matls->get(m);
      NCVariable<double> phi;
      new_dw->allocateAndPut(phi, d_phi_label, matl, patch);
      phi.initialize(0);

      for (Patch::FaceType face = Patch::startFace; face <= Patch::endFace;
           face                 = Patch::nextFace(face)) {

        if (patch->getBCType(face) == Patch::None) {
          int numChildren =
            patch->getBCDataArray(face)->getNumberChildren(matl);
          for (int child = 0; child < numChildren; child++) {
            Iterator nbound_ptr, nu;

            BoundCondBaseSP bcb =
              patch->getArrayBCValues(face, matl, "Phi", nu, nbound_ptr, child);

            BoundCond<double>::BoundCondP bc =
              std::dynamic_pointer_cast<BoundCond<double>>(bcb);
            double value = bc->getValue();
            for (nbound_ptr.reset(); !nbound_ptr.done(); nbound_ptr++) {
              phi[*nbound_ptr] = value;
            }
          }
        }
      }
#if 0        
      if(patch->getBCType(Patch::xminus) != Patch::Neighbor){
	IntVector l,h;
	patch->getFaceNodes(Patch::xminus, 0, l, h);
	for(NodeIterator iter(l,h); !iter.done(); iter++)
	  phi[*iter]=1;
      }
#endif
    }
  }
}

void
Poisson2::timeAdvance(const ProcessorGroup* pg,
                      const PatchSubset* patches,
                      const MaterialSubset* matls,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw,
                      LevelP level,
                      Scheduler* sched)
{
  SchedulerP subsched = sched->createSubScheduler();
  subsched->initialize();
  GridP grid = level->getGrid();

  // Create the tasks
  Task* task = scinew Task("iterate", this, &Poisson2::iterate);
  task->needs(Task::OldDW, d_phi_label, Ghost::AroundNodes, 1);
  task->computes(d_phi_label);
  task->computes(d_residual_label);
  subsched->addTask(
    task, level->eachPatch(), d_materialManager->allMaterials());

  // Compile the scheduler
  subsched->advanceDataWarehouse(grid);
  subsched->compile();

  int count = 0;
  double residual;
  subsched->get_dw(1)->transferFrom(old_dw, d_phi_label, patches, matls);
  // Iterate
  do {
    subsched->advanceDataWarehouse(grid);
    subsched->get_dw(0)->setScrubbing(DataWarehouse::ScrubComplete);
    subsched->get_dw(1)->setScrubbing(DataWarehouse::ScrubNonPermanent);
    subsched->execute();

    sum_vartype residual_var;
    subsched->get_dw(1)->get(residual_var, d_residual_label);
    residual = residual_var;

    if (pg->myRank() == 0) {
      std::cerr << "Iteration " << count++ << ", residual=" << residual << '\n';
    }
  } while (residual > d_maxresidual);

  new_dw->transferFrom(subsched->get_dw(1), d_phi_label, patches, matls);
}

void
Poisson2::iterate(const ProcessorGroup*,
                  const PatchSubset* patches,
                  const MaterialSubset* matls,
                  DataWarehouse* old_dw,
                  DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    for (int m = 0; m < matls->size(); m++) {
      int matl = matls->get(m);
      constNCVariable<double> phi;
      old_dw->get(phi, d_phi_label, matl, patch, Ghost::AroundNodes, 1);
      NCVariable<double> newphi;
      new_dw->allocateAndPut(newphi, d_phi_label, matl, patch);
      newphi.copyPatch(phi, newphi.getLow(), newphi.getHigh());
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
      for (NodeIterator iter(l, h); !iter.done(); iter++) {
        newphi[*iter] =
          (1. / 6) *
          (phi[*iter + IntVector(1, 0, 0)] + phi[*iter + IntVector(-1, 0, 0)] +
           phi[*iter + IntVector(0, 1, 0)] + phi[*iter + IntVector(0, -1, 0)] +
           phi[*iter + IntVector(0, 0, 1)] + phi[*iter + IntVector(0, 0, -1)]);
        double diff = newphi[*iter] - phi[*iter];
        residual += diff * diff;
      }
#if 0
      foreach boundary that exists on this patch {
	count face boundaries;
	switch(kind of bc for boundary){
	case Dirichlet:
	  set the value accordingly;
	  break;
	case Neumann:
	  Do the different derivative;
	  break;
	}
      }
      ASSERT(numFaceBoundaries == patch->getNumBoundaryFaces());
#endif
      new_dw->put(sum_vartype(residual), d_residual_label);
    }
  }
}
