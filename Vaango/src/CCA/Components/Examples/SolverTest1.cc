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

#include <CCA/Components/Examples/SolverTest1.h>

#include <CCA/Components/Examples/ExamplesLabel.h>

#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Scheduler.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/EmptyMaterial.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/SoleVariable.h>
#include <Core/Grid/Variables/Stencil7.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>

using namespace Uintah;

SolverTest1::SolverTest1(const ProcessorGroup* myworld,
                         const MaterialManagerP& mat_manager)
  : SimulationCommon(myworld, mat_manager)
{
  d_labels = std::make_unique<ExamplesLabel>();
}
//__________________________________
//
SolverTest1::~SolverTest1()
{
  delete d_solver_parameters;
}
//__________________________________
//
void
SolverTest1::problemSetup(const ProblemSpecP& prob_spec,
                          const ProblemSpecP& restart_prob_spec,
                          GridP& grid,
                          const std::string& input_ups_dir)
{
  d_solver = dynamic_cast<SolverInterface*>(getPort("solver"));
  if (!d_solver) {
    throw InternalError("ST1:couldn't get solver port", __FILE__, __LINE__);
  }

  ProblemSpecP st_ps = prob_spec->findBlock("SolverTest");
  d_solver->readParameters(st_ps, "implicitPressure");
  d_solver->getParameters()->setSolveOnExtraCells(false);

  st_ps->require("delt", d_delT);

  // whether or not to do laplacian in x,y,or z direction
  if (st_ps->findBlock("X_Laplacian")) {
    x_laplacian = true;
  } else {
    x_laplacian = false;
  }
  if (st_ps->findBlock("Y_Laplacian")) {
    y_laplacian = true;
  } else {
    y_laplacian = false;
  }
  if (st_ps->findBlock("Z_Laplacian")) {
    z_laplacian = true;
  } else {
    z_laplacian = false;
  }

  if (!x_laplacian && !y_laplacian && !z_laplacian) {
    throw ProblemSetupException("SolverTest: Must specify one of X_Laplacian, "
                                "Y_Laplacian, or Z_Laplacian",
                                __FILE__,
                                __LINE__);
  }
  d_mymat = std::make_shared<EmptyMaterial>();
  d_materialManager->registerEmptyMaterial(d_mymat);
}
//__________________________________
//
void
SolverTest1::scheduleInitialize(const LevelP& level, SchedulerP& sched)
{
  d_solver->scheduleInitialize(level, sched, d_materialManager->allMaterials());
}
//__________________________________
//
void
SolverTest1::scheduleComputeStableTimestep(const LevelP& level,
                                           SchedulerP& sched)
{
  Task* task = scinew Task(
    "computeStableTimestep", this, &SolverTest1::computeStableTimestep);
  task->computes(getDelTLabel(), level.get_rep());
  sched->addTask(task, level->eachPatch(), d_materialManager->allMaterials());
}
//__________________________________
//
void
SolverTest1::scheduleTimeAdvance(const LevelP& level, SchedulerP& sched)
{
  Task* task = scinew Task(
    "timeAdvance", this, &SolverTest1::timeAdvance, level, sched.get_rep());
  task->computes(d_labels->pressure_matrix);
  task->computes(d_labels->pressure_rhs);

  sched->addTask(task, level->eachPatch(), d_materialManager->allMaterials());

  d_solver->scheduleSolve(level,
                          sched,
                          d_materialManager->allMaterials(),
                          d_labels->pressure_matrix,
                          Task::NewDW,
                          d_labels->pressure,
                          false,
                          d_labels->pressure_rhs,
                          Task::NewDW,
                          0,
                          Task::OldDW,
                          false);
}
//__________________________________
//
void
SolverTest1::computeStableTimestep(const ProcessorGroup*,
                                   const PatchSubset* pss,
                                   const MaterialSubset*,
                                   DataWarehouse*,
                                   DataWarehouse* new_dw)
{
  new_dw->put(delt_vartype(d_delT), getDelTLabel(), getLevel(pss));
}
//__________________________________
//
void
SolverTest1::initialize(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset* matls,
                        DataWarehouse*,
                        DataWarehouse* new_dw)
{
}
//______________________________________________________________________
//
void
SolverTest1::timeAdvance(const ProcessorGroup* pg,
                         const PatchSubset* patches,
                         const MaterialSubset* matls,
                         DataWarehouse* old_dw,
                         DataWarehouse* new_dw,
                         LevelP level,
                         Scheduler* sched)
{
  int center = 0;
  int n = 0, s = 0, e = 0, w = 0, t = 0, b = 0;

  if (x_laplacian) {
    center += 2;
    e = -1;
    w = -1;
  }
  if (y_laplacian) {
    center += 2;
    n = -1;
    s = -1;
  }
  if (z_laplacian) {
    center += 2;
    t = -1;
    b = -1;
  }

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    for (int m = 0; m < matls->size(); m++) {
      int matl = matls->get(m);

      CCVariable<Stencil7> A;
      CCVariable<double> rhs;
      new_dw->allocateAndPut(A, d_labels->pressure_matrix, matl, patch);
      new_dw->allocateAndPut(rhs, d_labels->pressure_rhs, matl, patch);

      // bool first = true;
      for (CellIterator iter(patch->getExtraCellIterator()); !iter.done();
           iter++) {
        IntVector c     = *iter;
        Stencil7& A_tmp = A[c];
        A_tmp.p         = center;
        A_tmp.n         = n;
        A_tmp.s         = s;
        A_tmp.e         = e;
        A_tmp.w         = w;
        A_tmp.t         = t;
        A_tmp.b         = b;

        if (c == IntVector(0, 0, 0)) {
          rhs[c] = 1.0;
        } else {
          rhs[c] = 0;
        }
      }
    }
  }
}
