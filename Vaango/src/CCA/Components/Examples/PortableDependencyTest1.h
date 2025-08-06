/*
 * The MIT License
 *
 * Copyright (c) 1997-2024 The University of Utah
 * Copyright (c) 2024-2025 Biswajit Banerjee, Parresia Research Limited, NZ
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

#ifndef Packages_Uintah_CCA_Components_Examples_PortableDependencyTest1_h
#define Packages_Uintah_CCA_Components_Examples_PortableDependencyTest1_h

#include <CCA/Components/SimulationCommon/SimulationCommon.h>

#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/VarLabel.h>

#include <CCA/Components/Examples/ExamplesLabel.h>
#include <CCA/Components/Schedulers/DetailedTask.h>
#include <CCA/Components/Schedulers/OnDemandDataWarehouse.h>
#include <CCA/Ports/Scheduler.h>
#include <Core/Grid/BoundaryConditions/BCDataArray.h>
#include <Core/Grid/BoundaryConditions/BoundCond.h>
#include <Core/Grid/EmptyMaterial.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/KokkosViews.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Parallel/Portability.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>

using namespace std;
using namespace Uintah;

namespace Uintah {
class EmptyMaterial;

/**************************************

CLASS
   PortableDependencyTest1

   PortableDependencyTest1 simulation

GENERAL INFORMATION

   PortableDependencyTest1.h

   Steven G. Parker
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)


KEYWORDS
   PortableDependencyTest1

DESCRIPTION
   Long description...

WARNING

 ****************************************/
#define task_parameters                                                        \
  [[maybe_unused]] const PatchSubset *patches,                                 \
  [[maybe_unused]] const MaterialSubset *matls,                                \
  [[maybe_unused]] OnDemandDataWarehouse *old_dw,                              \
  [[maybe_unused]] OnDemandDataWarehouse *new_dw,                              \
  [[maybe_unused]] UintahParams &uintahParams,                                 \
  [[maybe_unused]] ExecutionObject<ExecSpace, MemSpace>&execObj

class PortableDependencyTest1 : public SimulationCommon
{
public:
  PortableDependencyTest1(const ProcessorGroup* myworld,
                          const MaterialManagerP materialManager)
    : SimulationCommon(myworld, materialManager)
  {
    phi_label =
      VarLabel::create("phi", NCVariable<double>::getTypeDescription());
  }

  virtual ~PortableDependencyTest1() { VarLabel::destroy(phi_label); }

  virtual void
  problemSetup(const ProblemSpecP& params,
               const ProblemSpecP& restart_prob_spec,
               GridP& grid,
               const std::string& input_ups_dir = "");

  void
  outputProblemSpec([[maybe_unused]] ProblemSpecP& ps)
  {
  }

  virtual void
  scheduleInitialize(const LevelP& level, SchedulerP& sched)
  {
    Task* task = scinew Task("PortableDependencyTest1::initialize",
                             this,
                             &PortableDependencyTest1::initialize);
    task->computes(phi_label);
    sched->addTask(task, level->eachPatch(), d_materialManager->allMaterials());
  }

  virtual void
  scheduleRestartInitialize([[maybe_unused]] const LevelP& level,
                            [[maybe_unused]] SchedulerP& sched)
  {
  }

  virtual void
  restartInitialize() {};

  virtual void
  scheduleComputeStableTimestep(const LevelP& level, SchedulerP& sched)
  {
    Task* task = scinew Task("PortableDependencyTest1::computeStableTimestep",
                             this,
                             &PortableDependencyTest1::computeStableTimestep);
    task->computes(getDelTLabel(), level.get_rep());
    sched->addTask(task, level->eachPatch(), d_materialManager->allMaterials());
  }

  // main tests start here. scheduleTimeAdvance will schedule tests as per
  // environment variables set.

  virtual void
  scheduleTimeAdvance(const LevelP& level, SchedulerP&);
  template<typename ExecSpace, typename MemSpace>
  void
  scheduleComputeTask(const LevelP& level, SchedulerP& sched);
  template<typename ExecSpace, typename MemSpace>
  void
  scheduleModifyTask(const LevelP& level, SchedulerP& sched);
  template<typename ExecSpace, typename MemSpace>
  void
  scheduleRequireTask(const LevelP& level, SchedulerP& sched);
  // the first compute task - creates phi in the new dw.
  template<typename ExecSpace, typename MemSpace>
  void computeTask(task_parameters);

  template<typename ExecSpace,
           typename MemSpace> // modifies phi after computeTask, if environment
                              // flag is set
                              void modifyTask(task_parameters);

  template<typename ExecSpace,
           typename MemSpace> // requires phi - verify values either after
                              // compute or after modify.
                              void requireTask(task_parameters);

private:
  void
  initialize(const ProcessorGroup*,
             const PatchSubset* patches,
             [[maybe_unused]] const MaterialSubset* matls,
             [[maybe_unused]] DataWarehouse* old_dw,
             DataWarehouse* new_dw)
  {
    for (int p = 0; p < patches->size(); p++) {
      const Patch* patch = patches->get(p);

      NCVariable<double> phi;
      new_dw->allocateAndPut(phi, phi_label, 0, patch);
      phi.initialize(0.);
    }
  }

  void
  computeStableTimestep(const ProcessorGroup*,
                        const PatchSubset* patches,
                        [[maybe_unused]] const MaterialSubset* matls,
                        [[maybe_unused]] DataWarehouse* old_dw,
                        DataWarehouse* new_dw)
  {
    new_dw->put(delt_vartype(delt_), getDelTLabel(), getLevel(patches));
  }

  double delt_;
  std::shared_ptr<EmptyMaterial> mymat_;
  const VarLabel* phi_label;
  std::string tasks, exespaces;

  PortableDependencyTest1(const PortableDependencyTest1&);
  PortableDependencyTest1&
  operator=(const PortableDependencyTest1&);
};
} // namespace Uintah

#endif
