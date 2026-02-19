/*
 * The MIT License
 *
 * Copyright (c) 2013-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#ifndef __VAANGO_CCA_COMPONENTS_MPM_DEM_TASKS_H__
#define __VAANGO_CCA_COMPONENTS_MPM_DEM_TASKS_H__

#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/SchedulerP.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {

class MPMFlags;
class MPMLabel;
class MPMMaterial;

class DEMTasks final
{
public:
  DEMTasks(ProblemSpecP& ps,
           const MaterialManagerP& mat_manager,
           const MPMLabel* mpm_labels,
           const MPMFlags* mpm_flags,
           int numGhostParticles,
           int numGhostNodes);

  ~DEMTasks() = default;

  void
  scheduleInitialize(SchedulerP& sched,
                     const PatchSet* patches,
                     const MaterialSet* matls);

  void
  actuallyInitialize(const Patch* patch,
                     const MPMMaterial* mpm_matl,
                     DataWarehouse* new_dw);

  void
  computeStableTimestep(const Patch* patch,
                        const MPMMaterial* mpm_matl,
                        DataWarehouse* old_dw,
                        DataWarehouse* new_dw,
                        double& dt_dem);

  void
  scheduleComputeDEMForces(SchedulerP& sched,
                           const PatchSet* patches,
                           const MaterialSet* matls);

  void
  computeDEMForces(const ProcessorGroup* pg,
                   const PatchSubset* patches,
                   const MaterialSubset* matls,
                   DataWarehouse* old_dw,
                   DataWarehouse* new_dw);

  void
  scheduleIntegrateDEMRotation(SchedulerP& sched,
                               const PatchSet* patches,
                               const MaterialSet* matls);

  void
  integrateDEMRotation(const ProcessorGroup* pg,
                       const PatchSubset* patches,
                       const MaterialSubset* matls,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw);

private:
  const MaterialManagerP d_mat_manager;
  const MPMLabel* d_mpm_labels;
  const MPMFlags* d_mpm_flags;
  int d_num_ghost_particles;
  int d_num_ghost_nodes;

  // Disallow copy and move
  DEMTasks(const DEMTasks&) = delete;
  DEMTasks(DEMTasks&&)      = delete;
  DEMTasks&
  operator=(const DEMTasks&) = delete;
  DEMTasks&
  operator=(DEMTasks&&) = delete;
};

} // end namespace Uintah

#endif //__VAANGO_CCA_COMPONENTS_MPM_DEM_TASKS_H__
