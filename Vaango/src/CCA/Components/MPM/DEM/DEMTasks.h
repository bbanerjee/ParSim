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

#include <CCA/Components/MPM/DEM/DEMCommon.h>
#include <CCA/Components/MPM/DEM/DEMRigidBodySpatialIndex.h>

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
class Task;
}
namespace Vaango {

class DEMTasks final
{
public:
  DEMTasks(Uintah::ProblemSpecP& ps,
           const Uintah::MaterialManagerP& mat_manager,
           const Uintah::MPMLabel* mpm_labels,
           const Uintah::MPMFlags* mpm_flags,
           int numGhostParticles,
           int numGhostNodes);

  ~DEMTasks() = default;

  void
  scheduleInitialize(Uintah::Task* task);

  void
  actuallyInitialize(const Uintah::Patch* patch,
                     const Uintah::MPMMaterial* mpm_matl,
                     Uintah::DataWarehouse* new_dw);

  void
  computeStableTimestep(const Uintah::Patch* patch,
                        const Uintah::MPMMaterial* mpm_matl,
                        Uintah::DataWarehouse* new_dw,
                        double& dt_dem);

  void
  scheduleComputeDEMForces(Uintah::SchedulerP& sched,
                           const Uintah::PatchSet* patches,
                           const Uintah::MaterialSet* matls);

  void
  computeDEMForces(const Uintah::ProcessorGroup* pg,
                   const Uintah::PatchSubset* patches,
                   const Uintah::MaterialSubset* matls,
                   Uintah::DataWarehouse* old_dw,
                   Uintah::DataWarehouse* new_dw);

  void
  scheduleIntegrateDEMRotation(Uintah::SchedulerP& sched,
                               const Uintah::PatchSet* patches,
                               const Uintah::MaterialSet* matls);

  void
  integrateDEMRotation(const Uintah::ProcessorGroup* pg,
                       const Uintah::PatchSubset* patches,
                       const Uintah::MaterialSubset* matls,
                       Uintah::DataWarehouse* old_dw,
                       Uintah::DataWarehouse* new_dw);

private:
  const Uintah::MaterialManagerP d_mat_manager;
  const Uintah::MPMLabel* d_mpm_labels;
  const Uintah::MPMFlags* d_mpm_flags;
  int d_num_ghost_particles;
  int d_num_ghost_nodes;

  // Contact properties
  DEMContactProps
  getContactProps(const Uintah::MPMMaterial* matl_i, const Uintah::MPMMaterial* matl_j) const;

  // Load particle data
  void
  getAndAllocateParticleData(const Uintah::Patch* patch,
                             Uintah::DataWarehouse* old_dw,
                             Uintah::DataWarehouse* new_dw,
                             DEMParticleSets& particle_sets,
                             DEMParticleInputData& input_data,
                             DEMParticleOutputData& output_data);

  // Master particle map
  DEMBodyIDToCurrentParticleIdxMap
  buildMasterParticleMap(int numMatls,
                         const DEMParticleSets& particle_sets,
                         const DEMParticleInputData& input_data) const;

  // Unique body selection
  DEMBodyIDToCurrentParticleIdxMap
  buildParticleInteractionMap(int m_i,
                    Uintah::particleIndex idx_i,
                    const DEMParticleInputData& input_data,
                    const DEMRigidBodySpatialIndex& spatial_index_j,
                    const Uintah::Vector& cell_size) const;

  // Per-case contact force routines
  // Case A: rigid SDF body (i) vs. MPM particle (j)
  DEMContactResult
  computeRigidVsParticle(Uintah::particleIndex idx_i,
                         int m_i,
                         Uintah::particleIndex idx_j,
                         int m_j,
                         const DEMContactProps& props,
                         const DEMParticleInputData& inputs,
                         const Uintah::MPMMaterial* matl_i) const;

  // Case B: MPM particle (i) vs. rigid SDF body (j)
  DEMContactResult
  computeParticleVsRigid(
    Uintah::particleIndex idx_i,
    int m_i,
    Uintah::particleIndex idx_j,
    int m_j,
    Uintah::long64 pDEMBodyID_j,
    const DEMContactProps& cp,
    const DEMParticleInputData& pd,
    const Uintah::MPMMaterial* matl_j,
    const DEMBodyIDToCurrentParticleIdxMap& master_particles) const;

  // Case C: rigid sphere (i) vs. rigid sphere (j)
  DEMContactResult
  computeRigidVsRigid(Uintah::particleIndex idx_i,
                      int m_i,
                      Uintah::particleIndex idx_j,
                      int m_j,
                      const DEMContactProps& props,
                      const DEMParticleInputData& inputs) const;

  // Force/torque application
  void
  applyContactForces(const DEMContactResult& contact,
                     Uintah::particleIndex idx_i,
                     int m_i,
                     Uintah::long64 pDEMBodyID_j,
                     int m_j,
                     const DEMParticleInputData& inputs,
                     const Uintah::MPMMaterial* matl_j,
                     const Uintah::Patch* patch,
                     const DEMBodyIDToCurrentParticleIdxMap& master_particles,
                     DEMParticleOutputData& outputs);

  // Disallow copy and move
  DEMTasks(const DEMTasks&) = delete;
  DEMTasks(DEMTasks&&)      = delete;
  DEMTasks&
  operator=(const DEMTasks&) = delete;
  DEMTasks&
  operator=(DEMTasks&&) = delete;
};

} // end namespace Vaango

#endif //__VAANGO_CCA_COMPONENTS_MPM_DEM_TASKS_H__
