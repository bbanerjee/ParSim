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
}
namespace Vaango {

using constParticleVarPoint = Uintah::constParticleVariable<Uintah::Point>;
using constParticleVarVector = Uintah::constParticleVariable<Uintah::Vector>;
using constParticleVarMatrix3 = Uintah::constParticleVariable<Uintah::Matrix3>;
using constParticleVarLong64 = Uintah::constParticleVariable<Uintah::long64>;
using constParticleVarDouble = Uintah::constParticleVariable<double>;
using ParticleVarVector = Uintah::ParticleVariable<Uintah::Vector>;
using ParticleVarLong64 = Uintah::ParticleVariable<Uintah::long64>;
using ParticleIDToCurrentIdxMap = std::map<Uintah::long64, int>;

// DEM force calculation data structures
struct DEMParticleInputData {
  explicit DEMParticleInputData(int num_mats)
    : pX_old(num_mats), pX0_old(num_mats), 
      pMass_old(num_mats), pRadius_old(num_mats), pSize_old(num_mats),
      pOrientation_old(num_mats), pVelocity_old(num_mats), pAngVel_old(num_mats),
      pRigidBodyID_old(num_mats) 
  {}

  std::vector<constParticleVarPoint>   pX_old, pX0_old;
  std::vector<constParticleVarDouble>  pMass_old, pRadius_old;
  std::vector<constParticleVarMatrix3> pSize_old, pOrientation_old;
  std::vector<constParticleVarVector>  pVelocity_old, pAngVel_old;
  std::vector<constParticleVarLong64>  pRigidBodyID_old;
};

struct DEMParticleSets {
  explicit DEMParticleSets(int num_mats)
    : psets_all(num_mats, nullptr), psets_real(num_mats, nullptr)
  {}

  std::vector<Uintah::ParticleSubset*>  psets_all;
  std::vector<Uintah::ParticleSubset*>  psets_real;
};

struct DEMParticleOutputData {
  explicit DEMParticleOutputData(int num_mats)
    : pExtForce_new(num_mats), pTorque_new(num_mats), pRigidBodyID_new(num_mats)
  {}

  std::vector<ParticleVarVector> pExtForce_new;
  std::vector<ParticleVarVector> pTorque_new;
  std::vector<ParticleVarLong64> pRigidBodyID_new;
};

struct DEMContactResult {
  Vector totalForce   { 0, 0, 0 };
  Vector arm_i        { 0, 0, 0 };
  Vector arm_j_center { 0, 0, 0 };
  bool   collision    { false };
};

// Contact property helpers 
struct DEMContactProps {
  double kn, kt, gamma, mu;
};

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
  scheduleInitialize(Uintah::SchedulerP& sched,
                     const Uintah::PatchSet* patches,
                     const Uintah::MaterialSet* matls);

  void
  actuallyInitialize(const Uintah::Patch* patch,
                     const Uintah::MPMMaterial* mpm_matl,
                     Uintah::DataWarehouse* new_dw);

  void
  computeStableTimestep(const Uintah::Patch* patch,
                        const Uintah::MPMMaterial* mpm_matl,
                        Uintah::DataWarehouse* old_dw,
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
  ParticleIDToCurrentIdxMap
  buildMasterParticleMap(int numMatls,
                         const DEMParticleSets& particle_sets,
                         const DEMParticleInputData& input_data) const;

  // Unique body selection
  ParticleIDToCurrentIdxMap
  buildUniqueBodies(int m_i,
                    Uintah::particleIndex idx_i,
                    int m_j,
                    const DEMParticleSets& particle_sets,
                    const DEMParticleInputData& input_data,
                    bool j_is_dem_material) const;

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
    Uintah::long64 rbID_j,
    const DEMContactProps& cp,
    const DEMParticleInputData& pd,
    const Uintah::MPMMaterial* matl_j,
    const ParticleIDToCurrentIdxMap& master_particles) const;

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
                     Uintah::long64 rbID_j,
                     int m_j,
                     const DEMParticleInputData& inputs,
                     const Uintah::MPMMaterial* matl_j,
                     const Uintah::Patch* patch,
                     const ParticleIDToCurrentIdxMap& master_particles,
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
