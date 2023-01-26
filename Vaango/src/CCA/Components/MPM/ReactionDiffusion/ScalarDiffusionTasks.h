/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2023 Biswajit Banerjee
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

#ifndef __VAANGO_CCA_COMPONENTS_MPM_SCALAR_DIFFUSION_TASKS_H__
#define __VAANGO_CCA_COMPONENTS_MPM_SCALAR_DIFFUSION_TASKS_H__

#include <CCA/Components/MPM/ReactionDiffusion/DiffusionInterfaces/SDInterfaceModel.h>
#include <CCA/Components/MPM/ReactionDiffusion/DiffusionModels/ScalarDiffusionModel.h>
#include <CCA/Components/MPM/ReactionDiffusion/SDInterfaceModelFactory.h>

#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/SchedulerP.h>

#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Parallel/ProcessorGroup.h>

#include <memory>

namespace Uintah {

class MPMFlags;
class MPMLabel;
class MaterialManager;
class FluxBCModel;

struct ScalarDiffusionGlobalConcData
{
  double concRate;
  double minPatchConc;
  double maxPatchConc;
  double totalConc;
};

struct ScalarDiffusionTaskData
{
  double maxEffectiveConc{ -999.0 };
  double minEffectiveConc{ -999.0 };
  constParticleVariable<double> pConcentration;
  constNCVariable<double> gConcentrationRate;
  ParticleVariable<double> pConcentrationNew;
  ParticleVariable<double> pConcPreviousNew;
};

class ScalarDiffusionTasks final
{
private:
  const MaterialManager* d_mat_manager;
  const MPMFlags* d_mpm_flags;
  const MPMLabel* d_mpm_labels;

public:
  std::unique_ptr<SDInterfaceModel> sdInterfaceModel;
  std::unique_ptr<FluxBCModel> fluxBC;

public:
  ScalarDiffusionTasks(ProblemSpecP& ps,
                       const MaterialManagerP& mat_manager,
                       const MPMLabel* mpm_labels,
                       const MPMFlags* mpm_flags);

  ~ScalarDiffusionTasks() = default;

  void
  outputProblemSpec(ProblemSpecP& ps);

  void
  scheduleInitialize(Task* task);

  void
  addInitialComputesAndRequires(Task* task,
                                const MPMMaterial* matl,
                                const PatchSet* patches);

  void
  actuallyInitialize(const Patch* patch,
                     const MPMMaterial* mpm_matl,
                     DataWarehouse* new_dw);

  void
  actuallyInitializeReductionVars(DataWarehouse* new_dw);

  void
  scheduleInitializeFluxBCs(const LevelP& level, SchedulerP& sched);

  void
  scheduleApplyExternalScalarFlux(SchedulerP& sched,
                                  const PatchSet* patches,
                                  const MaterialSet* matls);

  void
  scheduleConcInterpolated(SchedulerP& sched,
                           const PatchSet* patches,
                           const MaterialSet* matls);

  void
  scheduleCompute(SchedulerP& sched,
                  const PatchSet* patches,
                  const MaterialSet* matls);

  void
  scheduleComputeAMR(const LevelP& level,
                     SchedulerP& sched,
                     const MaterialSet* matls);

  void
  scheduleIntegrate(SchedulerP& sched,
                    const PatchSet* patches,
                    const MaterialSet* matls);

  void
  computeAndIntegrateDiffusion(const ProcessorGroup*,
                               const PatchSubset* patches,
                               const MaterialSubset*,
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw);

  void
  scheduleInterpolateParticlesToGrid(Task* task, int numGhostParticles);

  void
  scheduleInterpolateParticlesToGrid_CFI(Task* task, int numPaddingCells);

  void
  scheduleInterpolateToParticlesAndUpdate(Task* task, int numGhostNodes);

  void
  getAndAllocateForInterpolateToParticles(const Patch* patch,
                                          const MPMMaterial* mpm_matl,
                                          ParticleSubset* pset,
                                          DataWarehouse* old_dw,
                                          DataWarehouse* new_dw,
                                          int matl_dw_index,
                                          int num_ghost_particles,
                                          ScalarDiffusionTaskData& data);

  void
  interpolateToParticles(double delT,
                         particleIndex idx,
                         MPMMaterial* mpm_matl,
                         const std::vector<IntVector>& ni,
                         const std::vector<double>& S,
                         ScalarDiffusionTaskData& data,
                         ScalarDiffusionGlobalConcData& conc_data);

  void
  scheduleComputeFlux(SchedulerP& sched,
                      const PatchSet* patches,
                      const MaterialSet* matls);

  void
  computeFlux(const ProcessorGroup* procs,
              const PatchSubset* patch_subset,
              const MaterialSubset* matl_subset,
              DataWarehouse* old_dw,
              DataWarehouse* new_dw);

  void
  scheduleComputeDivergence(SchedulerP& sched,
                            const PatchSet* patches,
                            const MaterialSet*);

  void
  computeDivergence(const ProcessorGroup* procs,
                    const PatchSubset* patch_subset,
                    const MaterialSubset* matl_subset,
                    DataWarehouse* old_dw,
                    DataWarehouse* new_dw);

  void
  scheduleDiffusionInterfaceDiv(SchedulerP& sched,
                                const PatchSet* patches,
                                const MaterialSet* matls);

  void
  scheduleComputeDivergence_CFI(SchedulerP& sched,
                                const PatchSet* patches,
                                const MaterialSet* matls);
  void
  computeDivergence_CFI(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset* matls,
                        DataWarehouse* old_dw,
                        DataWarehouse* new_dw);

  // Disallow copy and move
  ScalarDiffusionTasks(const ScalarDiffusionTasks&) = delete;
  ScalarDiffusionTasks(ScalarDiffusionTasks&&)      = delete;

  ScalarDiffusionTasks&
  operator=(const ScalarDiffusionTasks&) = delete;
  ScalarDiffusionTasks&
  operator=(ScalarDiffusionTasks&&) = delete;
};

} // end namespace Uintah

#endif //__VAANGO_CCA_COMPONENTS_MPM_SCALAR_DIFFUSION_TASKS_H__
