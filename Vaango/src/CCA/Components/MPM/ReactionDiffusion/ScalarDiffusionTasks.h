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

#include <Core/Grid/MPMInterpolators/LinearInterpolator.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Parallel/ProcessorGroup.h>

#include <memory>

namespace Uintah {

class MPMFlags;
class MPMLabel;
class MPMBoundCond;
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
  constParticleVariable<double> pExternalScalarFlux;
  constParticleVariable<Vector> pConcentrationGrad;
  constParticleVariable<Matrix3> pStress;

  ParticleVariable<double> pConcentrationNew;
  ParticleVariable<double> pConcPreviousNew;

  constNCVariable<double> gConcentration;
  constNCVariable<double> gConcentrationRate;

  NCVariable<double> gConcentration_new;
  NCVariable<double> gExternalScalarFlux_new;

  NCVariable<double> gConcentrationNoBC_new;
  NCVariable<double> gHydrostaticStress_new;

  NCVariable<double> gConcentration_coarse;
  NCVariable<double> gExternalScalarFlux_coarse;
  NCVariable<double> gConcentrationRate_coarse;

  constNCVariable<double> gConcentration_fine;
  constNCVariable<double> gExternalScalarFlux_fine;
  constNCVariable<double> gConcentrationRate_fine;
};

class ScalarDiffusionTasks final
{
private:
  const MaterialManager* d_mat_manager;
  const MPMFlags* d_mpm_flags;
  const MPMLabel* d_mpm_labels;
  int d_num_ghost_particles;
  int d_num_ghost_nodes;

public:
  std::unique_ptr<SDInterfaceModel> sdInterfaceModel;
  std::unique_ptr<FluxBCModel> fluxBC;

public:
  ScalarDiffusionTasks(ProblemSpecP& ps,
                       const MaterialManagerP& mat_manager,
                       const MPMLabel* mpm_labels,
                       const MPMFlags* mpm_flags,
                       int numGhostParticles,
                       int numGhostNodes);

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
  scheduleInterpolateParticlesToGrid(Task* task);

  void
  scheduleInterpolateParticlesToGrid_CFI(Task* task, int numPaddingCells);

  void
  scheduleCoarsenNodalData_CFI(Task* task);

  void
  scheduleCoarsenNodalData_CFI2(Task* task);

  void
  scheduleNormalizeNodalConc(Task* task);

  void
  scheduleComputeAndIntegrateConcentration(Task* task);

  void
  scheduleComputeConcentrationGradient(SchedulerP& sched,
                                       const PatchSet* patches,
                                       const MaterialSet* matls);

  void
  scheduleComputeConcentrationGradientAMR(const LevelP& level_in,
                                          SchedulerP& sched,
                                          const MaterialSet* matls);

  void
  computeConcentrationGradient(const ProcessorGroup*,
                               const PatchSubset* patches,
                               const MaterialSubset*,
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw);

  void
  scheduleInterpolateToParticlesAndUpdate(Task* task);

  void
  scheduleAddParticles(Task* task);

  void
  scheduleRefine(Task* task);

  void
  getAndAllocateForParticlesToGrid(const Patch* patch,
                                   ParticleSubset* pset,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw,
                                   int matl_dw_index,
                                   ScalarDiffusionTaskData& data);

  void
  interpolateParticlesToGrid(const Patch* patch,
                             const std::vector<IntVector>& ni,
                             const std::vector<double>& S,
                             constParticleVariable<Point>& px,
                             constParticleVariable<double>& pMass,
                             particleIndex idx,
                             ScalarDiffusionTaskData& data);

  void
  getAndAllocateForParticlesToGrid_CFI(const Patch* patch,
                                       ParticleSubset* pset,
                                       DataWarehouse* old_dw,
                                       DataWarehouse* new_dw,
                                       int matl_dw_index,
                                       ScalarDiffusionTaskData& data);
  void
  getForParticlesToGrid_CFI(const Patch* coarsePatch,
                            ParticleSubset* pset,
                            DataWarehouse* old_dw,
                            DataWarehouse* new_dw,
                            int matl_dw_index,
                            ScalarDiffusionTaskData& data);

  void
  getModifiableForParticlesToGrid_CFI(const Patch* finePatch,
                                      DataWarehouse* old_dw,
                                      DataWarehouse* new_dw,
                                      int matl_dw_index,
                                      ScalarDiffusionTaskData& data);

  void
  interpolateParticlesToGrid_CFI(const Patch* patch,
                                 const std::vector<IntVector>& ni,
                                 const std::vector<double>& S,
                                 constParticleVariable<Point>& px,
                                 constParticleVariable<double>& pMass,
                                 particleIndex idx,
                                 ScalarDiffusionTaskData& data);

  void
  getModifiableCoarsenNodalData_CFI(const Patch* coarsePatch,
                                    DataWarehouse* new_dw,
                                    int dwi,
                                    ScalarDiffusionTaskData& data);

  void
  getRegionCoarsenNodalData_CFI(const Level* fineLevel,
                                const Patch* finePatch,
                                DataWarehouse* new_dw,
                                int dwi,
                                ScalarDiffusionTaskData& data);

  void
  coarsenNodalData_CFI(bool zero_flag,
                       ScalarDiffusionTaskData& coarse_data,
                       const ScalarDiffusionTaskData& fine_data,
                       const IntVector& coarse_node,
                       const IntVector& fine_node);

  void
  getModifiableCoarsenNodalData_CFI2(const Patch* coarsePatch,
                                     DataWarehouse* new_dw,
                                     int dwi,
                                     ScalarDiffusionTaskData& data);

  void
  getRegionCoarsenNodalData_CFI2(const Level* fineLevel,
                                 const Patch* finePatch,
                                 DataWarehouse* new_dw,
                                 int dwi,
                                 ScalarDiffusionTaskData& data);

  void
  coarsenNodalData_CFI2(bool zero_flag,
                        ScalarDiffusionTaskData& coarse_data,
                        const ScalarDiffusionTaskData& fine_data,
                        const IntVector& coarse_node,
                        const IntVector& fine_node);

  void
  getAndAllocateForNormalizeNodalConc(const Patch* patch,
                                      DataWarehouse* new_dw,
                                      int dwi,
                                      ScalarDiffusionTaskData& data);

  void
  normalizeNodalConc(const Patch* patch,
                     const NCVariable<double>& gMass,
                     MPMBoundCond& bc,
                     int dwi,
                     ScalarDiffusionTaskData& data);

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
  interpolateFluxBCsCBDI(LinearInterpolator* interpolator,
                         const Patch* patch,
                         ParticleSubset* pset,
                         constParticleVariable<Point>& pX,
                         constParticleVariable<int>& pLoadCurveID,
                         constParticleVariable<Matrix3>& pSize,
                         constParticleVariable<Matrix3>& pDefGrad,
                         ScalarDiffusionTaskData& data);

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

private:
  void
  getForParticlesToGrid(const Patch* patch,
                        ParticleSubset* pset,
                        DataWarehouse* old_dw,
                        DataWarehouse* new_dw,
                        int matl_dw_index,
                        ScalarDiffusionTaskData& data);

  void
  allocateForParticlesToGrid(const Patch* patch,
                             ParticleSubset* pset,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw,
                             int matl_dw_index,
                             ScalarDiffusionTaskData& data);

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
