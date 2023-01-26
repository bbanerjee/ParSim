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

#include <CCA/Components/MPM/PhysicalBC/FluxBCModelFactory.h>
#include <CCA/Components/MPM/ReactionDiffusion/ScalarDiffusionTasks.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>

#include <CCA/Components/MPM/Core/MPMBoundCond.h>
#include <CCA/Components/MPM/Core/MPMDiffusionLabel.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>

#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Ghost.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Task.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Util/DebugStream.h>

namespace Uintah {

static DebugStream cout_doing("MPM_SD", false);

ScalarDiffusionTasks::ScalarDiffusionTasks(ProblemSpecP& ps,
                                           const MaterialManagerP& mat_manager,
                                           const MPMLabel* mpm_labels,
                                           const MPMFlags* mpm_flags)
{
  d_mat_manager = mat_manager.get();
  d_mpm_flags   = mpm_flags;
  d_mpm_labels  = mpm_labels;

  if (d_mpm_flags->d_doScalarDiffusion) {
    sdInterfaceModel = SDInterfaceModelFactory::create(ps,
                                                       mat_manager.get(),
                                                       mpm_flags,
                                                       mpm_labels);
    fluxBC =
      FluxBCModelFactory::create(mat_manager.get(), mpm_labels, mpm_flags);
  }
}

void
ScalarDiffusionTasks::outputProblemSpec(ProblemSpecP& ps)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    sdInterfaceModel->outputProblemSpec(ps);
  }
}

void
ScalarDiffusionTasks::scheduleInitialize(Task* task)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    task->computes(d_mpm_labels->diffusion->pArea);
    task->computes(d_mpm_labels->diffusion->pConcentration);
    task->computes(d_mpm_labels->diffusion->pConcPrevious);
    task->computes(d_mpm_labels->diffusion->pGradConcentration);
    task->computes(d_mpm_labels->diffusion->pExternalScalarFlux);
  }
}

void
ScalarDiffusionTasks::addInitialComputesAndRequires(Task* task,
                                                    const MPMMaterial* mpm_matl,
                                                    const PatchSet* patches)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    ScalarDiffusionModel* sdm =
      static_cast<ScalarDiffusionModel*>(mpm_matl->getScalarDiffusionModel());
    sdm->addInitialComputesAndRequires(task, mpm_matl, patches);
  }
}

void
ScalarDiffusionTasks::actuallyInitialize(const Patch* patch,
                                         const MPMMaterial* mpm_matl,
                                         DataWarehouse* new_dw)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    ScalarDiffusionModel* sdm = mpm_matl->getScalarDiffusionModel();
    sdm->initializeTimeStep(patch, mpm_matl, new_dw);
    sdm->initializeSDMData(patch, mpm_matl, new_dw);
  }
}

void
ScalarDiffusionTasks::actuallyInitializeReductionVars(DataWarehouse* new_dw)
{
  if (d_mpm_flags->d_doAutoCycleBC && d_mpm_flags->d_doScalarDiffusion) {
    if (d_mpm_flags->d_autoCycleUseMinMax) {
      new_dw->put(min_vartype(5.0e11),
                  d_mpm_labels->diffusion->rMinConcentration);
      new_dw->put(max_vartype(-5.0e11),
                  d_mpm_labels->diffusion->rMaxConcentration);
    } else {
      new_dw->put(sum_vartype(0.0),
                  d_mpm_labels->diffusion->rTotalConcentration);
    }
  }
}

void
ScalarDiffusionTasks::scheduleInitializeFluxBCs(const LevelP& level,
                                                SchedulerP& sched)
{
  if (d_mpm_flags->d_useLoadCurves && d_mpm_flags->d_doScalarDiffusion) {
    // Schedule the initialization of scalar fluxBCs per particle
    fluxBC->scheduleInitializeScalarFluxBCs(level, sched);
  }
}

void
ScalarDiffusionTasks::scheduleApplyExternalScalarFlux(SchedulerP& sched,
                                                      const PatchSet* patches,
                                                      const MaterialSet* matls)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    fluxBC->scheduleApplyExternalScalarFlux(sched, patches, matls);
  }
}

void
ScalarDiffusionTasks::scheduleConcInterpolated(SchedulerP& sched,
                                               const PatchSet* patches,
                                               const MaterialSet* matls)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                   getLevel(patches)->getGrid()->numLevels())) {
      return;
    }

    printSchedule(patches, cout_doing, "MPM::scheduleConcInterpolated");

    sdInterfaceModel->addComputesAndRequiresInterpolated(sched, patches, matls);
  }
}

void
ScalarDiffusionTasks::scheduleCompute(SchedulerP& sched,
                                      const PatchSet* patches,
                                      const MaterialSet* matls)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    scheduleComputeFlux(sched, patches, matls);
    scheduleComputeDivergence(sched, patches, matls);
    scheduleDiffusionInterfaceDiv(sched, patches, matls);
  }
}

void
ScalarDiffusionTasks::scheduleComputeAMR(const LevelP& level,
                                         SchedulerP& sched,
                                         const MaterialSet* matls)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    GridP grid    = level->getGrid();
    int maxLevels = grid->numLevels();

    for (int l = 0; l < maxLevels; l++) {
      const LevelP& level     = grid->getLevel(l);
      const PatchSet* patches = level->eachPatch();
      scheduleComputeFlux(sched, patches, matls);
      scheduleComputeDivergence(sched, patches, matls);
      scheduleDiffusionInterfaceDiv(sched, patches, matls);
    }

    for (int l = 0; l < maxLevels; l++) {
      const LevelP& level     = grid->getLevel(l);
      const PatchSet* patches = level->eachPatch();
      scheduleComputeDivergence_CFI(sched, patches, matls);
    }
  }
}

void
ScalarDiffusionTasks::scheduleIntegrate(SchedulerP& sched,
                                        const PatchSet* patches,
                                        const MaterialSet* matls)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                   getLevel(patches)->getGrid()->numLevels())) {
      return;
    }
    printSchedule(patches,
                  cout_doing,
                  "ScalarDiffusionTasks::scheduleComputeAndIntegrateDiffusion");

    Task* t = scinew Task("ScalarDiffusionTasks::computeAndIntegrateDiffusion",
                          this,
                          &ScalarDiffusionTasks::computeAndIntegrateDiffusion);

    t->requires(Task::OldDW, d_mpm_labels->delTLabel);
    t->requires(Task::NewDW, d_mpm_labels->gMassLabel, Ghost::None);
    t->requires(Task::NewDW,
                d_mpm_labels->diffusion->gConcentrationNoBC,
                Ghost::None);
    t->requires(Task::NewDW,
                d_mpm_labels->diffusion->gConcentration,
                Ghost::None);
    t->requires(Task::NewDW,
                d_mpm_labels->diffusion->gExternalScalarFlux,
                Ghost::None);
    t->requires(Task::NewDW,
                sdInterfaceModel->getInterfaceFluxLabel(),
                Ghost::None);
    t->modifies(d_mpm_labels->diffusion->gConcentrationRate);
    t->computes(d_mpm_labels->diffusion->gConcentrationStar);

    sched->addTask(t, patches, matls);
  }
}

void
ScalarDiffusionTasks::computeAndIntegrateDiffusion(const ProcessorGroup*,
                                                   const PatchSubset* patches,
                                                   const MaterialSubset*,
                                                   DataWarehouse* old_dw,
                                                   DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); ++p) {
    const Patch* patch = patches->get(p);
    printTask(patches,
              patch,
              cout_doing,
              "Doing ScalarDiffusionTasks::computeAndIntegrateDiffusion");

    delt_vartype delT;
    old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

    Ghost::GhostType gnone = Ghost::None;
    for (size_t m = 0; m < d_mat_manager->getNumMaterials("MPM"); ++m) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      // Get required variables for this patch
      constNCVariable<double> mass;
      new_dw->get(mass, d_mpm_labels->gMassLabel, dwi, patch, gnone, 0);

      // Scalar Diffusion Related Variables -- JBH
      constNCVariable<double> gSD_IF_FluxRate;
      constNCVariable<double> gConcentration, gConcNoBC, gExternalScalarFlux;
      NCVariable<double> gConcRate, gConcStar;
      const VarLabel* SD_IF_FluxLabel =
        sdInterfaceModel->getInterfaceFluxLabel();
      new_dw->get(gSD_IF_FluxRate, SD_IF_FluxLabel, dwi, patch, gnone, 0);
      new_dw->get(gConcentration,
                  d_mpm_labels->diffusion->gConcentration,
                  dwi,
                  patch,
                  gnone,
                  0);
      new_dw->get(gConcNoBC,
                  d_mpm_labels->diffusion->gConcentrationNoBC,
                  dwi,
                  patch,
                  gnone,
                  0);
      new_dw->get(gExternalScalarFlux,
                  d_mpm_labels->diffusion->gExternalScalarFlux,
                  dwi,
                  patch,
                  gnone,
                  0);

      new_dw->getModifiable(gConcRate,
                            d_mpm_labels->diffusion->gConcentrationRate,
                            dwi,
                            patch);
      new_dw->allocateAndPut(gConcStar,
                             d_mpm_labels->diffusion->gConcentrationStar,
                             dwi,
                             patch);

      // JBH -- Variables associated with scalar diffusion
      for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
           ++iter) {
        IntVector node = *iter;
        //        gConcRate[node] /= mass[node];
        gConcStar[node] = gConcentration[node] + (gConcRate[node] / mass[node] +
                                                  gSD_IF_FluxRate[node]);
      }
      MPMBoundCond bc;
      bc.setBoundaryCondition(patch,
                              dwi,
                              "SD-Type",
                              gConcStar,
                              d_mpm_flags->d_interpolatorType);
      for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
           ++iter) {
        IntVector node  = *iter;
        gConcRate[node] = (gConcStar[node] - gConcNoBC[node]) / delT +
                          gExternalScalarFlux[node] / mass[node];
      }
    }
  }
}

void
ScalarDiffusionTasks::scheduleInterpolateParticlesToGrid(Task* task,
                                                         int numGhostParticles)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    task->requires(Task::OldDW,
                   d_mpm_labels->pStressLabel,
                   Ghost::AroundNodes,
                   numGhostParticles);
    task->requires(Task::OldDW,
                   d_mpm_labels->diffusion->pConcentration,
                   Ghost::AroundNodes,
                   numGhostParticles);
    if (d_mpm_flags->d_GEVelProj) {
      task->requires(Task::OldDW,
                     d_mpm_labels->diffusion->pGradConcentration,
                     Ghost::AroundNodes,
                     numGhostParticles);
    }
    task->requires(Task::NewDW,
                   d_mpm_labels->diffusion->pExternalScalarFlux_preReloc,
                   Ghost::AroundNodes,
                   numGhostParticles);
    task->computes(d_mpm_labels->diffusion->gConcentration);
    task->computes(d_mpm_labels->diffusion->gConcentrationNoBC);
    task->computes(d_mpm_labels->diffusion->gHydrostaticStress);
    task->computes(d_mpm_labels->diffusion->gExternalScalarFlux);

#ifdef CBDI_FLUXBCS
    if (d_mpm_flags->d_useLoadCurves) {
      task->requires(Task::OldDW,
                     d_mpm_labels->pLoadCurveIDLabel,
                     Ghost::AroundNodes,
                     numGhostParticles);
    }
#endif
  }
}

void
ScalarDiffusionTasks::scheduleInterpolateParticlesToGrid_CFI(
  Task* task,
  int numPaddingCells)
{
#define allPatches 0
#define allMatls 0

  if (d_mpm_flags->d_doScalarDiffusion) {
    task->requires(Task::OldDW,
                   d_mpm_labels->diffusion->pConcentration,
                   allPatches,
                   Task::CoarseLevel,
                   allMatls,
                   Task::NormalDomain,
                   Ghost::AroundCells,
                   numPaddingCells);
    task->requires(Task::OldDW,
                   d_mpm_labels->pStressLabel,
                   allPatches,
                   Task::CoarseLevel,
                   allMatls,
                   Task::NormalDomain,
                   Ghost::AroundCells,
                   numPaddingCells);
    task->requires(Task::NewDW,
                   d_mpm_labels->diffusion->pExternalScalarFlux_preReloc,
                   allPatches,
                   Task::CoarseLevel,
                   allMatls,
                   Task::NormalDomain,
                   Ghost::AroundCells,
                   numPaddingCells);

    task->modifies(d_mpm_labels->diffusion->gConcentration);
    task->modifies(d_mpm_labels->diffusion->gHydrostaticStress);
    task->modifies(d_mpm_labels->diffusion->gExternalScalarFlux);
  }
}

void
ScalarDiffusionTasks::scheduleInterpolateToParticlesAndUpdate(Task* task,
                                                              int numGhostNodes)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    task->requires(Task::OldDW,
                   d_mpm_labels->diffusion->pConcentration,
                   Ghost::None);
    task->requires(Task::NewDW,
                   d_mpm_labels->diffusion->gConcentrationRate,
                   Ghost::AroundCells,
                   numGhostNodes);

    task->computes(d_mpm_labels->diffusion->pConcentration_preReloc);
    task->computes(d_mpm_labels->diffusion->pConcPrevious_preReloc);

    if (d_mpm_flags->d_doAutoCycleBC) {
      if (d_mpm_flags->d_autoCycleUseMinMax) {
        task->computes(d_mpm_labels->diffusion->rMinConcentration);
        task->computes(d_mpm_labels->diffusion->rMaxConcentration);
      } else {
        task->computes(d_mpm_labels->diffusion->rTotalConcentration);
      }
    }
  }
}

void
ScalarDiffusionTasks::getAndAllocateForInterpolateToParticles(
  const Patch* patch,
  const MPMMaterial* mpm_matl,
  ParticleSubset* pset,
  DataWarehouse* old_dw,
  DataWarehouse* new_dw,
  int matl_dw_index,
  int num_ghost_particles,
  ScalarDiffusionTaskData& data)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    // Grab min/max concentration and conc. tolerance for particle loop.
    const ScalarDiffusionModel* sdm = mpm_matl->getScalarDiffusionModel();
    data.maxEffectiveConc =
      sdm->getMaxConcentration() - sdm->getConcentrationTolerance();
    data.minEffectiveConc =
      sdm->getMinConcentration() + sdm->getConcentrationTolerance();

    old_dw->get(data.pConcentration,
                d_mpm_labels->diffusion->pConcentration,
                pset);
    new_dw->get(data.gConcentrationRate,
                d_mpm_labels->diffusion->gConcentrationRate,
                matl_dw_index,
                patch,
                Ghost::AroundCells,
                num_ghost_particles);

    new_dw->allocateAndPut(data.pConcentrationNew,
                           d_mpm_labels->diffusion->pConcentration_preReloc,
                           pset);
    new_dw->allocateAndPut(data.pConcPreviousNew,
                           d_mpm_labels->diffusion->pConcPrevious_preReloc,
                           pset);
  }
}

void
ScalarDiffusionTasks::interpolateToParticles(
  double delT,
  particleIndex idx,
  MPMMaterial* mpm_matl,
  const std::vector<IntVector>& ni,
  const std::vector<double>& S,
  ScalarDiffusionTaskData& data,
  ScalarDiffusionGlobalConcData& conc_data)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    size_t NN = ni.size();
    for (size_t k = 0; k < NN; ++k) {
      IntVector node = ni[k];
      conc_data.concRate += data.gConcentrationRate[node] * S[k];
    }

    data.pConcentrationNew[idx] =
      data.pConcentration[idx] + conc_data.concRate * delT;
    if (data.pConcentrationNew[idx] < data.minEffectiveConc) {
      data.pConcentrationNew[idx] = data.minEffectiveConc;
    }
    if (data.pConcentrationNew[idx] > data.maxEffectiveConc) {
      data.pConcentrationNew[idx] = data.maxEffectiveConc;
    }

    data.pConcPreviousNew[idx] = data.pConcentration[idx];
    if (mpm_matl->doConcReduction()) {
      if (d_mpm_flags->d_autoCycleUseMinMax) {
        if (data.pConcentrationNew[idx] > conc_data.maxPatchConc) {
          conc_data.maxPatchConc = data.pConcentrationNew[idx];
        }
        if (data.pConcentrationNew[idx] < conc_data.minPatchConc) {
          conc_data.minPatchConc = data.pConcentrationNew[idx];
        }
      } else {
        conc_data.totalConc += data.pConcentration[idx];
      }
    }
  }
}

void
ScalarDiffusionTasks::scheduleComputeFlux(SchedulerP& sched,
                                          const PatchSet* patches,
                                          const MaterialSet* matls)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                   getLevel(patches)->getGrid()->numLevels())) {
      return;
    }

    printSchedule(patches, cout_doing, "MPM::scheduleComputeFlux");

    Task* t = scinew Task("ScalarDiifusionTasks::computeFlux",
                          this,
                          &ScalarDiffusionTasks::computeFlux);

    auto numMPM = d_mat_manager->getNumMaterials("MPM");
    for (size_t m = 0; m < numMPM; ++m) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      ScalarDiffusionModel* sdm = mpm_matl->getScalarDiffusionModel();
      sdm->scheduleComputeFlux(t, mpm_matl, patches);
    }

    sched->addTask(t, patches, matls);
  }
}

void
ScalarDiffusionTasks::computeFlux(const ProcessorGroup* procs,
                                  const PatchSubset* patches,
                                  const MaterialSubset* matls,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); ++p) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing MPM::computeFlux");

    auto numMatls = d_mat_manager->getNumMaterials("MPM");

    for (size_t m = 0; m < numMatls; ++m) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      ScalarDiffusionModel* sdm = mpm_matl->getScalarDiffusionModel();
      sdm->computeFlux(patch, mpm_matl, old_dw, new_dw);
    }
  }
}

void
ScalarDiffusionTasks::scheduleComputeDivergence(SchedulerP& sched,
                                                const PatchSet* patches,
                                                const MaterialSet* matls)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                   getLevel(patches)->getGrid()->numLevels())) {
      return;
    }

    printSchedule(patches, cout_doing, "MPM::scheduleComputeDivergence");

    Task* t = scinew Task("ScalarDiifusionTasks::computeDivergence",
                          this,
                          &ScalarDiffusionTasks::computeDivergence);

    auto numMPM = d_mat_manager->getNumMaterials("MPM");
    for (size_t m = 0; m < numMPM; ++m) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      ScalarDiffusionModel* sdm = mpm_matl->getScalarDiffusionModel();
      sdm->scheduleComputeDivergence(t, mpm_matl, patches);
    }

    sched->addTask(t, patches, matls);
  }
}

void
ScalarDiffusionTasks::computeDivergence(const ProcessorGroup* procs,
                                        const PatchSubset* patches,
                                        const MaterialSubset* matls,
                                        DataWarehouse* old_dw,
                                        DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); ++p) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing MPM::computeDivergence");

    auto numMatls = d_mat_manager->getNumMaterials("MPM");

    for (size_t m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      ScalarDiffusionModel* sdm = mpm_matl->getScalarDiffusionModel();
      sdm->computeDivergence(patch, mpm_matl, old_dw, new_dw);
    }
  }
}

void
ScalarDiffusionTasks::scheduleDiffusionInterfaceDiv(SchedulerP& sched,
                                                    const PatchSet* patches,
                                                    const MaterialSet* matls)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                   getLevel(patches)->getGrid()->numLevels())) {
      return;
    }

    printSchedule(patches, cout_doing, "MPM::scheduleDiffusionInterfaceDiv");

    sdInterfaceModel->addComputesAndRequiresDivergence(sched, patches, matls);
  }
}

void
ScalarDiffusionTasks::scheduleComputeDivergence_CFI(SchedulerP& sched,
                                                    const PatchSet* patches,
                                                    const MaterialSet* matls)
{
  const Level* fineLevel = getLevel(patches);
  int L_indx             = fineLevel->getIndex();

  if (L_indx > 0) {
    printSchedule(patches, cout_doing, "AMRMPM::scheduleComputeDivergence_CFI");

    Task* t = scinew Task("ScalarDiffusionTasks::computeDivergence_CFI",
                          this,
                          &ScalarDiffusionTasks::computeDivergence_CFI);

    int numMPM = d_mat_manager->getNumMaterials("MPM");
    for (int m = 0; m < numMPM; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      ScalarDiffusionModel* sdm = mpm_matl->getScalarDiffusionModel();
      sdm->scheduleComputeDivergence_CFI(t, mpm_matl, patches);
    }

    sched->addTask(t, patches, matls);
  }
}

void
ScalarDiffusionTasks::computeDivergence_CFI(const ProcessorGroup*,
                                            const PatchSubset* patches,
                                            const MaterialSubset* matls,
                                            DataWarehouse* old_dw,
                                            DataWarehouse* new_dw)
{
  int numMatls = d_mat_manager->getNumMaterials("MPM");

  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
    ScalarDiffusionModel* sdm = mpm_matl->getScalarDiffusionModel();
    sdm->computeDivergence_CFI(patches, mpm_matl, old_dw, new_dw);
  }
}

} // end namespace Uintah
