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

#define allPatches 0
#define allMatls 0

static DebugStream cout_doing("MPM_SD", false);

ScalarDiffusionTasks::ScalarDiffusionTasks(ProblemSpecP& ps,
                                           const MaterialManagerP& mat_manager,
                                           const MPMLabel* mpm_labels,
                                           const MPMFlags* mpm_flags,
                                           int num_ghost_particles,
                                           int num_ghost_nodes)
{
  d_mat_manager         = mat_manager.get();
  d_mpm_flags           = mpm_flags;
  d_mpm_labels          = mpm_labels;
  d_num_ghost_particles = num_ghost_particles;
  d_num_ghost_nodes     = num_ghost_nodes;

  if (d_mpm_flags->d_doScalarDiffusion) {
    sdInterfaceModel = SDInterfaceModelFactory::create(
      ps, mat_manager.get(), mpm_flags, mpm_labels);
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
    t->requires(
      Task::NewDW, d_mpm_labels->diffusion->gConcentrationNoBC, Ghost::None);
    t->requires(
      Task::NewDW, d_mpm_labels->diffusion->gConcentration, Ghost::None);
    t->requires(
      Task::NewDW, d_mpm_labels->diffusion->gExternalScalarFlux, Ghost::None);
    t->requires(
      Task::NewDW, sdInterfaceModel->getInterfaceFluxLabel(), Ghost::None);
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

      new_dw->getModifiable(
        gConcRate, d_mpm_labels->diffusion->gConcentrationRate, dwi, patch);
      new_dw->allocateAndPut(
        gConcStar, d_mpm_labels->diffusion->gConcentrationStar, dwi, patch);

      // JBH -- Variables associated with scalar diffusion
      for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
           ++iter) {
        IntVector node = *iter;
        //        gConcRate[node] /= mass[node];
        gConcStar[node] = gConcentration[node] + (gConcRate[node] / mass[node] +
                                                  gSD_IF_FluxRate[node]);
      }
      MPMBoundCond bc;
      bc.setBoundaryCondition(
        patch, dwi, "SD-Type", gConcStar, d_mpm_flags->d_interpolatorType);
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
ScalarDiffusionTasks::scheduleInterpolateParticlesToGrid(Task* task)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    task->requires(Task::OldDW,
                   d_mpm_labels->pStressLabel,
                   Ghost::AroundNodes,
                   d_num_ghost_particles);
    task->requires(Task::OldDW,
                   d_mpm_labels->diffusion->pConcentration,
                   Ghost::AroundNodes,
                   d_num_ghost_particles);
    if (d_mpm_flags->d_GEVelProj) {
      task->requires(Task::OldDW,
                     d_mpm_labels->diffusion->pGradConcentration,
                     Ghost::AroundNodes,
                     d_num_ghost_particles);
    }
    task->requires(Task::NewDW,
                   d_mpm_labels->diffusion->pExternalScalarFlux_preReloc,
                   Ghost::AroundNodes,
                   d_num_ghost_particles);
    task->computes(d_mpm_labels->diffusion->gConcentration);
    task->computes(d_mpm_labels->diffusion->gConcentrationNoBC);
    task->computes(d_mpm_labels->diffusion->gHydrostaticStress);
    task->computes(d_mpm_labels->diffusion->gExternalScalarFlux);

#ifdef CBDI_FLUXBCS
    if (d_mpm_flags->d_useLoadCurves) {
      task->requires(Task::OldDW,
                     d_mpm_labels->pLoadCurveIDLabel,
                     Ghost::AroundNodes,
                     d_num_ghost_particles);
    }
#endif
  }
}

void
ScalarDiffusionTasks::getAndAllocateForParticlesToGrid(
  const Patch* patch,
  ParticleSubset* pset,
  DataWarehouse* old_dw,
  DataWarehouse* new_dw,
  int matl_dw_index,
  ScalarDiffusionTaskData& data)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    getForParticlesToGrid(patch, pset, old_dw, new_dw, matl_dw_index, data);
    allocateForParticlesToGrid(
      patch, pset, old_dw, new_dw, matl_dw_index, data);
  }
}

void
ScalarDiffusionTasks::getForParticlesToGrid(const Patch* patch,
                                            ParticleSubset* pset,
                                            DataWarehouse* old_dw,
                                            DataWarehouse* new_dw,
                                            int matl_dw_index,
                                            ScalarDiffusionTaskData& data)
{
  new_dw->get(data.pExternalScalarFlux,
              d_mpm_labels->diffusion->pExternalScalarFlux_preReloc,
              pset);
  old_dw->get(
    data.pConcentration, d_mpm_labels->diffusion->pConcentration, pset);
  old_dw->get(data.pStress, d_mpm_labels->pStressLabel, pset);
  if (d_mpm_flags->d_GEVelProj) {
    old_dw->get(data.pConcentrationGrad,
                d_mpm_labels->diffusion->pGradConcentration,
                pset);
  }
}

void
ScalarDiffusionTasks::allocateForParticlesToGrid(const Patch* patch,
                                                 ParticleSubset* pset,
                                                 DataWarehouse* old_dw,
                                                 DataWarehouse* new_dw,
                                                 int matl_dw_index,
                                                 ScalarDiffusionTaskData& data)
{
  new_dw->allocateAndPut(data.gConcentration_new,
                         d_mpm_labels->diffusion->gConcentration,
                         matl_dw_index,
                         patch);
  new_dw->allocateAndPut(data.gConcentrationNoBC_new,
                         d_mpm_labels->diffusion->gConcentrationNoBC,
                         matl_dw_index,
                         patch);
  new_dw->allocateAndPut(data.gHydrostaticStress_new,
                         d_mpm_labels->diffusion->gHydrostaticStress,
                         matl_dw_index,
                         patch);
  new_dw->allocateAndPut(data.gExternalScalarFlux_new,
                         d_mpm_labels->diffusion->gExternalScalarFlux,
                         matl_dw_index,
                         patch);
  data.gConcentration_new.initialize(0);
  data.gConcentrationNoBC_new.initialize(0);
  data.gHydrostaticStress_new.initialize(0);
  data.gExternalScalarFlux_new.initialize(0);
}

void
ScalarDiffusionTasks::interpolateParticlesToGrid(
  const Patch* patch,
  const std::vector<IntVector>& ni,
  const std::vector<double>& S,
  constParticleVariable<Point>& px,
  constParticleVariable<double>& pMass,
  particleIndex idx,
  ScalarDiffusionTaskData& data)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    double one_third    = 1. / 3.;
    double pHydroStress = one_third * data.pStress[idx].Trace();
    double pConc_Ext    = data.pConcentration[idx];
    for (size_t k = 0; k < ni.size(); ++k) {
      auto node = ni[k];
      if (patch->containsNode(node)) {
        if (d_mpm_flags->d_GEVelProj) {
          Point gpos         = patch->getNodePosition(node);
          Vector pointOffset = px[idx] - gpos;
          pConc_Ext -= Dot(data.pConcentrationGrad[idx], pointOffset);
        }
        double massWeight = pMass[idx] * S[k];
        data.gHydrostaticStress_new[node] += pHydroStress * massWeight;
        data.gConcentration_new[node] += pConc_Ext * massWeight;
        data.gExternalScalarFlux_new[node] +=
          data.pExternalScalarFlux[idx] * massWeight;
      }
    }
  }
}

void
ScalarDiffusionTasks::interpolateFluxBCsCBDI(
  LinearInterpolator* interpolator,
  const Patch* patch,
  ParticleSubset* pset,
  constParticleVariable<Point>& pX,
  constParticleVariable<int>& pLoadCurveID,
  constParticleVariable<Matrix3>& pSize,
  constParticleVariable<Matrix3>& pDefGrad,
  ScalarDiffusionTaskData& data)
{
  std::vector<IntVector> ni_LPI(interpolator->size());
  std::vector<double> S_LPI(interpolator->size());

  Vector dx = patch->dCell();
  for (auto& idx : *pset) {

    Point flux_pos;
    if (pLoadCurveID[idx] == 1) {
      flux_pos = Point(pX[idx].x() - 0.5 * pSize[idx](0, 0) * dx.x(),
                       pX[idx].y(),
                       pX[idx].z());
    }
    if (pLoadCurveID[idx] == 2) {
      flux_pos = Point(pX[idx].x() + 0.5 * pSize[idx](0, 0) * dx.x(),
                       pX[idx].y(),
                       pX[idx].z());
    }
    if (pLoadCurveID[idx] == 3) {
      flux_pos = Point(pX[idx].x(),
                       pX[idx].y() + 0.5 * pSize[idx](1, 1) * dx.y(),
                       pX[idx].z());
    }
    interpolator->findCellAndWeights(
      flux_pos, ni_LPI, S_LPI, pSize[idx], pDefGrad[idx]);
    for (size_t k = 0; k < ni_LPI.size(); k++) {
      if (patch->containsNode(ni_LPI[k])) {
        data.gExternalScalarFlux_new[ni_LPI[k]] +=
          data.pExternalScalarFlux[idx] * S_LPI[k];
      }
    }
  }
}

void
ScalarDiffusionTasks::scheduleInterpolateParticlesToGrid_CFI(
  Task* task,
  int numPaddingCells)
{
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
ScalarDiffusionTasks::getForParticlesToGrid_CFI(const Patch* coarsePatch,
                                                ParticleSubset* pset,
                                                DataWarehouse* old_dw,
                                                DataWarehouse* new_dw,
                                                int matl_dw_index,
                                                ScalarDiffusionTaskData& data)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    getForParticlesToGrid(
      coarsePatch, pset, old_dw, new_dw, matl_dw_index, data);
  }
}

void
ScalarDiffusionTasks::getModifiableForParticlesToGrid_CFI(
  const Patch* finePatch,
  DataWarehouse* old_dw,
  DataWarehouse* new_dw,
  int matl_dw_index,
  ScalarDiffusionTaskData& data)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    new_dw->getModifiable(data.gConcentration_new,
                          d_mpm_labels->diffusion->gConcentration,
                          matl_dw_index,
                          finePatch);
    new_dw->getModifiable(data.gHydrostaticStress_new,
                          d_mpm_labels->diffusion->gHydrostaticStress,
                          matl_dw_index,
                          finePatch);
    new_dw->getModifiable(data.gExternalScalarFlux_new,
                          d_mpm_labels->diffusion->gExternalScalarFlux,
                          matl_dw_index,
                          finePatch);
  }
}

void
ScalarDiffusionTasks::interpolateParticlesToGrid_CFI(
  const Patch* patch,
  const std::vector<IntVector>& ni,
  const std::vector<double>& S,
  constParticleVariable<Point>& px,
  constParticleVariable<double>& pMass,
  particleIndex idx,
  ScalarDiffusionTaskData& data)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    double one_third    = 1. / 3.;
    double pConc_x_mass = data.pConcentration[idx] * pMass[idx];
    double pHydroStress_x_mass =
      one_third * data.pStress[idx].Trace() * pMass[idx];
    double pExternalFlux = data.pExternalScalarFlux[idx];

    for (size_t k = 0; k < ni.size(); ++k) {
      auto fineNode = ni[k];
      if (patch->containsNode(fineNode)) {
        data.gConcentration_new[fineNode] += pConc_x_mass * S[k];
        data.gHydrostaticStress_new[fineNode] += pHydroStress_x_mass * S[k];
        data.gExternalScalarFlux_new[fineNode] += pExternalFlux * S[k];
      }
    }
  }
}

void
ScalarDiffusionTasks::scheduleCoarsenNodalData_CFI(Task* task)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    task->requires(Task::NewDW,
                   d_mpm_labels->diffusion->gConcentration,
                   allPatches,
                   Task::FineLevel,
                   allMatls,
                   Task::NormalDomain,
                   Ghost::None,
                   0);
    task->modifies(d_mpm_labels->diffusion->gConcentration);
    task->requires(Task::NewDW,
                   d_mpm_labels->diffusion->gExternalScalarFlux,
                   allPatches,
                   Task::FineLevel,
                   allMatls,
                   Task::NormalDomain,
                   Ghost::None,
                   0);
    task->modifies(d_mpm_labels->diffusion->gExternalScalarFlux);
  }
}

void
ScalarDiffusionTasks::getModifiableCoarsenNodalData_CFI(
  const Patch* coarsePatch,
  DataWarehouse* new_dw,
  int dwi,
  ScalarDiffusionTaskData& data)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    new_dw->getModifiable(data.gConcentration_coarse,
                          d_mpm_labels->diffusion->gConcentration,
                          dwi,
                          coarsePatch);
    new_dw->getModifiable(data.gExternalScalarFlux_coarse,
                          d_mpm_labels->diffusion->gExternalScalarFlux,
                          dwi,
                          coarsePatch);
  }
}

void
ScalarDiffusionTasks::getRegionCoarsenNodalData_CFI(
  const Level* fineLevel,
  const Patch* finePatch,
  DataWarehouse* new_dw,
  int dwi,
  ScalarDiffusionTaskData& data)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    IntVector fl = finePatch->getNodeLowIndex();
    IntVector fh = finePatch->getNodeHighIndex();
    new_dw->getRegion(data.gConcentration_fine,
                      d_mpm_labels->diffusion->gConcentration,
                      dwi,
                      fineLevel,
                      fl,
                      fh);
    new_dw->getRegion(data.gExternalScalarFlux_fine,
                      d_mpm_labels->diffusion->gExternalScalarFlux,
                      dwi,
                      fineLevel,
                      fl,
                      fh);
  }
}

void
ScalarDiffusionTasks::coarsenNodalData_CFI(
  bool zero_flag,
  ScalarDiffusionTaskData& coarse_data,
  const ScalarDiffusionTaskData& fine_data,
  const IntVector& coarse_node,
  const IntVector& fine_node)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    coarse_data.gConcentration_coarse[coarse_node] =
      (zero_flag) ? 0.0 : fine_data.gConcentration_fine[fine_node];
    coarse_data.gExternalScalarFlux_coarse[coarse_node] =
      (zero_flag) ? 0.0 : fine_data.gExternalScalarFlux_fine[fine_node];
  }
}

void
ScalarDiffusionTasks::scheduleCoarsenNodalData_CFI2(Task* task)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    task->requires(Task::NewDW,
                   d_mpm_labels->diffusion->gConcentrationRate,
                   allPatches,
                   Task::FineLevel,
                   allMatls,
                   Task::NormalDomain,
                   Ghost::None,
                   0);
    task->modifies(d_mpm_labels->diffusion->gConcentrationRate);
  }
}

void
ScalarDiffusionTasks::getModifiableCoarsenNodalData_CFI2(
  const Patch* coarsePatch,
  DataWarehouse* new_dw,
  int dwi,
  ScalarDiffusionTaskData& data)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    new_dw->getModifiable(data.gConcentrationRate_coarse,
                          d_mpm_labels->diffusion->gConcentrationRate,
                          dwi,
                          coarsePatch);
  }
}

void
ScalarDiffusionTasks::getRegionCoarsenNodalData_CFI2(
  const Level* fineLevel,
  const Patch* finePatch,
  DataWarehouse* new_dw,
  int dwi,
  ScalarDiffusionTaskData& data)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    IntVector fl = finePatch->getNodeLowIndex();
    IntVector fh = finePatch->getNodeHighIndex();
    new_dw->getRegion(data.gConcentrationRate_fine,
                      d_mpm_labels->diffusion->gConcentrationRate,
                      dwi,
                      fineLevel,
                      fl,
                      fh);
  }
}

void
ScalarDiffusionTasks::coarsenNodalData_CFI2(
  bool zero_flag,
  ScalarDiffusionTaskData& coarse_data,
  const ScalarDiffusionTaskData& fine_data,
  const IntVector& coarse_node,
  const IntVector& fine_node)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    coarse_data.gConcentrationRate_coarse[coarse_node] =
      (zero_flag) ? 0.0 : fine_data.gConcentrationRate_fine[fine_node];
  }
}

void
ScalarDiffusionTasks::scheduleNormalizeNodalConc(Task* task)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    task->modifies(d_mpm_labels->diffusion->gConcentration);
    task->computes(d_mpm_labels->diffusion->gConcentrationNoBC);
    task->modifies(d_mpm_labels->diffusion->gHydrostaticStress);
  }
}

void
ScalarDiffusionTasks::getAndAllocateForNormalizeNodalConc(
  const Patch* patch,
  DataWarehouse* new_dw,
  int dwi,
  ScalarDiffusionTaskData& data)
{
  // NCVariable<double> gConcentration;
  // NCVariable<double> gConcentrationNoBC;
  // NCVariable<double> gHydroStress;
  if (d_mpm_flags->d_doScalarDiffusion) {
    new_dw->getModifiable(data.gConcentration_new,
                          d_mpm_labels->diffusion->gConcentration,
                          dwi,
                          patch,
                          Ghost::None,
                          0);
    new_dw->getModifiable(data.gHydrostaticStress_new,
                          d_mpm_labels->diffusion->gHydrostaticStress,
                          dwi,
                          patch,
                          Ghost::None,
                          0);
    new_dw->allocateAndPut(data.gConcentrationNoBC_new,
                           d_mpm_labels->diffusion->gConcentrationNoBC,
                           dwi,
                           patch);
  }
}

void
ScalarDiffusionTasks::normalizeNodalConc(const Patch* patch,
                                         const NCVariable<double>& gMass,
                                         MPMBoundCond& bc,
                                         int dwi,
                                         ScalarDiffusionTaskData& data)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
         iter++) {
      IntVector n = *iter;
      data.gConcentration_new[n] /= gMass[n];
      data.gHydrostaticStress_new[n] /= gMass[n];
      data.gConcentrationNoBC_new[n] = data.gConcentration_new[n];
    }

    // Apply boundary conditions (if symmetry)
    bc.setBoundaryCondition(patch,
                            dwi,
                            "SD-Type",
                            data.gConcentration_new,
                            d_mpm_flags->d_interpolatorType);
  }
}

void
ScalarDiffusionTasks::scheduleComputeAndIntegrateConcentration(Task* task)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    task->requires(
      Task::NewDW, sdInterfaceModel->getInterfaceFluxLabel(), Ghost::None);
    task->requires(
      Task::NewDW, d_mpm_labels->diffusion->gConcentration, Ghost::None);
    task->requires(
      Task::NewDW, d_mpm_labels->diffusion->gConcentrationNoBC, Ghost::None);
    task->requires(
      Task::NewDW, d_mpm_labels->diffusion->gExternalScalarFlux, Ghost::None);
    task->modifies(d_mpm_labels->diffusion->gConcentrationRate);
    task->computes(d_mpm_labels->diffusion->gConcentrationStar);
  }
}

void
ScalarDiffusionTasks::getAndAllocateForIntegrateConc(
  const Patch* patch,
  DataWarehouse* new_dw,
  int dwi,
  ScalarDiffusionTaskData& data)
{
  // constNCVariable<double> gConcentration, gConcNoBC, gExtScalarFlux;
  if (d_mpm_flags->d_doScalarDiffusion) {
    new_dw->get(data.gFluxRate,
                sdInterfaceModel->getInterfaceFluxLabel(),
                dwi,
                patch,
                Ghost::None,
                0);
    new_dw->get(data.gConcentration,
                d_mpm_labels->diffusion->gConcentration,
                dwi,
                patch,
                Ghost::None,
                0);
    new_dw->get(data.gConcentrationNoBC,
                d_mpm_labels->diffusion->gConcentrationNoBC,
                dwi,
                patch,
                Ghost::None,
                0);
    new_dw->get(data.gExternalScalarFlux,
                d_mpm_labels->diffusion->gExternalScalarFlux,
                dwi,
                patch,
                Ghost::None,
                0);

    new_dw->getModifiable(data.gConcentrationRate_new,
                          d_mpm_labels->diffusion->gConcentrationRate,
                          dwi,
                          patch);
    new_dw->allocateAndPut(data.gConcentrationStar_new,
                           d_mpm_labels->diffusion->gConcentrationStar,
                           dwi,
                           patch);
  }
}

void
ScalarDiffusionTasks::integrateConcentration(const Patch* patch,
                                             double delT,
                                             int dwi,
                                             constNCVariable<double>& gMass,
                                             ScalarDiffusionTaskData& data)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    for (auto iter = patch->getExtraNodeIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      data.gConcentrationRate_new[c] /= gMass[c];
      data.gConcentrationStar_new[c] =
        data.gConcentration[c] +
        (data.gConcentrationRate_new[c] + data.gFluxRate[c]) * delT;
    }

    MPMBoundCond bc;
    bc.setBoundaryCondition(patch,
                            dwi,
                            "SD-Type",
                            data.gConcentrationStar_new,
                            d_mpm_flags->d_interpolatorType);

    for (auto iter = patch->getExtraNodeIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      data.gConcentrationRate_new[c] =
        (data.gConcentrationStar_new[c] - data.gConcentrationNoBC[c]) / delT +
        data.gExternalScalarFlux[c] / gMass[c];
    }
  } // if doScalarDiffusion
}

void
ScalarDiffusionTasks::scheduleComputeConcentrationGradient(
  SchedulerP& sched,
  const PatchSet* patches,
  const MaterialSet* matls)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    const Level* level = getLevel(patches);
    if (!d_mpm_flags->doMPMOnLevel(level->getIndex(),
                                   level->getGrid()->numLevels())) {
      return;
    }

    printSchedule(patches,
                  cout_doing,
                  "ScalarDiffusionTasks::scheduleComputeConcentrationGradient");

    Task* t = scinew Task("ScalarDiffusionTasks::computeConcentrationGradient",
                          this,
                          &ScalarDiffusionTasks::computeConcentrationGradient);

    t->requires(Task::NewDW,
                d_mpm_labels->diffusion->gConcentrationStar,
                Ghost::AroundCells,
                d_num_ghost_nodes);
    t->requires(Task::OldDW, d_mpm_labels->diffusion->pArea, Ghost::None);
    t->computes(d_mpm_labels->diffusion->pGradConcentration_preReloc);
    t->computes(d_mpm_labels->diffusion->pArea_preReloc);

    sched->addTask(t, patches, matls);
  }
}

void
ScalarDiffusionTasks::scheduleComputeConcentrationGradientAMR(
  const LevelP& level_in,
  SchedulerP& sched,
  const MaterialSet* matls)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    GridP grid    = level_in->getGrid();
    int maxLevels = grid->numLevels();

    for (int l = 0; l < maxLevels; l++) {
      const LevelP& level     = grid->getLevel(l);
      const PatchSet* patches = level->eachPatch();
      scheduleComputeConcentrationGradient(sched, patches, matls);
    }
  }
}

void
ScalarDiffusionTasks::computeConcentrationGradient(const ProcessorGroup*,
                                                   const PatchSubset* patches,
                                                   const MaterialSubset*,
                                                   DataWarehouse* old_dw,
                                                   DataWarehouse* new_dw)
{
  if (!d_mpm_flags->d_doScalarDiffusion) {
    return;
  }

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches,
              patch,
              cout_doing,
              "Doing ScalarDiffusionTasks::computeConcGrad");

    auto interpolator = d_mpm_flags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());

    Vector dx      = patch->dCell();
    double oodx[3] = { 1. / dx.x(), 1. / dx.y(), 1. / dx.z() };

    size_t numMPMMatls = d_mat_manager->getNumMaterials("MPM");

    for (size_t m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      // Get the arrays of particle values to be changed
      constParticleVariable<Vector> pArea;
      ParticleVariable<Vector> pConcGradNew, pAreaNew;

      // Get the arrays of grid data on which the new particle values depend
      constNCVariable<double> gConcStar;

      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      old_dw->get(pArea, d_mpm_labels->diffusion->pArea, pset);
      new_dw->get(gConcStar,
                  d_mpm_labels->diffusion->gConcentrationStar,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_num_ghost_nodes);

      new_dw->allocateAndPut(
        pAreaNew, d_mpm_labels->diffusion->pArea_preReloc, pset);
      new_dw->allocateAndPut(
        pConcGradNew,
        d_mpm_labels->diffusion->pGradConcentration_preReloc,
        pset);

      for (auto& idx : *pset) {
        pConcGradNew[idx] = Vector(0.0, 0.0, 0.0);
        for (size_t k = 0; k < ni.size(); k++) {
          IntVector node = ni[k];
          for (int j = 0; j < 3; j++) {
            pConcGradNew[idx][j] += gConcStar[ni[k]] * d_S[k][j] * oodx[j];
          }
        }

        pAreaNew[idx] = pArea[idx];
      }
    }
  }
}

void
ScalarDiffusionTasks::scheduleInterpolateToParticlesAndUpdate(Task* task)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    task->requires(
      Task::OldDW, d_mpm_labels->diffusion->pConcentration, Ghost::None);
    task->requires(Task::NewDW,
                   d_mpm_labels->diffusion->gConcentrationRate,
                   Ghost::AroundCells,
                   d_num_ghost_nodes);

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

    old_dw->get(
      data.pConcentration, d_mpm_labels->diffusion->pConcentration, pset);
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
ScalarDiffusionTasks::updateReductionVars(
  DataWarehouse* new_dw,
  ScalarDiffusionGlobalConcData& conc_data)
{
  if (d_mpm_flags->d_doAutoCycleBC && d_mpm_flags->d_doScalarDiffusion) {
    if (d_mpm_flags->d_autoCycleUseMinMax) {
      new_dw->put(min_vartype(conc_data.minPatchConc),
                  d_mpm_labels->diffusion->rMinConcentration);
      new_dw->put(max_vartype(conc_data.maxPatchConc),
                  d_mpm_labels->diffusion->rMaxConcentration);
    } else {
      new_dw->put(sum_vartype(conc_data.totalConc),
                  d_mpm_labels->diffusion->rTotalConcentration);
    }
  }
}
void
ScalarDiffusionTasks::scheduleAddParticles(Task* task)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    task->modifies(d_mpm_labels->diffusion->pConcentration_preReloc);
    task->modifies(d_mpm_labels->diffusion->pConcPrevious_preReloc);
    task->modifies(d_mpm_labels->diffusion->pGradConcentration_preReloc);
    task->modifies(d_mpm_labels->diffusion->pExternalScalarFlux_preReloc);
    task->modifies(d_mpm_labels->diffusion->pArea_preReloc);
    task->modifies(d_mpm_labels->diffusion->pDiffusivity_preReloc);
  }
}

void
ScalarDiffusionTasks::scheduleAddSplitParticlesComputesAndRequires(
  Task* task,
  MPMMaterial* mpm_matl,
  const PatchSet* patches)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    ScalarDiffusionModel* sdm = mpm_matl->getScalarDiffusionModel();
    sdm->addSplitParticlesComputesAndRequires(task, mpm_matl, patches);
  }
}

void
ScalarDiffusionTasks::getModifiableForAddParticles(
  ParticleSubset* pset,
  DataWarehouse* new_dw,
  ScalarDiffusionTaskData& data)
{
  ParticleVariable<double> pESF, pConc, pConcpre, pConcgrad;
  if (d_mpm_flags->d_doScalarDiffusion) {
    new_dw->getModifiable(data.pConcentrationNew,
                          d_mpm_labels->diffusion->pConcentration_preReloc,
                          pset);
    new_dw->getModifiable(data.pConcPreviousNew,
                          d_mpm_labels->diffusion->pConcPrevious_preReloc,
                          pset);
    new_dw->getModifiable(data.pConcentrationGradNew,
                          d_mpm_labels->diffusion->pGradConcentration_preReloc,
                          pset);
    new_dw->getModifiable(data.pExternalScalarFluxNew,
                          d_mpm_labels->diffusion->pExternalScalarFlux_preReloc,
                          pset);
    new_dw->getModifiable(
      data.pAreaNew, d_mpm_labels->diffusion->pArea_preReloc, pset);
    new_dw->getModifiable(data.pDiffusivityNew,
                          d_mpm_labels->diffusion->pDiffusivity_preReloc,
                          pset);
  }
}

void
ScalarDiffusionTasks::copyTemporaryForAddParticles(
  size_t old_num_particles,
  ParticleSubset* pset,
  DataWarehouse* new_dw,
  const ScalarDiffusionTaskData& data,
  ScalarDiffusionTaskData& data_tmp)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    new_dw->allocateTemporary(data_tmp.pConcentrationNew, pset);
    new_dw->allocateTemporary(data_tmp.pConcPreviousNew, pset);
    new_dw->allocateTemporary(data_tmp.pConcentrationGradNew, pset);
    new_dw->allocateTemporary(data_tmp.pExternalScalarFluxNew, pset);
    new_dw->allocateTemporary(data_tmp.pAreaNew, pset);
    new_dw->allocateTemporary(data_tmp.pDiffusivityNew, pset);

    for (size_t pp = 0; pp < old_num_particles; ++pp) {
      data_tmp.pConcentrationNew[pp]      = data.pConcentrationNew[pp];
      data_tmp.pConcPreviousNew[pp]       = data.pConcPreviousNew[pp];
      data_tmp.pConcentrationGradNew[pp]  = data.pConcentrationGradNew[pp];
      data_tmp.pExternalScalarFluxNew[pp] = data.pExternalScalarFluxNew[pp];
      data_tmp.pAreaNew[pp]               = data.pAreaNew[pp];
      data_tmp.pDiffusivityNew[pp]        = data.pDiffusivityNew[pp];
    }
  }
}

void
ScalarDiffusionTasks::setDataForTmpAddParticles(
  double fourthOrEighth,
  int cell_idx,
  int idx,
  const ParticleVariable<IntVector>& pLoadCurveID,
  int new_idx,
  int last_idx,
  ParticleVariable<Point>& pX_tmp,
  ParticleVariable<IntVector>& pLoadCurveID_tmp,
  const ScalarDiffusionTaskData& data,
  ScalarDiffusionTaskData& data_tmp)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    data_tmp.pConcentrationNew[new_idx]      = data.pConcentrationNew[idx];
    data_tmp.pConcPreviousNew[new_idx]       = data.pConcPreviousNew[idx];
    data_tmp.pConcentrationGradNew[new_idx]  = data.pConcentrationGradNew[idx];
    data_tmp.pExternalScalarFluxNew[new_idx] = data.pExternalScalarFluxNew[idx];
    data_tmp.pDiffusivityNew[new_idx]        = data.pDiffusivityNew[idx];
    if ((std::abs(data.pAreaNew[idx].x()) > 0.0 &&
         std::abs(data.pAreaNew[idx].y()) > 0.0) ||
        (std::abs(data.pAreaNew[idx].x()) > 0.0 &&
         std::abs(data.pAreaNew[idx].z()) > 0.0) ||
        (std::abs(data.pAreaNew[idx].y()) > 0.0 &&
         std::abs(data.pAreaNew[idx].z()) > 0.0) ||
        (std::abs(data.pAreaNew[idx][0]) < 1.e-12)) {
      data_tmp.pAreaNew[new_idx] = fourthOrEighth * data.pAreaNew[idx];
    } else {
      if (cell_idx == 0) {
        data_tmp.pAreaNew[new_idx] = data.pAreaNew[idx];
      } else {
        if (pX_tmp[new_idx].asVector().length2() >
            pX_tmp[last_idx].asVector().length2()) {
          data_tmp.pAreaNew[last_idx] = 0.0;
          data_tmp.pAreaNew[new_idx]  = data.pAreaNew[idx];
          pLoadCurveID_tmp[last_idx]  = IntVector(0, 0, 0);
          pLoadCurveID_tmp[new_idx]   = pLoadCurveID[idx];
        } else {
          data_tmp.pAreaNew[new_idx] = 0.0;
          pLoadCurveID_tmp[new_idx]  = IntVector(0, 0, 0);
        } // if pxtmp
      }   // if i==0
    }     // if pArea
  }       // if diffusion
}

void
ScalarDiffusionTasks::splitCMSpecificParticleData(
  const Patch* patch,
  MPMMaterial* mpm_matl,
  int dwi,
  double fourOrEight,
  ParticleVariable<int>& pRefined_old,
  ParticleVariable<int>& pRefined,
  size_t oldNumPar,
  size_t numNewPartNeeded,
  DataWarehouse* old_dw,
  DataWarehouse* new_dw)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    ScalarDiffusionModel* sdm = mpm_matl->getScalarDiffusionModel();
    sdm->splitSDMSpecificParticleData(patch,
                                      dwi,
                                      fourOrEight,
                                      pRefined_old,
                                      pRefined,
                                      oldNumPar,
                                      numNewPartNeeded,
                                      old_dw,
                                      new_dw);
  }
}

void
ScalarDiffusionTasks::putTmpDataAddParticles(DataWarehouse* new_dw,
                                             ScalarDiffusionTaskData& data_tmp)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    new_dw->put(data_tmp.pConcentrationNew,
                d_mpm_labels->diffusion->pConcentration_preReloc,
                true);
    new_dw->put(data_tmp.pConcPreviousNew,
                d_mpm_labels->diffusion->pConcPrevious_preReloc,
                true);
    new_dw->put(data_tmp.pConcentrationGradNew,
                d_mpm_labels->diffusion->pGradConcentration_preReloc,
                true);
    new_dw->put(data_tmp.pExternalScalarFluxNew,
                d_mpm_labels->diffusion->pExternalScalarFlux,
                true);

    new_dw->put(
      data_tmp.pAreaNew, d_mpm_labels->diffusion->pArea_preReloc, true);
    new_dw->put(data_tmp.pDiffusivityNew,
                d_mpm_labels->diffusion->pDiffusivity_preReloc,
                true);
  }
}

void
ScalarDiffusionTasks::scheduleRefine(Task* task)
{
  if (d_mpm_flags->d_doScalarDiffusion) {
    task->computes(d_mpm_labels->diffusion->pConcentration);
    task->computes(d_mpm_labels->diffusion->pConcPrevious);
    task->computes(d_mpm_labels->diffusion->pGradConcentration);
    task->computes(d_mpm_labels->diffusion->pArea);
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
