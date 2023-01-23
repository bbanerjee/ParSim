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

#include <CCA/Components/MPM/ScalarDiffusionTasks.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>

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
                                           const MaterialManager* mat_manager,
                                           const MPMFlags* mpm_flags,
                                           const MPMLabel* mpm_labels)
{
  d_mat_manager = mat_manager;
  d_mpm_flags   = mpm_flags;
  d_mpm_labels  = mpm_labels;
  sdInterfaceModel =
    SDInterfaceModelFactory::create(ps, mat_manager, mpm_flags, mpm_labels);
}

void
ScarDiffusionTasks::scheduleInterpolateToParticlesAndUpdate(Task* task,
                                                            int numGhostNodes)
{
  t->requires(Task::OldDW,
              d_mpm_labels->diffusion->pConcentration,
              Ghost::None);
  t->requires(Task::NewDW,
              d_mpm_labels->diffusion->gConcentrationRate,
              Ghost::AroundCells,
              numGhostNodes);

  t->computes(d_mpm_labels->diffusion->pConcentration_preReloc);
  t->computes(d_mpm_labels->diffusion->pConcPrevious_preReloc);

  if (d_mpm_flags->d_doAutoCycleBC) {
    if (d_mpm_flags->d_autoCycleUseMinMax) {
      t->computes(d_mpm_labels->diffusion->rMinConcentration);
      t->computes(d_mpm_labels->diffusion->rMaxConcentration);
    } else {
      t->computes(d_mpm_labels->diffusion->rTotalConcentration);
    }
  }
}

void
ScalarDiffusionTasks::scheduleConcInterpolated(SchedulerP& sched,
                                               const PatchSet* patches,
                                               const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleConcInterpolated");

  sdInterfaceModel->addComputesAndRequiresInterpolated(sched, patches, matls);
}

void
ScalarDiffusionTasks::scheduleComputeFlux(SchedulerP& sched,
                                          const PatchSet* patches,
                                          const MaterialSet* matls)
{
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
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleDiffusionInterfaceDiv");

  sdInterfaceModel->addComputesAndRequiresDivergence(sched, patches, matls);
}

} // end namespace Uintah
