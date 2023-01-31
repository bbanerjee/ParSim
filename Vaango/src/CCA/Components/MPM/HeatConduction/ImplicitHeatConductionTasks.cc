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

#include <CCA/Components/MPM/HeatConduction/ImplicitHeatConductionTasks.h>

#include <CCA/Components/MPM/HeatConduction/ImplicitHeatConduction.h>

#include <CCA/Components/MPM/ThermalContact/ThermalContact.h>
#include <CCA/Components/MPM/ThermalContact/ThermalContactFactory.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/ImpMPMFlags.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPM/PhysicalBC/HeatFluxBC.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBCFactory.h>

#include <CCA/Ports/DataWarehouse.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Patch.h>

using namespace Uintah;

static DebugStream cout_doing("ImplicitHeatConduction_MPM", false);

ImplicitHeatConductionTasks::ImplicitHeatConductionTasks(
  const ProblemSpecP& ps,
  const MaterialManagerP& mat_manager,
  const MPMLabel* mpm_labels,
  const ImpMPMFlags* mpm_flags)
  : d_mat_manager(mat_manager)
  , d_mpm_labels(mpm_labels)
  , d_mpm_flags(mpm_flags)
{
  thermalContactModel =
    ThermalContactFactory::create(ps, mat_manager, mpm_labels, mpm_flags);

  heatConductionModel = std::make_unique<ImplicitHeatConduction>(
    mat_manager, mpm_labels, mpm_flags);

  heatConductionModel->problemSetup(d_mpm_flags->d_solverType);
}

void
ImplicitHeatConductionTasks::outputProblemSpec(ProblemSpecP& ps)
{
  thermalContactModel->outputProblemSpec(ps);
}

void
ImplicitHeatConductionTasks::scheduleInitializeHeatFluxBCs(const LevelP& level,
                                                           SchedulerP& sched)
{
  d_loadCurveIndex = std::make_shared<MaterialSubset>();
  d_loadCurveIndex->add(0);
  d_loadCurveIndex->addReference();

  int nofHeatFluxBCs = 0;
  for (auto particleBC : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    auto bcType = particleBC->getType();
    if (bcType == "HeatFlux") {
      d_loadCurveIndex->add(nofHeatFluxBCs++);
    }
  }

  if (nofHeatFluxBCs > 0) {

    // Create a task that calculates the total number of particles
    // associated with each load curve.
    Task* t = scinew Task(
      "ImplicitHeatConductionTasks::countMaterialPointsPerLoadCurve",
      this,
      &ImplicitHeatConductionTasks::countMaterialPointsPerLoadCurve);
    t->requires(Task::NewDW, d_mpm_labels->pLoadCurveIDLabel, Ghost::None);
    t->computes(d_mpm_labels->materialPointsPerLoadCurveLabel,
                d_loadCurveIndex.get(),
                Task::OutOfDomain);
    sched->addTask(t, level->eachPatch(), d_mat_manager->allMaterials("MPM"));

    // Create a task that calculates the heatflux to be associated with
    // each particle based on the HeatFluxBCs
    t = scinew Task("ImplicitHeatConductionTasks::initializeHeatFluxBC",
                    this,
                    &ImplicitHeatConductionTasks::initializeHeatFluxBC);
    t->requires(Task::NewDW, d_mpm_labels->pXLabel, Ghost::None);
    t->requires(Task::NewDW, d_mpm_labels->pLoadCurveIDLabel, Ghost::None);
    t->requires(Task::NewDW,
                d_mpm_labels->materialPointsPerLoadCurveLabel,
                d_loadCurveIndex.get(),
                Task::OutOfDomain,
                Ghost::None);
    t->modifies(d_mpm_labels->pExternalHeatFluxLabel);
    sched->addTask(t, level->eachPatch(), d_mat_manager->allMaterials("MPM"));
  }

  d_loadCurveIndex->removeReference();
}

void
ImplicitHeatConductionTasks::countMaterialPointsPerLoadCurve(
  const ProcessorGroup*,
  const PatchSubset* patches,
  const MaterialSubset*,
  DataWarehouse*,
  DataWarehouse* new_dw)
{
  // Find the number of pressure BCs in the problem
  int nofHeatFluxBCs = 0;
  for (auto particleBC : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    auto bcType = particleBC->getType();
    if (bcType == "HeatFlux") {
      nofHeatFluxBCs++;
      // std::cout << "nofHeatFluxBCs = " << nofHeatFluxBCs << std::endl;

      for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        size_t numMPMMatls = d_mat_manager->getNumMaterials("MPM");
        int numPts         = 0;
        for (size_t m = 0; m < numMPMMatls; m++) {
          MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
            static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m)));
          int dwi = mpm_matl->getDWIndex();

          ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
          constParticleVariable<int> pLoadCurveID;
          new_dw->get(pLoadCurveID, d_mpm_labels->pLoadCurveIDLabel, pset);

          for (auto idx : *pset) {
            if (pLoadCurveID[idx] == (nofHeatFluxBCs)) {
              ++numPts;
            }
          }
        } // matl loop
        // std::cout << "numPts found = " << numPts << std::endl;
        new_dw->put(sumlong_vartype(numPts),
                    d_mpm_labels->materialPointsPerLoadCurveLabel,
                    0,
                    nofHeatFluxBCs - 1);
      } // patch loop
    }   // end heat flux bc if
  }     // end particle BC loop
}

void
ImplicitHeatConductionTasks::initializeHeatFluxBC(const ProcessorGroup*,
                                                  const PatchSubset* patches,
                                                  const MaterialSubset*,
                                                  DataWarehouse*,
                                                  DataWarehouse* new_dw)
{
  // Get the current time
  double time = 0.0;

  // Calculate the heat flux at each particle
  int nofHeatFluxBCs = 0;
  for (auto particleBC : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    auto bcType = particleBC->getType();
    if (bcType == "HeatFlux") {

      // Get the material points per load curve
      sumlong_vartype numPart = 0;
      new_dw->get(numPart,
                  d_mpm_labels->materialPointsPerLoadCurveLabel,
                  0,
                  nofHeatFluxBCs++);

      double fluxPerPart = 0.;
      HeatFluxBC* phf    = nullptr;
      if (bcType == "HeatFlux") {
        phf = dynamic_cast<HeatFluxBC*>(particleBC.get());
        phf->numMaterialPoints(numPart);
        fluxPerPart = phf->fluxPerParticle(time);
        // std::cout << "numPart = " << numPart << std::endl;
        // std::cout << "fluxPerPart = " << fluxPerPart << std::endl;
      }

      // Loop through the patches and calculate the force vector
      // at each particle
      for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        int numMPMMatls    = d_mat_manager->getNumMaterials("MPM");
        for (int m = 0; m < numMPMMatls; m++) {
          MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
            static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m)));
          int matID = mpm_matl->getDWIndex();

          ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);
          constParticleVariable<Point> pX;
          constParticleVariable<int> pLoadCurveID;
          ParticleVariable<double> pExternalHeatFlux;
          new_dw->get(pX, d_mpm_labels->pXLabel, pset);
          new_dw->get(pLoadCurveID, d_mpm_labels->pLoadCurveIDLabel, pset);
          new_dw->getModifiable(
            pExternalHeatFlux, d_mpm_labels->pExternalHeatFluxLabel, pset);

          for (auto idx : *pset) {
            if (pLoadCurveID[idx] == nofHeatFluxBCs) {
              if (bcType == "HeatFlux") {
                pExternalHeatFlux[idx] = phf->getFlux(pX[idx], fluxPerPart);
              }
              // std::cout << "pExternalHeatFlux[idx] = " <<
              // pExternalHeatFlux[idx] << std::endl;
            }
          }
        } // matl loop
      }   // patch loop
    }
  }
}

void
ImplicitHeatConductionTasks::scheduleProjectHeatSource(
  SchedulerP& sched,
  const PatchSet* perproc_patches,
  const MaterialSubset* one_material,
  const MaterialSet* matls)
{
  if (d_mpm_flags->d_projectHeatSource) {
    scheduleComputeCCVolume(sched, perproc_patches, one_material, matls);
    scheduleProjectCCHeatSourceToNodes(
      sched, perproc_patches, one_material, matls);
  }
}

void
ImplicitHeatConductionTasks::scheduleComputeCCVolume(
  SchedulerP& sched,
  const PatchSet* patches,
  const MaterialSubset* one_matl,
  const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleComputeCCVolume");
  Task* t = scinew Task("ImplicitHeatConductionTasks::computeCCVolume",
                        this,
                        &ImplicitHeatConductionTasks::computeCCVolume);

  t->requires(Task::OldDW,
              d_mpm_labels->NC_CCweightLabel,
              one_matl,
              Ghost::AroundCells,
              1);
  t->requires(Task::NewDW, d_mpm_labels->gVolumeLabel, Ghost::AroundCells, 1);

  t->computes(d_mpm_labels->cVolumeLabel);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImplicitHeatConductionTasks::computeCCVolume(const ProcessorGroup*,
                                             const PatchSubset* patches,
                                             const MaterialSubset*,
                                             DataWarehouse* old_dw,
                                             DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing computeCCVolume");

    Ghost::GhostType gac = Ghost::AroundCells;

    constNCVariable<double> NC_CCweight;
    old_dw->get(NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gac, 1);

    size_t numMPMMatls = d_mat_manager->getNumMaterials("MPM");
    for (size_t m = 0; m < numMPMMatls; m++) {
      int dwi = d_mat_manager->getMaterial("MPM", m)->getDWIndex();

      constNCVariable<double> gVolume;
      CCVariable<double> cvolume;

      new_dw->get(gVolume, d_mpm_labels->gVolumeLabel, dwi, patch, gac, 1);
      new_dw->allocateAndPut(cvolume, d_mpm_labels->cVolumeLabel, dwi, patch);
      cvolume.initialize(1.e-20);

      for (CellIterator iter = patch->getExtraCellIterator(); !iter.done();
           iter++) {
        IntVector c = *iter;
        IntVector nodeIdx[8];
        patch->findNodesFromCell(c, nodeIdx);
        for (int in = 0; in < 8; in++) {
          cvolume[c] += NC_CCweight[nodeIdx[in]] * gVolume[nodeIdx[in]];
        }
      }
    }
  }
}

void
ImplicitHeatConductionTasks::scheduleProjectCCHeatSourceToNodes(
  SchedulerP& sched,
  const PatchSet* patches,
  const MaterialSubset* one_matl,
  const MaterialSet* matls)
{
  printSchedule(
    patches, cout_doing, "IMPM::scheduleProjectCCHeatSourceToNodes");
  Task* t =
    scinew Task("ImplicitHeatConductionTasks::projectCCHeatSourceToNodes",
                this,
                &ImplicitHeatConductionTasks::projectCCHeatSourceToNodes);

  t->requires(Task::OldDW,
              d_mpm_labels->NC_CCweightLabel,
              one_matl,
              Ghost::AroundCells,
              1);
  t->requires(Task::NewDW, d_mpm_labels->gVolumeLabel, Ghost::AroundCells, 1);
  t->requires(Task::NewDW, d_mpm_labels->cVolumeLabel, Ghost::AroundCells, 1);
  t->requires(
    Task::OldDW, d_mpm_labels->heatRate_CCLabel, Ghost::AroundCells, 1);

  t->computes(d_mpm_labels->heatRate_CCLabel);
  t->modifies(d_mpm_labels->gExternalHeatRateLabel);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImplicitHeatConductionTasks::projectCCHeatSourceToNodes(
  const ProcessorGroup*,
  const PatchSubset* patches,
  const MaterialSubset*,
  DataWarehouse* old_dw,
  DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches,
              patch,
              cout_doing,
              "Doing projectCCHeatSourceToNodes on patch\t\t");

    Ghost::GhostType gac = Ghost::AroundCells;

    constNCVariable<double> NC_CCweight;
    old_dw->get(NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gac, 1);

    size_t numMPMMatls = d_mat_manager->getNumMaterials("MPM");
    for (size_t m = 0; m < numMPMMatls; m++) {
      int dwi = d_mat_manager->getMaterial("MPM", m)->getDWIndex();

      constNCVariable<double> gVolume;
      NCVariable<double> gextHR;
      constCCVariable<double> CCheatrate, cvolume;
      CCVariable<double> CCheatrate_copy;

      new_dw->get(gVolume, d_mpm_labels->gVolumeLabel, dwi, patch, gac, 1);
      old_dw->get(
        CCheatrate, d_mpm_labels->heatRate_CCLabel, dwi, patch, gac, 1);
      new_dw->get(cvolume, d_mpm_labels->cVolumeLabel, dwi, patch, gac, 1);
      new_dw->getModifiable(
        gextHR, d_mpm_labels->gExternalHeatRateLabel, dwi, patch);
      new_dw->allocateAndPut(
        CCheatrate_copy, d_mpm_labels->heatRate_CCLabel, dwi, patch);

      // carry forward heat rate.
      CCheatrate_copy.copyData(CCheatrate);

      // Project  CC heat rate to nodes
      for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
        IntVector n = *iter;
        IntVector cIdx[8];
        patch->findCellsFromNode(n, cIdx);
        for (int ic = 0; ic < 8; ic++) {
          double solid_volume = cvolume[cIdx[ic]];
          gextHR[n] +=
            CCheatrate[cIdx[ic]] * (NC_CCweight[n] * gVolume[n]) / solid_volume;
        }
      }
    }
  }
}

void
ImplicitHeatConductionTasks::scheduleDestroyAndCreateMatrix(
  SchedulerP& sched,
  const PatchSet* patches,
  const MaterialSet* matls)
{
  scheduleDestroyMatrix(sched, patches, matls);
  scheduleCreateMatrix(sched, patches, matls);
}

void
ImplicitHeatConductionTasks::scheduleDestroyMatrix(SchedulerP& sched,
                                                   const PatchSet* patches,
                                                   const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleDestroyHCMatrix");
  heatConductionModel->scheduleDestroyHCMatrix(sched, patches, matls);
}

void
ImplicitHeatConductionTasks::scheduleCreateMatrix(SchedulerP& sched,
                                                  const PatchSet* patches,
                                                  const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleCreateHCMatrix");
  heatConductionModel->scheduleCreateHCMatrix(sched, patches, matls);
}

void
ImplicitHeatConductionTasks::scheduleApplyBoundaryConditions(
  SchedulerP& sched,
  const PatchSet* patches,
  const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleApplyHCBoundaryConditions");
  heatConductionModel->scheduleApplyHCBoundaryConditions(sched, patches, matls);
}

void
ImplicitHeatConductionTasks::scheduleSolve(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSet* matls)
{
  scheduleFindFixedHCDOF(sched, patches, matls);
  scheduleFormHCStiffnessMatrix(sched, patches, matls);
  scheduleFormHCQ(sched, patches, matls);
  scheduleAdjustHCQAndHCKForBCs(sched, patches, matls);
  scheduleSolveForTemp(sched, patches, matls);
  scheduleGetTemperatureIncrement(sched, patches, matls);
}

void
ImplicitHeatConductionTasks::scheduleFindFixedHCDOF(SchedulerP& sched,
                                                    const PatchSet* patches,
                                                    const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleFindFixedHCDOF");
  heatConductionModel->scheduleFindFixedHCDOF(sched, patches, matls);
}

void
ImplicitHeatConductionTasks::scheduleFormHCStiffnessMatrix(
  SchedulerP& sched,
  const PatchSet* patches,
  const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleFormHCStiffnessMatrix");
  heatConductionModel->scheduleFormHCStiffnessMatrix(sched, patches, matls);
}

void
ImplicitHeatConductionTasks::scheduleFormHCQ(SchedulerP& sched,
                                             const PatchSet* patches,
                                             const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleFormHCQ");
  heatConductionModel->scheduleFormHCQ(sched, patches, matls);
}

void
ImplicitHeatConductionTasks::scheduleAdjustHCQAndHCKForBCs(
  SchedulerP& sched,
  const PatchSet* patches,
  const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleFormHCQAndHCKForBCs");
  heatConductionModel->scheduleAdjustHCQAndHCKForBCs(sched, patches, matls);
}

void
ImplicitHeatConductionTasks::scheduleSolveForTemp(SchedulerP& sched,
                                                  const PatchSet* patches,
                                                  const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleSolveForTemp");
  heatConductionModel->scheduleSolveForTemp(sched, patches, matls);
}

void
ImplicitHeatConductionTasks::scheduleGetTemperatureIncrement(
  SchedulerP& sched,
  const PatchSet* patches,
  const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleGetTemperatureIncrement");
  heatConductionModel->scheduleGetTemperatureIncrement(sched, patches, matls);
}

void
ImplicitHeatConductionTasks::scheduleComputeHeatExchange(
  SchedulerP& sched,
  const PatchSet* patches,
  const MaterialSet* matls)
{
  /* computeHeatExchange
   *   in(G.MASS, G.TEMPERATURE, G.EXTERNAL_HEAT_RATE)
   *   operation(peform heat exchange which will cause each of
   *   velocity fields to exchange heat according to
   *   the temperature differences)
   *   out(G.EXTERNAL_HEAT_RATE) */

  printSchedule(patches, cout_doing, "IMPM::scheduleComputeHeatExchange");
  Task* t = scinew Task("ThermalContact::computeHeatExchange",
                        thermalContactModel.get(),
                        &ThermalContact::computeHeatExchange);

  thermalContactModel->addComputesAndRequires(t, patches, matls);
  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * scheduleComputeInternalHeatRate
 *-----------------------------------------------------------------------*/
void
ImplicitHeatConductionTasks::scheduleComputeInternalHeatRate(
  SchedulerP& sched,
  const PatchSet* patches,
  const MaterialSet* matls)
{
}

/*!----------------------------------------------------------------------
 * scheduleComputeNodalHeatFlux
 *-----------------------------------------------------------------------*/
void
ImplicitHeatConductionTasks::scheduleComputeNodalHeatFlux(
  SchedulerP& sched,
  const PatchSet* patches,
  const MaterialSet* matls)
{
}

/*!----------------------------------------------------------------------
 * scheduleSolveHeatEquations
 *-----------------------------------------------------------------------*/
void
ImplicitHeatConductionTasks::scheduleSolveHeatEquations(
  SchedulerP& sched,
  const PatchSet* patches,
  const MaterialSet* matls)
{
}

/*!----------------------------------------------------------------------
 * scheduleIntegrateTemperatureRate
 *-----------------------------------------------------------------------*/
void
ImplicitHeatConductionTasks::scheduleIntegrateTemperatureRate(
  SchedulerP& sched,
  const PatchSet* patches,
  const MaterialSet* matls)
{
}
