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

#include <CCA/Components/MPM/HeatConduction/HeatConductionTasks.h>

#include <CCA/Components/MPM/HeatConduction/HeatConduction.h>
#include <CCA/Components/MPM/ThermalContact/ThermalContact.h>
#include <CCA/Components/MPM/ThermalContact/ThermalContactFactory.h>

#include <CCA/Components/MPM/Core/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>

#include <CCA/Ports/DataWarehouse.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Patch.h>

using namespace Uintah;

static DebugStream cout_doing("HeatConduction_MPM", false);

HeatConductionTasks::HeatConductionTasks(const ProblemSpecP& ps,
                                         MaterialManagerP& mat_manager,
                                         const MPMLabel* mpm_labels,
                                         const MPMFlags* mpm_flags)
{
  d_mat_manager = mat_manager;
  d_mpm_labels  = mpm_labels;
  d_mpm_flags   = mpm_flags;

  thermalContactModel =
    ThermalContactFactory::create(ps, mat_manager, mpm_labels, mpm_flags);

  heatConductionModel =
    std::make_unique<HeatConduction>(mat_manager, mpm_labels, mpm_flags);
}

void
HeatConductionTasks::outputProblemSpec(ProblemSpecP& ps)
{
  thermalContactModel->outputProblemSpec(ps);
}

void
HeatConductionTasks::scheduleCompute(SchedulerP& sched,
                                     const PatchSet* patches,
                                     const MaterialSet* matls)
{
  if (d_mpm_flags->d_doExplicitHeatConduction) {
    scheduleComputeHeatExchange(sched, patches, matls);
    scheduleComputeInternalHeatRate(sched, patches, matls);
    scheduleComputeNodalHeatFlux(sched, patches, matls);
    scheduleSolveHeatEquations(sched, patches, matls);
    scheduleIntegrateTemperatureRate(sched, patches, matls);
  }
}

/*!----------------------------------------------------------------------
 * scheduleComputeHeatExchange
 *-----------------------------------------------------------------------*/
void
HeatConductionTasks::scheduleComputeHeatExchange(SchedulerP& sched,
                                                 const PatchSet* patches,
                                                 const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }
  /* computeHeatExchange
   *   in(G.MASS, G.TEMPERATURE, G.EXTERNAL_HEAT_RATE)
   *   operation(peform heat exchange which will cause each of
   *   velocity fields to exchange heat according to
   *   the temperature differences)
   *   out(G.EXTERNAL_HEAT_RATE) */

  printSchedule(patches, cout_doing, "MPM::scheduleComputeHeatExchange");

  Task* t = scinew Task("ThermalContact::computeHeatExchange",
                        thermalContactModel.get(),
                        &ThermalContact::computeHeatExchange);

  thermalContactModel->addComputesAndRequires(t, patches, matls);
  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * scheduleComputeInternalHeatRate
 *-----------------------------------------------------------------------*/
void
HeatConductionTasks::scheduleComputeInternalHeatRate(SchedulerP& sched,
                                                     const PatchSet* patches,
                                                     const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }
  printSchedule(patches, cout_doing, "MPM::scheduleComputeInternalHeatRate");
  heatConductionModel->scheduleComputeInternalHeatRate(sched, patches, matls);
}

/*!----------------------------------------------------------------------
 * scheduleComputeNodalHeatFlux
 *-----------------------------------------------------------------------*/
void
HeatConductionTasks::scheduleComputeNodalHeatFlux(SchedulerP& sched,
                                                  const PatchSet* patches,
                                                  const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }
  printSchedule(patches, cout_doing, "MPM::scheduleComputeNodalHeatFlux");
  heatConductionModel->scheduleComputeNodalHeatFlux(sched, patches, matls);
}

/*!----------------------------------------------------------------------
 * scheduleSolveHeatEquations
 *-----------------------------------------------------------------------*/
void
HeatConductionTasks::scheduleSolveHeatEquations(SchedulerP& sched,
                                                const PatchSet* patches,
                                                const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }
  printSchedule(patches, cout_doing, "MPM::scheduleSolveHeatEquations");
  heatConductionModel->scheduleSolveHeatEquations(sched, patches, matls);
}

/*!----------------------------------------------------------------------
 * scheduleIntegrateTemperatureRate
 *-----------------------------------------------------------------------*/
void
HeatConductionTasks::scheduleIntegrateTemperatureRate(SchedulerP& sched,
                                                      const PatchSet* patches,
                                                      const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }
  printSchedule(patches, cout_doing, "MPM::scheduleIntegrateTemperatureRate");
  heatConductionModel->scheduleIntegrateTemperatureRate(sched, patches, matls);
}
