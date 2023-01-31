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

#ifndef __CCA_COMPONENTS_MPM_HEATCONDUCTION_HEATCONDUCTIONTASKS_H__
#define __CCA_COMPONENTS_MPM_HEATCONDUCTION_HEATCONDUCTIONTASKS_H__

#include <CCA/Ports/Scheduler.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Task.h>

#include <map>
#include <vector>

namespace Uintah {

class Patch;
class DataWarehouse;
class MPMFlags;
class MPMLabel;
class ThermalContact;
class HeatConduction;

class HeatConductionTasks final
{
public:
  HeatConductionTasks(const ProblemSpecP& ps,
                      MaterialManagerP& ss,
                      const MPMLabel* mpm_labels,
                      const MPMFlags* mpm_flags);

  ~HeatConductionTasks() = default;

  void
  outputProblemSpec(ProblemSpecP& ps);

  void
  scheduleCompute(SchedulerP& sched,
                  const PatchSet* patches,
                  const MaterialSet* matls);

private:
  void
  scheduleComputeHeatExchange(SchedulerP& sched,
                              const PatchSet* patches,
                              const MaterialSet* matls);

  void
  scheduleComputeInternalHeatRate(SchedulerP& sched,
                                  const PatchSet* patches,
                                  const MaterialSet* matls);

  void
  scheduleComputeNodalHeatFlux(SchedulerP& sched,
                               const PatchSet* patches,
                               const MaterialSet* matls);

  void
  scheduleSolveHeatEquations(SchedulerP& sched,
                             const PatchSet* patches,
                             const MaterialSet* matls);

  void
  scheduleIntegrateTemperatureRate(SchedulerP& sched,
                                   const PatchSet* patches,
                                   const MaterialSet* matls);

private:
  const MPMLabel* d_mpm_labels;
  const MPMFlags* d_mpm_flags;
  MaterialManagerP d_mat_manager;

  std::unique_ptr<ThermalContact> thermalContactModel{ nullptr };
  std::unique_ptr<HeatConduction> heatConductionModel{ nullptr };
};

} // End of namespace Uintah

#endif // __CCA_COMPONENTS_MPM_HEATCONDUCTION_HEATCONDUCTIONTASKS_H__
