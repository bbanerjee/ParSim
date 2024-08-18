/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#ifndef UINTAH_HEAT_CONDUCTION_H
#define UINTAH_HEAT_CONDUCTION_H

#include <CCA/Ports/SchedulerP.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Variables/ComputeSet.h>

namespace Uintah {

class MPMLabel;
class MPMFlags;
class DataWarehouse;
class ProcessorGroup;

class HeatConduction
{
public:
  HeatConduction(const MaterialManagerP& mat_manager,
                 const MPMLabel* labels,
                 const MPMFlags* flags);

  ~HeatConduction() = default;

  void
  scheduleComputeInternalHeatRate(SchedulerP&,
                                  const PatchSet*,
                                  const MaterialSet*);

  void
  scheduleComputeNodalHeatFlux(SchedulerP&,
                               const PatchSet*,
                               const MaterialSet*);

  void
  scheduleSolveHeatEquations(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleIntegrateTemperatureRate(SchedulerP&,
                                   const PatchSet*,
                                   const MaterialSet*);

  void
  computeInternalHeatRate(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw);

  void
  computeNodalHeatFlux(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset*,
                       DataWarehouse* /*old_dw*/,
                       DataWarehouse* new_dw);

  void
  solveHeatEquations(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset* matls,
                     DataWarehouse* /*old_dw*/,
                     DataWarehouse* new_dw);

  void
  integrateTemperatureRate(const ProcessorGroup*,
                           const PatchSubset* patches,
                           const MaterialSubset* matls,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw);

private:
  const MaterialManagerP d_mat_manager{nullptr};
  const MPMLabel* d_mpm_labels{nullptr};
  const MPMFlags* d_mpm_flags{nullptr};
  int d_num_ghost_particles{2};
  int d_num_ghost_nodes{2};

  HeatConduction(const HeatConduction&) = delete;
  HeatConduction(HeatConduction&&) = delete;
  HeatConduction&
  operator=(const HeatConduction&) = delete;
  HeatConduction&
  operator=(HeatConduction&&) = delete;
};

} // end namespace Uintah
#endif
