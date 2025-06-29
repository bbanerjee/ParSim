/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2025 Biswajit Banerjee, Parresia Research Ltd, NZ
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

#ifndef UINTAH_CCA_COMPONENTS_MPM_PHYSICALBC_FLUXBCMODEL_H
#define UINTAH_CCA_COMPONENTS_MPM_PHYSICALBC_FLUXBCMODEL_H

#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/Scheduler.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Parallel/ProcessorGroup.h>

namespace Uintah {

class FluxBCModel
{
public:
  FluxBCModel(const MaterialManager* materialManager,
              const MPMLabel* mpm_labels,
              const MPMFlags* mpm_flags);

  virtual ~FluxBCModel() = default;

  virtual void
  scheduleInitializeScalarFluxBCs(const LevelP& level, SchedulerP& sched);

  virtual void
  scheduleApplyExternalScalarFlux(SchedulerP& sched,
                                  const PatchSet* patches,
                                  const MaterialSet* matls);

protected:
  virtual void
  initializeScalarFluxBC(const ProcessorGroup*,
                         const PatchSubset* patches,
                         const MaterialSubset*,
                         DataWarehouse* old_dw,
                         DataWarehouse* new_dw);

  virtual void
  applyExternalScalarFlux(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset*,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw);

  virtual void
  countMaterialPointsPerFluxLoadCurve(const ProcessorGroup*,
                                      const PatchSubset* patches,
                                      const MaterialSubset*,
                                      DataWarehouse* old_dw,
                                      DataWarehouse* new_dw);

  FluxBCModel(const FluxBCModel&) = delete;
  FluxBCModel&
  operator=(const FluxBCModel&) = delete;

  MaterialSubset* d_load_curve_index;
  const MaterialManager* d_materialManager;
  const MPMLabel* d_mpm_lb;
  const MPMFlags* d_mpm_flags;
};
} // namespace Uintah

#endif // UINTAH_CCA_COMPONENTS_MPM_PHYSICALBC_FLUXBCMODEL_H
