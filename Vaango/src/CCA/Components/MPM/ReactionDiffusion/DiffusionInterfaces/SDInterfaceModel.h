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

#ifndef __SDINTERFACEMODEL_H__
#define __SDINTERFACEMODEL_H__

#include <CCA/Components/MPM/Contact/ContactMaterialSpec.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/Scheduler.h>
#include <CCA/Ports/SchedulerP.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {

class SDInterfaceModel
{
public:
  SDInterfaceModel(ProblemSpecP& ps,
                   const MaterialManager* sS,
                   const MPMFlags* mpm_flags,
                   const MPMLabel* mpm_lb);

  virtual ~SDInterfaceModel();

  virtual void
  addComputesAndRequiresInterpolated(SchedulerP& sched,
                                     const PatchSet* patches,
                                     const MaterialSet* matls);

  virtual void
  sdInterfaceInterpolated(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw);

  virtual void
  addComputesAndRequiresDivergence(SchedulerP& sched,
                                   const PatchSet* patches,
                                   const MaterialSet* matls);

  virtual void
  sdInterfaceDivergence(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset* matls,
                        DataWarehouse* old_dw,
                        DataWarehouse* new_dw);

  virtual void
  outputProblemSpec(ProblemSpecP& ps);

  const VarLabel*
  getInterfaceFluxLabel() const;

  const VarLabel*
  getInterfaceFlagLabel() const;

protected:
  void
  setBaseComputesAndRequiresDivergence(Task* task, const MaterialSubset* matls);

  const MPMLabel* d_mpm_lb;
  const MaterialManager* d_materialManager;
  ContactMaterialSpec d_materials_list;
  const MPMFlags* d_mpm_flags;

  // Stores dC/dt at the interface points.
  VarLabel* sdInterfaceRate;
  VarLabel* sdInterfaceFlag; // True means interface at point

  SDInterfaceModel(const SDInterfaceModel&);
  SDInterfaceModel&
  operator=(const SDInterfaceModel&);
};

} // end namespace Uintah
#endif
