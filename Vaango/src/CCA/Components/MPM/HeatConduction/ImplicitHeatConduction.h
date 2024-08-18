/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef UINTAH_IMPLICIT_HEAT_CONDUCTION_H
#define UINTAH_IMPLICIT_HEAT_CONDUCTION_H

#include <sci_defs/petsc_defs.h>

#include <CCA/Ports/SchedulerP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Variables/Array3.h>
#include <Core/Grid/Variables/ComputeSet.h>

#include <cmath>
#include <string>
#include <vector>

namespace Uintah {

class MPMLabel;
class MPMFlags;
class DataWarehouse;
class ProcessorGroup;
class Solver;

class ImplicitHeatConduction
{
public:
  ImplicitHeatConduction(const MaterialManagerP& mat_manager,
                         const MPMLabel* labels,
                         const MPMFlags* flags);
  ~ImplicitHeatConduction();

  ImplicitHeatConduction(const ImplicitHeatConduction&) = delete;
  ImplicitHeatConduction(ImplicitHeatConduction&&)      = delete;
  ImplicitHeatConduction&
  operator=(const ImplicitHeatConduction&) = delete;
  ImplicitHeatConduction&
  operator=(ImplicitHeatConduction&&) = delete;

  void
  problemSetup(std::string solver_type);

  void
  scheduleFormHCStiffnessMatrix(SchedulerP&,
                                const PatchSet*,
                                const MaterialSet*);

  void
  scheduleFormHCQ(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleAdjustHCQAndHCKForBCs(SchedulerP&,
                                const PatchSet*,
                                const MaterialSet*);

  void
  scheduleDestroyHCMatrix(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleCreateHCMatrix(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleApplyHCBoundaryConditions(SchedulerP&,
                                    const PatchSet*,
                                    const MaterialSet*);

  void
  scheduleFindFixedHCDOF(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleSolveForTemp(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleGetTemperatureIncrement(SchedulerP&,
                                  const PatchSet*,
                                  const MaterialSet*);

  void
  destroyHCMatrix(const ProcessorGroup*,
                  const PatchSubset* patches,
                  const MaterialSubset* matls,
                  DataWarehouse* old_dw,
                  DataWarehouse* new_dw);

  void
  createHCMatrix(const ProcessorGroup*,
                 const PatchSubset* patches,
                 const MaterialSubset* matls,
                 DataWarehouse* old_dw,
                 DataWarehouse* new_dw);

  void
  applyHCBoundaryConditions(const ProcessorGroup*,
                            const PatchSubset* patches,
                            const MaterialSubset* matls,
                            DataWarehouse* old_dw,
                            DataWarehouse* new_dw);

  void
  findFixedHCDOF(const ProcessorGroup*,
                 const PatchSubset* patches,
                 const MaterialSubset* matls,
                 DataWarehouse* old_dw,
                 DataWarehouse* new_dw);

  void
  formHCStiffnessMatrix(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset* matls,
                        DataWarehouse* old_dw,
                        DataWarehouse* new_dw);

  void
  formHCQ(const ProcessorGroup*,
          const PatchSubset* patches,
          const MaterialSubset* matls,
          DataWarehouse* old_dw,
          DataWarehouse* new_dw);

  void
  adjustHCQAndHCKForBCs(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset* matls,
                        DataWarehouse* old_dw,
                        DataWarehouse* new_dw);

  void
  solveForTemp(const ProcessorGroup*,
               const PatchSubset* patches,
               const MaterialSubset* matls,
               DataWarehouse* old_dw,
               DataWarehouse* new_dw);

  void
  getTemperatureIncrement(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw);

  void
  fillgTemperatureRate(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset* matls,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw);

private:
  const MaterialManagerP d_mat_manager{ nullptr };
  const MPMLabel* d_mpm_labels{ nullptr };
  const MPMFlags* d_mpm_flags{ nullptr };
  const PatchSet* d_perproc_patches{ nullptr };

  MaterialSubset* d_one_matl{ nullptr };
  std::unique_ptr<Solver> d_HC_solver{ nullptr };

  bool do_IHC{ false };
  bool d_HC_transient{ false };

  int d_num_ghost_particles{ 1 };
  int d_num_ghost_nodes{ 1 };

  void
  findNeighbors(IntVector n, std::vector<int>& neigh, Array3<int>& l2g);

  inline bool
  compare(double num1, double num2)
  {
    double EPSILON = 1.e-16;
    return (std::abs(num1 - num2) <= EPSILON);
  };
};

} // end namespace Uintah
#endif
