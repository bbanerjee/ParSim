/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#ifndef Packages_Uintah_CCA_Components_Examples_Wave_h
#define Packages_Uintah_CCA_Components_Examples_Wave_h

#include <CCA/Components/SimulationCommon/SimulationCommon.h>

#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/ComputeSet.h>

namespace Uintah {

class EmptyMaterial;
class ExamplesLabel;
class VarLabel;

/**************************************
CLASS
   Wave
   Wave simulation
****************************************/

class Wave : public SimulationCommon
{
public:
  Wave(const ProcessorGroup* myworld, const MaterialManagerP& mat_manager);

  virtual ~Wave();

  Wave(const Wave&) = delete;
  Wave(Wave&&)      = delete;
  Wave&
  operator=(const Wave&) = delete;
  Wave&
  operator=(Wave&&) = delete;

  virtual void
  problemSetup(const ProblemSpecP& params,
               const ProblemSpecP& restart_prob_spec,
               GridP& grid,
               const std::string& input_ups_dir = "");

  virtual void
  scheduleInitialize(const LevelP& level, SchedulerP& sched);

  virtual void
  scheduleRestartInitialize([[maybe_unused]] const LevelP& level, [[maybe_unused]] SchedulerP& sched) override
  {
  }

  virtual void
  scheduleComputeStableTimestep(const LevelP& level, SchedulerP&);

  virtual void
  scheduleTimeAdvance(const LevelP& level, SchedulerP&);

  virtual void
  outputProblemSpec(Uintah::ProblemSpecP&) override
  {
  }

protected:
  struct Step
  {
    Task::WhichDW cur_dw;
    const VarLabel* curphi_label;
    const VarLabel* curpi_label;
    const VarLabel* newphi_label;
    const VarLabel* newpi_label;
    double stepweight;
    double totalweight;
  };

  inline double
  curl(constCCVariable<double>& phi,
       const IntVector& c,
       double sumdx2,
       Vector inv_dx2)
  {
    return sumdx2 * phi[c] +
           (phi[c + IntVector(1, 0, 0)] + phi[c - IntVector(1, 0, 0)]) *
             inv_dx2.x() +
           (phi[c + IntVector(0, 1, 0)] + phi[c - IntVector(0, 1, 0)]) *
             inv_dx2.y() +
           (phi[c + IntVector(0, 0, 1)] + phi[c - IntVector(0, 0, 1)]) *
             inv_dx2.z();
  }

  void
  initialize(const ProcessorGroup*,
             const PatchSubset* patches,
             const MaterialSubset* matls,
             DataWarehouse* old_dw,
             DataWarehouse* new_dw);

  void
  computeStableTimestep(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset* matls,
                        DataWarehouse* old_dw,
                        DataWarehouse* new_dw);
  void
  timeAdvanceEuler(const ProcessorGroup*,
                   const PatchSubset* patches,
                   const MaterialSubset* matls,
                   DataWarehouse* old_dw,
                   DataWarehouse* new_dw);

  void
  setupRK4(const ProcessorGroup*,
           const PatchSubset* patches,
           const MaterialSubset* matls,
           DataWarehouse* old_dw,
           DataWarehouse* new_dw);

  void
  timeAdvanceRK4(const ProcessorGroup*,
                 const PatchSubset* patches,
                 const MaterialSubset* matls,
                 DataWarehouse* old_dw,
                 DataWarehouse* new_dw,
                 Step* s);

  virtual void
  addRefineDependencies(Task* /*task*/,
                        const VarLabel* /*label*/,
                        [[maybe_unused]] bool needCoarseOld,
                        [[maybe_unused]] bool needCoarseNew)
  {
  }

  virtual void
  refineFaces(const Patch* /*finePatch*/,
              const Level* /*fineLevel*/,
              const Level* /*coarseLevel*/,
              CCVariable<double>& /*finevar*/,
              const VarLabel* /*label*/,
              int /*matl*/,
              DataWarehouse* /*coarse_old_dw*/,
              DataWarehouse* /*coarse_new_dw*/)
  {
  }

  const VarLabel* d_phi_label;
  const VarLabel* d_pi_label;
  double d_r0;
  string d_initial_condition;
  string d_integration;

  EmptyMaterial* d_mymat;
  Step d_rk4steps[4];
};

} // namespace Uintah

#endif
