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

#ifndef Packages_Uintah_CCA_Components_Examples_Poisson4_h
#define Packages_Uintah_CCA_Components_Examples_Poisson4_h

#include <CCA/Components/SimulationCommon/SimulationCommon.h>

#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/VarLabel.h>

namespace Uintah {

class EmptyMaterial;

/**************************************
CLASS
   Poisson4
   Poisson4 simulation
****************************************/

class Poisson4 : public SimulationCommon
{
public:
  Poisson4(const ProcessorGroup* myworld, const MaterialManagerP& mat_manager);
  virtual ~Poisson4();

  Poisson4(const Poisson4&) = delete;
  Poisson4(Poisson4&&) = delete;
  Poisson4&
  operator=(const Poisson4&) = delete;
  Poisson4&
  operator=(Poisson4&&) = delete;

  virtual void
  problemSetup(const ProblemSpecP& params,
               const ProblemSpecP& restart_prob_spec,
               GridP& grid,
               const std::string& input_ups_dir = "");

  virtual void
  scheduleInitialize(const LevelP& level, SchedulerP& sched);

  virtual void
  scheduleComputeStableTimestep(const LevelP& level, SchedulerP&);

  virtual void
  scheduleTimeAdvance(const LevelP& level, SchedulerP&);

  virtual void
  scheduleTimeAdvance0(SchedulerP&, const PatchSet*, const MaterialSet*);

  virtual void
  scheduleTimeAdvance1(SchedulerP&, const PatchSet*, const MaterialSet*);

  virtual void
  scheduleRestartInitialize([[maybe_unused]] const LevelP& level, [[maybe_unused]] SchedulerP& sched)
  {
  }

  virtual void
  outputProblemSpec(Uintah::ProblemSpecP&)
  {
  }

private:
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
  timeAdvance(const ProcessorGroup*,
              const PatchSubset* patches,
              const MaterialSubset* matls,
              DataWarehouse* old_dw,
              DataWarehouse* new_dw);

  void
  timeAdvance1(const ProcessorGroup*,
               const PatchSubset* patches,
               const MaterialSubset* matls,
               DataWarehouse* old_dw,
               DataWarehouse* new_dw);

  double d_delT;
  std::shared_ptr<EmptyMaterial> d_mymat;
  const VarLabel* d_phi_label;
  const VarLabel* d_residual_label;

};
}

#endif
