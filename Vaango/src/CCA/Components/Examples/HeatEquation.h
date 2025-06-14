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

#ifndef Packages_Uintah_CCA_Components_Examples_HeatEquation_h
#define Packages_Uintah_CCA_Components_Examples_HeatEquation_h

#include <CCA/Components/SimulationCommon/SimulationCommon.h>

#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/VarLabel.h>

namespace Uintah {
class EmptyMaterial;
class ExamplesLabel;

/**************************************

CLASS
   HeatEquation

   HeatEquation simulation

GENERAL INFORMATION

   HeatEquation.h

   John Schmidt
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)


KEYWORDS
   HeatEquation

DESCRIPTION
   Long description...

WARNING

****************************************/

class HeatEquation
  : public SimulationCommon
{
public:
  HeatEquation(const ProcessorGroup* myworld, const MaterialManagerP& mat_manager);

  virtual ~HeatEquation();

  virtual void
  problemSetup(const ProblemSpecP& params,
               const ProblemSpecP& restart_prob_spec,
               GridP& grid);

  virtual void
  scheduleInitialize(const LevelP& level, SchedulerP& sched);

  virtual void
  scheduleComputeStableTimestep(const LevelP& level, SchedulerP&);

  virtual void
  scheduleTimeAdvance(const LevelP& level, SchedulerP&);

  virtual void
  scheduleRestartInitialize([[maybe_unused]] const LevelP& level, [[maybe_unused]] SchedulerP& sched)
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
  iterate(const ProcessorGroup*,
          const PatchSubset* patches,
          const MaterialSubset* matls,
          DataWarehouse* old_dw,
          DataWarehouse* new_dw);

  const VarLabel* d_temperature_label;
  const VarLabel* d_residual_label;
  double d_delT;
  double d_maxresidual;
  std::shared_ptr<EmptyMaterial> d_mymat;

  HeatEquation(const HeatEquation&);
  HeatEquation&
  operator=(const HeatEquation&);
};
}

#endif
