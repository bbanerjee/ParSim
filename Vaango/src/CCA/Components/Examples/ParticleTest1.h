/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#ifndef Packages_Uintah_CCA_Components_Examples_ParticleTest1_h
#define Packages_Uintah_CCA_Components_Examples_ParticleTest1_h

#include <CCA/Components/SimulationCommon/SimulationCommon.h>

#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/ParticleVariable.h>

namespace Uintah {
class EmptyMaterial;
class ExamplesLabel;

/**************************************

CLASS
   ParticleTest1

   ParticleTest1 simulation

GENERAL INFORMATION

   ParticleTest1.h

   Steven G. Parker
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)


KEYWORDS
   ParticleTest1

DESCRIPTION
   Long description...

WARNING

****************************************/

class ParticleTest1
  : public SimulationCommon
{
public:
  ParticleTest1(const ProcessorGroup* myworld, const MaterialManagerP& mat_manager);
  virtual ~ParticleTest1();

  ParticleTest1(const ParticleTest1&) = delete;
  ParticleTest1(ParticleTest1&&) = delete;
  ParticleTest1&
  operator=(const ParticleTest1&) = delete;
  ParticleTest1&
  operator=(ParticleTest1&&) = delete;

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
  scheduleRestartInitialize(const LevelP& level, SchedulerP& sched)
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

  ExamplesLabel* d_labels;
  double d_delT;
  std::shared_ptr<EmptyMaterial> d_mymat;
  int d_doOutput;
  int d_doGhostCells;
};
}

#endif
