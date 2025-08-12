/*
 * The MIT License
 *
 * Copyright (c) 1997-2024 The University of Utah
 * Copyright (c) 2024-2025 Biswajit Banerjee, Parresia Research Limited, NZ
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

#ifndef Packages_Uintah_CCA_Components_Examples_Heat_hpp
#define Packages_Uintah_CCA_Components_Examples_Heat_hpp

#include <CCA/Components/SimulationCommon/SimulationCommon.h>

#include <Core/Grid/EmptyMaterial.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <CCA/Components/Examples/ExamplesLabel.h>

namespace Uintah{
  class Heat : public SimulationCommon {
  public:
    Heat(const ProcessorGroup* myworld,
	 const MaterialManagerP materialManager);
    
    virtual ~Heat();

    virtual void
    problemSetup(const ProblemSpecP& ) {};
    virtual void problemSetup(const ProblemSpecP&     ps,
                              const ProblemSpecP&     restart_ps,
                                    GridP&            grid,
                              const std::string& restart_file = "");

    virtual void scheduleInitialize(const LevelP&     level,
                                          SchedulerP& sched);

    virtual void scheduleRestartInitialize(const LevelP&     level,
                                                 SchedulerP& sched);

    virtual void scheduleComputeStableTimeStep(const LevelP&     level,
                                                     SchedulerP& sched);

    virtual void scheduleTimeAdvance(const LevelP&     level,
                                           SchedulerP& sched);

  protected:
    ExamplesLabel* d_lb;
    std::shared_ptr<EmptyMaterial> d_mat;
    double d_delt, d_alpha, d_r0, d_gamma;

  private:
    virtual void initialize(const ProcessorGroup* pg,
                            const PatchSubset*    patches,
                            const MaterialSubset* matls,
                                  DataWarehouse*  old_dw,
                                  DataWarehouse*  new_dw);

    virtual void computeStableTimeStep(const ProcessorGroup* pg,
                                       const PatchSubset*    patches,
                                       const MaterialSubset* matls,
                                             DataWarehouse*  old_dw,
                                             DataWarehouse*  new_dw);

    virtual void timeAdvance(const ProcessorGroup* pg,
                             const PatchSubset*    patches,
                             const MaterialSubset* matls,
                                   DataWarehouse*  old_dw,
                                   DataWarehouse*  new_dw);



    Heat(const Heat&);
    Heat& operator=(const Heat&);
  };
}

#endif // End Packages_Uintah_CCA_Components_Examples_Heat_hpp
