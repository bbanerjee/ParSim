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

#ifndef Packages_Uintah_CCA_Components_Examples_AMRHeat_hpp
#define Packages_Uintah_CCA_Components_Examples_AMRHeat_hpp

#include <CCA/Components/Examples/Heat.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Grid/Variables/NCVariable.h>

namespace Uintah {
  class AMRHeat : public Heat {
  public:
    AMRHeat(const ProcessorGroup* world,
	    const MaterialManagerP materialManager);
    
    virtual ~AMRHeat();

    virtual void
    problemSetup(const ProblemSpecP& ) {}
    virtual void problemSetup(const ProblemSpecP&     ps,
                              const ProblemSpecP&     restart_ps,
                                    GridP&            grid,
                              const std::string& name = "");

    virtual void scheduleRefineInterface(const LevelP&     fineLevel,
                                               SchedulerP& scheduler,
                                               bool        needCoarseOld,
                                               bool        needCoarseNew);

    virtual void scheduleCoarsen(const LevelP&     coarseLevel,
                                       SchedulerP& sched);

    virtual void scheduleRefine (const PatchSet*   patches,
                                       SchedulerP& sched);

    virtual void scheduleErrorEstimate(const LevelP&     coarseLevel,
                                             SchedulerP& sched);

    virtual void scheduleInitialErrorEstimate(const LevelP&     coarseLevel,
                                                    SchedulerP& sched);

  private:
    double d_refine_threshold;
    IntVector d_refinement_ratio;

    void errorEstimate(const ProcessorGroup* pg,
                       const PatchSubset*    patches,
                       const MaterialSubset* matls,
                             DataWarehouse*  old_dw,
                             DataWarehouse*  new_dw);

    double computeError(const IntVector& c,
                        const Patch* patch,
                              constNCVariable<double>& temp);

    void refine(const ProcessorGroup* pg,
                const PatchSubset*    patches,
                const MaterialSubset* matls,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw);

    void refineNode(NCVariable<double>& temp, constNCVariable<double>& coarse_temp,
                    IntVector fine_index,
                    const Level* fine_level, const Level* coarse_level);

    AMRHeat(const AMRHeat&);
    AMRHeat& operator=(const AMRHeat&);

  };
}
#endif // Packages_Uintah_CCA_Components_Examples_AMRHeat_hpp
