/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
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


#ifndef Packages_Uintah_CCA_Components_Examples_TestModel_h
#define Packages_Uintah_CCA_Components_Examples_TestModel_h

#include <CCA/Components/Models/FluidsBased/FluidsBasedModel.h>

#include <Core/Grid/Variables/ComputeSet.h>

namespace Uintah {
  class MPMICELabel;

/**************************************

CLASS
   TestModel
   
   TestModel simulation

GENERAL INFORMATION

   TestModel.h

   Steven G. Parker
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
  
   
KEYWORDS
   TestModel

DESCRIPTION
   Long description...
  
WARNING
  
****************************************/
  class ICELabel;
  
  class TestModel : public FluidsBasedModel {
  public:
    TestModel(const ProcessorGroup* myworld,
              const MaterialManagerP& materialManager,
              const ProblemSpecP& params);
    
    virtual ~TestModel();

    virtual void outputProblemSpec(ProblemSpecP& ps);

    virtual void problemSetup(GridP& grid,
                               const bool isRestart);
      
    virtual void scheduleInitialize(SchedulerP&,
                                    const LevelP& level);

    virtual void scheduleRestartInitialize(SchedulerP&,
                                           [[maybe_unused]] const LevelP& level){};
      
    virtual void scheduleComputeStableTimestep(SchedulerP&,
                                               const LevelP& level);
      
    virtual void scheduleComputeModelSources(SchedulerP&,
                                             const LevelP& level);
                                             
    virtual void scheduleModifyThermoTransportProperties(SchedulerP&,
                                                         const LevelP&,
                                                         const MaterialSet*);
                                               
    virtual void computeSpecificHeat(CCVariable<double>&,
                                    const Patch*,
                                    DataWarehouse*,
                                    const int);
                                    
   virtual void scheduleErrorEstimate(const LevelP& coarseLevel,
                                      SchedulerP& sched);
                                      
   virtual void scheduleTestConservation(SchedulerP&,
                                         const PatchSet* patches);

  private:    
    void computeModelSources(const ProcessorGroup*, 
                             const PatchSubset* patches,
                             const MaterialSubset* matls, 
                             DataWarehouse*, 
                             DataWarehouse* new_dw);

    TestModel(const TestModel&);
    TestModel& operator=(const TestModel&);

    ProblemSpecP d_params;
    const Material* matl0;
    const Material* matl1;
    ICELabel* Ilb;
    MPMICELabel* MIlb;
    MaterialSet* mymatls;
    Material* d_matl;
    double d_rate;
    double d_startTime;   // time to start converting m0->m1
    bool d_is_mpm_matl;  // Is matl 0 a mpm_matl?
    
    const VarLabel* totalMassXLabel;
    const VarLabel* totalIntEngXLabel;
  };
}

#endif
