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


#ifndef Packages_Uintah_CCA_Components_Examples_Simple_Burn_h
#define Packages_Uintah_CCA_Components_Examples_Simple_Burn_h

#include <CCA/Components/Models/HEChem/HEChemModel.h>

#include <Core/Grid/Variables/ComputeSet.h>

namespace Uintah {
  class ICELabel;
  class MPMLabel;
  class MPMICELabel;
/**************************************

CLASS
   Simple_Burn
  

GENERAL INFORMATION

   Simple_Burn.h

   Steven G. Parker
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
  
   
KEYWORDS
   Simple_Burn

DESCRIPTION
   Long description...
  
WARNING
  
****************************************/

  class Simple_Burn : public HEChemModel {
  public:
    Simple_Burn(const ProcessorGroup* myworld,
                const MaterialManagerP& materialManager,
                const ProblemSpecP& params,
                const ProblemSpecP& prob_spec);
    
    virtual ~Simple_Burn();

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
                                             
  private:    
    void computeModelSources(const ProcessorGroup*, 
                             const PatchSubset* patches,
                             const MaterialSubset* matls, 
                             DataWarehouse*, 
                             DataWarehouse* new_dw);
    
    Simple_Burn(const Simple_Burn&);
    Simple_Burn& operator=(const Simple_Burn&);

    const VarLabel* onSurfaceLabel;   // diagnostic labels
    const VarLabel* surfaceTempLabel;
    const VarLabel* totalMassBurnedLabel;
    const VarLabel* totalHeatReleasedLabel;

    ProblemSpecP d_params;
    ProblemSpecP d_prob_spec;
    const Material* matl0;
    const Material* matl1;
    
    MPMICELabel* MIlb;
    ICELabel* Ilb;
    MPMLabel* Mlb;
    MaterialSet* mymatls;
    
    double d_thresholdTemp;
    double d_thresholdPress;
    double d_Enthalpy;
    double d_BurnCoeff;
    double d_refPress;
          
    bool d_is_mpm_matl;  // Is matl 0 a mpm_matl?
    double d_cv_0;      //specific heat
    // flags for the conservation test
    
    struct saveConservedVars{
      bool onOff;
      bool mass;
      bool energy;
     };
     saveConservedVars* d_saveConservedVars;
  };
}

#endif
