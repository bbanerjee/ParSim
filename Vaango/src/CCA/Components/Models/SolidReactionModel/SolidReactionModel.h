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


#ifndef Packages_Uintah_CCA_Components_Models_SolidReactionModel_h
#define Packages_Uintah_CCA_Components_Models_SolidReactionModel_h

#include <CCA/Ports/ModelInterface.h>

namespace Uintah {

  class RateConstant;
  class RateModel;
  
  class ICELabel;

  /**************************************
     
     CLASS
     SolidReactionModel
     
     A generalized Rate Model class that very much borrows
     from the factory model idiom for rate expression composition.
     
     GENERAL INFORMATION
     
     SolidReactionModel.h
     
     Joseph R. Peterson
     Department of Chemistry
     University of Utah
     
     Center for the Simulation of Accidental Fires and Explosions (C-SAFE)     
          
     KEYWORDS
     SolidReactionModel
     
     DESCRIPTION
     Long description...
     
     WARNING
     
     ****************************************/
    
    class SolidReactionModel : public ModelInterface {
    public:
      SolidReactionModel(const ProcessorGroup* d_myworld,
                         const MaterialManagerP& materialManager,
                         const ProblemSpecP& params,
                         const ProblemSpecP& prob_spec);
      
      virtual ~SolidReactionModel();
      
      virtual void problemSetup(GridP& grid, const bool isRestart);
      
      virtual void outputProblemSpec(ProblemSpecP& ps);
        
      virtual void scheduleInitialize(SchedulerP&,
                                      const LevelP& level);

      virtual void scheduleRestartInitialize([[maybe_unused]] SchedulerP&,
                                             [[maybe_unused]] const LevelP& level){};
      
      virtual void scheduleComputeStableTimestep(SchedulerP& sched,
                                                 const LevelP& level);
      
      virtual void scheduleComputeModelSources(SchedulerP&,
                                                 const LevelP& level);
                
    private:

      void computeModelSources(const ProcessorGroup*,
                                 const PatchSubset* patches,
                                 const MaterialSubset* matls,
                                 DataWarehouse*,
                                 DataWarehouse* new_dw);

      // Functions
      SolidReactionModel(const SolidReactionModel&);
      SolidReactionModel& operator=(const SolidReactionModel&);
      
      // Innards
      RateConstant *rateConstant {nullptr};  // k(T)
      RateModel    *rateModel {nullptr};     // f(a)
      const Material* reactant {nullptr};
      const Material* product {nullptr};
      std::string fromMaterial;
      std::string toMaterial;
      double d_E0;                 // Enthalpy change for reaction in J/kg
       
      ICELabel *Ilb;               // Used to get handles on temperature, pressure, etc.
      MaterialSet *mymatls;        // All the materials referenced by this model

      // Variables used for tracking the Reaction
      const VarLabel* reactedFractionLabel;   // Fraction of reactant in cell
      const VarLabel* delFLabel;              // Change of fraction of reactant during timestep
      const VarLabel* totalMassBurnedLabel;  
      const VarLabel* totalHeatReleasedLabel;
 
      // flags for the conservation test
      struct saveConservedVars{
        bool onOff;
        bool mass;
        bool energy;
      };

      saveConservedVars* d_saveConservedVars;

      // Some Uintah Necessities
      ProblemSpecP d_params {nullptr};
      ProblemSpecP d_prob_spec {nullptr};
    };
}

#endif
