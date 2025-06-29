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


#ifndef Packages_Uintah_CCA_Components_Models_DDT0_h
#define Packages_Uintah_CCA_Components_Models_DDT0_h

#include <CCA/Components/Models/HEChem/HEChemModel.h>

#include <Core/Grid/Variables/ComputeSet.h>

namespace Uintah {
  class ICELabel;
  class MPMLabel;
  class MPMICELabel;
/**************************************

CLASS
   DDT0
  

GENERAL INFORMATION

   DDT0.h

   Joseph Peterson
   Department of Chemistry
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
  
   
KEYWORDS
   "DDT0" "JWL++" "Reactive Flow Model" "Simple Burn"

DESCRIPTION
   Model for deflagration to detonation transition using a simple Arrhenius rate
   equation for deflagration and the detonation model based on "JWL++:  A Simple
   Reactive Flow Code Package for Detonation", P.Clark Souers, Steve Anderson,
   James Mercer, Estella McGuire and Peter Vitello, Propellants,
   Explosives, Pyrotechnics, 25, 54-58, 2000.
  
WARNING

****************************************/

  class DDT0 : public HEChemModel {
  public:
    DDT0(const ProcessorGroup* myworld,
         const MaterialManagerP& materialManager,
         const ProblemSpecP& params,
         const ProblemSpecP& prob_spec);

    virtual ~DDT0();

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

    void initialize(const ProcessorGroup*,
                                const PatchSubset*,
                                const MaterialSubset*,
                                DataWarehouse*,
                                DataWarehouse*);

    void computeModelSources(const ProcessorGroup*, 
                             const PatchSubset* patches,
                             const MaterialSubset* matls,
                             DataWarehouse*,
                             DataWarehouse* new_dw);
      
    DDT0(const DDT0&);
    DDT0& operator=(const DDT0&);

      
    // Simple Burn
    const VarLabel* onSurfaceLabel; // diagnostic labels
    const VarLabel* surfaceTempLabel;
    const VarLabel* totalMassBurnedLabel;
    const VarLabel* totalHeatReleasedLabel;
    const VarLabel* burningLabel;   
    // JWL++
    const VarLabel* reactedFractionLabel;   // diagnostic labels
    const VarLabel* delFLabel;   // diagnostic labels
    const VarLabel* totalMassConvertedLabel;
    const VarLabel* detonatingLabel;

    const VarLabel* pCrackRadiusLabel {nullptr};
    
    ProblemSpecP d_params {nullptr};
    ProblemSpecP d_prob_spec {nullptr};

    const Material* d_matl0 {nullptr};
    const Material* d_matl1 {nullptr};

    ICELabel* Ilb {nullptr};
    MPMICELabel* MIlb {nullptr};
    MPMLabel* Mlb {nullptr};
    MaterialSet* d_mymatls {nullptr};
    MaterialSubset* d_one_matl {nullptr};   

    std::string fromMaterial, toMaterial;
    double d_G;
    double d_b;
    double d_E0;
    double d_threshold_press_JWL;    // JWL++
    double d_threshold_volFrac;
    double d_thresholdTemp;          // Simple Burn
    double d_thresholdPress_SB;      
    double d_Enthalpy;
    double d_BurnCoeff;
    double d_refPress;
    bool d_useCrackModel;
    double d_Gcrack;            // Crack Burning Rate Constant
    double d_nCrack;            // Crack Burning Pressure Exponent
    double d_crackVolThreshold; // for cracking

      
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
