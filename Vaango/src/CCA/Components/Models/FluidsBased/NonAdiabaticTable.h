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


#ifndef Packages_Uintah_CCA_Components_Examples_NonAdiabaticTable_h
#define Packages_Uintah_CCA_Components_Examples_NonAdiabaticTable_h

#include <CCA/Components/Models/FluidsBased/FluidsBasedModel.h>

#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/VarTypes.h>

#include <vector>

namespace Uintah {
  class ICELabel;
  class TableInterface;

/**************************************

CLASS
   NonAdiabaticTable
   
   NonAdiabaticTable simulation

GENERAL INFORMATION

   NonAdiabaticTable.h

   Todd Harman
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
  
   
KEYWORDS
   NonAdiabaticTable

DESCRIPTION
   Long description...
  
WARNING
  
****************************************/

  class NonAdiabaticTable :public FluidsBasedModel {
  public:
    NonAdiabaticTable(const ProcessorGroup* myworld,
                      const MaterialManagerP& materialManager,
                      const ProblemSpecP& params);
    
    virtual ~NonAdiabaticTable();
    
    virtual void problemSetup(GridP& grid,
                               const bool isRestart);
    
    virtual void scheduleInitialize(SchedulerP&,
                                    const LevelP& level);
      
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
    ICELabel* Ilb;
                                                
    void modifyThermoTransportProperties(const ProcessorGroup*, 
                                         const PatchSubset* patches,        
                                         const MaterialSubset*,             
                                         DataWarehouse*,                    
                                         DataWarehouse* new_dw);             
   
    void initialize(const ProcessorGroup*, 
                    const PatchSubset* patches,
                    const MaterialSubset* matls, 
                    DataWarehouse*, 
                    DataWarehouse* new_dw);
                                   
    void computeModelSources(const ProcessorGroup*, 
                             const PatchSubset* patches,
                             const MaterialSubset*,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw);
                             
    void testConservation(const ProcessorGroup*, 
                          const PatchSubset* patches,
                          const MaterialSubset*,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw);
    //__________________________________
    NonAdiabaticTable(const NonAdiabaticTable&);
    NonAdiabaticTable& operator=(const NonAdiabaticTable&);

    ProblemSpecP d_params {nullptr};

    const Material* d_matl;
    MaterialSet* d_matl_set;

    class Region {
    public:
      Region(GeometryPieceP piece, ProblemSpecP&);
      
      GeometryPieceP piece;
      double initialScalar;
    };

    class Scalar {
    public:
      int index;
      std::string name;
      // labels for this particular scalar
      VarLabel* scalar_CCLabel;
      VarLabel* scalar_src_CCLabel;
      VarLabel* diffusionCoeffLabel;
      VarLabel* varianceLabel;
      VarLabel* scaledVarianceLabel;
      VarLabel* sum_scalar_fLabel;
      
      std::vector<Region*> regions;
      double diff_coeff;
    };

    double oldProbeDumpTime;
    Scalar* d_scalar;
    // Release is J/kg
    VarLabel* cumulativeEnergyReleased_CCLabel;
    VarLabel* cumulativeEnergyReleased_src_CCLabel;
    
    std::vector<Vector> d_probePts;
    std::vector<std::string> d_probePtsNames;
    bool d_usingProbePts;
    double d_probeFreq;

    TableInterface* table;
    struct TableValue {
      std::string name;
      int index;
      VarLabel* label;
    };
    std::vector<TableValue*> tablevalues;

    int d_density_index;
    int d_gamma_index;
    int d_cv_index;
    int d_viscosity_index;
    int d_thermalcond_index;
    int d_temp_index;
    int d_ref_cv_index;
    int d_ref_gamma_index;
    int d_ref_temp_index;
    bool d_test_conservation;

    bool useVariance;
    double varianceScale;
    double varianceMax;
  };
}

#endif
