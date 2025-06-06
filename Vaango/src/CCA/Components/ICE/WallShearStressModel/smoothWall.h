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

#ifndef _SMOOTHWALL_H_
#define _SMOOTHWALL_H_

#include <CCA/Components/ICE/WallShearStressModel/WallShearStress.h>

namespace Uintah {

  class smoothwall : public WallShearStress{

  public:
  
    smoothwall(ProblemSpecP& ps, MaterialManagerP& materialManager);
    
    virtual ~smoothwall();
    
    virtual void sched_Initialize(SchedulerP& sched, 
                                  const LevelP& level,
                                  const MaterialSet* matls);
                                    
    virtual void sched_AddComputeRequires(Task* task, 
                                          const MaterialSubset* matls);

    virtual
    void computeWallShearStresses( DataWarehouse* old_dw,
                                   DataWarehouse* new_dw,
                                   const Patch* patch,
                                   const int indx,
                                   constCCVariable<double>& vol_frac_CC,  
                                   constCCVariable<Vector>& vel_CC,      
                                   const CCVariable<double>& viscosity,
                                   constCCVariable<double>& rho_CC,      
                                   SFCXVariable<Vector>& tau_X_FC,
                                   SFCYVariable<Vector>& tau_Y_FC,
                                   SFCZVariable<Vector>& tau_Z_FC );
    private:
    
    //__________________________________
    //
    template<class T> 
    void wallShearStresses(DataWarehouse* old_dw,
                           DataWarehouse* new_dw,
                           const Patch* patch,
                           const int indx,
                           const CCVariable<double>& viscosity,
                           constCCVariable<double>& vol_frac_CC,
                           constCCVariable<double>& rho_CC,
                           constCCVariable<Vector>& vel_CC,
                           T& Tau_FC);
                           
    void Initialize(const ProcessorGroup*, 
                    const PatchSubset* ,
                    const MaterialSubset* ,
                    DataWarehouse*, 
                    DataWarehouse*);
      
      MaterialManagerP d_materialManager;
      Patch::FaceType d_face;
      
      double d_invVonKarman;        // 1/VonKarman 
      double d_B_const;
      double d_convergence_uTau;    // convergence criteria
      double d_uTau_guess;         
      unsigned int d_max_iter;
      
    };// End class

}// End namespace Uintah
#endif
