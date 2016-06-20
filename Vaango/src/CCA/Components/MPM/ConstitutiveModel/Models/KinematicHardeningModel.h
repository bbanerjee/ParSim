/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2016 Parresia Research Limited, New Zealand
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

#ifndef __BB_KINEMATIC_HARDENING_MODEL_H__
#define __BB_KINEMATIC_HARDENING_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelStateBase.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/InternalVariableModel.h>

#include <Core/Math/Matrix3.h>
#include <vector>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/VarLabel.h>



namespace Vaango {

  ///////////////////////////////////////////////////////////////////////////
  /*!
    \class  KinematicHardeningModel
    \brief  Abstract Base class for kinematic hardening models 
    \author Biswajit Banerjee, \n
            C-SAFE and Department of Mechanical Engineering, \n
            University of Utah,\n
    \warn   Assumes vonMises yield condition and the associated flow rule for 
            all cases other than Gurson plasticity.
  */
  ///////////////////////////////////////////////////////////////////////////

  class KinematicHardeningModel {

  protected:

    InternalVariableModel* d_intvar;

  public:

    const Uintah::VarLabel* pBackStressLabel;
    const Uintah::VarLabel* pBackStressLabel_preReloc;

  public:
         
    KinematicHardeningModel();
    virtual ~KinematicHardeningModel();

    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps) = 0;
         
    /////////////////////////////////////////////////////////////////////////
    /*!
      \brief Get the model parameters
     */
    /////////////////////////////////////////////////////////////////////////
    virtual std::map<std::string, double> getParameters() const = 0 ;

    //////////
    /*! \brief Calculate the back stress */
    /* Note that df_dsigma_normal_new is the normalized value of df_dsigma */
    //////////
    virtual
    void computeBackStress(const ModelStateBase* state,
                           const double& delT,
                           const Uintah::particleIndex idx,
                           const double& delLambda,
                           const Uintah::Matrix3& df_dsigma_normal_new,
                           const Uintah::Matrix3& backStress_old,
                           Uintah::Matrix3& backStress_new) = 0;
 
    // virtual
    // void computeBackStress(const ModelStateBase* state,
    //                        Uintah::Matrix3& backStress_new) {}
    template <typename ModelState>
    void computeBackStress(const ModelState* state,
                           Uintah::Matrix3& backStress_new) {}
 
    /*! Compute the direction of back stress evolution (\f$h^beta\f$) 
        for the equation \f$ d/dt(\beta) = d/dt(\gamma) h^beta \f$ */
    virtual
    void eval_h_beta(const Uintah::Matrix3& df_dsigma,
                     const ModelStateBase* state,
                     Uintah::Matrix3& h_beta) = 0;

    /*! Template test */
    template <typename ModelState>
    void evalState(const ModelState* state);

    /*!  Data management apparatus */
    virtual void addInitialComputesAndRequires(Uintah::Task* task,
                                               const Uintah::MPMMaterial* matl,
                                               const Uintah::PatchSet* patches) const;

    virtual void addComputesAndRequires(Uintah::Task* task,
                                        const Uintah::MPMMaterial* matl,
                                        const Uintah::PatchSet* patches) const;

    virtual void addComputesAndRequires(Uintah::Task* task,
                                        const Uintah::MPMMaterial* matl,
                                        const Uintah::PatchSet* patches,
                                        bool recurse) const;


    virtual void allocateCMDataAddRequires(Uintah::Task* task, 
                                           const Uintah::MPMMaterial* matl,
                                           const Uintah::PatchSet* patch, 
                                           Uintah::MPMLabel* lb) const;

    virtual void allocateCMDataAdd(Uintah::DataWarehouse* new_dw,
                                   Uintah::ParticleSubset* addset,
                                   Uintah::map<const Uintah::VarLabel*, 
                                     Uintah::ParticleVariableBase*>* newState,
                                   Uintah::ParticleSubset* delset,
                                   Uintah::DataWarehouse* old_dw);

    virtual void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                                  std::vector<const Uintah::VarLabel*>& to);

    virtual void initializeBackStress(Uintah::ParticleSubset* pset,
                                      Uintah::DataWarehouse* new_dw);

    virtual void initializeLocalVariables(const Uintah::Patch* patch,
                                          Uintah::ParticleSubset* pset,
                                          Uintah::DataWarehouse* new_dw,
                                          Uintah::constParticleVariable<double>& pVolume) {};

    virtual void getBackStress(Uintah::ParticleSubset* pset,
                               Uintah::DataWarehouse* old_dw,
                               Uintah::constParticleVariable<Uintah::Matrix3>& pBackStress);

    virtual void allocateAndPutBackStress(Uintah::ParticleSubset* pset,
                                          Uintah::DataWarehouse* new_dw, 
                                          Uintah::ParticleVariable<Uintah::Matrix3>& pBackStress);

    virtual void allocateAndPutRigid(Uintah::ParticleSubset* pset,
                                     Uintah::DataWarehouse* new_dw); 

  };
} // End namespace Uintah
      


#endif  // __BB_KINEMATIC_HARDENING_MODEL_H__

