/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015 Parresia Research Limited, New Zealand
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

#ifndef __INTERNAL_VARIABLE_MODEL_H__
#define __INTERNAL_VARIABLE_MODEL_H__


#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelStateBase.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Task.h>
#include <Core/Math/Matrix3.h>
#include <vector>
#include <map>


namespace Vaango {

  ///////////////////////////////////////////////////////////////////////////
  /*!
    \class  InternalVariableModel
    \brief  Abstract base class for the evolution of internal variables
            for plasticity models
  */
  ///////////////////////////////////////////////////////////////////////////

  class InternalVariableModel {

  private:

  public:
         
    InternalVariableModel();
    virtual ~InternalVariableModel();

    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps) = 0;
         
    /////////////////////////////////////////////////////////////////////////
    /*!
      \brief Get the model parameters
     */
    /////////////////////////////////////////////////////////////////////////
    virtual std::map<std::string, double> getParameters() const = 0 ;

    // Initial computes and requires for internal evolution variables
    virtual void addInitialComputesAndRequires(Uintah::Task* task,
                                               const Uintah::MPMMaterial* matl,
                                               const Uintah::PatchSet* patches) {};

    // Actual initialization
    virtual void initializeInternalVariable(Uintah::ParticleSubset* pset,
                                            Uintah::DataWarehouse* new_dw){};

    // Sometimes extra parameters from the YieldCondition/Modulus models may
    // need to be passed
    virtual void initializeInternalVariable(Uintah::ParticleSubset* pset,
                                            Uintah::DataWarehouse* new_dw,
                                            std::map<std::string, double>& params) {}

    // Computes and requires for internal evolution variables
    virtual void addComputesAndRequires(Uintah::Task* task,
                                        const Uintah::MPMMaterial* matl,
                                        const Uintah::PatchSet* patches) {};

    virtual void allocateCMDataAddRequires(Uintah::Task* task, 
                                           const Uintah::MPMMaterial* matl,
                                           const Uintah::PatchSet* patch, 
                                           Uintah::MPMLabel* lb){};

    virtual void allocateCMDataAdd(Uintah::DataWarehouse* new_dw,
                                   Uintah::ParticleSubset* addset,
                                   Uintah::ParticleLabelVariableMap* newstate,
                                   Uintah::ParticleSubset* delset,
                                   Uintah::DataWarehouse* old_dw){};

    virtual void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                                  std::vector<const Uintah::VarLabel*>& to){};

    /* If there is only one internal variable */
    virtual void getInternalVariable(Uintah::ParticleSubset* pset,
                                     Uintah::DataWarehouse* old_dw,
                                     Uintah::constParticleVariableBase& intvar){};

    /* If there is only one internal variable */
    virtual void allocateAndPutInternalVariable(Uintah::ParticleSubset* pset,
                                                Uintah::DataWarehouse* new_dw,
                                                Uintah::ParticleVariableBase& intvar){}; 

    /* If there is only one internal variable */
    virtual void allocateAndPutRigid(Uintah::ParticleSubset* pset,
                                     Uintah::DataWarehouse* new_dw,
                                     Uintah::constParticleVariableBase& intvar){}; 

    /* For more than one internal variable */
    virtual void getInternalVariable(Uintah::ParticleSubset* pset,
                                     Uintah::DataWarehouse* old_dw,
                                     Uintah::constParticleLabelVariableMap& vars){};

    /* For more than one internal variable */
    virtual void allocateAndPutInternalVariable(Uintah::ParticleSubset* pset,
                                                Uintah::DataWarehouse* new_dw,
                                                Uintah::ParticleLabelVariableMap& vars){}; 

    /* For more than one internal variable */
    virtual void allocateAndPutRigid(Uintah::ParticleSubset* pset,
                                     Uintah::DataWarehouse* new_dw,
                                     Uintah::constParticleLabelVariableMap& intvars){}; 

    ///////////////////////////////////////////////////////////////////////////
    /*! \brief Compute the internal variable and return new value  */
    virtual
    double computeInternalVariable(const ModelStateBase* state) const = 0;

    ///////////////////////////////////////////////////////////////////////////
    // Compute derivative of internal variable with respect to volumetric
    // elastic strain
    virtual
    double computeVolStrainDerivOfInternalVariable(const ModelStateBase* state) const = 0;

  };
} // End namespace Uintah
      


#endif  // __INTERNAL_VARIABLE_MODEL_H__

