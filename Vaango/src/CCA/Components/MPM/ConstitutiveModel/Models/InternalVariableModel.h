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

  using ParameterDict = std::map<std::string, double>;

  class ElasticModuliModel;
  class ShearModulusModel;

  ///////////////////////////////////////////////////////////////////////////
  /*!
    \class  InternalVariableModel
    \brief  Abstract base class for the evolution of internal variables
            for plasticity models
  */
  ///////////////////////////////////////////////////////////////////////////

  class InternalVariableModel {

  protected:

    ElasticModuliModel* d_elastic;
    ShearModulusModel* d_shear;

  public:
         
    InternalVariableModel();
    virtual ~InternalVariableModel();

    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps) = 0;
         
    /*!
      \brief Get the local <double> particle variables
             For more than one internal variable 
     */
    virtual
    std::vector<Uintah::constParticleVariable<double> >
         getInternalVariables(Uintah::ParticleSubset* pset,
                              Uintah::DataWarehouse* old_dw,
                              const double& dummy) {
      Uintah::constParticleVariable<double> pNull;
      std::vector<Uintah::constParticleVariable<double> > pIntVars;
      pIntVars.emplace_back(pNull);
      return pIntVars;
    } 

    /*!
      \brief Get the local <Matrix3> particle variables
             For more than one internal variable 
     */
    virtual
    std::vector<Uintah::constParticleVariable<Uintah::Matrix3> >
         getInternalVariables(Uintah::ParticleSubset* pset,
                              Uintah::DataWarehouse* old_dw,
                              const Uintah::Matrix3& dummy) {
      Uintah::constParticleVariable<Uintah::Matrix3> pNull;
      std::vector<Uintah::constParticleVariable<Uintah::Matrix3> > pIntVars;
      pIntVars.emplace_back(pNull);
      return pIntVars;
    } 

    /*!
      \brief Allocate and put the local <double> particle variables
      For more than one internal variable */
    typedef std::vector<Uintah::ParticleVariable<double>* > vectorParticleDoubleP;
    virtual
    void allocateAndPutInternalVariable(Uintah::ParticleSubset* pset,
                                        Uintah::DataWarehouse* new_dw,
                                        vectorParticleDoubleP& pVars){}

    /*!
      \brief Allocate and put the local <Matrix3> particle variables
      For more than one internal variable */
    typedef std::vector<Uintah::ParticleVariable<Uintah::Matrix3>* > vectorParticleMatrix3P;
    virtual
    void allocateAndPutInternalVariable(Uintah::ParticleSubset* pset,
                                        Uintah::DataWarehouse* new_dw,
                                        vectorParticleMatrix3P& pVars){}

    /*!
      \brief Get the internal variable labels
     */
    virtual std::vector<const Uintah::VarLabel*> getLabels() const  = 0;

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
    virtual void initializeInternalVariable(const Uintah::Patch* patch,
                                            const Uintah::MPMMaterial* matl,
                                            Uintah::ParticleSubset* pset,
                                            Uintah::DataWarehouse* new_dw,
                                            Uintah::MPMLabel* lb,
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

