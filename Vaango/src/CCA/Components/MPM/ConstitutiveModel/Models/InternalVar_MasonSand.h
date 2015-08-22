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

#ifndef __MASON_SAND_INT_VAR_MODEL_H__
#define __MASON_SAND_INT_VAR_MODEL_H__


#include <CCA/Components/MPM/ConstitutiveModel/Models/InternalVariableModel.h>    
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_MasonSand.h>    
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

  ////////////////////////////////////////////////////////////////////////////
  /*! 
    \class InternalVar_MasonSand
    \brief The evolution of the kappa, X, porosity, and saturation internal variables 
           in the partially saturated Arenisca model

    Reference: Arenisca manual.

    The rate of change of kappa is given by
 
     dkappa/deps_v = F(kappa) - G(eps_v) - H(eps_v)

    where

     eps_v = volumetric plastic strain

    and

     F(kappa) = 1.0/(p1*p3)*exp(-p1*kappa - p0)
     G(eps_v) = B1 exp(p3 + p4 + eps_v)/[exp(p3 + p4 + eps_v) - 1]^2
     H(eps_v) = B1 exp(p3 + eps_v)/[exp(p3 + eps_v) - 1]^2

    with 
     B1 = 3 B0 [exp(p3 + p4) - 1]
      

    The incremental update of the consolidation pressure is given by

       kappa_{n+1} = kappa_n + [F(kappa_{n+1}) - G(eps_v) - H(eps_v)] Delta eps_v
  */
  ////////////////////////////////////////////////////////////////////////////

  class InternalVar_MasonSand : public InternalVariableModel {

  public:

    // Internal variables
    const Uintah::VarLabel* pKappaLabel; 
    const Uintah::VarLabel* pKappaLabel_preReloc; 

  private:

    // Crush Curve Model parameters
    struct CrushParameters {
      double p0;
      double p1;
      double p2;
      double p3;
    };

    CrushParameters d_crushParam;

    // Prevent copying of this class
    // copy constructor
    //InternalVar_MasonSand(const InternalVar_MasonSand &cm);
    InternalVar_MasonSand& operator=(const InternalVar_MasonSand &cm);

  public:
    // constructors
    InternalVar_MasonSand(Uintah::ProblemSpecP& ps);
    InternalVar_MasonSand(const InternalVar_MasonSand* cm);
         
    // destructor 
    virtual ~InternalVar_MasonSand();

    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps);
         
    // Computes and requires for internal evolution variables
    virtual void addInitialComputesAndRequires(Uintah::Task* task,
                                               const Uintah::MPMMaterial* matl,
                                               const Uintah::PatchSet* patches);

    virtual void addComputesAndRequires(Uintah::Task* task,
                                        const Uintah::MPMMaterial* matl,
                                        const Uintah::PatchSet* patches);

    virtual void allocateCMDataAddRequires(Uintah::Task* task, 
                                           const Uintah::MPMMaterial* matl,
                                           const Uintah::PatchSet* patch, 
                                           Uintah::MPMLabel* lb);

    virtual void allocateCMDataAdd(Uintah::DataWarehouse* new_dw,
                                   Uintah::ParticleSubset* addset,
                                   std::map<const Uintah::VarLabel*, 
                                     Uintah::ParticleVariableBase*>* newState,
                                   Uintah::ParticleSubset* delset,
                                   Uintah::DataWarehouse* old_dw);

    virtual void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                                  std::vector<const Uintah::VarLabel*>& to);

    virtual void initializeInternalVariable(Uintah::ParticleSubset* pset,
                                            Uintah::DataWarehouse* new_dw);

    virtual void getInternalVariable(Uintah::ParticleSubset* pset,
                                     Uintah::DataWarehouse* old_dw,
                                     Uintah::constParticleVariableBase& intvar);

    virtual void allocateAndPutInternalVariable(Uintah::ParticleSubset* pset,
                                                Uintah::DataWarehouse* new_dw,
                                                Uintah::ParticleVariableBase& intvar); 

    virtual void allocateAndPutRigid(Uintah::ParticleSubset* pset, 
                                     Uintah::DataWarehouse* new_dw,
                                     Uintah::constParticleVariableBase& intvar);

    ///////////////////////////////////////////////////////////////////////////
    /*! \brief Compute the internal variable */
    double computeInternalVariable(const ModelStateBase* state) const;

    ///////////////////////////////////////////////////////////////////////////
    // Compute derivative of internal variable with respect to volumetric
    // elastic strain
    double computeVolStrainDerivOfInternalVariable(const ModelStateBase*) const
    {
      return 0.0;
    }

 };

} // End namespace Uintah

#endif  // __MASON_SAND_INT_VAR_MODEL_H__ 
