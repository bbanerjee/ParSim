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
  */
  ////////////////////////////////////////////////////////////////////////////

  class InternalVar_MasonSand : public InternalVariableModel {

  public:

    // Internal variables
    const Uintah::VarLabel* pKappaLabel;  // Branch point
    const Uintah::VarLabel* pKappaLabel_preReloc; 

    const Uintah::VarLabel* pCapXLabel;   // Hydrostatic strength
    const Uintah::VarLabel* pCapXLabel_preReloc; 

    const Uintah::VarLabel* pPorosityLabel;     // Porosity
    const Uintah::VarLabel* pPorosityLabel_preReloc; 

    const Uintah::VarLabel* pSaturationLabel;     // Porosity
    const Uintah::VarLabel* pSaturationLabel_preReloc; 

  private:

    // Crush Curve Model parameters
    struct CrushParameters {
      double p0;
      double p1;
      double p2;
      double p3;
    };

    // Initial porosity and saturation parameters
    struct FluidEffectParameters {
      double phi0;  // initial porosity
      double Sw0;   // initial water saturation
    };

    CrushParameters d_crushParam;
    FluidEffectParameters d_fluidParam;

    // Prevent copying of this class
    // copy constructor
    //InternalVar_MasonSand(const InternalVar_MasonSand &cm);
    InternalVar_MasonSand& operator=(const InternalVar_MasonSand &cm);

  public:
    // constructors
    InternalVar_MasonSand(Uintah::ProblemSpecP& ps,
                          ElasticModuliModel* elastic);
    InternalVar_MasonSand(const InternalVar_MasonSand* cm);
         
    // destructor 
    virtual ~InternalVar_MasonSand();

    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps);
         
    /*! Get parameters */
    ParameterDict getParameters() const {
      ParameterDict params;
      params["p0"] = d_crushParam.p0;
      params["p1"] = d_crushParam.p1;
      params["p2"] = d_crushParam.p2;
      params["p3"] = d_crushParam.p3;
      params["phi0"] = d_fluidParam.phi0;
      params["Sw0"] = d_fluidParam.Sw0;
      return params;
    }

    // Computes and requires for internal evolution variables
    virtual void addInitialComputesAndRequires(Uintah::Task* task,
                                               const Uintah::MPMMaterial* matl,
                                               const Uintah::PatchSet* patches);

    virtual void initializeInternalVariable(Uintah::ParticleSubset* pset,
                                            Uintah::DataWarehouse* new_dw) {}

    virtual void initializeInternalVariable(Uintah::ParticleSubset* pset,
                                            Uintah::DataWarehouse* new_dw,
                                            ParameterDict& params);

    virtual void addComputesAndRequires(Uintah::Task* task,
                                        const Uintah::MPMMaterial* matl,
                                        const Uintah::PatchSet* patches);

    virtual void allocateCMDataAddRequires(Uintah::Task* task, 
                                           const Uintah::MPMMaterial* matl,
                                           const Uintah::PatchSet* patch, 
                                           Uintah::MPMLabel* lb);

    virtual void allocateCMDataAdd(Uintah::DataWarehouse* new_dw,
                                   Uintah::ParticleSubset* addset,
                                   Uintah::ParticleLabelVariableMap* newState,
                                   Uintah::ParticleSubset* delset,
                                   Uintah::DataWarehouse* old_dw);

    virtual void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                                  std::vector<const Uintah::VarLabel*>& to);

    virtual void getInternalVariable(Uintah::ParticleSubset* pset,
                                     Uintah::DataWarehouse* old_dw,
                                     Uintah::constParticleLabelVariableMap& intvar);

    virtual void allocateAndPutInternalVariable(Uintah::ParticleSubset* pset,
                                                Uintah::DataWarehouse* new_dw,
                                                Uintah::ParticleLabelVariableMap& intvar); 

    virtual void allocateAndPutRigid(Uintah::ParticleSubset* pset, 
                                     Uintah::DataWarehouse* new_dw,
                                     Uintah::constParticleLabelVariableMap& intvar);

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

 private:

    double computeX(const double& ev_p, 
                    const double& I1,
                    const double& phi,
                    const double& Sw,
                    ParameterDict& params);
    double elasticVolStrainYield(const double& ev_p_bar,
                                 ParameterDict& params);
    double crushCurveDrainedSandX(const double& ev_p_bar) ;
 };

} // End namespace Uintah

#endif  // __MASON_SAND_INT_VAR_MODEL_H__ 
