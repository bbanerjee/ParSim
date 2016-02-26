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
    const Uintah::VarLabel* pKappaLabel;                         // Branch point
    const Uintah::VarLabel* pKappaLabel_preReloc; 

    const Uintah::VarLabel* pCapXLabel;                          // Hydrostatic strength
    const Uintah::VarLabel* pCapXLabel_preReloc; 

    const Uintah::VarLabel* pPlasticStrainLabel;                 // Plastic Strain
    const Uintah::VarLabel* pPlasticStrainLabel_preReloc;

    const Uintah::VarLabel* pPlasticVolStrainLabel;              // Plastic Volumetric Strain
    const Uintah::VarLabel* pPlasticVolStrainLabel_preReloc;
    
    const Uintah::VarLabel* pP3Label;                            // Evolution of parameter P3
    const Uintah::VarLabel* pP3Label_preReloc;

    // Return the internal variable labels
    std::vector<const Uintah::VarLabel*> getLabels() const
    {
       std::vector<const Uintah::VarLabel*> labels;
       labels.push_back(pKappaLabel);                         // Branch point
       labels.push_back(pKappaLabel_preReloc); 

       labels.push_back(pCapXLabel);                          // Hydrostatic strength
       labels.push_back(pCapXLabel_preReloc); 

       labels.push_back(pPlasticStrainLabel);                 // Plastic Strain
       labels.push_back(pPlasticStrainLabel_preReloc);

       labels.push_back(pPlasticVolStrainLabel);              // Plastic Volumetric Strain
       labels.push_back(pPlasticVolStrainLabel_preReloc);
    
       labels.push_back(pP3Label);                            // Evolution of parameter P3
       labels.push_back(pP3Label_preReloc);

       return labels;
    }

  private:

    // Crush Curve Model parameters
    struct CrushParameters {
      double p0;
      double p1;
      double p2;
      double p3;
    };

    CrushParameters d_crushParam;
    bool d_use_disaggregation_algorithm;

    // Prevent copying of this class
    // copy constructor
    //InternalVar_MasonSand(const InternalVar_MasonSand &cm);
    InternalVar_MasonSand& operator=(const InternalVar_MasonSand &cm);

    // Initialize local VarLabels
    void initializeLocalMPMLabels() 
    {
      pKappaLabel = Uintah::VarLabel::create("p.kappa",
        Uintah::ParticleVariable<double>::getTypeDescription());
      pKappaLabel_preReloc = Uintah::VarLabel::create("p.kappa+",
        Uintah::ParticleVariable<double>::getTypeDescription());

      pCapXLabel = Uintah::VarLabel::create("p.capX",
        Uintah::ParticleVariable<double>::getTypeDescription());
      pCapXLabel_preReloc = Uintah::VarLabel::create("p.capX+",
        Uintah::ParticleVariable<double>::getTypeDescription());

      pPlasticStrainLabel = Uintah::VarLabel::create("p.plasticStrain",
        Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());
      pPlasticStrainLabel_preReloc = Uintah::VarLabel::create("p.plasticStrain+",
        Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());

      pPlasticVolStrainLabel = Uintah::VarLabel::create("p.plasticVolStrain",
        Uintah::ParticleVariable<double>::getTypeDescription());
      pPlasticVolStrainLabel_preReloc = Uintah::VarLabel::create("p.plasticVolStrain+",
        Uintah::ParticleVariable<double>::getTypeDescription());

      pP3Label = Uintah::VarLabel::create("p.p3",
        Uintah::ParticleVariable<double>::getTypeDescription());
      pP3Label_preReloc = Uintah::VarLabel::create("p.p3+",
        Uintah::ParticleVariable<double>::getTypeDescription());
    }

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
      return params;
    }

    // Computes and requires for internal evolution variables
    virtual void addInitialComputesAndRequires(Uintah::Task* task,
                                               const Uintah::MPMMaterial* matl,
                                               const Uintah::PatchSet* patches);

    virtual void initializeInternalVariable(Uintah::ParticleSubset* pset,
                                            Uintah::DataWarehouse* new_dw) {}

    virtual void initializeInternalVariable(const Uintah::Patch* patch,
                                            const Uintah::MPMMaterial* matl,
                                            Uintah::ParticleSubset* pset,
                                            Uintah::DataWarehouse* new_dw,
                                            Uintah::MPMLabel* lb,
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
                    ParameterDict& params,
                    const double& pP3);
    double elasticVolStrainYield(const double& ev_p_bar,
                                 ParameterDict& params);
    double crushCurveDrainedSandX(const double& ev_p_bar) ;
 };

} // End namespace Uintah

#endif  // __MASON_SAND_INT_VAR_MODEL_H__ 
