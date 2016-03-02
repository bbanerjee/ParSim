/*
 * The MIT License
 *
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

    // Get the internal variables
    std::vector<Uintah::constParticleVariable<double> >
         getInternalVariables(Uintah::ParticleSubset* pset,
                              Uintah::DataWarehouse* old_dw,
                              const double& dummy)
    {
      Uintah::constParticleVariable<double> pKappa, pCapX, pPlasticVolStrain, pP3; 
      old_dw->get(pKappa,            pKappaLabel,            pset);
      old_dw->get(pCapX,             pCapXLabel,             pset);
      old_dw->get(pPlasticVolStrain, pPlasticVolStrainLabel, pset);
      old_dw->get(pP3,               pP3Label,               pset);

      std::vector<Uintah::constParticleVariable<double> > pIntVars;
      pIntVars.emplace_back(pKappa);
      pIntVars.emplace_back(pCapX);
      pIntVars.emplace_back(pPlasticVolStrain);
      pIntVars.emplace_back(pP3);
    
      return pIntVars;
    }

    std::vector<Uintah::constParticleVariable<Uintah::Matrix3> >
         getInternalVariables(Uintah::ParticleSubset* pset,
                              Uintah::DataWarehouse* old_dw,
                              const Uintah::Matrix3& dummy)
    {
      Uintah::constParticleVariable<Uintah::Matrix3> pPlasticStrain; 
      old_dw->get(pPlasticStrain, pPlasticStrainLabel, pset);

      std::vector<Uintah::constParticleVariable<Uintah::Matrix3> > pIntVars;
      pIntVars.emplace_back(pPlasticStrain);
    
      return pIntVars;
    }

    // Allocate and put the local particle internal variables
    void allocateAndPutInternalVariable(Uintah::ParticleSubset* pset,
                                        Uintah::DataWarehouse* new_dw,
                                        vectorParticleDoubleP& pVars) {
      new_dw->allocateAndPut(*pVars[0], pKappaLabel_preReloc,            pset);
      new_dw->allocateAndPut(*pVars[1], pCapXLabel_preReloc,             pset);
      new_dw->allocateAndPut(*pVars[2], pPlasticVolStrainLabel_preReloc, pset);
      new_dw->allocateAndPut(*pVars[3], pP3Label_preReloc,               pset);
    }

    // Allocate and put the local <Matrix3> particle variables
    void allocateAndPutInternalVariable(Uintah::ParticleSubset* pset,
                                        Uintah::DataWarehouse* new_dw,
                                        vectorParticleMatrix3P& pVars) {
      new_dw->allocateAndPut(*pVars[0], pPlasticStrainLabel_preReloc, pset);
    }

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

    /**
     * Function: computeDrainedHydrostaticStrength
     *
     * Purpose:
     *   Compute the drained hydrostatic strength (X) using the crush curve
     *   and the initial porosity
     *
     * Inputs:
     *   ep_v_bar = -tr(ep) = volumetric part of plastic strain tensor
     *   phi0     = initial porosity
     *
     * Returns:
     *   Xbar     = hydrostatic compressive strength     
     */
    double computeDrainedHydrostaticStrength(const double& ep_v_bar,
                                             const double& phi0) const;

    /**
     * Function: computeP3
     *
     * Purpose:
     *   Compute the crush curve parameter P3 from the initial porosity
     *
     * Inputs:
     *   phi0 = initial porosity
     *
     * Returns:
     *   p3 = crush curve parameter
     */
    inline double computeP3(const double& phi0) const {
      double p3 = -std::log(1.0 - phi0);
      return p3;
    };

    /**
     * Function: computePorosity
     *
     * Purpose:
     *   Compute the porosity from crush curve parameter P3
     *
     * Inputs:
     *   ep_v = volumetric plastic strain
     *   p3 = crush
     *
     * Returns:
     *   porosity = porosity from crush curve parameter
     */
    inline double computePorosity(const double& ep_v, const double& p3) const {
      double porosity = 1.0 - std::exp(-p3 + ep_v);
      return porosity;
    }

    /**
     * Function: computeElasticVolStrainAtYield
     *
     * Purpose:
     *   Compute the elastic volumetric strain at yield
     *
     * Inputs:
     *   ep_v_bar = -tr(ep) = volumetric part of plastic strain tensor
     *   phi0     = initial porosity
     *
     * Returns:
     *   ev_e_yield = elastic volumetric strain at yield
     */
    double computeElasticVolStrainAtYield(const double& ep_v_bar,
                                          const double& phi0) const;

    /**
     * Function: computePartSatHydrostaticStrength
     *
     * Purpose:
     *   Compute the partially saturated hydrostatic strength (X_sat)
     *
     * Inputs:
     *   I1_bar   = -tr(sigma) = isotropic part of stress tensor
     *   ep_v_bar = -tr(ep) = volumetric part of plastic strain tensor
     *   phi      = porosity
     *   Sw       = saturation
     *   phi0     = initial porosity
     *
     * Returns:
     *   Xbar_sat = hydrostatic compressive strength     
     */
    double computePartSatHydrostaticStrength(const double& I1_bar, 
                                             const double& ep_v_bar,
                                             const double& phi,
                                             const double& Sw,
                                             const double& phi0) const;
 };

} // End namespace Uintah

#endif  // __MASON_SAND_INT_VAR_MODEL_H__ 
