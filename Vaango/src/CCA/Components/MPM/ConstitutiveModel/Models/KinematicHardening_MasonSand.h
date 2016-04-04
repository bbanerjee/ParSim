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

#ifndef __MASONSAND_POREPRESSURE_KINEMATIC_HARDENING_MODEL_H__
#define __MASONSAND_POREPRESSURE_KINEMATIC_HARDENING_MODEL_H__



#include <CCA/Components/MPM/ConstitutiveModel/Models/KinematicHardeningModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/Pressure_Air.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/Pressure_Water.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/Pressure_Granite.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelStateBase.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

  /////////////////////////////////////////////////////////////////////////////
  /*!
    \class KinematicHardening_MasonSand
    \brief Backstress model for the pore pressure (Mason sand)
    \author Biswajit Banerjee, 
   
    See ONR-MURI Reports from August 2015 through February 2016.
  */
  /////////////////////////////////////////////////////////////////////////////

  class KinematicHardening_MasonSand : public KinematicHardeningModel {

    friend class InternalVar_MasonSand;

  private:

    static const Uintah::Matrix3 Identity; 

    struct CMData {
      double fluid_pressure_initial;
    };

    CMData d_cm;

    /* Tangent bulk modulus models for air, water, granite */
    Pressure_Air     d_air;
    Pressure_Water   d_water;
    Pressure_Granite d_granite;

    // Prevent copying of this class
    //KinematicHardening_MasonSand(const KinematicHardening_MasonSand &cm);
    KinematicHardening_MasonSand& operator=(const KinematicHardening_MasonSand &cm);

  public:
    // constructors
    KinematicHardening_MasonSand(Uintah::ProblemSpecP& ps,
                                 InternalVariableModel* intvar);
    KinematicHardening_MasonSand(const KinematicHardening_MasonSand* cm);
         
    // destructor 
    virtual ~KinematicHardening_MasonSand();

    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps);
         
    /*! Get parameters */
    std::map<std::string, double> getParameters() const {
      std::map<std::string, double> params;
      params["Pf0"] = d_cm.fluid_pressure_initial;
      return params;
    }

    //////////
    /*! \brief Calculate the back stress */
    //////////
    void computeBackStress(const ModelStateBase* state,
                           const double& delT,
                           const Uintah::particleIndex idx,
                           const double& delLambda,
                           const Uintah::Matrix3& df_dsigma_new,
                           const Uintah::Matrix3& backStress_old,
                           Uintah::Matrix3& backStress_new) {}

    void computeBackStress(const ModelStateBase* state,
                           Uintah::Matrix3& backStress_new);

    void eval_h_beta(const Uintah::Matrix3& df_dsigma,
                     const ModelStateBase* state,
                     Uintah::Matrix3& h_beta) {}

  public:

    // We use the Matrix3 pBackStressLabel instead of pZetaLabel.
    // pBackStressLabel is defined in the base class.
    //const Uintah::VarLabel*   pZetaLabel;  
    //const Uintah::VarLabel*   pZetaLabel_preReloc;

    // Add particle state for these labels
    void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                          std::vector<const Uintah::VarLabel*>& to) 
    {
      from.push_back(pBackStressLabel);
      to.push_back(pBackStressLabel_preReloc);
    }

    /**
     * Initialize pore pressure label
     */
    void initializeLocalMPMLabels() 
    {
      pBackStressLabel            = Uintah::VarLabel::create("p.porePressure",
        Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());
      pBackStressLabel_preReloc   = Uintah::VarLabel::create("p.porePressure+",
        Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());
    }

    /**
     * Set up task graph for initialization
     */
    void addInitialComputesAndRequires(Uintah::Task* task,
                                       const Uintah::MPMMaterial* matl,
                                       const Uintah::PatchSet* patch) const 
    {
      const Uintah::MaterialSubset* matlset = matl->thisMaterial(); 
      task->computes(pBackStressLabel,   matlset);
    }

    /**
     *  Actually initialize the pore pressure
     */
    void initializeLocalVariables(const Uintah::Patch* patch,
                                  Uintah::ParticleSubset* pset,
                                  Uintah::DataWarehouse* new_dw,
                                  Uintah::constParticleVariable<double>& pVolume)
    {
      Uintah::ParticleVariable<Uintah::Matrix3> pBackStress;
      new_dw->allocateAndPut(pBackStress,      pBackStressLabel,      pset);

      for (auto iter = pset->begin(); iter != pset->end(); iter++) {
        Uintah::particleIndex idx = *iter;
        pBackStress[idx] = -d_cm.fluid_pressure_initial*Identity;
      }
    }

    /**
     * Set up task graph for parameter copying to new datawarehouse
     */
    void addComputesAndRequires(Uintah::Task* task,
                                const Uintah::MPMMaterial* matl,
                                const Uintah::PatchSet* patches) const 
    {
      const Uintah::MaterialSubset* matlset = matl->thisMaterial(); 
      task->requires(Uintah::Task::OldDW, pBackStressLabel,   matlset, Uintah::Ghost::None);
      task->computes(pBackStressLabel_preReloc,   matlset);
    }

    /**
     *  Update the pore pressure
     */
    void computeBackStress(Uintah::ParticleSubset* pset,
                           Uintah::DataWarehouse* old_dw,
                           Uintah::DataWarehouse* new_dw) 
    {
    }

    void 
    getBackStress(Uintah::ParticleSubset* pset,
                  Uintah::DataWarehouse* old_dw,
                  Uintah::constParticleVariable<Uintah::Matrix3>& pBackStress)
    {
      old_dw->get(pBackStress, pBackStressLabel, pset);
    }

    void 
    allocateAndPutBackStress(Uintah::ParticleSubset* pset,
                             Uintah::DataWarehouse* new_dw,
                             Uintah::ParticleVariable<Uintah::Matrix3>& pBackStress_new)
    {
      new_dw->allocateAndPut(pBackStress_new, pBackStressLabel_preReloc, pset); 
    } 

  };

} // End namespace Uintah

#endif  // __MASONSAND_POREPRESSURE_KINEMATIC_HARDENING_MODEL_H__ 
