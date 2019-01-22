/*
 * The MIT License
 *
 * Copyright (c) 2015-2017 Parresia Research Limited, New Zealand
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

#ifndef __MPM_CONSTITUTIVEMODEL_TABULAR_PLASTICITY_H__
#define __MPM_CONSTITUTIVEMODEL_TABULAR_PLASTICITY_H__

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ElasticModuliModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_Tabular.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/TabularData.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/YieldCondition.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <cmath>

namespace Uintah {
class MPMLabel;
class MPMFlags;
class Matrix3;
class VarLabel;
class Patch;
class Task;
class MPMMaterial;
}

namespace Vaango {

class TabularPlasticity : public Uintah::ConstitutiveModel
{

public:

  enum class Status 
  {
    SUCCESS,
    UNREASONABLE_SUBSTEPS,
    TOO_SMALL_TIMESTEP,
    TOO_LARGE_PLASTIC_STRAIN,
    TOO_LARGE_YIELD_NORMAL_CHANGE,
    TOO_MANY_CONSISTENCY_ITERATIONS,
    UNREASONABLE_INTERNAL_VARIABLE_VALUE
  };

  // Create datatype for storing model parameters
  struct CMData
  {
    double yield_scale_fac;
    double subcycling_characteristic_number;
  };

  const Uintah::VarLabel* pElasticStrainLabel;
  const Uintah::VarLabel* pElasticStrainLabel_preReloc;
  const Uintah::VarLabel* pPlasticStrainLabel;
  const Uintah::VarLabel* pPlasticStrainLabel_preReloc;
  const Uintah::VarLabel* pPlasticCumEqStrainLabel; // Equivalent plastic strain
  const Uintah::VarLabel* pPlasticCumEqStrainLabel_preReloc; 
  const Uintah::VarLabel* pElasticVolStrainLabel; // Elastic Volumetric Strain
  const Uintah::VarLabel* pElasticVolStrainLabel_preReloc;
  const Uintah::VarLabel* pPlasticVolStrainLabel; // Plastic Volumetric Strain
  const Uintah::VarLabel* pPlasticVolStrainLabel_preReloc;
  const Uintah::VarLabel* pBulkModulusLabel; 
  const Uintah::VarLabel* pBulkModulusLabel_preReloc;

  TabularPlasticity(Uintah::ProblemSpecP& ps, Uintah::MPMFlags* flag);
  TabularPlasticity(const TabularPlasticity* cm);
  TabularPlasticity(const TabularPlasticity& cm);
  ~TabularPlasticity() override;
  TabularPlasticity& operator=(const TabularPlasticity& cm) = delete;


  ModelType modelType() const override
  {
    return ModelType::INCREMENTAL;
  }

  void outputProblemSpec(Uintah::ProblemSpecP& ps,
                         bool output_cm_tag = true) override;

  // clone
  TabularPlasticity* clone() override;

  /*! Get parameters */
  ParameterDict getParameters() const
  {
    ParameterDict params;
    return params;
  }

  // compute stable timestep for this patch
  virtual void computeStableTimestep(const Uintah::Patch* patch,
                                     const Uintah::MPMMaterial* matl,
                                     Uintah::DataWarehouse* new_dw);

  // compute stress at each particle in the patch
  void computeStressTensor(const Uintah::PatchSubset* patches,
                           const Uintah::MPMMaterial* matl,
                           Uintah::DataWarehouse* old_dw,
                           Uintah::DataWarehouse* new_dw) override;

  ////////////////////////////////////////////////////////////////////////
  /* Make the value for pLocalized computed locally available outside of the
   * model. */
  ////////////////////////////////////////////////////////////////////////
  void addRequiresDamageParameter(
    Uintah::Task* task, const Uintah::MPMMaterial* matl,
    const Uintah::PatchSet* patches) const override;

  ////////////////////////////////////////////////////////////////////////
  /* Make the value for pLocalized computed locally available outside of the
   * model */
  ////////////////////////////////////////////////////////////////////////
  void getDamageParameter(const Uintah::Patch* patch,
                          Uintah::ParticleVariable<int>& damage, int dwi,
                          Uintah::DataWarehouse* old_dw,
                          Uintah::DataWarehouse* new_dw) override;

  // carry forward CM data for RigidMPM
  void carryForward(const Uintah::PatchSubset* patches,
                    const Uintah::MPMMaterial* matl,
                    Uintah::DataWarehouse* old_dw,
                    Uintah::DataWarehouse* new_dw) override;

  // initialize  each particle's constitutive model data
  void initializeCMData(const Uintah::Patch* patch,
                        const Uintah::MPMMaterial* matl,
                        Uintah::DataWarehouse* new_dw) override;

  void addInitialComputesAndRequires(
    Uintah::Task* task, const Uintah::MPMMaterial* matl,
    const Uintah::PatchSet* patches) const override;

  void addComputesAndRequires(Uintah::Task* task,
                              const Uintah::MPMMaterial* matl,
                              const Uintah::PatchSet* patches) const override;

  void addComputesAndRequires(Uintah::Task* task,
                              const Uintah::MPMMaterial* matl,
                              const Uintah::PatchSet* patches,
                              const bool recursion,
                              const bool dummy) const override;

  void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                        std::vector<const Uintah::VarLabel*>& to) override;

  double computeRhoMicroCM(double pressure, const double p_ref,
                           const Uintah::MPMMaterial* matl, double temperature,
                           double rho_guess) override;

  void computePressEOSCM(double rho_m, double& press_eos, double p_ref,
                         double& dp_drho, double& ss_new,
                         const Uintah::MPMMaterial* matl,
                         double temperature) override;

  double getCompressibility() override;

  /*! This is for adding/deleting particles when a particle is switched
      from one material to another */
  void allocateCMDataAdd(Uintah::DataWarehouse* new_dw,
                         Uintah::ParticleSubset* addset,
                         Uintah::ParticleLabelVariableMap* newState,
                         Uintah::ParticleSubset* delset,
                         Uintah::DataWarehouse* old_dw) override;

protected:

  ElasticModuliModel* d_elastic;
  YieldCondition* d_yield;
  TabularData d_hydrostat;

  CMData d_cm;

  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: computeElasticProperties
   * Purpose:
   *   Compute the bulk and shear mdoulus at a given state
   * Inputs:
   *   state = state
   * Modifies
   *   state.bulkModulus
   *   state.shearModulus
   */
  //////////////////////////////////////////////////////////////////////////
  void computeElasticProperties(ModelState_Tabular& state) const;

  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: computeTrialStress
   * Purpose:
   *   Compute the trial stress for some increment in strain assuming linear
   *   elasticity over the step.
   * Inputs:
   *   state_old = state at t = t_n
   *   strain_inc = strain increment (D * delT)
   *     D = rate of deformation
   *     delT = time step
   * Returns
   *   stress_trial
   */
  //////////////////////////////////////////////////////////////////////////
  Uintah::Matrix3 computeTrialStress(const ModelState_Tabular& state_old,
                                     const Uintah::Matrix3& strain_inc);

  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: computeZMatrix
   * Purpose:
   *   Compute the elastic-plastic coupling Z matrix
   * Inputs:
   *   state 
   * Returns
   *   Z matrix
   */
  //////////////////////////////////////////////////////////////////////////
  Matrix3 computeZMatrix(const ElasticModuli& moduli,
                         const ElasticModuli& derivs,
                         double p,
                         const Matrix3& s,
                         const Matrix3& M) const;

private:

  void initializeLocalMPMLabels();
  void checkInputParameters();

  void initializeInternalVariables(const Patch* patch, const MPMMaterial* matl,
                                   ParticleSubset* pset, DataWarehouse* new_dw,
                                   ParameterDict& params);

  //////////////////////////////////////////////////////////////////////////
  /**
   * Function:
   *   rateIndependentPlasticUpdate
   * Purpose:
   *   Rate-independent plasticity:
   *   Divides the strain increment into substeps, and calls substep function
   *   to copute the stress and internal variables
   * Inputs:
   *   D          = rate of deformation
   *   delT       = time increment
   *   idx        = Particle index
   *   particleID = long64 ID
   *   state_old  = state at t = t_n
   * Outputs:
   *   state_new  = state at t = t_(n+1)
   * Returns: True for success; False for failure
   */
  //////////////////////////////////////////////////////////////////////////
  Status rateIndependentPlasticUpdate(const Uintah::Matrix3& D,
                                    const double& delT,
                                    Uintah::particleIndex idx,
                                    Uintah::long64 pParticleID,
                                    const ModelState_Tabular& state_old,
                                    ModelState_Tabular& state_new);

  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: computeStepDivisions
   * Purpose:
   *   Compute the number of step divisions (substeps) based on a comparison
   *   of the trial stress relative to the size of the yield surface, as well
   *   as change in elastic properties between sigma_n and sigma_trial.
   * Inputs:
   *   idx            = particle index
   *   particleID     = long64 ID
   *   state_substep  = state at the current subsetp
   *   state_trial    = the trial state
   * Returns
   *   nSub = the number of substeps
   */
  //////////////////////////////////////////////////////////////////////////
  int computeStepDivisions(Uintah::particleIndex idx, Uintah::long64 particleID,
                           const ModelState_Tabular& state_substep,
                           const ModelState_Tabular& state_trial);

  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: computeSubstep
   * Purpose:
   *   Computes the updated stress state for a substep that may be either
   *   elastic, plastic, or partially elastic.
   * Inputs:
   *   D              = rate of deformation tensor
   *   dt             = substep time increment
   *   state_old      = state at start of substep
   * Outputs:
   *   state_new      = updated state at end of substep
   * Returns:
   *   True for success; false for failure
   */
  //////////////////////////////////////////////////////////////////////////
  Status computeSubstep(const Uintah::Matrix3& D, const double& dt,
                      const ModelState_Tabular& state_old,
                      ModelState_Tabular& state_new);

  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: nonHardeningReturn
   * Purpose:
   *   Computes a non-hardening return to the yield surface in the meridional
   *   profile (constant Lode angle) based on the current values of the 
   *   internal state variables and elastic properties.  Returns the updated 
   *   stress and  the increment in plastic strain corresponding to this 
   *   return.
   *   NOTE: all values of r and z in this function are transformed!
   * Inputs:
   *   strain_inc     = total strain icremenet = D*dt
   *     D              = rate of deformation tensor
   *     dt             = substep time increment
   *   state_old      = state at start of substep
   *   state_trial    = trial state at start of substep
   *   params         = yield condition parameters
   * Outputs:
   *   sig_new                 = updated stress at end of substep
   *   elasticStrain_inc_new   = updated elastic strain increment at end of
   * substep
   *   plasticStrain_inc_new   = updated plastic strain increment at end of
   * substep
   * Returns:
   *   true  = success
   *   false = failure
   */
  //////////////////////////////////////////////////////////////////////////
  Status nonHardeningReturn(const Uintah::Matrix3& strain_inc,
                          const ModelState_Tabular& state_old,
                          const ModelState_Tabular& state_trial,
                          Uintah::Matrix3& sig_new,
                          Uintah::Matrix3& elasticStrain_inc_new,
                          Uintah::Matrix3& plasticStrain_inc_new);

  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: consistencyBisectionSimplified
   * Purpose:
   *   Find the updated stress for hardening plasticity using the consistency
   *   bisectionalgorithm
   *   Returns whether the procedure is sucessful orhas failed
   * Inputs:
   *   deltaEps_new = strain increment for the substep
   *   state_old    = state at the beginning of the substep
   *   state_trial  = trial state
   *   deltaEps_e_0 = elastic strain increment at the beginning of substep
   *   deltaEps_p_0 = plastic strain increment at the beginning of substep
   *   sig_0        = stress at the beginning of substep
   *   params       = yield condition parameters
   * Outputs:
   *   state_new    = state at the end of the substep
   * Returns:
   *   isSuccess    = true if success, else false
   */
  //////////////////////////////////////////////////////////////////////////
  Status consistencyBisectionSimplified(const Matrix3& deltaEps_new,
                                      const ModelState_Tabular& state_old,
                                      const ModelState_Tabular& state_trial,
                                      const Matrix3& deltaEps_e_0,
                                      const Matrix3& deltaEps_p_0,
                                      const Matrix3& sig_0,
                                      ModelState_Tabular& state_new);

  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: computeInternalVariables
   * Purpose:
   *   Update an old state with new values of internal variables given the old
   *   state and an increment in volumetric plastic strain
   * Inputs:
   *   state         - Old state
   *   delta_eps_p_v - negative in comression volumetruc plastic strain
   * Outputs:
   *   state         - Modified state
   * Returns:  true if success
   *           false if failure
   */
  //////////////////////////////////////////////////////////////////////////
  Status computeInternalVariables(ModelState_Tabular& state,
                                const double& delta_eps_p_v);

};
} // End namespace Vaango

#endif // __MPM_CONSTITUTIVEMODEL_TABULAR_PLASTICITY_H__
