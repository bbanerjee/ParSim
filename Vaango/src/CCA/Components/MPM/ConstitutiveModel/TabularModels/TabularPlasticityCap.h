/*
 * The MIT License
 *
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#ifndef __MPM_CONSTITUTIVEMODEL_TABULAR_PLASTICITY_CAP_H__
#define __MPM_CONSTITUTIVEMODEL_TABULAR_PLASTICITY_CAP_H__

#include <CCA/Components/MPM/ConstitutiveModel/TabularModels/TabularPlasticity.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuliModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_TabularCap.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCondition.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/SearchUtils.h>
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

class TabularPlasticityCap : public TabularPlasticity
{

public:

  TabularPlasticityCap(Uintah::ProblemSpecP& ps, Uintah::MPMFlags* flag);
  TabularPlasticityCap(const TabularPlasticityCap* cm);
  TabularPlasticityCap(const TabularPlasticityCap& cm);
  ~TabularPlasticityCap() override;
  TabularPlasticityCap& operator=(const TabularPlasticityCap& cm) = delete;


  ModelType modelType() const override
  {
    return ModelType::INCREMENTAL;
  }

  void outputProblemSpec(Uintah::ProblemSpecP& ps,
                         bool output_cm_tag = true) override;

  // clone
  std::unique_ptr<ConstitutiveModel> clone() override;

  /*! Get parameters */
  ParameterDict getParameters() const
  {
    ParameterDict params;
    return params;
  }

  // compute stable timestep for this patch
  void computeStableTimestep(const Uintah::Patch* patch,
                             const Uintah::MPMMaterial* matl,
                             Uintah::DataWarehouse* new_dw) override;

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

  /*! This is for adding/deleting particles when a particle is switched
      from one material to another */
  void allocateCMDataAdd(Uintah::DataWarehouse* new_dw,
                         Uintah::ParticleSubset* addset,
                         Uintah::ParticleLabelVariableMap* newState,
                         Uintah::ParticleSubset* delset,
                         Uintah::DataWarehouse* old_dw) override;

private:

  double d_consistency_bisection_tolerance;
  double d_max_bisection_iterations;
  bool d_decrease_substep;

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
                                    const ModelState_TabularCap& state_old,
                                    ModelState_TabularCap& state_new);

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
  void computeElasticProperties(ModelState_TabularCap& state) const;
  void computeElasticProperties(ModelState_TabularCap& state,
                                double elasticVolStrainInc,
                                double plasticVolStrainInc) const;
  std::tuple<double, double>
    computeElasticProperties(const ModelState_TabularCap& state,
                             const Matrix3& elasticStrain,
                             const Matrix3& plasticStrain) const;
  std::tuple<double, double>
    computeElasticProperties(const ModelState_TabularCap& state) const;

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
  Uintah::Matrix3 computeTrialStress(const ModelState_TabularCap& state_old,
                                     const Uintah::Matrix3& strain_inc);

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
                           const ModelState_TabularCap& state_substep,
                           const ModelState_TabularCap& state_trial);

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
                      const ModelState_TabularCap& state_old,
                      ModelState_TabularCap& state_new);

  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: Compute the stress state at the end of a purely elastic
   *         substep.
   * Purpose:
   *   Find the updated stress before an elastic-plastic non-hardening
   *   stress update step
   *
   * Returns:
   *   Status  :  whether the procedure is sucessful or has failed
   *   double  :  elastic delta T
   *    
   */
  //////////////////////////////////////////////////////////////////////////
  std::tuple<Status, double>
  computePurelyElasticSubstep(double deltaT, const Matrix3& D,
                              const ModelState_TabularCap& state_k_old,
                              const ModelState_TabularCap& state_k_trial,
                              ModelState_TabularCap& state_k_new);

  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: nonHardeningReturnElasticPlastic
   * Purpose:
   *   Computes a non-hardening return to the yield surface in the meridional
   *   profile (constant Lode angle) based on the current values of the 
   *   internal state variables and elastic properties.  Returns the updated 
   *   stress.
   *   NOTE: all values of r and z in this function are transformed!
   * Inputs:
   *   elasticplastic_dt = timestep size after purely elastic update
   *   D                 = rate of deformation 
   *   state_k_elastic   = state at the end of the purely elastic substep
   * Outputs:
   *   state_k_nonhardening = state at the end of the non-hardening 
   *                          elastic-plastic substep
   * Returns:
   *   Status: true  = success
   *           false = failure
   *   Matrix3: plastic strain increment
   */
  //////////////////////////////////////////////////////////////////////////
  std::tuple<Status, Matrix3>
  nonHardeningReturnElasticPlastic(double elasticplastic_dt,
                          const Matrix3& D,
                          const ModelState_TabularCap& state_k_elastic,
                          ModelState_TabularCap& state_k_nonhardening);

  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: consistencyBisectionHardeningSoftening
   * Purpose:
   *   Find the updated stress for hardening/softening plasticity using the consistency
   *   bisection algorithm
   *
   *   Returns - whether the procedure is sucessful or has failed
   */
  //////////////////////////////////////////////////////////////////////////
  Status consistencyBisectionHardeningSoftening(double elastic_plastic_dt, 
                                      const Matrix3& D,              
                                      const ModelState_TabularCap& state_k_elastic,
                                      const ModelState_TabularCap& state_k_nonhardening,
                                      const Matrix3& deltaEps_p_nonhardening,
                                      ModelState_TabularCap& state_new);

  
  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: computeSigmaFixed
   * Purpose:
   *   Compute sigma_F (stress state on the yield surface with updated normal)
   * Inputs:
   *   state_old      = state at start of substep
   *   state_trial    = trial state at start of substep
   * Outputs:
   *   Tuple containing the sequence
   *     sig_F   = projected stress on fixed field surface
   *     P_F     = projection direction
   *     M_F     = unit normal to th eyield surface
   *     H_F     = hardening modulus
   *     Gamma_F = plastic parameter increment
   */
  //////////////////////////////////////////////////////////////////////////
  std::tuple<Matrix3, Matrix3, Matrix3, double, double>
  computeSigmaFixed(const ModelState_TabularCap& state_old, 
                    const ModelState_TabularCap& state_trial) const;

  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: computeSigmaHardening
   * Purpose:
   *   Make correction to sigma_F by accounting for hardening
   * Inputs:
   *   state_trial    = trial state at start of substep
   *   Gamma_F = plastic parameter increment
   *   P_F     = projection direction
   *   N_F_norm = unit normal to th eyield surface
   * Outputs:
   *   sig_H = projected stress on hardening yield surface
   *   Gamma_H = corrected plastic parameter increment
   */
  //////////////////////////////////////////////////////////////////////////
  std::tuple<Uintah::Matrix3 , double>
  computeSigmaHardening(const ModelState_TabularCap& state_trial, 
                        double Gamma_F, 
                        const Matrix3& P_F, 
                        const Matrix3& N_F_norm, 
                        double H_F) const;

  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: closestPointInZRSpace
   * Purpose:
   *   Find closest point from the trial stress to the fixed yield surface in z-r space
   * Inputs:
   *   state_old             = state at start of substep
   *   state_trial           = trial state at start of substep
   *   z_r_index (optional)  = precompute KDtree index for yield surface polyline
   * Outputs:
   *   sig_closest           = stress at closest point
   *   closest               = closest point in pbar-sqrtJ2 coordinates
   *   tangent               = tangent vector in pbar-sqrtJ2 coordinates
   */
  //////////////////////////////////////////////////////////////////////////
  std::tuple<Uintah::Matrix3, Uintah::Point, Uintah::Vector> 
  closestPointInZRSpace(const ModelState_TabularCap& state_old,
                        const ModelState_TabularCap& state_trial) const;

  std::tuple<Uintah::Matrix3, Uintah::Point, Uintah::Vector> 
  closestPointInZRSpace(const ModelState_TabularCap& state_old,
                        const ModelState_TabularCap& state_trial,
                        const Polyline& z_r_table, 
                        const Util::PolylineKDTree& z_r_index) const;

  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: computeYieldSurfaceNormal
   * Purpose:
   *   Find yield surface normal at the closest point
   * Inputs:
   *   state_old      = state at start of substep
   *   sig_closest    = stress at closest point
   * Outputs:
   *   df_dsigma      = normal to yield surface
   */
  //////////////////////////////////////////////////////////////////////////
  Uintah::Matrix3 
  computeYieldSurfaceNormal(const ModelState_TabularCap& state_old,
                            const Matrix3& sig_closest) const;

  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: computeProjectionTensor
   * Purpose:
   *   Compute projection tensor (P = C:M + Z) at the closest point
   * Inputs:
   *   state_old       = state at start of substep
   *   sig_closest     = stress at closest point
   *   df_dsig_closest = unit normal to yield surface at closest point
   * Outputs:
   *   P               = projection tensor
   */
  //////////////////////////////////////////////////////////////////////////
  Uintah::Matrix3 
  computeProjectionTensor(const ModelState_TabularCap& state_old, 
                          const Matrix3& sig_closest,
                          const Matrix3& df_dsig_closest) const;

  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: computeGammaClosest
   * Purpose:
   *   Compute plastic strain increment (Gamma = lambda_{n+1} - \lambda_n)
   * Inputs:
   *   sig_closest     = stress at closest point
   *   sig_trial       = trial stress
   *   df_dsig_closest = unit normal to yield surface at closest point
   *   P_closest       = projection tensor
   * Outputs:
   *   Gamma_F         = value of gamma at the closest point
   */
  //////////////////////////////////////////////////////////////////////////
  double 
  computeGammaClosest(const Matrix3& sig_closest, 
                      const Matrix3& sig_trial, 
                      const Matrix3& df_dsig_closest, 
                      const Matrix3& P_closest) const;

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
   *   state_old      = state at start of substep
   *   state_trial    = trial state at start of substep
   * Outputs:
   *   sig_new                 = updated stress at end of substep
   *   elasticStrain_inc_new   = updated elastic strain increment at end of
   *                             substep
   *   plasticStrain_inc_new   = updated plastic strain increment at end of
   *                             substep
   * Returns:
   *   true  = success
   *   false = failure
   */
  //////////////////////////////////////////////////////////////////////////
  Status nonHardeningReturn(const Uintah::Matrix3& strain_inc,
                          const ModelState_TabularCap& state_old,
                          const ModelState_TabularCap& state_trial,
                          Uintah::Matrix3& sig_new,
                          Uintah::Matrix3& elasticStrain_inc_new,
                          Uintah::Matrix3& plasticStrain_inc_new);

  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: firstOrderHardeningUpdate
   * Purpose:
   *   Find the updated stress for hardening plasticity using a first-order
   *   update based on the velocity of the yield surface
   *
   *   Returns whether the procedure is sucessful or has failed
   */
  //////////////////////////////////////////////////////////////////////////
  Status firstOrderHardeningUpdate(const Matrix3& deltaEps_new,
                                   const ModelState_TabularCap& state_k_old,
                                   const ModelState_TabularCap& state_k_trial,
                                   const Matrix3& deltaEps_e_fixed,
                                   const Matrix3& deltaEps_p_fixed,
                                   const Matrix3& sig_fixed,
                                   ModelState_TabularCap& state_k_new);

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
                                      const ModelState_TabularCap& state_old,
                                      const ModelState_TabularCap& state_trial,
                                      const Matrix3& deltaEps_e_0,
                                      const Matrix3& deltaEps_p_0,
                                      const Matrix3& sig_0,
                                      ModelState_TabularCap& state_new);

  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: computeInternalVariables
   * Purpose:
   *   Update an old state with new values of internal variables given the old
   *   state and an increment in volumetric plastic strain
   * Inputs:
   *   state         - Old state
   *   delta_eps_e_v - negative in comression volumetric elastic strain
   *   delta_eps_p_v - negative in comression volumetric plastic strain
   * Outputs:
   *   state         - Modified state
   * Returns:  true if success
   *           false if failure
   */
  //////////////////////////////////////////////////////////////////////////
  Status computeInternalVariables(ModelState_TabularCap& state,
                                const double& delta_eps_e_v,
                                const double& delta_eps_p_v);

  double computeInternalVariable(const ModelState_TabularCap& state,
                                 const Matrix3& elasticStrain,
                                 const Matrix3& plasticStrain) const;
  double computeInternalVariable(const ModelState_TabularCap& state) const;

  //////////////////////////////////////////////////////////////////////////
  /**
   * Method: computeElasticDeltaT
   * Purpose:
   *   Compute the fraction of the time step that is elastic
   * Inputs:
   *   double       - Total delta t
   *   ModelState   - Old state
   *   ModelState   - Trial state
   * Returns:
   *   bool         - true if success
   *                  false if failure
   *   double       - t value of intersection of line segment joining
   *                  old state and trial state with yield surface
   *                  in pbar-sqrtJ2 space
   *   double       - elastic delta t
   *   Point        - point of intersection in pbar-sqrtJ2 space
   */
  //////////////////////////////////////////////////////////////////////////
  std::tuple<bool, double, double, Uintah::Point>
    computeElasticDeltaT(double dt,
                         const ModelState_TabularCap& state_old,
                         const ModelState_TabularCap& state_trial);
};
} // End namespace Uintah

#endif // __MPM_CONSTITUTIVEMODEL_TABULAR_PLASTICITY_CAP_H__
