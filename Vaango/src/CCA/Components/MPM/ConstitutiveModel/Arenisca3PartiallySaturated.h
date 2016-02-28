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



#ifndef __ARENISCA3_PARTIALLY_SATURATED__
#define __ARENISCA3_PARTIALLY_SATURATED__


#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_MasonSand.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ElasticModuliModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/YieldCondition.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/InternalVariableModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/KinematicHardeningModel.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <CCA/Ports/DataWarehouseP.h>

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

  class Arenisca3PartiallySaturated : public Uintah::ConstitutiveModel {

  public:
    static const double one_third;
    static const double two_third;
    static const double four_third;
    static const double sqrt_two;
    static const double one_sqrt_two;
    static const double sqrt_three;
    static const double one_sqrt_three;
    static const double one_sixth;
    static const double one_ninth;
    static const double pi;
    static const double pi_fourth;
    static const double pi_half;
    static const Uintah::Matrix3 Identity;

    static const int NMAX;
    static const std::vector<double> sinV;
    static const std::vector<double> cosV;

    // Create datatype for storing model parameters
    struct CMData {
      double subcycling_characteristic_number;
      bool   use_disaggregation_algorithm;
      double K0_Murnaghan_EOS;
      double n_Murnaghan_EOS;
    };

    // Initial porosity and saturation parameters
    struct FluidEffectParameters {
      double phi0;  // initial porosity
      double Sw0;   // initial water saturation
    };

    const Uintah::VarLabel* pPorosityLabel;                      // Porosity
    const Uintah::VarLabel* pPorosityLabel_preReloc; 

    const Uintah::VarLabel* pSaturationLabel;                    // Saturation
    const Uintah::VarLabel* pSaturationLabel_preReloc; 

    const Uintah::VarLabel* pLocalizedLabel;
    const Uintah::VarLabel* pLocalizedLabel_preReloc;
    const Uintah::VarLabel* pElasticVolStrainLabel;              //Elastic Volumetric Strain
    const Uintah::VarLabel* pElasticVolStrainLabel_preReloc;
    const Uintah::VarLabel* pStressQSLabel;
    const Uintah::VarLabel* pStressQSLabel_preReloc;

  protected:

    ElasticModuliModel*      d_elastic;
    YieldCondition*          d_yield;
    InternalVariableModel*   d_intvar;
    KinematicHardeningModel* d_backstress;

  private:
    double small_number;
    double big_number;

    double d_Kf,
           d_Km,
           d_phi_i,
           d_ev0,
           d_C1;

    CMData d_cm;
    FluidEffectParameters d_fluidParam;

    // Prevent copying of this class
    // copy constructor

    Arenisca3PartiallySaturated& operator=(const Arenisca3PartiallySaturated &cm);

    void initializeLocalMPMLabels();
    void checkInputParameters();

  public:
    // constructor
    Arenisca3PartiallySaturated(Uintah::ProblemSpecP& ps, Uintah::MPMFlags* flag);
    Arenisca3PartiallySaturated(const Arenisca3PartiallySaturated* cm);

    // destructor
    virtual ~Arenisca3PartiallySaturated();

    virtual void outputProblemSpec(SCIRun::ProblemSpecP& ps,bool output_cm_tag = true);

    // clone
    Arenisca3PartiallySaturated* clone();

    /*! Get parameters */
    ParameterDict getParameters() const {
      ParameterDict params;
      params["phi0"] = d_fluidParam.phi0;
      params["Sw0"] = d_fluidParam.Sw0;
      return params;
    }

    // compute stable timestep for this patch
    virtual void computeStableTimestep(const Uintah::Patch* patch,
                                       const Uintah::MPMMaterial* matl,
                                       Uintah::DataWarehouse* new_dw);

    // compute stress at each particle in the patch
    virtual void computeStressTensor(const Uintah::PatchSubset* patches,
                                     const Uintah::MPMMaterial* matl,
                                     Uintah::DataWarehouse* old_dw,
                                     Uintah::DataWarehouse* new_dw);

  private: //non-uintah mpm constitutive specific functions

    //////////////////////////////////////////////////////////////////////////
    /**
     * Function: 
     *   rateIndependentPlasticUpdate
     *
     * Purpose:
     *   Rate-independent plasticity:
     *   Divides the strain increment into substeps, and calls substep function
     *   to copute the stress and internal variables
     *
     * Inputs:
     *   D          = rate of deformation
     *   delT       = time increment
     *   yieldParam = parameters of the yield condition model (mean values)
     *   idx        = Particle index
     *   particleID = long64 ID
     *   state_old  = state at t = t_n
     *
     * Outputs:
     *
     *   state_new  = state at t = t_(n+1)
     *
     * Returns: True for success; False for failure
     */
    //////////////////////////////////////////////////////////////////////////
    bool rateIndependentPlasticUpdate(const Uintah::Matrix3& D,
                                      const double& delT, 
                                      const ParameterDict& yieldParam,
                                      Uintah::particleIndex idx, 
                                      Uintah::long64 pParticleID, 
                                      const ModelState_MasonSand& state_old,
                                      ModelState_MasonSand& state_new);

    //////////////////////////////////////////////////////////////////////////
    /** 
     * Function: rateDependentPlasticUpdate
     *
     * Purpose:
     *   Rate-dependent plastic step
     *   Compute the new dynamic stress from the old dynamic stress and the new and old QS stress
     *   using Duvaut-Lions rate dependence, as described in "Elements of Phenomenological Plasticity",
     *   by RM Brannon.
     *
     * Inputs:
     *   D                = rate of deformation
     *   delT             = time increment
     *   yieldParams      = parameters of the yield condition model (mean values)
     *   stateStatic_old  = Quasistatic stress state at t = t_n
     *   stateStatic_new  = Quasistatic stress state at t = t_n+1
     *   stateDynamic_old = Dynamic stress state at t = t_n
     *
     * Outputs:
     *
     *   pStress_new  = stress state at t = t_(n+1)
     *
     * Returns: True for rate-dependent plasticity; False for rate-independent plasticity
     */
    bool rateDependentPlasticUpdate(const Uintah::Matrix3& D,
                                    const double& delT,
                                    const ParameterDict& yieldParams,
                                    const ModelState_MasonSand& stateStatic_old,
                                    const ModelState_MasonSand& stateStatic_new,
                                    const ModelState_MasonSand& stateDynamic_old,
                                    Uintah::Matrix3& pStress_new);

    //////////////////////////////////////////////////////////////////////////
    /** 
     * Method: computeElasticProperties
     *
     * Purpose: 
     *   Compute the bulk and shear mdoulus at a given state
     *
     * Inputs:
     *   state = state
     * 
     * Modifies
     *   state.bulkModulus
     *   state.shearModulus
     */
     //////////////////////////////////////////////////////////////////////////
    void computeElasticProperties(ModelState_MasonSand& state);

    //////////////////////////////////////////////////////////////////////////
    /** 
     * Method: computeTrialStress
     *
     * Purpose: 
     *   Compute the trial stress for some increment in strain assuming linear elasticity
     *   over the step.
     *
     * Inputs:
     *   state_old = state at t = t_n 
     *   strain_inc = strain increment (D * delT)
     *     D = rate of deformation
     *     delT = time step 
     * 
     * Returns
     *   stress_trial
     */
     //////////////////////////////////////////////////////////////////////////
    Uintah::Matrix3 computeTrialStress(const ModelState_MasonSand& state_old,
                                       const Uintah::Matrix3& strain_inc);

    //////////////////////////////////////////////////////////////////////////
    /** 
     * Method: computeStepDivisions
     *
     * Purpose: 
     *   Compute the number of step divisions (substeps) based on a comparison
     *   of the trial stress relative to the size of the yield surface, as well
     *   as change in elastic properties between sigma_n and sigma_trial.
     *
     * Inputs:
     *   idx            = particle index
     *   particleID     = long64 ID
     *   state_substep  = state at the current subsetp 
     *   state_trial    = the trial state
     *   yieldParam     = the mean parameters of the yield surface model
     * 
     * Returns
     *   nSub = the number of substeps
     */
     //////////////////////////////////////////////////////////////////////////
    int computeStepDivisions(Uintah::particleIndex idx,
                             Uintah::long64 particleID,
                             const ModelState_MasonSand& state_substep,
                             const ModelState_MasonSand& state_trial,
                             const ParameterDict& yieldParam);

    //////////////////////////////////////////////////////////////////////////
    /** 
     * Method: computeSubstep
     *
     * Purpose: 
     *   Computes the updated stress state for a substep that may be either 
     *   elastic, plastic, or partially elastic.   
     *
     * Inputs:
     *   D              = rate of deformation tensor
     *   dt             = substep time increment
     *   yieldParams    = yield parameter dictionary
     *   state_old      = state at start of substep
     * 
     * Outputs:
     *   state_new      = updated state at end of substep
     *
     * Returns:
     *   True for success; false for failure
     */
    //////////////////////////////////////////////////////////////////////////
    bool computeSubstep(const Uintah::Matrix3& D,
                        const double& dt,
                        const ParameterDict& yieldParams,
                        const ModelState_MasonSand& state_old,
                        ModelState_MasonSand& state_new);

    //////////////////////////////////////////////////////////////////////////
    /** 
     * Method: nonHardeningReturn
     *
     * Purpose: 
     *   Computes a non-hardening return to the yield surface in the meridional profile
     *   (constant Lode angle) based on the current values of the internal state variables
     *   and elastic properties.  Returns the updated stress and  the increment in plastic
     *   strain corresponding to this return.
     *
     *   NOTE: all values of r and z in this function are transformed!
     *
     * Inputs:
     *   strain_inc     = total strain icremenet = D*dt
     *     D              = rate of deformation tensor
     *     dt             = substep time increment
     *   state_old      = state at start of substep
     *   state_trial    = trial state at start of substep
     *   params         = yield condition parameters
     * 
     * Outputs:
     *   sig_new                 = updated stress at end of substep
     *   plasticStrain_inc_new   = updated plastic strain incremente at end of substep
     *
     * Returns:
     *   0 for success; not 0 for failure
     */
    //////////////////////////////////////////////////////////////////////////
    int nonHardeningReturn(const Uintah::Matrix3& strain_inc,
                           const ModelState_MasonSand& state_old,
                           const ModelState_MasonSand& state_trial,
                           const ParameterDict& params,
                           Uintah::Matrix3& sig_new,
                           Uintah::Matrix3& plasticStrain_inc_new);

    //////////////////////////////////////////////////////////////////////////
    /**
     * Method: applyBisectionAlgorithm
     *
     * Purpose: 
     *   Uses bisection to find the intersction of a loading path with the yield surface
     *   Returns location of intersection point in transformed stress space
     *
     * Inputs:
     *   z_0, rprime_0         : Transformed Lode coordinates of stress at the start of substep
     *   z_trial, rprime_trial : Transformed Lode coordinates of trial stress 
     *   state_old             : State at the bginning of the timestep
     *   params                : Yield condition parameters
     *
     * Outputs:
     *   z_new, rprime_new     : Transformed Lode coordinates of stress at the end of substep
     */
    //////////////////////////////////////////////////////////////////////////
    void applyBisectionAlgorithm(const double& z_0, const double& rprime_0, 
                                 const double& z_trial, const double& rprime_trial, 
                                 const ModelState_MasonSand& state_old, 
                                 const ParameterDict& params,
                                 double &z_new, double &rprime_new);

    //////////////////////////////////////////////////////////////////////////
    /**
     * Method: findNewInternalPoint
     * Purpose: 
     *   Apply rotation algorithm to rotate the stress around the trial state to find
     *   a new internal point
     *   Returns location of the internal point and the angle
     *
     * Inputs:
     *   z_trial, rprime_trial : Transformed Lode coordinates of trial stress 
     *   z_new, rprime_new     : Transformed Lode coordinates of stress at the end of substep
     *   theta_old             : Angle at the beginning of the substep
     *   state_old             : State at the beginning of the timestep
     *   params                : Yield condition parameters
     *
     * Outputs:
     *   z_rot, rprime_rot     : Transformed Lode coordinates of internal point
     *
     * Returns:
     *   theta                 : Rotation angle
     */
    //////////////////////////////////////////////////////////////////////////
    double findNewInternalPoint(const double& z_trial, const double& rprime_trial, 
                                const double& z_new, const double& rprime_new, 
                                const double& theta_old,
                                const ModelState_MasonSand& state_old, 
                                const ParameterDict& params,
                                double& z_rot, double& rprime_rot);

    //////////////////////////////////////////////////////////////////////////
    /**
     * Method: evalYieldCondition
     *
     * Purpose: 
     *   Evaluate the yield condition in transformed Lode coordinates
     *   Returns whether the stress is elastic or not
     *
     * Inputs:
     *   z_stress, rprime_stress : Transformed Lode coordinates of stress point
     *   state_old               : State at the beginning of the timestep
     *   params                  : Yield condition parameters
     *
     * Returns:
     *   isElastic               : Whether state is inside yield surface or not
     */
    //////////////////////////////////////////////////////////////////////////
    bool evalYieldCondition(const double& z_stress, const double& rprime_stress, 
                            const ModelState_MasonSand& state_old, 
                            const ParameterDict& params);

    //////////////////////////////////////////////////////////////////////////
    /**
     * Method: consistencyBisection
     *
     * Purpose: 
     *   Find the updated stress for hardening plasticity using the consistency bisection 
     *   algorithm
     *   Returns whether the procedure is sucessful orhas failed
     *
     * Inputs:
     *   deltaEps_new = strain increment for the substep
     *   state_old    = state at the beginning of the substep 
     *   state_trial  = trial state
     *   deltaEps_p_0 = plastic strain increment at the beginning of substep
     *   sig_0        = stress at the beginning of substep
     *   params       = yield condition parameters
     *
     * Outputs:
     *   state_old    = state at the end of the substep 
     *
     * Returns:
     *   isSuccess    = true if success, else false
     */
    //////////////////////////////////////////////////////////////////////////
    bool consistencyBisection(const Matrix3& deltaEps_new,
                              const ModelState_MasonSand& state_old, 
                              const ModelState_MasonSand& state_trial,
                              const Matrix3& deltaEps_p_0, 
                              const Matrix3& sig_0, 
                              const ParameterDict& params, 
                              ModelState_MasonSand& state_new);

    //////////////////////////////////////////////////////////////////////////
    /**
     * Method: computeHydrostaticStrength
     *
     * Purpose: 
     *   Compute state variable X, the Hydrostatic Compressive strength (cap position)
     *   X is the value of (I1 - Zeta) at which the cap function crosses
     *   the hydrostat. 
     *   In tension, M. Homel's piecewsie formulation is used.
     *
     * Inputs:
     *   state - Model state containing updated values of
     *     evp - volumetric plastic strain
     *     P3  - Disaggregation strain P3
     *
     * Returns:
     *   double scalar value
     */
     //////////////////////////////////////////////////////////////////////////
    double computeHydrostaticStrength(const ModelState_MasonSand& state);

    double computeDerivativeOfBackstress(double Zeta,
                            double evp);

    double computePorePressure(const double ev);




  public: //Uintah MPM constitutive model specific functions
    ////////////////////////////////////////////////////////////////////////
    /* Make the value for pLocalized computed locally available outside of the model. */
    ////////////////////////////////////////////////////////////////////////
    virtual void addRequiresDamageParameter(Uintah::Task* task,
                                            const Uintah::MPMMaterial* matl,
                                            const Uintah::PatchSet* patches) const;


    ////////////////////////////////////////////////////////////////////////
    /* Make the value for pLocalized computed locally available outside of the model */
    ////////////////////////////////////////////////////////////////////////
    virtual void getDamageParameter(const Uintah::Patch* patch,
                                    Uintah::ParticleVariable<int>& damage, int dwi,
                                    Uintah::DataWarehouse* old_dw,
                                    Uintah::DataWarehouse* new_dw);


    // carry forward CM data for RigidMPM
    virtual void carryForward(const Uintah::PatchSubset* patches,
                              const Uintah::MPMMaterial* matl,
                              Uintah::DataWarehouse* old_dw,
                              Uintah::DataWarehouse* new_dw);


    // initialize  each particle's constitutive model data
    virtual void initializeCMData(const Uintah::Patch* patch,
                                  const Uintah::MPMMaterial* matl,
                                  Uintah::DataWarehouse* new_dw);


    virtual void addInitialComputesAndRequires(Uintah::Task* task,
                                               const Uintah::MPMMaterial* matl,
                                               const Uintah::PatchSet* patches) const;

    virtual void addComputesAndRequires(Uintah::Task* task,
                                        const Uintah::MPMMaterial* matl,
                                        const Uintah::PatchSet* patches) const;

    virtual void addComputesAndRequires(Uintah::Task* task,
                                        const Uintah::MPMMaterial* matl,
                                        const Uintah::PatchSet* patches,
                                        const bool recursion,
                                        const bool dummy) const;

    virtual void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                                  std::vector<const Uintah::VarLabel*>& to);

    virtual double computeRhoMicroCM(double pressure,
                                     const double p_ref,
                                     const Uintah::MPMMaterial* matl,
                                     double temperature,
                                     double rho_guess);

    virtual void computePressEOSCM(double rho_m, double& press_eos,
                                   double p_ref,
                                   double& dp_drho, double& ss_new,
                                   const Uintah::MPMMaterial* matl,
                                   double temperature);

    virtual double getCompressibility();

    /*! This is for adding/deleting particles when a particle is switched
        from one material to another */
    virtual void allocateCMDataAdd(Uintah::DataWarehouse* new_dw,
                                   Uintah::ParticleSubset* addset,
                                   Uintah::ParticleLabelVariableMap* newState,
                                   Uintah::ParticleSubset* delset,
                                   Uintah::DataWarehouse* old_dw);

  };
} // End namespace Uintah


#endif  // __ARENISCA3_PARTIALLY_SATURATED__
