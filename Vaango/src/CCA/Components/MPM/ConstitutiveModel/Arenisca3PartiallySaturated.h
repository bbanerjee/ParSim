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

    // Create datatype for storing model parameters
    struct CMData {
      double subcycling_characteristic_number;
      bool   use_disaggregation_algorithm;
      double K0_Murnaghan_EOS;
      double n_Murnaghan_EOS;
    };

    const Uintah::VarLabel* pLocalizedLabel;
    const Uintah::VarLabel* pLocalizedLabel_preReloc;
    const Uintah::VarLabel* pAreniscaFlagLabel;          //0: ok, 1: pevp<-p3
    const Uintah::VarLabel* pAreniscaFlagLabel_preReloc;
    const Uintah::VarLabel* pScratchDouble1Label;
    const Uintah::VarLabel* pScratchDouble1Label_preReloc;
    const Uintah::VarLabel* pScratchDouble2Label;
    const Uintah::VarLabel* pScratchDouble2Label_preReloc;
    const Uintah::VarLabel* pElasticVolStrainLabel;              //Elastic Volumetric Strain
    const Uintah::VarLabel* pElasticVolStrainLabel_preReloc;
    const Uintah::VarLabel* pStressQSLabel;
    const Uintah::VarLabel* pStressQSLabel_preReloc;
    const Uintah::VarLabel* pScratchMatrixLabel;
    const Uintah::VarLabel* pScratchMatrixLabel_preReloc;

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

    // Prevent copying of this class
    // copy constructor

    Arenisca3PartiallySaturated& operator=(const Arenisca3PartiallySaturated &cm);

    void initializeLocalMPMLabels();

  public:
    // constructor
    Arenisca3PartiallySaturated(Uintah::ProblemSpecP& ps, Uintah::MPMFlags* flag);
    Arenisca3PartiallySaturated(const Arenisca3PartiallySaturated* cm);

    // destructor
    virtual ~Arenisca3PartiallySaturated();

    virtual void outputProblemSpec(SCIRun::ProblemSpecP& ps,bool output_cm_tag = true);

    // clone
    Arenisca3PartiallySaturated* clone();

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

    bool computeStep(Uintah::particleIndex idx,
                     Uintah::long64 particleID, 
                     const Uintah::Matrix3& D,
                     const double & dt,
                     const ModelState_MasonSand& state_n,
                     const double & coher,
                     const double & P3,
                     ModelState_MasonSand& state_p,
                     Uintah::long64 ParticleID);

    void computeElasticProperties(double & bulk,
                                  double & shear);

    void computeElasticProperties(const ModelState_MasonSand& state,
                                  const double& P3,
                                  double & bulk,
                                  double & shear
                                 );

    void computeElasticProperties(const Uintah::Matrix3& sigma,
                                  const Uintah::Matrix3& ep,
                                  const double& P3,
                                  double & bulk,
                                  double & shear
                                 );

    Uintah::Matrix3 computeTrialStress(const Uintah::Matrix3& sigma_old,  // old stress
                               const Uintah::Matrix3& d_e,        // Strain increment
                               const double& bulk,        // bulk modulus
                               const double& shear);      // shear modulus

    int computeStepDivisions(Uintah::particleIndex idx,
                             Uintah::long64 particleID,
                             const ModelState_MasonSand& state,
                             const double& P3,
                             const Uintah::Matrix3& sigma_trial);

    bool computeSubstep(Uintah::particleIndex idx,
                        Uintah::long64 particleID,
                        const Uintah::Matrix3& d_e,             // total strain increment for substep
                        const ModelState_MasonSand& state_old, // state at start of substep
                        const double & coher,           // scalar valued coher
                        const double & P3,              // initial disaggregation strain
                        ModelState_MasonSand& state_new        // state at end of substep
                       );

    double computeX(const double& evp, const double& P3);

    double computedZetadevp(double Zeta,
                            double evp);

    double computePorePressure(const double ev);

    int nonHardeningReturn(const ModelState_MasonSand& invar_trial,
                           const ModelState_MasonSand& invar_old,
                           const Uintah::Matrix3& d_e,
                           const ModelState_MasonSand& state,
                           const double& coher,
                           const double& bulk,
                           const double& shear,
                           ModelState_MasonSand& invar_new,
                           Uintah::Matrix3& d_ep_new,
                           double& kappa);

    void transformedBisection(double& z_0,
                              double& r_0,
                              const double& z_trial,
                              const double& r_trial,
                              const ModelState_MasonSand& state,
                              const double& coher,
                              const double  limitParameters[4], // XXX
                              const double& r_to_rJ2,
                              double& kappa
                             );

    int transformedYieldFunction(const double& z,
                                 const double& r,
                                 const ModelState_MasonSand& state,
                                 const double& coher,
                                 const double  limitParameters[4], // XXX
                                 const double& r_to_rJ2,
                                 double& kappa
                                );

    int computeYieldFunction(const ModelState_MasonSand& invariants,
                             const ModelState_MasonSand& state,
                             const double& coher,
                             const double  limitParameters[4],
                             double& kappa // XXX
                            );

    int computeYieldFunction(const double& I1,
                             const double& rJ2,
                             const ModelState_MasonSand& state,
                             const double& coher,
                             const double  limitParameters[4],
                             double& kappa // XXX
                            );

    void computeLimitParameters(double *limitParameters,
                                const double& coher //XXX
                               );
    void checkInputParameters();



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
