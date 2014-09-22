#ifndef __VAANGO_MPM_POLAR_ORTHOTROPIC_HYPOELASTIC_MODEL_H__
#define __VAANGO_MPM_POLAR_ORTHOTROPIC_HYPOELASTIC_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>

#include <Core/Grid/Variables/ComputeSet.h>    // For PatchSubset
#include <Core/Math/SymmMatrix6.h>

/*!  A polar orthotropic linear elastic material model which uses a Green Naghdi stress rate */

namespace Uintah {

  class PolarOrthotropicHypoElastic : public ConstitutiveModel {

  public:

    struct CMData {
      SCIRun::Point top;
      SCIRun::Point bottom;
      double Er;
      double Etheta;
      double Ez;
      double nuthetar;
      double nuzr;
      double nuztheta;
      double Gthetaz;
      double Gzr;
      double Grtheta;
      Uintah::SymmMatrix6 stiffnessMatrix;
    };

    const VarLabel* pRCoordLabel;
    const VarLabel* pThetaCoordLabel;
    const VarLabel* pZCoordLabel;
    const VarLabel* pRCoordLabel_preReloc;
    const VarLabel* pThetaCoordLabel_preReloc;
    const VarLabel* pZCoordLabel_preReloc;

  private:

    CMData d_cm;  // Constitutive model data

    // Prevent assignment
    PolarOrthotropicHypoElastic& operator=(const PolarOrthotropicHypoElastic& cm);

  public:
         
    // Default constructor
    PolarOrthotropicHypoElastic(Uintah::ProblemSpecP& ps,
                                MPMFlags* flags);

    // Copy constructor
    PolarOrthotropicHypoElastic(const PolarOrthotropicHypoElastic* cm);

    // Make a clone of the constitutive model
    PolarOrthotropicHypoElastic* clone();

    // Destroy
    virtual ~PolarOrthotropicHypoElastic();

    /*!  Output the problem spec for restart */
    void outputProblemSpec(Uintah::ProblemSpecP& ps,
                           bool output_cm_tag = true);

    /*! Initial computes and requires for the constitutive model */
    void addInitialComputesAndRequires(Uintah::Task* task,
                                       const MPMMaterial* matl,
                                       const Uintah::PatchSet* patches) const;

    /*! Initialize the variables used in the CM */
    void initializeCMData(const Uintah::Patch* patch,
                          const MPMMaterial* matl,
                          Uintah::DataWarehouse* new_dw);

    /*! Compute a stable initial timestep */
    void computeStableTimestep(const Uintah::Patch* patch,
                               const MPMMaterial* matl,
                               Uintah::DataWarehouse* new_dw);

    /*! Set up the computes and requires for the task */
    void addComputesAndRequires(Uintah::Task* task,
                                const MPMMaterial* matl,
                                const Uintah::PatchSet* patches) const;

    /*! Compute the stress tensor */
    void computeStressTensor(const Uintah::PatchSubset* patches,
                             const MPMMaterial* matl,
                             Uintah::DataWarehouse* old_dw,
                             Uintah::DataWarehouse* new_dw);

    /*! Add particle state for variables that have to be relocated
        when particles move to another patch at the end of a time step */
    void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                          std::vector<const Uintah::VarLabel*>& to);

    /*! Set up computes and requires for implicit time integration.
        @todo:  This task has not been implemented yet. 
     */
    void addComputesAndRequires(Task* task,
                                const MPMMaterial* matl,
                                const PatchSet* patches,
                                const bool recursion,
                                const bool schedPar = true) const;

    /*! Compute grid cell microscopic density for MPMICE calculations.
        @todo:  This task has not been implemented yet. 
     */
    double computeRhoMicroCM(double pressure,
                             const double p_ref,
                             const MPMMaterial* matl, 
                             double temperature,
                             double rho_guess);

    /*! Compute grid cell pressure using an equation of state 
        for MPMICE calculations.
        @todo:  This task has not been implemented yet. 
     */
    void computePressEOSCM(double rho_m, 
                           double& press_eos,
                           double p_ref,
                           double& dp_drho, double& ss_new,
                           const MPMMaterial* matl, 
                           double temperature);

    /*! Get the compressiblity (inverse of bulk modulus) 
        for MPMICE calculations.
        @todo:  This task has not been implemented yet. 
     */
    double getCompressibility();

    /*! Carry forward CM data for RigidMPM */
    void carryForward(const PatchSubset* patches,
                      const MPMMaterial* matl,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw);

    /*! Set up task variables for situations where particles are moved to
        another material type */
    void allocateCMDataAddRequires(Task* task, 
                                   const MPMMaterial* matl,
                                   const PatchSet* patch, 
                                   MPMLabel* lb) const;

    /*! Allocate variables for situations where particles are to be 
        transformed into a different type of material */
    void allocateCMDataAdd(DataWarehouse* new_dw,
                           ParticleSubset* subset,
                           map<const VarLabel*, ParticleVariableBase*>* newState,
                           ParticleSubset* delset,
                           DataWarehouse* old_dw);

  };
} // End namespace Uintah
      

#endif  // __VAANGO_MPM_POLAR_ORTHOTROPIC_HYPOELASTIC_MODEL_H__

