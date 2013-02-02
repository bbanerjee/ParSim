#ifndef __VAANGO_BASIC_DAMAGE_MODEL_H__
#define __VAANGO_BASIC_DAMAGE_MODEL_H__

#include <Core/Grid/Variables/ComputeSet.h>
#include <vector>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>
#include <Core/Containers/StaticArray.h>
#include <Core/Grid/SimulationStateP.h>
#include <Core/Grid/Variables/Array3.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/LinearInterpolator.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/Math/FastMatrix.h>
#include <CCA/Components/MPM/MPMFlags.h>


namespace Vaango {

  class Uintah::Task;
  class Uintah::Patch;
  class Uintah::VarLabel;
  class Uintah::MPMLabel;
  class Uintah::MPMFlags;
  class Uintah::MPMMaterial;
  class Uintah::DataWarehouse;
  class Uintah::ParticleSubset;
  class Uintah::ParticleVariableBase;

  //////////////////////////////////////////////////////////////////////////
  /*!
    \class BasicDamageModel
   
    \brief Base class for the default damage models.

  */
  //////////////////////////////////////////////////////////////////////////

  class BasicDamageModel {
  public:
         
    BasicDamageModel(MPMFlags* MFlag);
    BasicDamageModel(const BasicDamageModel* cm);
    virtual ~BasicDamageModel();

    virtual void getDamageModelData(Uintah::ProblemSpecP& ps);
    virtual void setDamageModelData(BasicDamageModel* cm);

    virtual void getFailureStressOrStrainData(Uintah::ProblemSpecP& ps);
    virtual void getBrittleDamageData(Uintah::ProblemSpecP& ps);
    virtual void outputProblemSpecDamage(Uintah::ProblemSpecP& cm_ps);

    virtual void setFailureStressOrStrainData(const BasicDamageModel* cm);
    virtual void setBrittleDamageData(const BasicDamageModel* cm);
    virtual void setErosionAlgorithm();
    virtual void setErosionAlgorithm(const BasicDamageModel* cm);

    virtual void initializeDamageVarLabels();
    virtual void deleteDamageVarLabels();

    virtual void initializeDamageData(const Uintah::Patch* patch,
                                      const Uintah::MPMMaterial* matl,
                                      Uintah::DataWarehouse* new_dw);

    virtual void copyDamageDataFromDeletedToAddedParticle(Uintah::DataWarehouse* new_dw,
                                                          Uintah::ParticleSubset* addset,
                                                          map<const Uintah::VarLabel*,
                                                            Uintah::ParticleVariableBase*>* newState,
                                                          Uintah::ParticleSubset* delset,
                                                          Uintah::DataWarehouse* old_dw );

    virtual void allocateDamageDataAddRequires(Uintah::Task* task,
                                               const Uintah::MPMMaterial* matl,
                                               const Uintah::PatchSet* patches,
                                               Uintah::MPMLabel*lb ) const;

    virtual void carryForwardDamageData(Uintah::ParticleSubset* pset,
                                        Uintah::DataWarehouse*  old_dw,
                                        Uintah::DataWarehouse*  new_dw,
                                        const Uintah::MPMMaterial* matl);

    virtual void addComputesAndRequiresForDamage(Uintah::Task* task,
                                                 const Uintah::MPMMaterial* matl,
                                                 const Uintah::PatchSet* patches) const;
        
    virtual void addInitialComputesAndRequiresForDamage(Uintah::Task* task,
                                                        const Uintah::MPMMaterial* matl,
                                                        const Uintah::PatchSet* patches) const;

    virtual void addRequiresDamageParameterDefault(Uintah::Task* task,
                                                   const Uintah::MPMMaterial* matl,
                                                   const Uintah::PatchSet* patches);

    virtual void addParticleStateDamage(std::vector<const Uintah::VarLabel*>& from,
                                        std::vector<const Uintah::VarLabel*>& to);

    virtual void updateDamageAndModifyStress(const Uintah::Matrix3& defGrad, 
                                             const double& pFailureStrain, 
	        			     double& pFailureStrain_new,
                                             const double& pVolume, 
                                             const double& pDamage,
                                             double& pDamage_new, 
                                             Uintah::Matrix3& pStress,
	        			     const long64 particleID);

    virtual void updateFailedParticlesAndModifyStress(const Uintah::Matrix3& defGrad,
                                                      const double& pFailureStr,
                                                      const int& pLocalized,
                                                      int& pLocalized_new,
                                                      const double& pTimeOfLoc,
                                                      double& pTimeOfLoc_new,
                                                      Uintah::Matrix3& pStress,
                                                      const long64 particleID,
                                                      double time);

  protected:

    MPMFlags* flag;

    // Damage Requirements //
    /////////////////////////
    // Create datatype for failure strains
    struct FailureStressOrStrainData {
      double mean;         /* Mean failure stress, strain or cohesion */
      double std;          /* Standard deviation of failure strain */
                           /* or Weibull modulus */
      double exponent;     /* Exponent used in volume scaling of failure crit */
      double refVol;       /* Reference volume for scaling failure criteria */
      std::string scaling; /* Volume scaling method: "none" or "kayenta" */
      std::string dist;    /* Failure distro: "constant", "gauss" or "weibull"*/
      int seed;            /* seed for random number distribution generator */
      double t_char;       /* characteristic time for damage to occur */
    };

    //Create datatype for brittle damage
    struct BrittleDamageData {
      double r0b;          /* Initial energy threshold (\sqrt{Pa}) */
      double Gf;           /* Fracture energy (J/m^3) */
      double constant_D;   /* Shape factor in softening function */
      double maxDamageInc; /* Maximum damage increment in a time step */
      bool allowRecovery;  /* Recovery of stiffness allowed */
      double recoveryCoeff;  /* Fraction of stiffness to be recovered */
      bool printDamage;    /* Flag to print damage */
    };

    const VarLabel* pFailureStressOrStrainLabel;
    const VarLabel* pLocalizedLabel;
    const VarLabel* pDamageLabel;
    const VarLabel* pTimeOfLocLabel;
    const VarLabel* pFailureStressOrStrainLabel_preReloc;
    const VarLabel* pLocalizedLabel_preReloc;
    const VarLabel* pDamageLabel_preReloc;
    const VarLabel* pTimeOfLocLabel_preReloc;
    
    FailureStressOrStrainData d_epsf;
    BrittleDamageData d_brittle_damage;

    // Erosion algorithms
    bool d_setStressToZero; /* set stress tensor to zero*/
    bool d_allowNoTension;  /* retain compressive mean stress after failue*/
    bool d_allowNoShear;    /* retain mean stress after failure - no deviatoric stress */
                            /* i.e., no deviatoric stress */
    bool d_brittleDamage;   /* use brittle damage with mesh size control*/

    std::string d_failure_criteria; /* Options are:  "MaximumPrincipalStrain" */
                                    /* "MaximumPrincipalStress", "MohrColoumb"*/

    // These three are for the MohrColoumb option
    double d_friction_angle;  // Assumed to come in degrees
    double d_tensile_cutoff;  // Fraction of the cohesion at which 
                              // tensile failure occurs

  };
} // End namespace Uintah
      


#endif  // __VAANGO_BASIC_DAMAGE_MODEL_H__

