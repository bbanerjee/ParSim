#ifndef __VAANGO_PERIDYNAMICS_SPHERICAL_STRAIN_ENERGY_DAMAGE_MODEL_H__
#define __VAANGO_PERIDYNAMICS_SPHERICAL_STRAIN_ENERGY_DAMAGE_MODEL_H__

#include <CCA/Components/Peridynamics/DamageModels/PeridynamicsDamageModel.h>

namespace Vaango {

  class SphericalStrainEnergyDamageModel : public PeridynamicsDamageModel {

  public:
         
    /*! Default constructor */
    SphericalStrainEnergyDamageModel(Uintah::ProblemSpecP& ps,
                                     PeridynamicsLabel* labels, 
                                     PeridynamicsFlags* flags);

    /*! Copy constructor */
    SphericalStrainEnergyDamageModel(const SphericalStrainEnergyDamageModel* cm);

    /*! Make a clone of the constitutive model */
    SphericalStrainEnergyDamageModel* clone();

    /*! Destructor */
    virtual ~SphericalStrainEnergyDamageModel();

    /*! Write out the input data for restart purposes */
    void outputProblemSpec(Uintah::ProblemSpecP& ps,
                           bool output_cm_tag = true);

    /*! Set up initial computes and requires for the damage model compute task */
    void addInitialComputesAndRequires(Uintah::Task* task,
                                       const PeridynamicsMaterial* matl,
                                       const Uintah::PatchSet* patches) const;

    /*! Actually initialize the variables used in the damage models */
    void initialize(const Uintah::Patch* patch,
                    const PeridynamicsMaterial* matl,
                    Uintah::DataWarehouse* new_dw);

    /*! Set up the computes and requires for the task that computes the
      damage tensor */
    void addComputesAndRequires(Uintah::Task* task,
                                const PeridynamicsMaterial* matl,
                                const Uintah::PatchSet* patches) const;

    /*! Actually compute damage tensor */
    void computeDamageTensor(const Uintah::PatchSubset* patches,
                             const PeridynamicsMaterial* matl,
                             Uintah::DataWarehouse* old_dw,
                             Uintah::DataWarehouse* new_dw);
                                     
    /*! Tag variables that are to be relocated */
    void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                          std::vector<const Uintah::VarLabel*>& to);


  private:

    double d_GIc;  // Critical mode I strain energy release rate

    // Prevent assignment
    SphericalStrainEnergyDamageModel& operator=(const SphericalStrainEnergyDamageModel& cm);

  };

} // End namespace Vaango
      

#endif  // __VAANGO_PERIDYNAMICS_SPHERICAL_STRAIN_ENERGY_DAMAGE_MODEL_H__

