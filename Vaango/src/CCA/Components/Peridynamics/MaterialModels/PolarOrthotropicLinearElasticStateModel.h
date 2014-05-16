#ifndef __VAANGO_PERIDYNAMICS_POLAR_ORTHOTROPIC_LINEAR_ELASTIC_STATE_MODEL_H__
#define __VAANGO_PERIDYNAMICS_POLAR_ORTHOTROPIC_LINEAR_ELASTIC_STATE_MODEL_H__

#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModel.h>

#include <Core/Grid/Variables/ComputeSet.h>    // For PatchSubset

/*!  A polar orthotropic linear elastic material model which uses a Green Naghdi stress rate */

namespace Vaango {

  class PolarOrthotropicLinearElasticStateModel : public PeridynamicsMaterialModel {

  public:

    struct CMData {
      double bulkModulus;
      double shearModulus;
    };

  private:

    CMData d_cm;  // Constitutive model data

    // Prevent assignment
    PolarOrthotropicLinearElasticStateModel& operator=(const PolarOrthotropicLinearElasticStateModel& cm);

  public:
         
    // Default constructor
    PolarOrthotropicLinearElasticStateModel(Uintah::ProblemSpecP& ps,
                                            PeridynamicsFlags* flags);

    // Copy constructor
    PolarOrthotropicLinearElasticStateModel(const PolarOrthotropicLinearElasticStateModel* cm);

    // Make a clone of the constitutive model
    PolarOrthotropicLinearElasticStateModel* clone();

    virtual ~PolarOrthotropicLinearElasticStateModel();

    /*!  Output the problem spec for restart */
    void outputProblemSpec(Uintah::ProblemSpecP& ps,
                           bool output_cm_tag = true);

    /*! Initial computes and requires for the constitutive model */
    void addInitialComputesAndRequires(Uintah::Task* task,
                                       const PeridynamicsMaterial* matl,
                                       const Uintah::PatchSet* patches) const;

    /*! Initialize the variables used in the CM */
    void initialize(const Uintah::Patch* patch,
                    const PeridynamicsMaterial* matl,
                    Uintah::DataWarehouse* new_dw);

    /*! Set up the computes and requires for the task */
    void addComputesAndRequires(Uintah::Task* task,
                                const PeridynamicsMaterial* matl,
                                const Uintah::PatchSet* patches) const;

    /*! Compute the stress tensor */
    void computeStress(const Uintah::PatchSubset* patches,
                       const PeridynamicsMaterial* matl,
                       Uintah::DataWarehouse* old_dw,
                       Uintah::DataWarehouse* new_dw);

    void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                          std::vector<const Uintah::VarLabel*>& to);

  };
} // End namespace Vaango
      

#endif  // __VAANGO_PERIDYNAMICS_POLAR_ORTHOTROPIC_LINEAR_ELASTIC_STATE_MODEL_H__

