#ifndef __VAANGO_PERIDYNAMICS_LINEAR_ELASTIC_STATE_MODEL_H__
#define __VAANGO_PERIDYNAMICS_LINEAR_ELASTIC_STATE_MODEL_H__

#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModel.h>

#include <Core/Grid/Variables/ComputeSet.h>    // For PatchSubset

/*!  The is essentially an isotropic compressible Neo-Hookean material model */

namespace Vaango {

  class LinearElasticStateModel : public PeridynamicsMaterialModel {

  public:

    struct CMData {
      double bulkModulus;
      double shearModulus;
    };

  private:

    CMData d_cm;  // Constitutive model data

    // May be needed at a later stage (BB: 5/15/14)
    // friend const TypeDescription* fun_getTypeDescription(CMData* cm);

    // Prevent assignment
    LinearElasticStateModel& operator=(const LinearElasticStateModel& cm);

  public:
         
    // Default constructor
    LinearElasticStateModel(Uintah::ProblemSpecP& ps,
                           PeridynamicsFlags* flags);

    // Copy constructor
    LinearElasticStateModel(const LinearElasticStateModel* cm);

    // Make a clone of the constitutive model
    LinearElasticStateModel* clone();

    virtual ~LinearElasticStateModel();

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
      

#endif  // __VAANGO_PERIDYNAMICS_LINEAR_ELASTIC_STATE_MODEL_H__
