#ifndef __VAANGO_PERIDYNAMICS_LINEAR_ELASTIC_BOND_MODEL_H__
#define __VAANGO_PERIDYNAMICS_LINEAR_ELASTIC_BOND_MODEL_H__

#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModel.h>

namespace Vaango {

  class LinearElasticBondModel : public PeridynamicsMaterialModel {

  public:
         
    LinearElasticBondModel(Uintah::ProblemSpecP& ps,
                           PeridynamicsFlags* flags);
    LinearElasticBondModel(const LinearElasticBondModel* cm);
    virtual ~LinearElasticBondModel();

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

    void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                          std::vector<const Uintah::VarLabel*>& to);

    // Make a clone of the constitutive model
    LinearElasticBondModel* clone();
  };

} // End namespace Vaango
      

#endif  // __VAANGO_PERIDYNAMICS_LINEAR_ELASTIC_BOND_MODEL_H__

