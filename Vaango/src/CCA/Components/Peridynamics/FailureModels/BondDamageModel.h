#ifndef __VAANGO_PERIDYNAMICS_BOND_DAMAGE_MODEL_H__
#define __VAANGO_PERIDYNAMICS_BOND_DAMAGE__MODELH__

#include <CCA/Components/Peridynamics/FailureModels/PeridynamicsFailureModel.h>

#include <vector>

namespace Vaango {

  class BondDamageModel: public PeridynamicsFailureModel {

  public:
         
    BondDamageModel(Uintah::ProblemSpecP& ps,
                    PeridynamicsFlags* flags);
    BondDamageModel(const BondDamageModel* cm);
    virtual ~BondDamageModel();

    void outputProblemSpec(Uintah::ProblemSpecP& ps,
                           bool output_cm_tag = true);
         
    /*! Initial computes and requires for the failure model */
    void addInitialComputesAndRequires(Uintah::Task* task,
                                       const PeridynamicsMaterial* matl,
                                       const Uintah::PatchSet* patches) const;

    /*! Initialize the variables used in the CM */
    void initialize(const Uintah::Patch* patch,
                    const PeridynamicsMaterial* matl,
                    Uintah::DataWarehouse* new_dw);

    /*! Set up the computes and requires for the task  */
    void addComputesAndRequires(Uintah::Task* task,
                                const PeridynamicsMaterial* matl,
                                const Uintah::PatchSet* patches) const;

    // Basic failure model calculations
    void updateDamage(const Uintah::PatchSubset* patches,
                      const PeridynamicsMaterial* matl,
                      Uintah::DataWarehouse* old_dw,
                      Uintah::DataWarehouse* new_dw);
                                     
    void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                          std::vector<const Uintah::VarLabel*>& to);

    BondDamageModel* clone();

  };

} // End namespace Vaango
      

#endif  // __VAANGO_PERIDYNAMICS_BOND_DAMAGE_MODEL_H__

