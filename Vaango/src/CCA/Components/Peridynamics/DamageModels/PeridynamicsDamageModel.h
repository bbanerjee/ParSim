#ifndef __VAANGO_PERIDYNAMICS_DAMAGE_MODEL_H__
#define __VAANGO_PERIDYNAMICS_DAMAGE_MODEL_H__

#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <Core/Grid/SimulationStateP.h>
#include <Core/Grid/Variables/ComputeSet.h>

namespace Uintah {
  class Task;
  class Patch;
  class VarLabel;
  class DataWarehouse;
  class ParticleSubset;
  class ParticleVariableBase;
}

namespace Vaango {

  class PeridynamicsLabel;
  class PeridynamicsFlags;
  class PeridynamicsMaterial;

  class PeridynamicsDamageModel {

  public:
         
    /*! Default constructor */
    PeridynamicsDamageModel(PeridynamicsLabel* labels, PeridynamicsFlags* flags);

    /*! Copy constructor */
    PeridynamicsDamageModel(const PeridynamicsDamageModel* cm);

    /*! Make a clone of the constitutive model */
    virtual PeridynamicsDamageModel* clone() = 0;

    /*! Destructor */
    virtual ~PeridynamicsDamageModel();

    /*! Write out the input data for restart purposes */
    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps,
                                   bool output_cm_tag = true) = 0;

    /*! Set up initial computes and requires for the damage model compute task */
    virtual void addInitialComputesAndRequires(Uintah::Task* task,
                                               const PeridynamicsMaterial* matl,
                                               const Uintah::PatchSet* patches) const;

    /*! Actually initialize the variables used in the damage models */
    virtual void initialize(const Uintah::Patch* patch,
                            const PeridynamicsMaterial* matl,
                            Uintah::DataWarehouse* new_dw) = 0;

    /*! Set up the computes and requires for the task that computes the
      damage tensor */
    virtual void addComputesAndRequires(Uintah::Task* task,
                                        const PeridynamicsMaterial* matl,
                                        const Uintah::PatchSet* patches) const;

    /*! Actually compute damage tensor */
    virtual void computeDamageTensor(const Uintah::PatchSubset* patches,
                                     const PeridynamicsMaterial* matl,
                                     Uintah::DataWarehouse* old_dw,
                                     Uintah::DataWarehouse* new_dw) = 0;
                                     
    /*! Tag variables that are to be relocated */
    virtual void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                                  std::vector<const Uintah::VarLabel*>& to) = 0;


  protected:

    PeridynamicsLabel* d_label;
    PeridynamicsFlags* d_flag;

  };

} // End namespace Vaango
      

#endif  // __VAANGO_PERIDYNAMICS_DAMAGE_MODEL_H__

