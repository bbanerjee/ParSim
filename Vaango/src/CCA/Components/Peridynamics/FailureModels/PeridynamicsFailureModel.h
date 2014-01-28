#ifndef __VAANGO_PERIDYNAMICS_FAILURE_MODEL_H__
#define __VAANGO_PERIDYNAMICS_FAILURE_MODEL_H__

#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsSimulationStateP.h>

#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Containers/StaticArray.h>
#include <Core/Grid/Variables/Array3.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Parallel/ProcessorGroup.h>

#include <vector>

namespace Vaango {

  class PeridynamicsLabel;
  class PeridynamicsFlags;
  class PeridynamicsMaterial;

  class Uintah::Task;
  class Uintah::Patch;
  class Uintah::VarLabel;
  class Uintah::DataWarehouse;
  class Uintah::ParticleSubset;
  class Uintah::ParticleVariableBase;

  class PeridynamicsFailureModel {

  public:
         
    PeridynamicsFailureModel(PeridynamicsFlags* flags);
    PeridynamicsFailureModel(const PeridynamicsFailureModel* cm);
    virtual ~PeridynamicsFailureModel();

    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps,
                                   bool output_cm_tag = true) = 0;

         
    /*! Initial computes and requires for the failure model */
    virtual void addInitialComputesAndRequires(Uintah::Task* task,
                                               const PeridynamicsMaterial* matl,
                                               const Uintah::PatchSet* patches) const;

    /*! Initialize the variables used in the CM */
    virtual void initializeDamage(const Uintah::Patch* patch,
                                  const PeridynamicsMaterial* matl,
                                  Uintah::DataWarehouse* new_dw) = 0;

    /*! Set up the computes and requires for the task  */
    virtual void addComputesAndRequires(Uintah::Task* task,
                                        const PeridynamicsMaterial* matl,
                                        const Uintah::PatchSet* patches) const;

    // Basic failure model calculations
    virtual void updateDamage(const Uintah::PatchSubset* patches,
                               const PeridynamicsMaterial* matl,
                               Uintah::DataWarehouse* old_dw,
                               Uintah::DataWarehouse* new_dw);
                                     
    virtual void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                                  std::vector<const Uintah::VarLabel*>& to) = 0;

    // Make a clone of the failure model
    virtual PeridynamicsFailureModel* clone() = 0;

  protected:

    PeridynamicsLabel* d_varLabel;
    PeridynamicsFlags* d_flag;

  };

} // End namespace Vaango
      

#endif  // __VAANGO_PERIDYNAMICS_FAILURE_MODEL_H__

