#ifndef __VAANGO_PERIDYNAMICS_MATERIAL_MODEL_H__
#define __VAANGO_PERIDYNAMICS_MATERIAL_MODEL_H__

#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <Core/Grid/SimulationStateP.h>

#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Containers/StaticArray.h>
#include <Core/Grid/Variables/Array3.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Parallel/ProcessorGroup.h>

#include <vector>

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

  class PeridynamicsMaterialModel {

  public:
         
    PeridynamicsMaterialModel(PeridynamicsFlags* flags);
    PeridynamicsMaterialModel(const PeridynamicsMaterialModel* cm);
    virtual ~PeridynamicsMaterialModel();

    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps,
                                   bool output_cm_tag = true) = 0;

         
    // Basic constitutive model calculations
    virtual void computeStressTensor(const Uintah::PatchSubset* patches,
                                     const PeridynamicsMaterial* matl,
                                     Uintah::DataWarehouse* old_dw,
                                     Uintah::DataWarehouse* new_dw);
                                     
    // Basic constitutive model calculations
    virtual void computeInternalForce(const Uintah::PatchSubset* patches,
                                      const PeridynamicsMaterial* matl,
                                      Uintah::DataWarehouse* old_dw,
                                      Uintah::DataWarehouse* new_dw);
                                     
    /*! Initial computes and requires for the constitutive model */
    virtual void addInitialComputesAndRequires(Uintah::Task* task,
                                               const PeridynamicsMaterial* matl,
                                               const Uintah::PatchSet* patches) const;

    /*! Initialize the variables used in the CM */
    virtual void initialize(const Uintah::Patch* patch,
                            const PeridynamicsMaterial* matl,
                            Uintah::DataWarehouse* new_dw) = 0;

    /*! Set up the computes and requires for the task that computes the
        stress tensor and associated kinematic and thermal quantities */
    virtual void addComputesAndRequires(Uintah::Task* task,
                                        const PeridynamicsMaterial* matl,
                                        const Uintah::PatchSet* patches) const;

    virtual void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                                  std::vector<const Uintah::VarLabel*>& to) = 0;


    inline void setWorld(const Uintah::ProcessorGroup* myworld)
      {
        d_world = myworld;
      }

    inline void setSharedState(Uintah::SimulationState* sharedState)
      {
        d_sharedState = sharedState;
      }

    // Make a clone of the constitutive model
    virtual PeridynamicsMaterialModel* clone() = 0;

  protected:

    PeridynamicsLabel* d_label;
    PeridynamicsFlags* d_flag;
    int d_numGhostNodes;

    const Uintah::ProcessorGroup* d_world;

    // don't store SimulationStateP or it will add a reference 
    // that will never be removed
    Uintah::SimulationState* d_sharedState;

  };

} // End namespace Vaango
      

#endif  // __VAANGO_PERIDYNAMICS_MATERIAL_MODEL_H__

