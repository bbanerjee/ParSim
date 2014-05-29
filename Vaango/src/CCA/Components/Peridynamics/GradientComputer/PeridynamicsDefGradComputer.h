#ifndef __VAANGO_PERIDYNAMICS_DEFORMATION_GRADIENT_COMPUTER_H__
#define __VAANGO_PERIDYNAMICS_DEFORMATION_GRADIENT_COMPUTER_H__

#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Task.h>

#include <vector>

namespace Uintah {
  class Patch;
  class DataWarehouse;
  class VarLabel;
}

namespace Vaango {

  class PeridynamicsFlags;
  class PeridynamicsMaterial;
  class PeridynamicsLabel;

  class PeridynamicsDefGradComputer {

  public:

    PeridynamicsDefGradComputer(PeridynamicsFlags* flags, 
                                PeridynamicsLabel* labels);
    virtual ~PeridynamicsDefGradComputer();

    void addInitialComputesAndRequires(Uintah::Task* task,
                                       const PeridynamicsMaterial* matl,
                                       const Uintah::PatchSet*);

    void addComputesAndRequires(Uintah::Task* task,
                                const PeridynamicsMaterial* matl,
                                const Uintah::PatchSet*);

    void initialize(const Uintah::Patch* patch,    
                    PeridynamicsMaterial* matl, 
                    Uintah::DataWarehouse* new_dw);

    void computeDeformationGradient(const Uintah::Patch* patch,
                                    const PeridynamicsMaterial* matl,
                                    Uintah::DataWarehouse* old_dw,
                                    Uintah::DataWarehouse* new_dw);

   void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                         std::vector<const Uintah::VarLabel*>& to);

  private:

    PeridynamicsLabel* d_labels;
    PeridynamicsFlags* d_flags;
  };

} // end namespace Vaango

#endif
