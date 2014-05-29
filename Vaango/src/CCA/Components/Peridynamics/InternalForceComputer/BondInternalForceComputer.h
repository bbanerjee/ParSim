#ifndef __VAANGO_BOND_INTERNAL_FORCE_COMPUTER_H__
#define __VAANGO_BOND_INTERNAL_FORCE_COMPUTER_H__

#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/ComputeSet.h>

namespace Uintah {
  class Patch;
  class DataWarehouse;
}

namespace Vaango {

  class PeridynamicsFlags;
  class PeridynamicsMaterial;
  class PeridynamicsLabel;

  class BondInternalForceComputer {

  public:
    
    BondInternalForceComputer(PeridynamicsFlags* flags, PeridynamicsLabel* labels);
    virtual ~BondInternalForceComputer();

    /*! Computes and requires for the internal force computer */
    void addComputesAndRequires(Uintah::Task* task,
                                const PeridynamicsMaterial* matl,
                                const Uintah::PatchSet* patches) const;

    /*! Actually compute the internal force */
    void computeInternalForce(const Uintah::PatchSubset* patches,
                              const PeridynamicsMaterial* matl,
                              Uintah::DataWarehouse* old_dw,
                              Uintah::DataWarehouse* new_dw);

  private:

    PeridynamicsLabel* d_label;
    PeridynamicsFlags* d_flags;
    
  };

} // End of namespace Vaango

#endif // __VAANGO_INTERNAL_FORCE_COMPUTER_H__
