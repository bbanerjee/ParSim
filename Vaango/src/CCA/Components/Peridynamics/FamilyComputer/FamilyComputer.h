#ifndef __VAANGO_FAMILY_COMPUTER_H__
#define __VAANGO_FAMILY_COMPUTER_H__

#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Task.h>

#include <Core/Thread/CrowdMonitor.h>

#include <vector>
#include <map>

namespace Uintah {
  class Patch;
  class DataWarehouse;
  class ParticleSubset;
  class VarLabel;
}

namespace Vaango {

  typedef int particleIndex;
  typedef int particleId;

  class PeridynamicsFlags;
  class PeridynamicsMaterial;
  class PeridynamicsLabel;

  class FamilyComputer {

  public:
    
    FamilyComputer(PeridynamicsFlags* flags, PeridynamicsLabel* labels);
    virtual ~FamilyComputer();

    void createNeighborList(PeridynamicsMaterial* matl,
                            const Uintah::Patch* patch,
                            Uintah::DataWarehouse* new_dw);

  protected:

    PeridynamicsLabel* d_varLabel;
    PeridynamicsFlags* d_flags;
    
    mutable Uintah::CrowdMonitor d_lock;
  };

} // End of namespace Vaango

#endif // __VAANGO_FAMILY_COMPUTER_H__
