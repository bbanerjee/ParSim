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

    /*! Initial computes and requires for the family computer */
    void addInitialComputesAndRequires(Uintah::Task* task,
                                       const PeridynamicsMaterial* matl,
                                       const Uintah::PatchSet* patches) const;

    void createNeighborList(PeridynamicsMaterial* matl,
                            const Uintah::Patch* patch,
                            Uintah::DataWarehouse* new_dw);

  protected:

    void findCellsInHorizon(const Uintah::Patch* patch,
                            const SCIRun::Point& pos,
                            const double& horizon,
                            SCIRun::IntVector& cellLow,
                            SCIRun::IntVector& cellHigh);

  private:

    PeridynamicsLabel* d_label;
    PeridynamicsFlags* d_flags;
    
  };

} // End of namespace Vaango

#endif // __VAANGO_FAMILY_COMPUTER_H__
