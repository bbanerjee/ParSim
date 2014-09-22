#ifndef __VAANGO_NULL_CONTACT_H__
#define __VAANGO_NULL_CONTACT_H__

#include <CCA/Components/Peridynamics/ContactModels/ContactModelBase.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Grid/SimulationState.h>
#include <Core/Grid/SimulationStateP.h>

namespace Vaango {

  /////////////////////////////////////////////////////////////////////////////
  /*!
    \class   NullContact
    \brief   The default scenario where contact detection is inactive.
    \author  
    \warning 
  */
  /////////////////////////////////////////////////////////////////////////////

  class NullContact : public ContactModelBase {

  public:

    NullContact(const Uintah::ProcessorGroup* myworld,
                Uintah::SimulationStateP& ss, 
                PeridynamicsLabel* labels,
                PeridynamicsFlags* flags);
      
    virtual ~NullContact();

    void outputProblemSpec(Uintah::ProblemSpecP& ps);

    void exchangeMomentumInterpolated(const Uintah::ProcessorGroup*,
                                      const Uintah::PatchSubset* patches,
                                      const Uintah::MaterialSubset* matls,
                                      Uintah::DataWarehouse* old_dw,
                                      Uintah::DataWarehouse* new_dw);

    void exchangeMomentumIntegrated(const Uintah::ProcessorGroup*,
                                    const Uintah::PatchSubset* patches,
                                    const Uintah::MaterialSubset* matls,
                                    Uintah::DataWarehouse* old_dw,
                                    Uintah::DataWarehouse* new_dw);
      
    void addComputesAndRequiresInterpolated(Uintah::SchedulerP & sched,
                                            const Uintah::PatchSet* patches,
                                            const Uintah::MaterialSet* matls);

    void addComputesAndRequiresIntegrated(Uintah::SchedulerP & sched,
                                          const Uintah::PatchSet* patches,
                                          const Uintah::MaterialSet* matls);

  private:
      
    Uintah::SimulationStateP d_sharedState;

    // Prevent copying 
    NullContact(const NullContact &con);
    NullContact& operator=(const NullContact &con);
      
  };
} // End namespace Vaango
    


#endif /* __VAANGO_NULL_CONTACT_H__ */

