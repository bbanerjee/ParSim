#ifndef __VAANGO_SINGLE_VELOCITY_CONTACT_H__
#define __VAANGO_SINGLE_VELOCITY_CONTACT_H__

#include <CCA/Components/Peridynamics/ContactModels/ContactModelBase.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/SimulationStateP.h>
#include <Core/Grid/Task.h>

namespace Vaango {

  /////////////////////////////////////////////////////////////////////////////
  /*!
    \class   SingleVelocityContact
    \brief   Recaptures single velocity field behavior from objects belonging to 
             multiple velocity fields.  Ensures that one can get the same answer 
             using prescribed contact as can be gotten using "automatic" 
             MPM-type contact on the background grid.
    \author  Steven G. Parker
    \warning 
  */
  /////////////////////////////////////////////////////////////////////////////

  class SingleVelocityContact : public ContactModelBase {

  public:

    SingleVelocityContact(const Uintah::ProcessorGroup* myworld,
                          Uintah::ProblemSpecP& ps,
                          Uintah::SimulationStateP& d_sS,
                          PeridynamicsLabel* labels,
                          PeridynamicsFlags* flags);
         
    virtual ~SingleVelocityContact();

    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps);

    virtual void exchangeMomentumInterpolated(const Uintah::ProcessorGroup*,
                                              const Uintah::PatchSubset* patches,
                                              const Uintah::MaterialSubset* matls,
                                              Uintah::DataWarehouse* old_dw,
                                              Uintah::DataWarehouse* new_dw);
         
    virtual void exchangeMomentumIntegrated(const Uintah::ProcessorGroup*,
                                            const Uintah::PatchSubset* patches,
                                            const Uintah::MaterialSubset* matls,
                                            Uintah::DataWarehouse* old_dw,
                                            Uintah::DataWarehouse* new_dw);

    virtual void addComputesAndRequiresInterpolated(Uintah::SchedulerP & sched,
                                                    const Uintah::PatchSet* patches,
                                                    const Uintah::MaterialSet* matls);

    virtual void addComputesAndRequiresIntegrated(Uintah::SchedulerP & sched,
                                                  const Uintah::PatchSet* patches,
                                                  const Uintah::MaterialSet* matls);
  protected:

    Uintah::SimulationStateP    d_sharedState;

  private:
         
    // Prevent copying of this class
    SingleVelocityContact(const SingleVelocityContact &con);
    SingleVelocityContact& operator=(const SingleVelocityContact &con);

  };
} // End namespace Vaango
      

#endif /* __VAANGO_SINGLE_VELOCITY_CONTACT_H__ */

