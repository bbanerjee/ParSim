#ifndef _VAANGO_CONTACT_MODEL_FACTORY_H_
#define _VAANGO_CONTACT_MODEL_FACTORY_H_

#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Grid/SimulationStateP.h>

namespace Vaango {

  class ContactModelBase;
  class PeridynamicsLabel;
  class PeridynamicsFlags;

  class ContactModelFactory
  {
  public:
        
    static ContactModelBase* create(const Uintah::ProcessorGroup* myworld,
                                    const Uintah::ProblemSpecP& ps,
                                    Uintah::SimulationStateP& ss,
                                    PeridynamicsLabel* labels, 
                                    PeridynamicsFlags* flags);
  };
} // End namespace Vaango
  
#endif /* _VAANGO_CONTACT_MODEL_FACTORY_H_ */

