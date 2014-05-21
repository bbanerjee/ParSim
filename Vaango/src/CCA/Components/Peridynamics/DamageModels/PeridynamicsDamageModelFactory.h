#ifndef _VAANGO_PERIDYNAMICS_DAMAGE_MODEL_FACTORY_H_
#define _VAANGO_PERIDYNAMICS_DAMAGE_MODEL_FACTORY_H_

#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

  class PeridynamicsDamageModel;
  class PeridynamicsLabel;
  class PeridynamicsFlags;

  class PeridynamicsDamageModelFactory
  {
  public:
    static PeridynamicsDamageModel* create(Uintah::ProblemSpecP& ps, 
                                           PeridynamicsLabel* labels,
                                           PeridynamicsFlags* flags);

  };
} // End namespace Vaango
      
#endif /* _VAANGO_PERIDYNAMICS_DAMAGE_MODEL_FACTORY_H_ */
