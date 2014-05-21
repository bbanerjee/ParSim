#ifndef _VAANGO_PERIDYNAMICS_FAILUREMODELFACTORY_H_
#define _VAANGO_PERIDYNAMICS_FAILUREMODELFACTORY_H_

#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

  class PeridynamicsFailureModel;
  class PeridynamicsFlags;

  class PeridynamicsFailureModelFactory
  {
  public:
    static PeridynamicsFailureModel* create(Uintah::ProblemSpecP& ps, PeridynamicsFlags* flags);

  };
} // End namespace Vaango
      
#endif /* _VAANGO_PERIDYNAMICS_FAILUREMODELFACTORY_H_ */
