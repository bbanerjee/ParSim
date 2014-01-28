#ifndef _VAANGO_PERIDYNAMICS_FAILUREMODELFACTORY_H_
#define _VAANGO_PERIDYNAMICS_FAILUREMODELFACTORY_H_

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <string>

namespace Vaango {

  class PeridynamicsFailureModel;
  class PeridynamicsLabel;
  class PeridynamicsFlags;

  class PeridynamicsFailureModelFactory
  {
  public:
    static PeridynamicsFailureModel* create(Uintah::ProblemSpecP& ps, PeridynamicsFlags* flags);

  };
} // End namespace Vaango
      
#endif /* _VAANGO_PERIDYNAMICS_FAILUREMODELFACTORY_H_ */
