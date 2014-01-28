#ifndef _VAANGO_PERIDYNAMICS_MATERIALMODELFACTORY_H_
#define _VAANGO_PERIDYNAMICS_MATERIALMODELFACTORY_H_

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <string>

namespace Vaango {

  class PeridynamicsMaterialModel;
  class PeridynamicsLabel;
  class PeridynamicsFlags;

  class PeridynamicsMaterialModelFactory
  {
  public:
    static PeridynamicsMaterialModel* create(Uintah::ProblemSpecP& ps, PeridynamicsFlags* flags);

  };
} // End namespace Vaango
      
#endif /* _VAANGO_PERIDYNAMICS_MATERIALMODELFACTORY_H_ */
