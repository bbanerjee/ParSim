#ifndef __MATITI_DAMAGE_MODEL_FACTORY_H__
#define __MATITI_DAMAGE_MODEL_FACTORY_H__

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Pointers/DamageModelSP.h>

namespace Matiti {

  class DamageModelFactory
  {
  public:

    static DamageModelSP create(Uintah::ProblemSpecP& ps);

  }; // end class

}  // end namespace

#endif
