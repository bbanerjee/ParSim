#ifndef __MATITI_LOAD_BC_FACTORY_H__
#define __MATITI_LOAD_BC_FACTORY_H__

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Pointers/LoadBCSP.h>

namespace Matiti {

  class LoadBCFactory
  {
  public:

    static LoadBCSP create(Uintah::ProblemSpecP& ps);

  }; // end class

}  // end namespace

#endif
