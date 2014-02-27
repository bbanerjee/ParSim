#ifndef MATITI_LOAD_BC_SP_H
#define MATITI_LOAD_BC_SP_H

#include <memory>

namespace Matiti {
  
  // Forward declaration.  Make sure <LoadBC.h> is included before using LoadBCP.
  // using stdlib shared_ptr 
  class LoadBC;
  typedef std::shared_ptr<LoadBC> LoadBCSP;
}

#endif
