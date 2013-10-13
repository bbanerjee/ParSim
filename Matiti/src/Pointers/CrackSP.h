#ifndef MATITI_CRACK_SP_H
#define MATITI_CRACK_SP_H

#include <memory>

namespace Matiti {
  
  // Forward declaration.  Make sure <Crack.h> is included before using CrackP.
  // using stdlib shared_ptr 
  class Crack;
  typedef std::shared_ptr<Crack> CrackSP;
}

#endif
