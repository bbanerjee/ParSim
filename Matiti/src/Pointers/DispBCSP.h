#ifndef MATITI_DISP_BC_SP_H
#define MATITI_DISP_BC_SP_H

#include <memory>

namespace Matiti {
  
  // Forward declaration.  Make sure <DisplacementBC.h> is included before using DispBCP.
  // using stdlib shared_ptr 
  class DisplacementBC;
  typedef std::shared_ptr<DisplacementBC> DispBCSP;
}

#endif
