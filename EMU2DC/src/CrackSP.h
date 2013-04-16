#ifndef EMU2DC_CRACK_SP_H
#define EMU2DC_CRACK_SP_H

#include <memory>

namespace Emu2DC {
  
  // Forward declaration.  Make sure <Crack.h> is included before using CrackP.
  // using stdlib shared_ptr 
  class Crack;
  typedef std::shared_ptr<Crack> CrackSP;
}

#endif
