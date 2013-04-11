#ifndef EMU2DC_BODY_SP_H
#define EMU2DC_BODY_SP_H

#include <memory>

namespace Emu2DC {
  
  // Forward declaration.  Make sure <Body.h> is included before using BodyP.
  // using stdlib shared_ptr 
  class Body;
  typedef std::shared_ptr<Body> BodySP;
}

#endif
