#ifndef EMU2DC_ELEMENTP_H
#define EMU2DC_ELEMENTP_H

#include <memory>

namespace Emu2DC {
  
  // Forward declaration.  Make sure <Node.h> is included before using NodeP.
  // using stdlib shared_ptr instead of SCIRun::Handle
  class Element;
  typedef std::shared_ptr<Element> ElementP;
}

#endif
