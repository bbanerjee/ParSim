#ifndef MATITI_ELEMENTP_H
#define MATITI_ELEMENTP_H

#include <memory>

namespace Matiti {
  
  // Forward declaration.  Make sure <Node.h> is included before using NodeP.
  // using stdlib shared_ptr instead of SCIRun::Handle
  class Element;
  typedef std::shared_ptr<Element> ElementP;
}

#endif
