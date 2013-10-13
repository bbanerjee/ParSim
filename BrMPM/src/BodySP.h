#ifndef MATITI_BODY_SP_H
#define MATITI_BODY_SP_H

#include <memory>

namespace BrMPM {
  
  // Forward declaration.  Make sure <Body.h> is included before using BodyP.
  // using stdlib shared_ptr 
  class Body;
  typedef std::shared_ptr<Body> BodySP;
}

#endif
