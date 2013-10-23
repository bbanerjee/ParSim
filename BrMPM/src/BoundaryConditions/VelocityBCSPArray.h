#ifndef MATITI_VELOCITYBC_SP_ARRAY_H
#define MATITI_VELOCITYBC_SP_ARRAY_H

#include <VelocityBCSP.h>
#include <vector>

namespace BrMPM {
  
  typedef std::vector<VelocityBCSP> VelocityBCSPArray;
  typedef std::vector<VelocityBCSP>::iterator VelocityBCSPIterator;
  typedef std::vector<VelocityBCSP>::const_iterator constVelocityBCSPIterator;
}

#endif
