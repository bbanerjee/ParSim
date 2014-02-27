#ifndef MATITI_LOAD_BC_SP_ARRAY_H
#define MATITI_LOAD_BC_SP_ARRAY_H

#include <Pointers/LoadBCSP.h>
#include <vector>

namespace Matiti {
  
  typedef std::vector<LoadBCSP> LoadBCSPArray;
  typedef std::vector<LoadBCSP>::iterator LoadBCSPIterator;
  typedef std::vector<LoadBCSP>::const_iterator constLoadBCSPIterator;
}

#endif
