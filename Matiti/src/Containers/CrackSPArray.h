#ifndef MATITI_CRACK_SP_ARRAY_H
#define MATITI_CRACK_SP_ARRAY_H

#include <Pointers/CrackSP.h>
#include <vector>

namespace Matiti {
  
  typedef std::vector<CrackSP> CrackSPArray;
  typedef std::vector<CrackSP>::iterator CrackSPIterator;
  typedef std::vector<CrackSP>::const_iterator constCrackSPIterator;
}

#endif
