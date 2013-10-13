#ifndef MATITI_BONDPARRAY_H
#define MATITI_BONDPARRAY_H

#include <Pointers/BondP.h>
#include <vector>

namespace Matiti {
  
  typedef std::vector<BondP> BondPArray;
  typedef std::vector<BondP>::iterator BondPIterator;
  typedef std::vector<BondP>::const_iterator constBondPIterator;
}

#endif
