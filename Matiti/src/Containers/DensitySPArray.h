#ifndef MATITI_DENSITY_SP_ARRAY_H
#define MATITI_DENSITY_SP_ARRAY_H

#include <Pointers/DensitySP.h>
#include <vector>

namespace Matiti {
  
  typedef std::vector<DensitySP> DensitySPArray;
  typedef std::vector<DensitySP>::iterator DensitySPIterator;
  typedef std::vector<DensitySP>::const_iterator constDensitySPIterator;
}

#endif
