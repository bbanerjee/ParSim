#ifndef MATITI_WOODS_SP_ARRAY_H
#define MATITI_WOODS_SP_ARRAY_H

#include <Pointers/WoodSP.h>
#include <vector>

namespace Matiti {
  
  typedef std::vector<WoodSP> WoodSPArray;
  typedef std::vector<WoodSP>::iterator WoodSPIterator;
  typedef std::vector<WoodSP>::const_iterator constWoodSPIterator;
}

#endif
