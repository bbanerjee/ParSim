#ifndef MATITI_ELEMENTPARRAY_H
#define MATITI_ELEMENTPARRAY_H

#include <Pointers/ElementP.h>
#include <vector>

namespace Matiti {
  
  typedef std::vector<ElementP> ElementPArray;
  typedef std::vector<ElementP>::iterator ElementPIterator;
  typedef std::vector<ElementP>::const_iterator constElementPIterator;
}

#endif
