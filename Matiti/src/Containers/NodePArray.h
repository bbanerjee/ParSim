#ifndef MATITI_NODEPARRAY_H
#define MATITI_NODEPARRAY_H

#include <Pointers/NodeP.h>
#include <vector>

namespace Matiti {
  
  typedef std::vector<NodeP> NodePArray;
  typedef std::vector<NodeP>::iterator NodePIterator;
  typedef std::vector<NodeP>::const_iterator constNodePIterator;
}

#endif
