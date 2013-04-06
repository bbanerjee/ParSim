#ifndef EMU2DC_NODEPARRAY_H
#define EMU2DC_NODEPARRAY_H

#include <NodeP.h>
#include <vector>

namespace Emu2DC {
  
  typedef std::vector<NodeP> NodePArray;
  typedef std::vector<NodeP>::iterator NodePIterator;
  typedef std::vector<NodeP>::const_iterator constNodePIterator;
}

#endif
