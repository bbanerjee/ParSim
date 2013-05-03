#ifndef EMU2DC_ELEMENTPARRAY_H
#define EMU2DC_ELEMENTPARRAY_H

#include <ElementP.h>
#include <vector>

namespace Emu2DC {
  
  typedef std::vector<ElementP> ElementPArray;
  typedef std::vector<ElementP>::iterator ElementPIterator;
  typedef std::vector<ElementP>::const_iterator constElementPIterator;
}

#endif
