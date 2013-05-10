#ifndef EMU2DC_BONDPARRAY_H
#define EMU2DC_BONDPARRAY_H

#include <BondP.h>
#include <vector>

namespace Emu2DC {
  
  typedef std::vector<BondP> BondPArray;
  typedef std::vector<BondP>::iterator BondPIterator;
  typedef std::vector<BondP>::const_iterator constBondPIterator;
}

#endif
