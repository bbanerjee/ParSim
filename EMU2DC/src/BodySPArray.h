#ifndef EMU2DC_BODY_SP_ARRAY_H
#define EMU2DC_BODY_SP_ARRAY_H

#include <BodySP.h>
#include <vector>

namespace Emu2DC {
  
  typedef std::vector<BodySP> BodySPArray;
  typedef std::vector<BodySP>::iterator BodySPIterator;
  typedef std::vector<BodySP>::const_iterator constBodySPIterator;
}

#endif
