#ifndef EMU2DC_CRACK_SP_ARRAY_H
#define EMU2DC_CRACK_SP_ARRAY_H

#include <CrackSP.h>
#include <vector>

namespace Emu2DC {
  
  typedef std::vector<CrackSP> CrackSPArray;
  typedef std::vector<CrackSP>::iterator CrackSPIterator;
  typedef std::vector<CrackSP>::const_iterator constCrackSPIterator;
}

#endif
