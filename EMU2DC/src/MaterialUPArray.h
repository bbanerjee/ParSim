#ifndef EMU2DC_MATERIAL_UP_ARRAY_H
#define EMU2DC_MATERIAL_UP_ARRAY_H

#include <MaterialUP.h>
#include <vector>

namespace Emu2DC {
  
  typedef std::vector<MaterialUP> MaterialUPArray;
  typedef std::vector<MaterialUP>::iterator MaterialUPIterator;
  typedef std::vector<MaterialUP>::const_iterator constMaterialUPIterator;
}

#endif
