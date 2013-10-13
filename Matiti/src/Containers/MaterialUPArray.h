#ifndef MATITI_MATERIAL_UP_ARRAY_H
#define MATITI_MATERIAL_UP_ARRAY_H

#include <Pointers/MaterialUP.h>
#include <vector>

namespace Matiti {
  
  typedef std::vector<MaterialUP> MaterialUPArray;
  typedef std::vector<MaterialUP>::iterator MaterialUPIterator;
  typedef std::vector<MaterialUP>::const_iterator constMaterialUPIterator;
}

#endif
