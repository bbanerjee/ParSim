#ifndef MATITI_MATERIAL_SP_ARRAY_H
#define MATITI_MATERIAL_SP_ARRAY_H

#include <Pointers/MaterialSP.h>
#include <vector>

namespace Matiti {
  
  typedef std::vector<MaterialSP> MaterialSPArray;
  typedef std::vector<MaterialSP>::iterator MaterialSPIterator;
  typedef std::vector<MaterialSP>::const_iterator constMaterialSPIterator;
}

#endif
