#ifndef MATITI_MATERIAL_SP_ARRAY_H
#define MATITI_MATERIAL_SP_ARRAY_H

#include <MaterialSP.h>
#include <vector>

namespace BrMPM {
  
  typedef std::vector<MaterialSP> MaterialSPArray;
  typedef std::vector<MaterialSP>::iterator MaterialSPIterator;
  typedef std::vector<MaterialSP>::const_iterator constMaterialSPIterator;
}

#endif
