#ifndef EMU2DC_CELLNODEPMAP_H
#define EMU2DC_CELLNODEPMAP_H

#include <NodeP.h>
#include <Hash64.h>
#include <tr1/unordered_map>
#include <cstdint>

namespace Emu2DC {
  
  typedef std::tr1::unordered_multimap<long64, NodeP, Hash64> CellNodePMap;
  typedef CellNodePMap::iterator CellNodePMapIterator;
  typedef CellNodePMap::const_iterator constCellNodePMapIterator;
  typedef std::pair<long64, NodeP> CellNodePPair; 
  typedef std::pair<CellNodePMapIterator, CellNodePMapIterator> CellNodePPairIterator; 
  typedef std::pair<constCellNodePMapIterator, constCellNodePMapIterator> constCellNodePPairIterator; 

  // alternative function object class for Hashing with cell id
  /*
  typedef std::tr1::unordered_multimap<long64, NodeP, HashCell> CellNodePMap;
  struct HashCell {
    std::size_t operator() (const long64& cellID) const {
      return std::hash<int>()(cellID >> 16) ^ std::hash<int>()(cellID >> 32) ^ std::hash<int>()(cellID >> 48);
    }
  };
  */

}

#endif
