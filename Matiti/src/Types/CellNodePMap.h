/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef MATITI_CELLNODEPMAP_H
#define MATITI_CELLNODEPMAP_H

#include <Pointers/NodeP.h>
#include <Core/Hash64.h>
#include <tr1/unordered_map>
#include <cstdint>

namespace Matiti {
  
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
