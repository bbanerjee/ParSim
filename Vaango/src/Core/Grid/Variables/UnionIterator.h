/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#ifndef __CORE_GRID_VARIABLES_UnionIterator_H__
#define __CORE_GRID_VARIABLES_UnionIterator_H__

#include <Core/Geometry/IntVector.h>

#include <Core/Grid/Variables/BaseIterator.h>
#include <Core/Grid/Variables/Iterator.h>
#include <Core/Grid/Variables/ListOfCellsIterator.h>

namespace Uintah {

using Uintah::IntVector;

/**************************************

  CLASS
  UnionIterator

  This iterator will iterator over the union between two iterators

  GENERAL INFORMATION

  UnionIterator.h

  Justin Luitjens
  Department of Computer Science
  University of Utah

  Center for the Simulation of Accidental Fires and Explosions (C-SAFE)


  KEYWORDS
  UnionIterator

  DESCRIPTION
  This iterator will iterator over the union between two iterators

  WARNING

 ****************************************/

class UnionIterator : public ListOfCellsIterator
{
  friend std::ostream&
  operator<<(std::ostream& out, const Uintah::UnionIterator& b);

public:
  UnionIterator(Iterator iter1, Iterator iter2);

  UnionIterator() {}

  std::ostream&
  put(std::ostream& out) const
  {
    out << *this;
    return out;
  }

}; // end class UnionIterator

} // End namespace Uintah

#endif //__CORE_GRID_VARIABLES_UnionIterator_H__
