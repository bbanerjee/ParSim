/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef __CORE_GRID_VARIABLES_ListOfCellsIterator_H__
#define __CORE_GRID_VARIABLES_ListOfCellsIterator_H__

#include <Core/Grid/Variables/BaseIterator.h>
#include <Core/Grid/Variables/Iterator.h>

#include <Core/Geometry/IntVector.h>
#include <Core/Malloc/Allocator.h>

#include <climits>
#include <vector>

namespace Uintah {

class ListOfCellsIterator : public BaseIterator
{
  friend std::ostream&
  operator<<(std::ostream& out, const Uintah::ListOfCellsIterator& b);

public:

  ListOfCellsIterator(int size)
    : d_size{ 0 }
    , d_index{ 0 }
    , d_list_of_cells(size + 1)
  {
    d_list_of_cells[d_size] = IntVector(INT_MAX, INT_MAX, INT_MAX);
  }

  ListOfCellsIterator(const ListOfCellsIterator& copy)
    : d_list_of_cells(copy.d_list_of_cells)
  {
    reset();
  }

  ListOfCellsIterator&
  operator=(const ListOfCellsIterator& copy)
  {
    d_list_of_cells = copy.d_list_of_cells;
    reset();
    return *this;
  }

  /**
   * prefix operator to move the iterator forward
   */
  ListOfCellsIterator&
  operator++()
  {
    d_index++;
    return *this;
  }

  /**
   * postfix operator to move the iterator forward
   * does not return the iterator because of performance issues
   */
  void
  operator++(int)
  {
    d_index++;
  }

  /**
   * returns true if the iterator is done
   */
  bool
  done() const
  {
    return d_index == d_size;
  }

  /**
   * returns the IntVector that the current iterator is pointing at
   */
  IntVector
  operator*() const
  {
    ASSERT(d_index < d_list_of_cells.size());
    return d_list_of_cells[d_index];
  }

  /**
   * Assignment operator - this is expensive as we have to allocate new memory
   */
  inline ListOfCellsIterator&
  operator=(Iterator& copy)
  {
    // delete old iterator
    int i = 0;

    // copy iterator into portable container
    for (copy.reset(); !copy.done(); copy++) {
      d_list_of_cells[i] = (*copy);
      i++;
    }

    d_size = i;
    return *this;
  }

  /**
   * Return the first element of the iterator
   */
  inline IntVector
  begin() const
  {
    return d_list_of_cells[0];
  }

  /**
   * Return one past the last element of the iterator
   */
  inline IntVector
  end() const
  {
    return d_list_of_cells[d_size];
  }

  /**
   * Return the number of cells in the iterator
   */
  inline unsigned int
  size() const
  {
    return d_size;
  };

  /**
   * adds a cell to the list of cells
   */
  inline void
  add(const IntVector& c)
  {
    // place at back of list
    d_list_of_cells[d_size] = c;

    // read sentinel to list
    d_list_of_cells[d_size] = IntVector(INT_MAX, INT_MAX, INT_MAX);
  }

  /**
   * resets the iterator
   */
  inline void
  reset()
  {
    d_index = 0;
  }

  inline std::vector<IntVector>&
  get_ref_to_iterator()
  {
    return d_list_of_cells;
  }

protected:

  ListOfCellsIterator() {}

  /**
   * Returns a pointer to a deep copy of the virtual class
   * this should be used only by the Iterator class
   */
  ListOfCellsIterator*
  clone() const
  {
    return scinew ListOfCellsIterator(*this);
  };

  virtual std::ostream&
  put(std::ostream& out) const
  {
    out << *this;
    return out;
  }

  virtual std::ostream&
  limits(std::ostream& out) const
  {
    out << begin() << " " << end();
    return out;
  }

  size_t d_size{ 0 };

  // index into the iterator
  size_t d_index{ 0 };

  // vector to store cells
  std::vector<IntVector> d_list_of_cells;

}; // end class ListOfCellsIterator

} // End namespace Uintah

#endif //__CORE_GRID_VARIABLES_ListOfCellsIterator_H__
