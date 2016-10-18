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

/*
 * BaseMarcher.cc
 *
 *  Created on: 4/12/2013
 *      Author: banerjee
 */

#include "Contact/FastMarching/BaseMarcher.h"

using namespace BrMPM;

BaseMarcher::BaseMarcher(Double3D* phi, const Vector3D& dx, Int3D* flag,
                         Double3D* distance, int ndim, const Double3DSizeType* shape,
                         bool self_test, int order)
  : d_distance(distance), d_phi(phi), d_dx(dx), d_flag(flag), d_error(1), d_dim(ndim),
    d_size(1), d_idx2(ndim, 0.0), d_order(order), d_heap(0), d_self_test(self_test)
{
  for (int ii = 0; ii < ndim; ii++) {
    d_shape[ii] = shape[ii];
    d_size *= shape[ii];
    d_idx2[ii] = 1.0/(dx[ii]*dx[ii]);
  }
  for (int ii = 0; ii < ndim; ii++) {
    int prod = 1;
    for (int jj = ii+1; jj < ndim; jj++) prod *= d_shape[jj];
    d_shift[ii] = prod;
  }
}

BaseMarcher::~BaseMarcher()
{
  delete d_heap;
}

void
BaseMarcher::march()
{
  initializeFrozen();

  int maxHeap = 0;
  for (int ii = 0; ii < d_size; ii++) {
    if (d_flag[ii] == BrMPM::Far) maxHeap++;
  }

  d_heap = new Heap(maxHeap, d_self_test);
  d_heap_ptr.resize(d_size);

  initializeNarrow();
  solve();
  cleanUp();
}

// assume c order.
// for a given point find its neighbor in the given dimension
// and direction. Return -1 if not possible.
// consult d_shape information
int
BaseMarcher::getN(int current, int dim, int dir, int flag)
{
  std::vector<int> coord(MaximumDimension, 0);
  getIndex(current, coord);
  unsigned int newc = coord[dim]+dir;
  if (newc >= d_shape[dim] || newc < 0) return -1;
  int newa = current + dir*d_shift[dim];
  if (d_flag[newa]==flag)  return -1;
  getIndex(newa, coord);
  return newa;
}

// For each point in the far field check if neighbor is frozen
// if so calculate distance and insert into heap
void
BaseMarcher::initializeNarrow()
{
  for (int ii = 0; ii < d_size; ii++) {
    if (d_flag[ii] == BrMPM::Far) {
      for (int jj = -1; jj < 2; jj += 2) {
        int naddr = getN(ii, d_dim, jj, BrMPM::Mask);
        if (naddr != 1 && d_flag[naddr] == BrMPM::Frozen) {
          if (d_flag[ii] == BrMPM::Far) {
            d_flag[ii] = BrMPM::Narrow;
            double dd;
            if (d_order == 2) {
              dd = updatePointSecondOrder(ii);
            } else {
              dd = updatePointFirstOrder(ii);
            }
            d_distance[ii] = dd;
            d_heap_ptr[ii] = d_heap->push(ii, std::abs(dd));
          } // end if far and frozen
        } // end if frozen
      } // each direction
    } // each dimension
  } // each far field point
}

// Fast marching algorithm main loop
// (1) take the smallest narrow band element and
//     freeze it.
// (2) for each neighbor of the frozen point calculate distance based
// on frozen elements
//     - mark each neighbor as narrow band and stick it
//       into the heap
//     - if the neighbor is already in the heap update the
//       distance value.

void BaseMarcher::solve()
{
  // Update count of frozen points
  int frozenCount = 0;
  for (int ii = 0; ii < d_size; ii++) {
    if (d_flag[ii] == BrMPM::Frozen) frozenCount++;
  }
  if (!frozenCount) {
    d_error = 2;
    return;
  }

  int ii = 0;
  while (!d_heap->empty()) {
    ii++;
    double value = 0.0;
    int addr = 0;

    d_heap->pop(addr, value);
    d_flag[addr] = BrMPM::Frozen;
    finalizePoint(addr, value);

    for (int dim = 0; dim < d_dim; dim++) {
      for (int jj = -1; jj < 2; jj += 2) {
        int naddr = getN(addr, d_dim, jj, BrMPM::Frozen);
        if (naddr != 1 && d_flag[naddr] == BrMPM::Frozen) {
          if (d_flag[naddr] == BrMPM::Narrow) {
            double dd;
            if (d_order == 2) {
              dd = updatePointSecondOrder(naddr);
            } else {
              dd = updatePointFirstOrder(naddr);
            }
            if (dd) {
              d_heap->set(d_heap_ptr[naddr], std::abs(dd));
              d_distance[naddr] = dd;
            }
          } // end if narrow
          else if (d_flag[naddr] == BrMPM::Far) {
            double dd;
            if (d_order == 2) {
              dd = updatePointSecondOrder(naddr);
            } else {
              dd = updatePointFirstOrder(naddr);
            }
            if (dd) {
              d_distance[naddr] = dd;
              d_flag[naddr] = BrMPM::Narrow;
              d_heap_ptr[naddr] = d_heap->push(naddr, std::abs(dd));
            }
          }
        } // end if frozen
        //==========================================================
        // update the far point in the second order stencil
        // "jump" over a Frozen point if needed
        if (d_order == 2) {
          int local_naddr = getN(addr, dim, jj, BrMPM::Mask);
          if (local_naddr!=-1 && d_flag[local_naddr]==BrMPM::Frozen) {
            int naddr2 = getN(addr, dim, jj*2, BrMPM::Frozen);
            if (naddr2 !=-1 && d_flag[naddr2]== BrMPM::Narrow) {
              double dd = updatePointSecondOrder(naddr2);
              if (dd) {
                d_heap->set(d_heap_ptr[naddr2], std::abs(dd));
                d_distance[naddr2]=dd;
              }
            }
          }
        }
        //==========================================================
      } // end for each direction
    } // end for each dimension
  } // end of main loop

  // Add back mask (needs to be preprocessed by calling routine)
  for (int ii = 0; ii < d_size; ii++) {
    if (d_flag[ii] == BrMPM::Mask) d_distance[ii] = BrMPM::MaxDouble;
    if (d_flag[ii] == BrMPM::Far) d_distance[ii] = BrMPM::MaxDouble;
  }
  d_error = 0;
  return;
}

void
BaseMarcher::getIndex(int current, std::vector<int>& coord)
{
  int rem = current;
  for (int ii = 0; ii < d_dim; ii++) {
    coord[ii] = rem/d_shift[ii];
    rem -= coord[ii]*d_shift[ii];
  }
}

