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
 * Heap.h
 *
 *  Created on: 4/12/2013
 *      Author: banerjee
 *      Original version: scikit-fmm/skfmm/heap.h
 *
 * This class implements a binary min heap data structure to support the
 * fast marching algorithm. A min heap is a list which has the property
 * that the smallest element is always the first element.
 *
 * The fast marching method uses this data structure to track elements in
 * the solution narrow band. The fast marching method needs to know which
 * element in the narrow band is nearest the zero level-set at each
 * iteration.
 *
 * When a new point enters the solution narrow band the element is added
 * to the heap.
 *
 * New elements are added to the heap with the push() method. The address
 * passed to pop is the (flat) address of the element in the grid and the
 * value passed is the distance. The push() method returns an integer (a
 * heap index) which is used later to refer to the element.
 *
 * As the solution evolves the distance of points already in the heap can
 * be updated via the set() method. The set() method takes the heap index
 * returned by the push() along with a new distance value.
 *
 * The narrow band element nearest the zero-level set is taken off the
 * heap with the pop() method. The grid address of the top element is
 * returned to the caller.
 *
 * The constructor for heap needs to know the number of elements that
 * will enter the narrow band. See Heap.cc for implementation details.
 */

#ifndef HEAP_H_
#define HEAP_H_

#include <vector>

namespace BrMPM
{

  class Heap
  {
  public:
    Heap(int depth, bool selfTest=false);
    virtual ~Heap();

    int push(int address, double value);
    void pop(int& address, double& value);
    void set(int index, double value);
    bool empty() const;

  private:
    void test() const;
    inline void siftUp(int pos);
    inline void siftDown(int startPos, int pos);

  private:

    int d_maxLength;
    int d_heapLength;
    int d_listLength;

    std::vector<double> d_distance;
    std::vector<int> d_heap;
    std::vector<int> d_address;
    std::vector<int> d_backPointer;

    bool d_selfTest;

  };

} /* namespace BrMPM */
#endif /* HEAP_H_ */
