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
 * BaseMarcher.h
 *
 *  Created on: 4/12/2013
 *      Author: banerjee
 *      Original version: scikit-fmm/skfmm/base_marcher.h
 */

#ifndef BASEMARCHER_H_
#define BASEMARCHER_H_

#include <Contact/FastMarching/Heap.h>
#include <GeometryMath/Vector3D.h>

#include <boost/multi_array.hpp>
#include <vector>
#include <limits>

namespace BrMPM
{
  typedef boost::multi_array<double, 3> Double3DArray;
  typedef boost::multi_array<int, 3> Int3DArray;
  typedef std::vector<int> Int1DArray;
  typedef Double3DArray::index DoubleIndex;
  typedef Int3DArray::index IntIndex;
  typedef Double3DArray::element Double3D;
  typedef Int3DArray::element Int3D;
  typedef Double3DArray::size_type Double3DSizeType;

  const unsigned int MaximumDimension = 12;
  const char Far = 0;
  const char Narrow = 0;
  const char Frozen = 0;
  const char Mask = 0;

  const double DoubleEpsilon = std::numeric_limits<double>::epsilon();
  const double MaxDouble = std::numeric_limits<double>::max();

  class BaseMarcher
  {
  public:
    BaseMarcher(Double3D* phi, const Vector3D& dx, Int3D* flag,
                Double3D* distance, int ndim, const Double3DSizeType* shape,
                bool self_test, int order);

    virtual ~BaseMarcher();

    void march();
    int getError() const {return d_error;}

  protected:

    virtual void initializeFrozen() = 0;
    virtual double updatePointSecondOrder(int ii) = 0;
    virtual double updatePointFirstOrder(int ii) = 0;

    virtual void cleanUp() {}
    virtual void finalizePoint(int ii, double phi_i) {}

    int getN(int current, int dim, int dir, int flag);

  private:

    void initializeNarrow();
    void solve();
    void getIndex(int current, std::vector<int>& coord);

  protected:

    Double3D* d_distance;
    Double3D* d_phi;
    Vector3D d_dx;
    Int3D* d_flag;
    int d_error;
    int d_dim;
    int d_size;
    std::vector<double> d_idx2;

  private:

    int d_order;
    std::vector<int> d_heap_ptr;
    Heap* d_heap;
    Double3DSizeType* d_shape;
    std::vector<int> d_shift;
    bool d_self_test;

  };

} /* namespace BrMPM */

#endif /* BASEMARCHER_H_ */
