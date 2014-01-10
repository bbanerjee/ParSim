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
    BaseMarcher(Double3D* phi, Vector3D& dx, Int3D* flag,
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
    Double3D* d_flag;
    Double3DSizeType* d_shape;
    int d_error;
    int d_dim;
    int d_size;
    std::vector<double> d_idx2;

  private:

    int d_order;
    std::vector<int> d_heap_ptr;
    Heap* d_heap;
    std::vector<int> d_shift;
    bool d_self_test;

  };

} /* namespace BrMPM */

#endif /* BASEMARCHER_H_ */
