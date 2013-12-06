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

#include <vector>
#include <limits>

namespace BrMPM
{
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
    BaseMarcher(double* phi, double* dx, long* flag,
                double* distance, int ndim, const std::vector<int>& shape,
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

    double* d_distance;
    double* d_phi;
    double* d_dx;
    double* d_flag;
    int d_error;
    int d_dim;
    int d_size;
    std::vector<double> d_idx2;

  private:

    int d_order;
    std::vector<int> d_heap_ptr;
    Heap* d_heap;
    std::vector<int> d_shape;
    std::vector<int> d_shift;
    bool d_self_test;

  };

} /* namespace BrMPM */

#endif /* BASEMARCHER_H_ */
