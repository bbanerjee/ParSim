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
 * DistanceMarcher.h
 *
 *  Created on: 5/12/2013
 *      Author: banerjee
 *      Original version: scikit-fmm/skfmm/distance_marcher.h
 */

#ifndef DISTANCEMARCHER_H_
#define DISTANCEMARCHER_H_

#include <Contact/FastMarching/BaseMarcher.h>

namespace BrMPM
{

  class DistanceMarcher: public BaseMarcher
  {
  public:
    DistanceMarcher(Double3D* phi, const Vector3D& dx, Int3D* flag,
                    Double3D* distance, int ndim, const Double3DSizeType* shape,
                    bool self_test, int order);
    virtual ~DistanceMarcher();

  protected:

    virtual double solveQuadratic(int ii, const double& aa, const double& bb, double& cc);

    virtual void initializeFrozen();
    virtual double updatePointSecondOrder(int ii);
    virtual double updatePointFirstOrder(int ii);
  };

} /* namespace BrMPM */
#endif /* DISTANCEMARCHER_H_ */
