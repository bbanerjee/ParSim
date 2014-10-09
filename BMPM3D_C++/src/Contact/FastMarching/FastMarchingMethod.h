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
 * FastMarchingMethod.h
 *
 *  Created on: 5/12/2013
 *      Author: banerjee
 *      Original: scikit-fmm/skfmm/__init__.py
 *                scikit-fmm/skfmm/pfmm.py
 *                scikit-fmm/skfmm/fmm.cpp
 *
 *  The fast marching method is used to model the evolution of boundaries
 *  and interfaces in a variety of application areas. More specifically,
 *  the fast marching method is a numerical technique for finding
 *  approximate solutions to boundary value problems of the Eikonal
 *  equation:
 *
 *  F(x) | grad T(x) | = 1.
 *
 * Typically, such a problem describes the evolution of a closed curve as
 * a function of time T with speed F(x)>0 in the normal direction at a
 * point x on the curve. The speed function is specified, and the time at
 * which the contour crosses a point x is obtained by solving the
 * equation.
 *
 * The FastMarching classes provide functions to calculate
 * the signed distance and travel time to an interface described by the
 * zero contour of the input array phi.
 */

#ifndef FASTMARCHINGMETHOD_H_
#define FASTMARCHINGMETHOD_H_

#include <Contact/FastMarching/DistanceMarcher.h>
#include <GeometryMath/Vector3D.h>

namespace BrMPM
{
  class FastMarchingMethod
  {
  public:

    FastMarchingMethod();
    virtual ~FastMarchingMethod();

    /*  Return the distance from the zero contour of the array phi.
     *
     *    Parameters
     *    ----------
     *    phi : array-like
     *          the zero contour of this array is the boundary location for
     *          the distance calculation. Phi can of 1,2,3 or higher
     *          dimension and can be a masked array.
     *
     *    dx  : float or an array-like of shape len(phi), optional
     *          the cell length in each dimension.
     *
     *    self_test : bool, optional
     *                if True consistency checks are made on the binary min
     *                heap during the calculation. This is used in testing and
     *                results in a slower calculation.
     *
     *    order : int, optional
     *            order of computational stencil to use in updating points during
     *            the fast marching method. Must be 1 or 2, the default is 2.
     *
     *    Returns
     *    -------
     *    d : an array the same shape as phi
     *        contains the distance from the zero contour (zero level set)
     *        of phi to each point in the array.
     *
     */
    void distance(Double3DArray& phi, const Vector3D& dx,
                  bool self_test, int order,
                  Double3DArray& distance);

  };

} /* namespace BrMPM */

#endif /* FASTMARCHINGMETHOD_H_ */
