/*
 * FastMarchingMethod.h
 *
 *  Created on: 5/12/2013
 *      Author: banerjee
 *      Original: scikit-fmm/skfmm/__init__.py
 *                scikit-fmm/skfmm/pfmm.py
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

#include <vector>

namespace BrMPM
{

  class FastMarchingMethod
  {
  public:
    FastMarchingMethod();
    virtual ~FastMarchingMethod();

    /* Get and check input data
     * Returns: flag
     */
    void preprocess(std::vector<double>& phi, double* dx, std::vector<int>& ext_mask,
                    std::vector<int>& flag);

    void postprocess(double* result);

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
    void distance(double* phi, double* dx, bool self_test=false, int order=2);

    /*  Return the travel from the zero contour of the array phi given the
     *  scalar velocity field speed.
     *
     *  Parameters
     *  ----------
     *  phi : array-like
     *        the zero contour of this array is the boundary location for
     *        the travel time calculation. Phi can of 1,2,3 or higher
     *        dimension and can be a masked array.
     *
     *  speed : array-like, the same shape as phi
     *          contains the speed of interface propagation at each point
     *          in the domain.
     *
     *  dx  : float or an array-like of shape len(phi), optional
     *        the cell length in each dimension.
     *
     *  self_test : bool, optional
     *              if True consistency checks are made on the binary min
     *              heap during the calculation. This is used in testing and
     *              results in a slower calculation.
     *
     *  order : int, optional
     *          order of computational stencil to use in updating points during
     *          the fast marching method. Must be 1 or 2, the default is 2.
     *
     *  Returns
     *  -------
     *  t : an array the same shape as phi
     *      contains the travel time from the zero contour (zero level
     *      set) of phi to each point in the array given the scalar
     *      velocity field speed. If the input array speed has values less
     *      than or equal to zero the return value will be a masked array.
     */
    void travel_time(double* phi, double* speed, double* dx, bool self_test=false, int order=2);

    /*  Extend the velocities defined at the zero contour of phi to the
     *  rest of the domain. Extend the velocities such that
     *  grad f_ext dot grad d = 0 where where f_ext is the
     *  extension velocity and d is the signed distance function.
     *
     *  Parameters
     *  ----------
     *  phi : array-like
     *        the zero contour of this array is the boundary location for
     *        the travel time calculation. Phi can of 1,2,3 or higher
     *        dimension and can be a masked array.
     *
     *  speed : array-like, the same shape as phi
     *          contains the speed of interface propagation at each point
     *          in the domain.
     *
     *  dx  : float or an array-like of shape len(phi), optional
     *        the cell length in each dimension.
     *
     *  self_test : bool, optional
     *              if True consistency checks are made on the binary min
     *              heap during the calculation. This is used in testing and
     *              results in a slower calculation.
     *
     *  order : int, optional
     *          order of computational stencil to use in updating points during
     *          the fast marching method. Must be 1 or 2, the default is 2.
     *
     *  ext_mask : array-like, the same shape as phi, optional
     *             enables initial front values to be eliminated when
     *             calculating the value at the interface before the
     *             values are extended away from the interface.
     *
     *  Returns
     *  -------
     *  (d, f_ext) : tuple
     *      a tuple containing the signed distance function d and the
     *      extension velocities f_ext.
     *
    */
    void extension_velocities(std::vector<double>& phi, std::vector<double>& speed,
                              double* dx,
                              bool self_test=false, int order=2,
                              std::vector<int>& ext_mask,
                              std::vector<double>& dd, std::vector<double>& f_ext);

  private:

    void distance_method();
  };

} /* namespace BrMPM */

#endif /* FASTMARCHINGMETHOD_H_ */
