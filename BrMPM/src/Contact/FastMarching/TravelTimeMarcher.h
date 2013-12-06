/*
 * TravelTimeMarcher.h
 *
 *  Created on: 5/12/2013
 *      Author: banerjee
 *      Original version: scikit-fmm/skfmm/travel_time_marcher.h
 */

#ifndef TRAVELTIMEMARCHER_H_
#define TRAVELTIMEMARCHER_H_

#include <Contact/FastMarching/DistanceMarcher.h>

namespace BrMPM
{

  class TravelTimeMarcher: public DistanceMarcher
  {
  public:
    TravelTimeMarcher(double* phi, double* dx, long* flag,
                      double* distance, int ndim, const std::vector<int>& shape,
                      bool self_test, int order, double* speed);
    virtual ~TravelTimeMarcher();

  protected:

    virtual double solveQuadratic(int ii, const double& aa, const double& bb, double& cc);

    virtual void initializeFrozen();
    virtual double updatePointSecondOrder(int ii);

  private:

    double* d_speed;
  };

} /* namespace BrMPM */
#endif /* TRAVELTIMEMARCHER_H_ */
