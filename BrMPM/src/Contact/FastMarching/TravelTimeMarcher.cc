/*
 * TravelTimeMarcher.cpp
 *
 *  Created on: 5/12/2013
 *      Author: banerjee
 */

#include "Contact/FastMarching/TravelTimeMarcher.h"
#include <cmath>

using namespace BrMPM;

TravelTimeMarcher::TravelTimeMarcher(double* phi, double* dx, long * flag,
    double* distance, int ndim, const std::vector<int>& shape, bool self_test,
    int order, double* speed)
  : DistanceMarcher(phi, dx, flag, distance, ndim, shape, self_test, order),
    d_speed(speed)
{
  for (int ii =0; ii < d_size; ii++) {
    if (d_speed[ii] < BrMPM::DoubleEpsilon) d_flag[ii] = BrMPM::Mask;
  }
}

TravelTimeMarcher::~TravelTimeMarcher()
{
}

double
TravelTimeMarcher::solveQuadratic(int ii, const double& aa,
                                  const double& bb, double& cc)
{
  cc -= 1.0/(d_speed[ii]*d_speed[ii]);
  double r0 = 0.0;
  double det = bb*bb - 4.0*aa*cc;
  if (det > 0.0) {
    r0 = (-bb + std::sqrt(det))/(2.0*aa);
  } else {
    return 0.0;
  }
  return r0;
}

void
TravelTimeMarcher::initializeFrozen()
{
  DistanceMarcher::initializeFrozen();
  for (int ii =0; ii < d_size; ii++) {
    if (d_flag[ii] == BrMPM::Frozen) {
      d_distance[ii] = std::abs(d_distance[ii]/d_speed[ii]); // convert distance to time
    }
  }
}

const double sqrThreeByTwo = 9.0/4.0;
const double oneThird = 1.0/3.0;
double
TravelTimeMarcher::updatePointSecondOrder(int ii)
{
  double aa = 0.0, bb = 0.0, cc = 0.0;
  int naddr = 0;

  for (int dim = 0; dim < d_dim; dim++) {
    double value1 = BrMPM::MaxDouble;
    double value2 = BrMPM::MaxDouble;

    for (int jj = -1; jj < 2; jj += 2) {  // each direction
       naddr = getN(ii, dim, jj, BrMPM::Mask);
       if (naddr != -1 && d_flag[naddr] == BrMPM::Frozen) {
         if (d_distance[naddr] < value1) {
           value1 = d_distance[naddr];
           int naddr2 = getN(ii, dim, jj*2, BrMPM::Mask);
           if (naddr2 !=-1 && d_flag[naddr2] == BrMPM::Frozen &&
                ((d_distance[naddr2] <= value1 && value1 >= 0.0) ||
                 (d_distance[naddr2] >= value1 && value1 <= 0.0))) {
             value2 = d_distance[naddr2];
             if (d_phi[naddr2]*d_phi[naddr] < 0.0 || d_phi[naddr2](d_phi[ii] < 0.0)) {
               value2 *= -1.0;
             }
           }
         }
       }
    } // end for each direction

    if (value2 < BrMPM::MaxDouble) {
      double tp = oneThird*(4.0*value1-value2);
      aa += d_idx2[dim]*sqrThreeByTwo;
      bb -= d_idx2[dim]*2.0*sqrThreeByTwo*tp;
      cc += d_idx2[dim]*sqrThreeByTwo*tp*tp;
    } else if (value1 < BrMPM::MaxDouble) {
      aa += d_idx2[dim];
      bb -= d_idx2[dim]*2.0*value1;
      cc += d_idx2[dim]*value1*value1;
    }

  } // end for through phi values

  return solveQuadratic(ii, aa, bb, cc);

}

