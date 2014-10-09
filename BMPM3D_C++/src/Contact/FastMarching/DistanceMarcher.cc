/*
 * DistanceMarcher.cc
 *
 *  Created on: 5/12/2013
 *      Author: banerjee
 *      Original version: scikit-fmm/skfmm/distance_marcher.cpp
 */

#include <Contact/FastMarching/DistanceMarcher.h>
#include <cmath>

using namespace BrMPM;

DistanceMarcher::DistanceMarcher(Double3D* phi, const Vector3D& dx, Int3D* flag,
                                 Double3D* distance, int ndim, const Double3DSizeType* shape,
                                 bool self_test, int order)
  : BaseMarcher(phi, dx, flag, distance, ndim, shape, self_test, order)
{
}

DistanceMarcher::~DistanceMarcher()
{
}

double
DistanceMarcher::solveQuadratic(int ii, const double& aa,
                                const double& bb, double& cc)
{
  cc -= 1;
  double r0 = 0.0;
  double r1 = 0.0;
  double det = bb*bb - 4.0*aa*cc;
  if (det > 0.0) {
    r0 = (-bb + std::sqrt(det))/(2.0*aa);
    r1 = (-bb - std::sqrt(det))/(2.0*aa);
  } else {
    return 0.0;
  }
  if (d_phi[ii] > BrMPM::DoubleEpsilon) return r0;
  return r1;
}

void
DistanceMarcher::initializeFrozen()
{
  // loop over phi to find zero values and mark them as frozen
  for (int ii = 0; ii < d_size; ii++) {
    if (d_flag[ii] != BrMPM::Mask && d_phi[ii] == 0.0) {
      d_flag[ii] = BrMPM::Frozen;
      d_distance[ii] = 0.0;
    }
  }

  // loop over all of phi and for each point check each direction
  //  to see if we cross the zero level set
  //     if so calculate the minimum distance to the zero level set
  //     mark as frozen.
  for (int ii = 0; ii < d_size; ii++) {
    if (d_flag[ii] == BrMPM::Far) {
      std::vector<double> ldistance(MaximumDimension, 0.0);
      bool borders=false;
      for (int dim = 0; dim < d_dim; dim++) {
        ldistance[dim] = 0.0;
        for (int jj = -1; jj < 2; jj += 2) { // each direction
          int naddr = getN(ii, dim, jj, BrMPM::Mask);
          if (naddr !=-1 && d_phi[ii] * d_phi[naddr] < 0) {
            // this cell and neighbor span the zero level set.
            borders=true;
            //calculate the distance to the zero level set.
            double dd = d_dx[dim]*d_phi[ii]/(d_phi[ii]-d_phi[naddr]);
            if (ldistance[dim] == 0.0 || ldistance[dim] > dd) {
              ldistance[dim] = dd;
            }
          }
        } // for each direction
      } // for each dimension
      if (borders) {
        double dsum = 0.0;
        for (int dim = 0; dim < d_dim; dim++) {
          if (ldistance[dim] > 0) dsum += 1.0/(ldistance[dim]*ldistance[dim]);
        }
        if (d_phi[ii] < 0.0) {
          d_distance[ii] = -sqrt(1.0/dsum);
        } else {
          d_distance[ii] = sqrt(1.0/dsum);
        }
        d_flag[ii] = BrMPM::Frozen;
      }
    }// for each point in the far field
  } // end loop through phi values
}

// second order point update
// update the distance from the frozen points
const double sqrThreeByTwo = 9.0/4.0;
const double oneThird = 1.0/3.0;
double
DistanceMarcher::updatePointSecondOrder(int ii)
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

double
DistanceMarcher::updatePointFirstOrder(int ii)
{
  double aa = 0.0, bb = 0.0, cc = 0.0;
  int naddr = 0;

  for (int dim = 0; dim < d_dim; dim++) {
    double value = BrMPM::MaxDouble;

    for (int jj = -1; jj < 2; jj += 2) {  // each direction
       naddr = getN(ii, dim, jj, BrMPM::Mask);
       if (naddr != -1 && d_flag[naddr] == BrMPM::Frozen) {
         if (d_distance[naddr] < value) {
           value = d_distance[naddr];
         }
       }
    } // end for each direction

    if (value < BrMPM::MaxDouble) {
      aa += d_idx2[dim];
      bb -= d_idx2[dim]*2.0*value;
      cc += d_idx2[dim]*value*value;
    }

  } // end for through phi values

  return solveQuadratic(ii, aa, bb, cc);
}

