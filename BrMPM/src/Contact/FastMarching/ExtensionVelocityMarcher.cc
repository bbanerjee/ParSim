/*
 * ExtensionVelocityMarcher.cc
 *
 *  Created on: 5/12/2013
 *      Author: banerjee
 *      Original: scikit-fmm/skfmm/extension_velocity_marcher.cpp
 */

#include <Contact/FastMarching/ExtensionVelocityMarcher.h>
#include <iostream>
#include <cmath>
#include <stdexcept>

using namespace BrMPM;

ExtensionVelocityMarcher::ExtensionVelocityMarcher(double* phi, double* dx,
    long * flag, double* distance, int ndim, const std::vector<int>& shape,
    bool self_test, int order, double* speed, double* f_ext, long* ext_mask)
  : DistanceMarcher(phi, dx, flag, distance, ndim, shape, self_test, order),
    d_speed(speed), d_f_ext(f_ext), d_ext_mask(ext_mask)
{
}

ExtensionVelocityMarcher::~ExtensionVelocityMarcher()
{
}

void
ExtensionVelocityMarcher::initializeFrozen()
{
  // loop over phi to find zero values and mark them as frozen
  for (int ii = 0; ii < d_size; ii++) {
    if (d_flag[ii] != BrMPM::Mask && d_phi[ii] == 0.0) {
      d_flag[ii] = BrMPM::Frozen;
      d_distance[ii] = 0.0;
      d_f_ext[ii] = d_speed[ii];
    }
  }

  // loop over all of phi and for each point check each direction
  //  to see if we cross the zero level set
  //     if so calculate the minimum distance to the zero level set
  //     mark as frozen.
  for (int ii = 0; ii < d_size; ii++) {
    if (d_flag[ii] == BrMPM::Far) {
      std::vector<double> ldistance(MaximumDimension, 0.0);
      std::vector<double> lspeed(MaximumDimension, 0.0);
      bool borders=false;
      for (int dim = 0; dim < d_dim; dim++) {
        ldistance[dim] = 0.0;
        lspeed[dim] = 0.0;
        for (int jj = -1; jj < 2; jj += 2) { // each direction
          int naddr = getN(ii, dim, jj, BrMPM::Mask);
          if (naddr !=-1 && d_phi[ii] * d_phi[naddr] < 0) {
            // this cell and neighbor span the zero level set.
            borders=true;
            //calculate the distance to the zero level set.
            double dd = d_dx[dim]*d_phi[ii]/(d_phi[ii]-d_phi[naddr]);
            if (ldistance[dim] == 0.0 || ldistance[dim] > dd) {
              ldistance[dim] = dd;
              if (d_ext_mask[ii]) {
                lspeed[dim] = d_speed[naddr];
              } else if (d_ext_mask[naddr]) {
                lspeed[dim] = d_speed[ii];
              } else {
                lspeed[dim] = d_speed[ii] + dd/d_dx[dim]*(d_speed[naddr]-d_speed[ii]);
              }
            }
          }
        } // for each direction
      } // for each dimension
      if (borders) {
        double numerator = 0.0;
        double denominator = 0.0;
        for (int dim = 0; dim < d_dim; dim++) {
          if (ldistance[dim] != 0.0) {
            numerator += lspeed[dim]/(ldistance[dim]*ldistance[dim]);
            denominator += 1.0/(ldistance[dim]*ldistance[dim]);
          }
        }
        if (denominator != 0.0) {
          d_f_ext[ii] = numerator/denominator;
        } else {
          throw std::runtime_error("MPMContact::ExtensionVelocityMarcher:initialize Division by zero.");
        }

        double dsum = 0.0;
        for (int dim = 0; dim < d_dim; dim++) {
          if (ldistance[dim] > 0.0) {
            dsum += 1.0/(ldistance[dim]*ldistance[dim]);
          }
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

// Set the extension velocity of this point
// find f_ext where grad f_ext . grad phi = 0
// as described in Adalsteinsson and Sethian

// Technically we do not need to calculate this extension velocity
// until the point is frozen.
void
ExtensionVelocityMarcher::finalizePoint(int ii, double phi_i)
{
  std::vector<double> ldistance(MaximumDimension, 0.0);
  std::vector<double> lspeed(MaximumDimension, 0.0);

  double numerator = 0.0;
  double denominator = 0.0;

  for (int dim = 0; dim < d_dim; dim++) {
    ldistance[dim] = 0.0;
    lspeed[dim] = 0.0;
    for (int jj = -1; jj < 2; jj += 2) { // each direction
      int naddr = getN(ii, dim, jj, BrMPM::Mask);
      if (naddr !=-1 && d_flag[naddr] == BrMPM::Frozen) {
        // Determine which direction, in this dimension, is nearest to
        // the front. Calculate the distance to front in this direction
        double dd = d_distance[ii] - d_distance[naddr];
        if (ldistance[dim] == 0.0 || ldistance[dim] > dd) {
          ldistance[dim] = dd;
          lspeed[dim] = d_f_ext[naddr];
        }
      }
    } // for each direction
  } // for each dimension

  for (int dim = 0; dim < d_dim; dim++) {
    numerator += std::abs(ldistance[dim])*lspeed[dim]*d_idx2[dim];
    denominator += std::abs(ldistance[dim])*d_idx2[dim];
  }
  if (denominator != 0.0) {
    d_f_ext[ii] = numerator/denominator;
  } else {
    throw std::runtime_error("MPMContact::ExtensionVelocityMarcher:finalize Division by zero.");
  }
}

void ExtensionVelocityMarcher::cleanUp()
{
  for (int ii =0; ii < d_size; ii++) {
    if (d_flag[ii] == BrMPM::Mask) d_f_ext[ii] = BrMPM::MaxDouble;
    if (d_flag[ii] == BrMPM::Far) d_f_ext[ii] = BrMPM::MaxDouble;
  }
}

