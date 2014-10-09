/*
 * MPMContact.cc
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#include <Contact/MPMContact.h>
#include <Contact/FastMarching/FastMarchingMethod.h>
#include <MPMDatawarehouse.h>
#include <MPMPatch.h>
#include <MPMDataTypes.h>

using namespace BrMPM;

// Initialize default contact algorithm
MPMContact::MPMContact()
  :d_mtol(1.0e-10)
{
}

// Initialize default contact algorithm
MPMContact::MPMContact(std::vector<int>& dwis,
                       MPMPatchP& patch)
  : d_dwis(dwis), d_patch(patch), d_mtol(1.0e-10)
{
}

MPMContact::~MPMContact()
{
}

// Find whether two bodies intersect using a fast marching algorithm
void
MPMContact::findIntersection(MPMDatawarehouseP& dw)
{
  // Calculate a distance measure
  IntVector3D ppe = d_patch->ppe();
  double zeroLevel = 1.0 - std::sqrt(2.0)/(double) ppe.max();
  double tol0 = 0.5*(d_patch->dX()).max();
  double tol = 0.0;

  // Get the current distances from the DW
  DoubleNodeData gDistance_mat_1, gDistance_mat_2;
  dw->get("gDist", d_dwis[0], gDistance_mat_1);
  dw->get("gDist", d_dwis[1], gDistance_mat_2);

  // Update the distances
  for (auto iter = gDistance_mat_1.begin(); iter != gDistance_mat_1.end(); iter++) {
    int node = iter - gDistance_mat_1.begin();
    gDistance_mat_1[node] = zeroLevel - gDistance_mat_1[node];
  }

  for (auto iter = gDistance_mat_2.begin(); iter != gDistance_mat_2.end(); iter++) {
    int node = iter - gDistance_mat_2.begin();
    gDistance_mat_2[node] = zeroLevel - gDistance_mat_2[node];
  }

  // Reshape the distance vector into a three-dimensional array
  Double3DArray::extent_gen extents;
  int nx = d_patch->nC()[0];
  int ny = d_patch->nC()[1];
  int nz = d_patch->nC()[2];
  Double3DArray phi0(extents[nx][ny][nz]);
  phi0.assign(gDistance_mat_1.data(), gDistance_mat_1.data()+gDistance_mat_1.size());
  Double3DArray phi1(extents[nx][ny][nz]);
  phi1.assign(gDistance_mat_2.data(), gDistance_mat_1.data()+gDistance_mat_1.size());

  // Find distance using fast marching method
  FastMarchingMethod fmm;
  Double3DArray dist0(extents[nx][ny][nz]);
  std::fill(dist0.origin(), dist0.origin() + dist0.size(), 0.0);
  fmm.distance(phi0, d_patch->dX(), false, 2, dist0);

  Double3DArray dist1(extents[nx][ny][nz]);
  std::fill(dist1.origin(), dist1.origin() + dist1.size(), 0.0);
  fmm.distance(phi1, d_patch->dX(), false, 2, dist1);

  // Compute mask
  Double3DArray gmask(extents[nx][ny][nz]);
  for (DoubleIndex ii = 0; ii != nx; ++ii) {
    for (DoubleIndex jj = 0; jj != ny; ++jj)  {
      for (DoubleIndex kk = 0; kk != nz; ++kk)  {
        gmask[ii][jj][kk] = (double)(dist0[ii][jj][kk] < tol0)*(dist1[ii][jj][kk] < tol0);
      }
    }
  }
  for (DoubleIndex ii = 0; ii != nx; ++ii) {
    for (DoubleIndex jj = 0; jj != ny; ++jj)  {
      for (DoubleIndex kk = 0; kk != nz; ++kk)  {
        gmask[ii][jj][kk] *= (double)((dist0[ii][jj][kk] + dist1[ii][jj][kk]) < tol);
      }
    }
  }

  // The required nodes are where gmask is true
  int numNodes = nx*ny*nz;
  boost::array<DoubleIndex, 3> dims = {{numNodes,1,1}};
  gmask.reshape(dims);
  int nodeID = 0;
  for (DoubleIndex ii = 0; ii < numNodes; ++ii) {
    if (gmask[ii][0][0] == 1.0) {
      d_nodes.emplace_back(nodeID);
    }
    ++nodeID;
  }

}

// Find intersection nodes between two bodies that share a common grid
// using grid masses
void
MPMContact::findIntersectionSimple(MPMDatawarehouseP& dw)
{
  DoubleNodeData gm_mat_0, gm_mat_1;
  dw->get("gm", d_dwis[0], gm_mat_0);
  dw->get("gm", d_dwis[1], gm_mat_1);

  // Iterate through nodes on the common grid
  for (auto iter = gm_mat_0.begin(); iter != gm_mat_0.end(); iter++) {

    // Get node id
    int node = iter - gm_mat_0.begin();

    // Check intersection and add intersection nodes to list
    int intersect = (gm_mat_0[node] > d_mtol)*(gm_mat_1[node] > d_mtol);
    if (intersect) {
      d_nodes.emplace_back(node);
    }
  }
}

void
MPMContact::exchMomentumInterpolated(MPMDatawarehouseP& dw)
{
}

void
MPMContact::exchForceInterpolated(MPMDatawarehouseP& dw)
{
}

void
MPMContact::exchMomentumIntegrated(MPMDatawarehouseP& dw)
{
}
