/*
 * MPMContact.cc
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#include <Contact/MPMContact.h>
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
  IntVector3D ppe = d_patch->ppe();
  double zeroLevel = 1.0 - std::sqrt(2.0)/(double) ppe.max();
  double tol0 = 0.5*(d_patch->dX()).max();
  double tol = 0.0;

  DoubleNodeData gDistance_mat_1, gDistance_mat_2;
  dw->get("gDist", d_dwis[0], gDistance_mat_1);
  dw->get("gDist", d_dwis[1], gDistance_mat_2);

  for (auto iter = gDistance_mat_1.begin(); iter != gDistance_mat_1.end(); iter++) {
    int node = iter - gDistance_mat_1.begin();
    gDistance_mat_1[node] = zeroLevel - gDistance_mat_1[node];
  }

  for (auto iter = gDistance_mat_2.begin(); iter != gDistance_mat_2.end(); iter++) {
    int node = iter - gDistance_mat_2.begin();
    gDistance_mat_2[node] = zeroLevel - gDistance_mat_2[node];
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
