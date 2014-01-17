/*
 * MPMFreeContact.cc
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#include <Contact/MPMFreeContact.h>
#include <MPMDatawarehouse.h>
#include <MPMDataTypes.h>

using namespace BrMPM;

MPMFreeContact::MPMFreeContact(std::vector<int>& dwis, MPMPatchP& patch)
   : MPMContact(dwis, patch)
{
}

MPMFreeContact::~MPMFreeContact() {
}

void
MPMFreeContact::exchMomentumInterpolated(MPMDatawarehouseP& dw) {

  findIntersection(dw);
  if (d_nodes.size() > 0) {
    DoubleNodeData gm_mat_0, gm_mat_1;
    DoubleNodeData gw_mat_0, gw_mat_1;
    dw->get("gm", d_dwis[0], gm_mat_0);
    dw->get("gw", d_dwis[0], gw_mat_0);
    dw->get("gm", d_dwis[1], gm_mat_1);
    dw->get("gw", d_dwis[1], gw_mat_1);

    for (unsigned int ii = 0; ii < d_nodes.size(); ii++) {
      gm_mat_0[ii] += gm_mat_1[ii];
      gw_mat_0[ii] += gw_mat_1[ii];
    }
    for (unsigned int ii = 0; ii < d_nodes.size(); ii++) {
      gm_mat_1[ii] = gm_mat_0[ii];
      gw_mat_1[ii] = gw_mat_0[ii];
    }
  }

}

void
MPMFreeContact::exchForceInterpolated(MPMDatawarehouseP& dw) {

  if (d_nodes.size() > 0) {
    exchVals("gfi", dw);
    Vector3DNodeData gfi_mat_0, gfi_mat_1;

    dw->get("gfi", d_dwis[0], gfi_mat_0);
    dw->get("gfi", d_dwis[1], gfi_mat_1);
    for (unsigned int ii = 0; ii < d_nodes.size(); ii++) {
      gfi_mat_0[ii] += gfi_mat_1[ii];
    }
    for (unsigned int ii = 0; ii < d_nodes.size(); ii++) {
      gfi_mat_1[ii] = gfi_mat_0[ii];
    }
  }
}

