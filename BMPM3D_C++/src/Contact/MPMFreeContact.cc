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
    dw->put("gm", d_dwis[0], gm_mat_0);
    dw->put("gm", d_dwis[1], gm_mat_1);

    dw->put("gw", d_dwis[0], gw_mat_0);
    dw->put("gw", d_dwis[1], gw_mat_1);
  }

}

void
MPMFreeContact::exchForceInterpolated(MPMDatawarehouseP& dw) {

  if (d_nodes.size() > 0) {
    Vector3DNodeData gfi_mat_0, gfi_mat_1;

    dw->get("gfi", d_dwis[0], gfi_mat_0);
    dw->get("gfi", d_dwis[1], gfi_mat_1);
    for (unsigned int ii = 0; ii < d_nodes.size(); ii++) {
      gfi_mat_0[ii] += gfi_mat_1[ii];
    }
    for (unsigned int ii = 0; ii < d_nodes.size(); ii++) {
      gfi_mat_1[ii] = gfi_mat_0[ii];
    }

    dw->put("gfi", d_dwis[0], gfi_mat_0);
    dw->put("gfi", d_dwis[1], gfi_mat_1);
  }
}

