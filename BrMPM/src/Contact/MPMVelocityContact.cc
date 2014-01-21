/*
 * MPMVelocityContact.cc
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#include <Contact/MPMVelocityContact.h>
#include <MPMDatawarehouse.h>
#include <MPMDataTypes.h>
#include <MPMUtils.h>

using namespace BrMPM;

MPMVelocityContact::MPMVelocityContact(std::vector<int>& dwis,
                                       MPMPatchP& patch)
   : MPMFrictionlessContact(dwis, patch)
{
}

MPMVelocityContact::~MPMVelocityContact()
{
}

void
MPMVelocityContact::exchMomentumIntegrated(MPMDatawarehouseP& dw)
{
  DoubleNodeData mr, ms;
  dw->get("gm", d_dwis[0], mr);
  dw->get("gm", d_dwis[1], ms);

  Vector3DNodeData Vr, Vs;
  dw->get("gv", d_dwis[0], Vr);
  dw->get("gv", d_dwis[1], Vs);

  DoubleNodeData gfc;
  dw->get("gfc", d_dwis[1], gfc);

  Vector3DNodeData gmr, gms;
  dw->get("gGm", d_dwis[0], gmr);
  dw->get("gGm", d_dwis[1], gms);

  for (auto iter = d_nodes.begin(); iter != d_nodes.end(); ++iter) {
    int ii = iter - d_nodes.begin();

    double mrVal = mr[ii];
    double msVal = ms[ii];

    Vector3D Pr = Vr[ii]*mrVal;
    Vector3D Ps = Vs[ii]*msVal;

    Vector3D nn = (gmr[ii].normalized()- gms[ii].normalized())/2.0;

    Vector3D dp0 = (Pr*msVal - Ps*mrVal)/(mrVal+msVal);
    double dp = dp0.dot(nn);
    dp = (dp > 0.0) ? dp : 0.0;

    gfc[ii] = 1.0;
    Vr[ii] -= nn*(dp/mrVal);
    Vs[ii] += nn*(dp/msVal);
  }

  dw->put("gv", d_dwis[0], Vr);
  dw->put("gv", d_dwis[1], Vs);
  dw->put("gfc", d_dwis[1], gfc);
}
