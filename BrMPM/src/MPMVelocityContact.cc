/*
 * MPMVelocityContact.cc
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#include <MPMVelocityContact.h>
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
MPMVelocityContact::exchForceInterpolated(MPMDatawarehouseP& dw)
{
  std::cout << "No exchange force interpolated implemented in MPMVelocityContact" << std::endl;
}

void
MPMVelocityContact::exchMomentumIntegrated(MPMDatawarehouseP& dw)
{
  findIntersection(dw);

  DoubleNodeData mr, ms;
  dw->get("gm", d_dwis[0], mr);
  dw->get("gm", d_dwis[1], ms);

  Vector3DNodeData Pr, Ps;
  dw->get("gw", d_dwis[0], Pr);
  dw->get("gw", d_dwis[1], Ps);

  Vector3DNodeData Vr, Vs;
  DoubleNodeData gfc;
  dw->get("gv", d_dwis[0], Vr);
  dw->get("gv", d_dwis[1], Vs);
  dw->get("gfc", d_dwis[1], gfc);

  Vector3DNodeData gmr, gms;
  dw->get("gGm", d_dwis[0], gmr);
  dw->get("gGm", d_dwis[1], gms);

  for (auto iter = d_nodes.begin(); iter != d_nodes.end(); ++iter) {
    int ii = iter - d_nodes.begin();
    Vector3D nn = (gmr[ii].normalized()- gms[ii].normalized())/2.0;

    double mrVal = mr[ii];
    double msVal = ms[ii];
    Vector3D dp0 = (Pr[ii]*msVal - Ps[ii]*mrVal)/(mrVal+msVal);
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
