/*
 * MPMFrictionContact.cc
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#include <Contact/MPMFrictionContact.h>
#include <MPMPatch.h>
#include <MPMDatawarehouse.h>

using namespace BrMPM;

MPMFrictionContact::MPMFrictionContact(std::vector<int>& dwis, MPMPatchP& patch,
                                       double mu)
   : MPMFrictionlessContact(dwis, patch), d_mu(mu)
{
  d_dt = patch->dt();
}

MPMFrictionContact::~MPMFrictionContact() {
}

void MPMFrictionContact::exchForceInterpolated(MPMDatawarehouseP& dw)
{
  // Mass
  DoubleNodeData mr, ms;
  dw->get("gm", d_dwis[0], mr);
  dw->get("gm", d_dwis[1], ms);

  // Momentum
  Vector3DNodeData Pr, Ps;
  dw->get("gw", d_dwis[0], Pr);
  dw->get("gw", d_dwis[1], Ps);

  // External force
  Vector3DNodeData fr, fs;
  dw->get("gfe", d_dwis[0], fr);
  dw->get("gfe", d_dwis[1], fs);

  // Internal force
  Vector3DNodeData fir, fis;
  dw->get("gfi", d_dwis[0], fir);
  dw->get("gfi", d_dwis[1], fis);

  // Grid gradient*mass
  Vector3DNodeData gmr, gms;
  dw->get("gGm", d_dwis[0], gmr);
  dw->get("gGm", d_dwis[1], gms);

  for (auto iter = d_nodes.begin(); iter != d_nodes.end(); ++iter) {
    int ii = iter - d_nodes.begin();

    // Add small amount to mass
    double mrVal = mr[ii] + d_mtol;
    double msVal = ms[ii] + d_mtol;

    // Compute grid velocities from momentum (vr - vs)
    Vector3D vr_vs = Pr[ii]/mrVal - Ps[ii]/msVal;

    // Find normal
    Vector3D nn = (gmr[ii].normalized()- gms[ii].normalized())/2.0;

    // Find velocity tangent
    Vector3D tt0 = vr_vs - nn*(vr_vs.dot(nn));
    Vector3D tt = tt0/(tt0.length() + d_mtol);

    // Find magnitude mass weighted internal force resultant along normal
    Vector3D mf = fir[ii]*msVal - fis[ii]*mrVal;
    double psi = mf.dot(nn);
    psi = (psi > 0.0) ? psi : 0.0;

    // Find normal force
    Vector3D fnor = nn*(psi/(mrVal+msVal));

    // Update normal external forces
    fr[ii] -= (nn*(fr[ii].dot(nn)) + fnor);
    fs[ii] -= (nn*(fs[ii].dot(nn)) - fnor);

    // Update tangential external forces
    Vector3D mftan = Pr[ii]*msVal - Ps[ii]*mrVal + mf*d_dt;
    double ftan = mftan.dot(tt)/((mrVal+msVal)*d_dt);

    // Compute friction corrections
    double fr_tan = fr[ii].dot(tt) + fnor.length()*d_mu;
    double fs_tan = fs[ii].dot(tt) + fnor.length()*d_mu;
    Vector3D fr_inc = tt*std::min(fr_tan, ftan);
    Vector3D fs_inc = tt*std::min(fs_tan, ftan);

    // Update external forces
    fr[ii] -= fr_inc;
    fs[ii] += fs_inc;
  }

  dw->put("gfe", d_dwis[0], fr);
  dw->put("gfe", d_dwis[1], fs);
}
