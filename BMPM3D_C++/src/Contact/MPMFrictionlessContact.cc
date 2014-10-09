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
 * MPMFrictionlessContact.cc
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#include <Contact/MPMFrictionlessContact.h>
#include <MPMDatawarehouse.h>
#include <MPMUtils.h>

using namespace BrMPM;

MPMFrictionlessContact::MPMFrictionlessContact(std::vector<int>& dwis,
                                               MPMPatchP& patch)
  : MPMContact(dwis, patch)
{
}

MPMFrictionlessContact::~MPMFrictionlessContact()
{
}

void
MPMFrictionlessContact::findIntersection(MPMDatawarehouseP& dw)
{
  for (auto iter = d_dwis.begin(); iter != d_dwis.end(); ++iter) {
    int dwi = *iter;

    // Get the particle interpolation information
    VectorIntParticleData cIdx;
    VectorDoubleParticleData cGradx;
    VectorDoubleParticleData cGrady;
    VectorDoubleParticleData cGradz;
    dw->get("cIdx", dwi, cIdx);
    dw->get("cGradx", dwi, cGradx);
    dw->get("cGrady", dwi, cGrady);
    dw->get("cGradz", dwi, cGradz);

    // Get the particle state information
    DoubleParticleData pm;
    DoubleParticleData pVol;
    dw->get("pm", dwi, pm);
    dw->get("pVol", dwi, pVol);

    // Get the node data
    Vector3DNodeData gGm;
    dw->get("gGm", dwi, gGm);

    // Compute gradient
    MPMUtils::gradscalar(cIdx, cGradx, cGrady, cGradz, pm, gGm);
  }

  // Call Contact findIntersection
  MPMContact::findIntersection(dw);
}

void
MPMFrictionlessContact::exchMomentumInterpolated(MPMDatawarehouseP& dw)
{
  findIntersection(dw);

  DoubleNodeData mr, ms;
  Vector3DNodeData Pr, Ps;
  Vector3DNodeData Vr, Vs;
  DoubleNodeData gfc;
  dw->get("gm", d_dwis[0], mr);
  dw->get("gm", d_dwis[1], ms);
  dw->get("gw", d_dwis[0], Pr);
  dw->get("gw", d_dwis[1], Ps);
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
    Pr[ii] -= nn*dp;
    Ps[ii] += nn*dp;
  }

  dw->put("gw", d_dwis[0], Pr);
  dw->put("gw", d_dwis[1], Ps);
  dw->put("gfc", d_dwis[1], gfc);
}

void
MPMFrictionlessContact::exchForceInterpolated(MPMDatawarehouseP& dw)
{
  DoubleNodeData mr, ms;
  dw->get("gm", d_dwis[0], mr);
  dw->get("gm", d_dwis[1], ms);

  Vector3DNodeData fr, fs;
  dw->get("gfe", d_dwis[0], fr);
  dw->get("gfe", d_dwis[1], fs);

  Vector3DNodeData fir, fis;
  dw->get("gfi", d_dwis[0], fir);
  dw->get("gfi", d_dwis[1], fis);

  Vector3DNodeData gmr, gms;
  dw->get("gGm", d_dwis[0], gmr);
  dw->get("gGm", d_dwis[1], gms);

  for (auto iter = d_nodes.begin(); iter != d_nodes.end(); ++iter) {
    int ii = iter - d_nodes.begin();
    Vector3D nn = (gmr[ii].normalized()- gms[ii].normalized())/2.0;
    double mrVal = mr[ii];
    double msVal = ms[ii];
    Vector3D mf = fir[ii]*msVal - fis[ii]*mrVal;
    double psi = mf.dot(nn);
    psi = (psi > 0.0) ? psi : 0.0;

    Vector3D fnor = nn*(psi/(mrVal+msVal));
    fr[ii] -= (nn*(fr[ii].dot(nn)) + fnor);
    fs[ii] -= (nn*(fs[ii].dot(nn)) - fnor);
  }

  dw->put("gfe", d_dwis[0], fr);
  dw->put("gfe", d_dwis[1], fs);

}
