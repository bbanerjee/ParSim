/*
 * MPMMaterial.cc
 *
 *  Created on: 11/10/2013
 *      Author: banerjee
 */

#include <MPMMaterial.h>
#include <MPMUtils.h>
#include <MPMDatawarehouse.h>
#include <MPMConstitutiveModel.h>
#include <MPMPatch.h>
#include <ShapeFunctions/MPMShapeFunction.h>

using namespace BrMPM;

// Material - holds update functions - default is deformable
// overridden by RigidMaterial for rigid materials
MPMMaterial::MPMMaterial(const int dwi,
		                 const MPMConstitutiveModelP& model,
		                 const MPMShapeFunctionP& shape)
  : d_dwi(dwi), d_model(model), d_shape(shape)
{
}

MPMMaterial::~MPMMaterial() {
}


void
MPMMaterial::updateContributions(MPMDatawarehouseP& dw,
		                             MPMPatchP& patch)
{
  dw->zeroGrid(d_dwi);
  d_shape->updateContribList(dw, patch, d_dwi);
}

void
MPMMaterial::setVelocity(MPMDatawarehouseP& dw,
		                     Vector3DParticleData& vel)
{
  DoubleParticleData pm;
  Point3DParticleData px;
  dw->get("pm", d_dwi, pm);
  dw->get("px", d_dwi, px);

  Vector3DParticleData pw;
  pw.resize(pm.size());
  for (auto iter = pm.begin(); iter != pm.end(); ++ iter) {
    int ii = iter - pm.begin();
	pw[ii] = vel[ii]*pm[ii];
  }
  dw->put("pw", d_dwi, pw);
}

void
MPMMaterial::setExternalLoad(MPMDatawarehouseP& dw,
                             Vector3DParticleData& fe)
{
  Vector3DParticleData pfe;
  dw->get("pfe", d_dwi, pfe);
  for (auto iter = pfe.begin(); iter != pfe.end(); ++ iter) {
    int ii = iter - pfe.begin();
    pfe[ii] = fe[ii];
  }
  dw->put("pfe", d_dwi, pfe);
}

void
MPMMaterial::setExternalAcceleration(MPMDatawarehouseP& dw,
                                     Vector3DParticleData& acc)
{
  DoubleParticleData pm;
  Vector3DParticleData pfe;
  dw->get("pm", d_dwi, pm);
  dw->get("pfe", d_dwi, pfe);
  for (auto iter = pm.begin(); iter != pm.end(); ++ iter) {
    int ii = iter - pm.begin();
	pfe[ii] = acc[ii]*pm[ii];
  }
  dw->put("pfe", d_dwi, pfe);
}

// Apply external loads to each material
void
MPMMaterial::applyExternalLoads(MPMDatawarehouseP& dw,
                                MPMPatchP& patch)
{
  VectorIntParticleData cIdx;
  VectorDoubleParticleData cW;
  dw->get("cIdx", d_dwi, cIdx);
  dw->get("cW", d_dwi, cW);

  Vector3DParticleData pp;
  Vector3DNodeData gg;
  dw->get("pfe", d_dwi, pp);
  dw->get("gfe", d_dwi, gg);

  MPMUtils::integrate(cIdx, cW, pp, gg);
  dw->put("gfe", d_dwi, gg);
}

// Interpolate particle mass and momentum to the grid
void
MPMMaterial::interpolateParticlesToGrid(MPMDatawarehouseP& dw,
                                        MPMPatchP& patch)
{
  VectorIntParticleData cIdx;
  VectorDoubleParticleData cW;
  dw->get("cIdx", d_dwi, cIdx);
  dw->get("cW", d_dwi, cW);

  // Mass
  DoubleParticleData pp;
  DoubleNodeData gg;
  dw->get("pm", d_dwi, pp);
  dw->get("gm", d_dwi, gg);
  MPMUtils::integrate(cIdx, cW, pp, gg);
  dw->put("gm", d_dwi, gg);

  // Momentum
  dw->get("pw", d_dwi, pp);
  dw->get("gw", d_dwi, gg);
  MPMUtils::integrate( cIdx, cW, pp, gg);
  dw->put("gw", d_dwi, gg);
}

void
MPMMaterial::computeStressTensor(MPMDatawarehouseP& dw,
                                 MPMPatchP& patch)
{
  Matrix3DParticleData pf, pvs;
  DoubleParticleData pv;
  Vector3DParticleData pn;

  // Deformation Gradient
  dw->get("pF", d_dwi, pf);
  // Volume * Stress
  dw->get("pVS", d_dwi, pvs);
  // Volume
  dw->get("pVol", d_dwi, pv);
  // Volume
  dw->get("pn", d_dwi, pn);

  for (auto iter = pv.begin(); iter != pv.end(); ++iter) {
    int ii = iter - pv.begin();
    Matrix3D S;
    double Ja;
    d_model->getStress(pf[ii], S, Ja);
    pn[ii] = Ja;
    pvs[ii] = S*(pv[ii]*Ja);
  }
}

// Compute internal body forces - integrate divergence of stress to grid
void
MPMMaterial::computeInternalForce(MPMDatawarehouseP& dw,
                                  MPMPatchP& patch)
{
  VectorIntParticleData cIdx;
  VectorDoubleParticleData cGradx;
  VectorDoubleParticleData cGrady;
  VectorDoubleParticleData cGradz;
  dw->get("cIdx", d_dwi, cIdx);
  dw->get("cGradx", d_dwi, cGradx);
  dw->get("cGrady", d_dwi, cGrady);
  dw->get("cGradz", d_dwi, cGradz);

  Matrix3DParticleData pp;
  Vector3DNodeData gg;
  dw->get("pVS", d_dwi, pp);
  dw->get("gfi", d_dwi, gg);

  MPMUtils::divergence(cIdx, cGradx, cGrady, cGradz, pp, gg);
  dw->put("gfi", d_dwi, gg);
}

// Integrate grid acceleration
void
MPMMaterial::computeAndIntegrateAcceleration(MPMDatawarehouseP& dw,
                                             MPMPatchP& patch,
                                             double& tol)
{
  // Initializes leap-frog
  double a_leap = 1.0 - (patch->it()==0) * 0.5;

  DoubleParticleData pm;
  dw->get("pm", d_dwi, pm);

  DoubleNodeData gm, gwc;
  Vector3DNodeData gw, gfi, gfe;

  // Mass
  dw->get("gm", d_dwi, gm);
  // Momentum
  dw->get("gw", d_dwi, gw);
  dw->get("gwc", d_dwi, gwc);
  // Internal force
  dw->get("gfi", d_dwi, gfi);
  // External force
  dw->get("gfe", d_dwi, gfe);

  VectorIntParticleData cIdx;
  VectorDoubleParticleData cW;
  VectorDoubleParticleData cGradx;
  VectorDoubleParticleData cGrady;
  VectorDoubleParticleData cGradz;
  dw->get("cIdx", d_dwi, cIdx);
  dw->get("cW", d_dwi, cW);
  dw->get("cGradx", d_dwi, cGradx);
  dw->get("cGrady", d_dwi, cGrady);
  dw->get("cGradz", d_dwi, cGradz);

  Vector3DParticleData pfi;
  pfi.resize(pm.size());
  MPMUtils::interpolate(cIdx, cW, pfi, gfi);
  dw->put("pfi", d_dwi, pfi);

  Vector3DParticleData pfc;
  pfc.resize(pm.size());
  MPMUtils::interpolate(cIdx, cW, pfc, gfe);
  dw->put("pfc", d_dwi, pfc);

  // Velocity
  Vector3DNodeData gv, ga;
  dw->get("gv", d_dwi, gv);
  dw->get("ga", d_dwi, ga);

  for (auto iter = gm.begin(); iter != gm.end(); ++iter) {
    int ii = iter - gm.begin();
    gm[ii] += tol;
    gv[ii] = (gw[ii] + gwc[ii])/gm[ii];
    ga[ii] = ((gfe[ii] + gfi[ii])/gm[ii])*a_leap;
    gv[ii] += ga[ii]*patch->dt();
  }

  dw->put("gv", d_dwi, gv);
  dw->put("ga", d_dwi, ga);
}

void
MPMMaterial::interpolateToParticlesAndUpdate(MPMDatawarehouseP& dw,
                                             MPMPatchP& patch)
{
  VectorIntParticleData cIdx;
  VectorDoubleParticleData cW;
  VectorDoubleParticleData cGradx;
  VectorDoubleParticleData cGrady;
  VectorDoubleParticleData cGradz;
  dw->get("cIdx", d_dwi, cIdx);
  dw->get("cW", d_dwi, cW);
  dw->get("cGradx", d_dwi, cGradx);
  dw->get("cGrady", d_dwi, cGrady);
  dw->get("cGradz", d_dwi, cGradz);

  Vector3DParticleData pvI;
  Vector3DParticleData pxI;
  Matrix3DParticleData pGv;
  Vector3DNodeData ga, gv;
  dw->get("pvI", d_dwi, pvI);
  dw->get("pxI", d_dwi, pxI);
  dw->get("pGv", d_dwi, pGv);
  dw->get("ga", d_dwi, ga);
  dw->get("gv", d_dwi, gv);

  MPMUtils::interpolate(cIdx, cW, pvI, ga);
  MPMUtils::interpolate(cIdx, cW, pxI, gv);
  MPMUtils::gradient(cIdx, cGradx, cGrady, cGradz, pGv, gv);

  Point3DParticleData px;
  Vector3DParticleData pw;
  DoubleParticleData pm;
  Matrix3DParticleData pF;
  dw->get("px", d_dwi, px);
  dw->get("pw", d_dwi, pw);
  dw->get("pm", d_dwi, pm);
  dw->get("pF", d_dwi, pF);

  for (auto iter = pm.begin(); iter != pm.end(); ++iter) {
    int ii = iter - pm.begin();
    pw[ii] = pxI[ii] * pm[ii];
    px[ii] += pxI[ii] * patch->dt();
    // pF += (pGv*dt).pF
    pF[ii] += (pGv[ii]*patch->dt())*pF[ii];
  }

  dw->put("pGv", d_dwi, pGv);
  dw->put("pF", d_dwi, pF);
  dw->put("pw", d_dwi, pw);
  dw->put("pvI", d_dwi, pvI);
  dw->put("pxI", d_dwi, pxI);
}

