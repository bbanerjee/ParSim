/*
 * GIMPShapeFunction.cc
 *
 *  Created on: 3/12/2013
 *      Author: banerjee
 */

#include <ShapeFunctions/GIMPShapeFunction.h>
#include <MPMDatawarehouse.h>
#include <MPMPatch.h>

using namespace BrMPM;

GIMPShapeFunction::GIMPShapeFunction()
{
}

GIMPShapeFunction::~GIMPShapeFunction()
{
}

// Update the nodes that contribute to each particle
void
GIMPShapeFunction::updateContribList(MPMDatawarehouseP& dw,
                                     const MPMPatchP& patch,
                                     int dwi)
{
  int nx = (patch->nC())[0];
  int ny = patch->nC()[1];
  Vector3D hh = patch->dX();
  std::vector<int> indices = {0, 1, 2, nx, nx+1, nx+2, nx*ny, nx*ny+1, nx*ny+2,
                              nx*ny+nx, nx*ny+nx+1, nx*ny+nx+2};
  Vector3D weights(0.0), derivatives(0.0);

  // Get the particle interpolation information
  VectorIntParticleData cIdx;
  VectorDoubleParticleData cW;
  VectorDoubleParticleData cGradx;
  VectorDoubleParticleData cGrady;
  VectorDoubleParticleData cGradz;
  dw->get("cIdx", dwi, cIdx);
  dw->get("cW", dwi, cW);
  dw->get("cGradx", dwi, cGradx);
  dw->get("cGrady", dwi, cGrady);
  dw->get("cGradz", dwi, cGradz);

  // Get the positions
  Point3DParticleData px;
  dw->get("px", dwi, px);

  Point3DNodeData gx;
  dw->get("gx", dwi, gx);

  // Get the particle volume and deformation gradient
  DoubleParticleData pVol;
  Matrix3DParticleData pDefGrad;
  dw->get("pVol", dwi, pVol);
  dw->get("pF", dwi, pDefGrad);

  // Get particle distance information for contact computation
  DoubleParticleData gDist;
  dw->get("gDist", dwi, gDist);

  // Loop thru particles
  for (auto piter = px.begin(); piter != px.end(); ++piter) {

    // Get the particle index
    int ii = piter - px.begin();

    // Compute the particle size based on deformation gradient (TODO)
    // For now just use 1/dx
    Vector3D pSize = hh.invDirection();

    // Get the cell index
    int cc = getCell(patch, px[ii]);

    // Loop through the node indices that contribute the the GIMP shape function
    for (auto niter = indices.begin(); niter != indices.end(); ++niter) {

      // Get the node index from the neighbor list
      int jj = niter - indices.begin();

      // Compute the node index in the 1-D full node array
      int idx = indices[jj] + cc;

      // Compute vector from particle to the node
      Vector3D rr = px[ii] - gx[idx];
      double dist = rr.length();

      // Compute the distance measure for the node-particle pair
      gDist[idx] = std::max(0.0, std::max(gDist[idx], (1.0 - dist/hh.min())));

      // Get the node weights and shape derivatives
      for (int kk = 0; kk < 3; ++kk) {
        getWeightAndDerivative(rr[kk], hh[kk], pSize[kk], weights[kk], derivatives[kk]);
      }

      // Update the data arrays
      cIdx[ii][jj] = idx;
      cW[ii][jj] = weights[0]*weights[1]*weights[2];
      cGradx[ii][jj] = weights[1]*weights[2]*derivatives[0];
      cGrady[ii][jj] = weights[0]*weights[2]*derivatives[1];
      cGradz[ii][jj] = weights[0]*weights[1]*derivatives[2];
    }
  }

  // Save the data
  dw->put("cIdx", dwi, cIdx);
  dw->put("cW", dwi, cW);
  dw->put("cGradx", dwi, cGradx);
  dw->put("cGrady", dwi, cGrady);
  dw->put("cGradz", dwi, cGradz);
  dw->put("gDist", dwi, gDist);

}

// Get the cell index for a particle at position "pos"
int
GIMPShapeFunction::getCell(const MPMPatchP& patch, const Point3D& pos)
{
  // Find the cell index based on the location wrt to the grid origin and cell size
  Vector3D x_sc = (pos - patch->x0())/patch->dX() + patch->nGhost();

  // Take the floor to get the cell index
  IntVector3D idx = x_sc.floor();

  // If the computed cell index and the floor cell index differ by 0.5 then decrease
  // cell index by 1
  for (int cc = 0; cc < 3; ++cc) {
    if ((x_sc[cc] - (double) idx[cc]) < 0.5) idx[cc] -= 1;
  }

  // Compute the cell index location in the 1D array
  int ii = idx[0];
  int jj = idx[1];
  int kk = idx[2];
  return (kk*patch->nC()[1]*patch->nC()[0]+jj*patch->nC()[0]+ii);
}

// Get the shape functions weights and derivatives
void
GIMPShapeFunction::getWeightAndDerivative(const double& xx,
                                          const double& hh,
                                          const double& ll,
                                          double& weight,
                                          double& gradient)
{
  double rr = std::abs(xx);
  double signx = std::copysign(1.0, xx);

  if (rr < ll) {
    weight = 1.0 - (rr*rr + ll*ll)/(2.0*hh*ll);
    gradient = -xx/(hh*ll);
  } else if (rr < hh - ll) {
    weight = 1.0 - rr/hh;
    gradient = -signx/hh;
  } else if (rr < hh + ll) {
    weight = (hh + ll - rr)*(hh + ll - rr)/(4.0 *hh*ll);
    gradient = (hh + ll -rr)/(-2.0*signx*hh*ll);
  } else {
    weight = 0.0;
    gradient = 0.0;
  }
}
