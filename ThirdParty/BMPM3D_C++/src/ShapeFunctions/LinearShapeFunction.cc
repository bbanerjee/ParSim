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
 * LinearShapeFunction.cc
 *
 *  Created on: 22/10/2013
 *      Author: banerjee
 */

#include <ShapeFunctions/LinearShapeFunction.h>
#include <MPMPatch.h>
#include <MPMDatawarehouse.h>
#include <vector>
#include <cmath>

using namespace BrMPM;

LinearShapeFunction::LinearShapeFunction()
{
}

LinearShapeFunction::~LinearShapeFunction()
{
}

// Update node contribution list
void
LinearShapeFunction::updateContribList(MPMDatawarehouseP& dw,
                                       const MPMPatchP& patch, int dwi)
{
  int nx = (patch->nC())[0];
  int ny = patch->nC()[1];
  Vector3D hh = patch->dX();
  std::vector<int> indices = {0, 1, nx, nx+1, nx*ny, nx*ny+1, nx*ny+nx, nx*ny+nx+1};
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

  // Loop thru particles
  for (auto piter = px.begin(); piter != px.end(); ++piter) {
    int ii = piter - px.begin();
    int cc = getCell(patch, px[ii]);
    for (auto niter = indices.begin(); niter != indices.end(); ++niter) {
      int jj = niter - indices.begin();
      int idx = indices[jj] + cc;
      Vector3D rr = px[ii] - gx[idx];
      for (int kk = 0; kk < 3; ++kk) {
        getWeightAndDerivative(rr[kk], hh[kk], weights[kk], derivatives[kk]);
      }
      cIdx[ii][jj] = idx;
      cW[ii][jj] = weights[0]*weights[1]*weights[2];
      cGradx[ii][jj] = weights[1]*weights[2]*derivatives[0];
      cGrady[ii][jj] = weights[0]*weights[2]*derivatives[1];
      cGradz[ii][jj] = weights[0]*weights[1]*derivatives[2];
    }
  }

  dw->put("cIdx", dwi, cIdx);
  dw->put("cW", dwi, cW);
  dw->put("cGradx", dwi, cGradx);
  dw->put("cGrady", dwi, cGrady);
  dw->put("cGradz", dwi, cGradz);

}

// Gets lower left node of 8-cell block
int
LinearShapeFunction::getCell(const MPMPatchP& patch,
                             const Point3D& pos)
{
  Vector3D x_sc = (pos - patch->x0())/patch->dX() + patch->nGhost();
  IntVector3D idx = x_sc.floor();
  int ii = idx[0];
  int jj = idx[1];
  int kk = idx[2];
  return (kk*patch->nC()[1]*patch->nC()[0]+jj*patch->nC()[0]+ii);
}

void
LinearShapeFunction::getWeightAndDerivative(const double& xx, const double& hh,
                                            double& weight, double& gradient)
{
  double rr = std::abs(xx);
  double signx = std::copysign(1.0, xx);

  if (rr < hh) {
    weight = 1.0 - rr/hh;
    gradient = -signx/hh;
  } else {
    weight = 0.0;
    gradient = 0.0;
  }
}


