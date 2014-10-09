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
 * MPMUtils.cc
 *
 *  Created on: 15/10/2013
 *      Author: banerjee
 */

#include <MPMUtils.h>

template<typename T1, typename T2>
void MPMUtils::integrate(const VectorIntParticleData& cIdx,
                         const VectorDoubleParticleData& cW,
                         const T1& pp,
                         T2& gg)
{
  unsigned int nParts = pp.size();
  auto cIdxIter = cIdx.begin();
  unsigned int nContrib = (*cIdxIter).size();
  for (unsigned int ii = 0; ii < nParts; ++ii) {
    for (unsigned int jj = 0; jj < nContrib; ++jj) {
      int ixc = cIdx[ii][jj];
      double weight = cW[ii][jj];
      gg[ixc] += (pp[ii]*weight);
    }
  }
}

void MPMUtils::interpolate(const VectorIntParticleData& cIdx,
                           const VectorDoubleParticleData& cW,
                           Vector3DParticleData& pp,
                           const Vector3DNodeData& gg)
{
  unsigned int nParts = pp.size();
  auto cIdxIter = cIdx.begin();
  unsigned int nContrib = (*cIdxIter).size();
  for (unsigned int ii = 0; ii < nParts; ++ii) {
    pp[ii].set(0.0);
    for (unsigned int jj = 0; jj < nContrib; ++jj) {
      int ixc = cIdx[ii][jj];
      double weight = cW[ii][jj];
      pp[ii] += (gg[ixc]*weight);
    }
  }
}

void MPMUtils::gradient(const VectorIntParticleData& cIdx,
                        const VectorDoubleParticleData& cGradx,
                        const VectorDoubleParticleData& cGrady,
                        const VectorDoubleParticleData& cGradz,
                        Matrix3DParticleData& pp,
                        const Vector3DNodeData& gg)
{
  unsigned int nParts = pp.size();
  auto cIdxIter = cIdx.begin();
  unsigned int nContrib = (*cIdxIter).size();
  for (unsigned int ii = 0; ii < nParts; ++ii) {
    pp[ii].set(0.0);
    for (unsigned int jj = 0; jj < nContrib; ++jj) {
      int ixc = cIdx[ii][jj];
      double gradx = cGradx[ii][jj];
      double grady = cGrady[ii][jj];
      double gradz = cGradz[ii][jj];
      Vector3D grad(gradx, grady, gradz);
      pp[ii] += dyadicProduct(gg[ixc], grad);
    }
  }
}

void MPMUtils::divergence(const VectorIntParticleData& cIdx,
                          const VectorDoubleParticleData& cGradx,
                          const VectorDoubleParticleData& cGrady,
                          const VectorDoubleParticleData& cGradz,
                          const Matrix3DParticleData& pp,
                          Vector3DNodeData& gg)
{
  unsigned int nParts = pp.size();
  auto cIdxIter = cIdx.begin();
  unsigned int nContrib = (*cIdxIter).size();
  for (unsigned int ii = 0; ii < nParts; ++ii) {
    for (unsigned int jj = 0; jj < nContrib; ++jj) {
      int ixc = cIdx[ii][jj];
      double gradx = cGradx[ii][jj];
      double grady = cGrady[ii][jj];
      double gradz = cGradz[ii][jj];
      Vector3D grad(gradx, grady, gradz);
      gg[ixc] -= pp[ii]*grad;
    }
  }
}

void MPMUtils::gradscalar(const VectorIntParticleData& cIdx,
                          const VectorDoubleParticleData& cGradx,
                          const VectorDoubleParticleData& cGrady,
                          const VectorDoubleParticleData& cGradz,
                          const DoubleParticleData& pp,
                          Vector3DNodeData& gg)
{
  unsigned int nParts = pp.size();
  auto cIdxIter = cIdx.begin();
  unsigned int nContrib = (*cIdxIter).size();
  for (unsigned int ii = 0; ii < nParts; ++ii) {
    for (unsigned int jj = 0; jj < nContrib; ++jj) {
      int ixc = cIdx[ii][jj];
      double gradx = cGradx[ii][jj];
      double grady = cGrady[ii][jj];
      double gradz = cGradz[ii][jj];
      Vector3D grad(gradx, grady, gradz);
      gg[ixc] += (grad*pp[ii]);
    }
  }
}

void MPMUtils::dotAdd(Matrix3DParticleData& pp,
                      const Matrix3DParticleData& qq)
{
  auto pIter = pp.begin();
  auto qIter = qq.begin();
  for (; pIter != pp.end(); ++pIter, ++qIter){
    Matrix3D pMat = *pIter;
    Matrix3D qMat = *qIter;
    pMat += qMat*pMat;
    pp[pIter - pp.begin()] = pMat;
  }
}






