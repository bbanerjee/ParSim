/*
 * MPMUtils.cc
 *
 *  Created on: 15/10/2013
 *      Author: banerjee
 */

#include <MPMUtils.h>

void MPMUtils::integrate(const VectorIntParticleData& cIdx,
                         const VectorDoubleParticleData& cW,
                         const Vector3DParticleData& pp,
                         Vector3DNodeData& gg)
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






