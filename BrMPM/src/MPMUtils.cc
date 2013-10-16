/*
 * MPMUtils.cc
 *
 *  Created on: 15/10/2013
 *      Author: banerjee
 */

#include <MPMUtils.h>

void MPMUtils::integrate(const VectorIntParticleData& cIdx,
                         const VectorDoubleParticleData& cW,
                         const DoubleParticleData& ppx,
                         const DoubleParticleData& ppy,
                         const DoubleParticleData& ppz,
                         DoubleNodeData& ggx,
                         DoubleNodeData& ggy,
                         DoubleNodeData& ggz)
{
  unsigned int nParts = ppx.size();
  auto cIdxIter = cIdx.begin();
  unsigned int nContrib = (*cIdxIter).size();
  for (unsigned int ii = 0; ii < nParts; ++ii) {
    for (unsigned int jj = 0; jj < nContrib; ++jj) {
      int ixc = cIdx[ii][jj];
      double weight = cW[ii][jj];
      ggx[ixc] += ppx[ii]*weight;
      ggy[ixc] += ppy[ii]*weight;
      ggz[ixc] += ppz[ii]*weight;
    }
  }
}

void MPMUtils::interpolate(const VectorIntParticleData& cIdx,
                           const VectorDoubleParticleData& cW,
                           DoubleParticleData& ppx,
                           DoubleParticleData& ppy,
                           DoubleParticleData& ppz,
                           const DoubleNodeData& ggx,
                           const DoubleNodeData& ggy,
                           const DoubleNodeData& ggz)
{
  unsigned int nParts = ppx.size();
  auto cIdxIter = cIdx.begin();
  unsigned int nContrib = (*cIdxIter).size();
  for (unsigned int ii = 0; ii < nParts; ++ii) {
    ppx[ii] = 0.0;
    ppy[ii] = 0.0;
    ppz[ii] = 0.0;
    for (unsigned int jj = 0; jj < nContrib; ++jj) {
      int ixc = cIdx[ii][jj];
      double weight = cW[ii][jj];
      ppx[ii] += ggx[ixc]*weight;
      ppy[ii] += ggx[ixc]*weight;
      ppz[ii] += ggx[ixc]*weight;
    }
  }
}


void MPMUtils::gradient(const VectorIntParticleData& cIdx,
                        const VectorDoubleParticleData& cGradx,
                        const VectorDoubleParticleData& cGrady,
                        const VectorDoubleParticleData& cGradz,
                        DoubleParticleData& ppxx,
                        DoubleParticleData& ppxy,
                        DoubleParticleData& ppxz,
                        DoubleParticleData& ppyx,
                        DoubleParticleData& ppyy,
                        DoubleParticleData& ppyz,
                        DoubleParticleData& ppzx,
                        DoubleParticleData& ppzy,
                        DoubleParticleData& ppzz,
                        const DoubleNodeData& ggx,
                        const DoubleNodeData& ggy,
                        const DoubleNodeData& ggz)
{
  unsigned int nParts = ppxx.size();
  auto cIdxIter = cIdx.begin();
  unsigned int nContrib = (*cIdxIter).size();
  for (unsigned int ii = 0; ii < nParts; ++ii) {
    ppxx[ii] = 0.0;
    ppxy[ii] = 0.0;
    ppxz[ii] = 0.0;
    ppyx[ii] = 0.0;
    ppyy[ii] = 0.0;
    ppyz[ii] = 0.0;
    ppzx[ii] = 0.0;
    ppzy[ii] = 0.0;
    ppzz[ii] = 0.0;
    for (unsigned int jj = 0; jj < nContrib; ++jj) {
      int ixc = cIdx[ii][jj];
      double gradx = cGradx[ii][jj];
      double grady = cGrady[ii][jj];
      double gradz = cGradz[ii][jj];
      ppxx[ii] += ggx[ixc]*gradx;
      ppxy[ii] += ggx[ixc]*grady;
      ppxz[ii] += ggx[ixc]*gradz;
      ppyx[ii] += ggy[ixc]*gradx;
      ppyy[ii] += ggy[ixc]*grady;
      ppyz[ii] += ggy[ixc]*gradz;
      ppzx[ii] += ggz[ixc]*gradx;
      ppzy[ii] += ggz[ixc]*grady;
      ppzz[ii] += ggz[ixc]*gradz;
    }
  }
}

void MPMUtils::divergence(const VectorIntParticleData& cIdx,
                          const VectorDoubleParticleData& cGradx,
                          const VectorDoubleParticleData& cGrady,
                          const VectorDoubleParticleData& cGradz,
                          const DoubleParticleData& ppxx,
                          const DoubleParticleData& ppxy,
                          const DoubleParticleData& ppxz,
                          const DoubleParticleData& ppyx,
                          const DoubleParticleData& ppyy,
                          const DoubleParticleData& ppyz,
                          const DoubleParticleData& ppzx,
                          const DoubleParticleData& ppzy,
                          const DoubleParticleData& ppzz,
                          DoubleNodeData& ggx,
                          DoubleNodeData& ggy,
                          DoubleNodeData& ggz)
{
  unsigned int nParts = ppxx.size();
  auto cIdxIter = cIdx.begin();
  unsigned int nContrib = (*cIdxIter).size();
  for (unsigned int ii = 0; ii < nParts; ++ii) {
    for (unsigned int jj = 0; jj < nContrib; ++jj) {
      int ixc = cIdx[ii][jj];
      double gradx = cGradx[ii][jj];
      double grady = cGrady[ii][jj];
      double gradz = cGradz[ii][jj];
      ggx[ixc] -= ppxx[ii]*gradx + ppxy[ii]*grady + ppxz[ii]*gradz;
      ggy[ixc] -= ppyx[ii]*gradx + ppyy[ii]*grady + ppyz[ii]*gradz;
      ggz[ixc] -= ppzx[ii]*gradx + ppzy[ii]*grady + ppzz[ii]*gradz;
    }
  }
}

void MPMUtils::gradscalar(const VectorIntParticleData& cIdx,
                          const VectorDoubleParticleData& cGradx,
                          const VectorDoubleParticleData& cGrady,
                          const VectorDoubleParticleData& cGradz,
                          const DoubleParticleData& pp,
                          DoubleNodeData& ggx,
                          DoubleNodeData& ggy,
                          DoubleNodeData& ggz)
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
      ggx[ixc] += pp[ii]*gradx;
      ggy[ixc] += pp[ii]*grady;
      ggz[ixc] += pp[ii]*gradz;
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






