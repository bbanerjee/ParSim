/*
 * MPMUtils.cc
 *
 *  Created on: 15/10/2013
 *      Author: banerjee
 */

#include <MPMUtils.h>

void MPMUtils::integrate(VectorIntParticleData& cIdx,
                         VectorDoubleParticleData& cW,
                         DoubleParticleData& ppx,
                         DoubleParticleData& ppy,
                         DoubleParticleData& ppz,
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

void MPMUtils::interpolate(VectorIntParticleData& cIdx,
                           VectorDoubleParticleData& cW,
                           DoubleParticleData& ppx,
                           DoubleParticleData& ppy,
                           DoubleParticleData& ppz,
                           DoubleNodeData& ggx,
                           DoubleNodeData& ggy,
                           DoubleNodeData& ggz)
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

void MPMUtils::gradient(NodeIndexParticleData& cIdx,
    NodeWeightParticleData& cGrad, VectorParticleData& pp,
    VectorParticleData& gg)
{
}

void MPMUtils::divergence(NodeIndexParticleData& cIdx,
    NodeWeightParticleData& cGrad, VectorParticleData& pp,
    VectorParticleData& gg)
{
}

void MPMUtils::gradscalar(NodeIndexParticleData& cIdx,
    NodeWeightParticleData& cGrad, VectorParticleData& pp,
    VectorParticleData& gg)
{
}

void MPMUtils::dotAdd(VectorParticleData& pp, VectorParticleData& gg)
{
}


