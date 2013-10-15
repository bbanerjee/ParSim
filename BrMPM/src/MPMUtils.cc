/*
 * MPMUtils.cc
 *
 *  Created on: 15/10/2013
 *      Author: banerjee
 */

#include <MPMUtils.h>

void MPMUtils::integrate(std::vector<int>& cIdx,
                         std::vector<double>& cW,
                         VectorParticleData& pp,
                         VectorGridData& gg)
{
  for (auto pIter = pp.begin(); pIter != pp.end(); ++pIter) {
    Vector3D pVec = *pIter;
  }
}

void MPMUtils::interpolate(std::vector<int>& cIdx,
    std::vector<double>& cW, std::vector<double>& pp, std::vector<double>& gg)
{
}

void MPMUtils::gradient(std::vector<int>& cIdx,
    std::vector<double>& cGrad, std::vector<double>& pp,
    std::vector<double>& gg)
{
}

void MPMUtils::divergence(std::vector<int>& cIdx,
    std::vector<double>& cGrad, std::vector<double>& pp,
    std::vector<double>& gg)
{
}

void MPMUtils::gradscalar(std::vector<int>& cIdx,
    std::vector<double>& cGrad, std::vector<double>& pp,
    std::vector<double>& gg)
{
}

void MPMUtils::dotAdd(std::vector<double>& pp, std::vector<double>& gg)
{
}


