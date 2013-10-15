/*
 * MPMUtils.h
 *
 *  Created on: 15/10/2013
 *      Author: banerjee
 */

#ifndef MPMUTILS_H_
#define MPMUTILS_H_

#include <MPMParticleData.h>
#include <vector>

namespace MPMUtils
{
  // Integrate particle values and move to grid
  void integrate(std::vector<int>& cIdx,
                 std::vector<double>& cW,
                 VectorParticleData& pp,
                 std::vector<double>& gg);

  void interpolate(std::vector<int>& cIdx,
                   std::vector<double>& cW,
                   VectorParticleData& pp,
                   std::vector<double>& gg);

  void gradient(std::vector<int>& cIdx,
                std::vector<double>& cGrad,
                VectorParticleData& pp,
                std::vector<double>& gg);

  void divergence(std::vector<int>& cIdx,
                  std::vector<double>& cGrad,
                  VectorParticleData& pp,
                  std::vector<double>& gg);

  void gradscalar(std::vector<int>& cIdx,
                  std::vector<double>& cGrad,
                  VectorParticleData& pp,
                  std::vector<double>& gg);

  void dotAdd(VectorParticleData& pp,
              std::vector<double>& gg);

} // end namespace


#endif /* MPMUTILS_H_ */
