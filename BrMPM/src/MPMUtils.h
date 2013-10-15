/*
 * MPMUtils.h
 *
 *  Created on: 15/10/2013
 *      Author: banerjee
 */

#ifndef MPMUTILS_H_
#define MPMUTILS_H_

#include <MPMData.h>

using namespace BrMPM;

namespace MPMUtils
{
  // Integrate particle values and move to grid
  void integrate(VectorIntParticleData& cIdx,
                 VectorDoubleParticleData& cW,
                 DoubleParticleData& ppx,
                 DoubleParticleData& ppy,
                 DoubleParticleData& ppz,
                 DoubleParticleData& ggx,
                 DoubleParticleData& ggy,
                 DoubleParticleData& ggz);

  void interpolate(VectorIntParticleData& cIdx,
                   VectorDoubleParticleData& cW,
                   DoubleParticleData& ppx,
                   DoubleParticleData& ppy,
                   DoubleParticleData& ppz,
                   DoubleParticleData& ggx,
                   DoubleParticleData& ggy,
                   DoubleParticleData& ggz);

  void gradient(VectorIntParticleData& cIdx,
                VectorDoubleParticleData& cGradx,
                VectorDoubleParticleData& cGrady,
                VectorDoubleParticleData& cGradz,
                VectorDoubleParticleData& ppx,
                VectorDoubleParticleData& ppy,
                VectorDoubleParticleData& ppz,
                VectorDoubleParticleData& ggx,
                VectorDoubleParticleData& ggy,
                VectorDoubleParticleData& ggz);

  void divergence(VectorIntParticleData& cIdx,
                  VectorDoubleParticleData& cGradx,
                  VectorDoubleParticleData& cGrady,
                  VectorDoubleParticleData& cGradz,
                  VectorDoubleParticleData& ppx,
                  VectorDoubleParticleData& ppy,
                  VectorDoubleParticleData& ppz,
                  VectorDoubleParticleData& ggx,
                  VectorDoubleParticleData& ggy,
                  VectorDoubleParticleData& ggz);

  void gradscalar(VectorIntParticleData& cIdx,
                  VectorDoubleParticleData& cGradx,
                  VectorDoubleParticleData& cGrady,
                  VectorDoubleParticleData& cGradz,
                  VectorDoubleParticleData& ppx,
                  VectorDoubleParticleData& ppy,
                  VectorDoubleParticleData& ppz,
                  VectorDoubleParticleData& ggx,
                  VectorDoubleParticleData& ggy,
                  VectorDoubleParticleData& ggz);

  void dotAdd(VectorDoubleParticleData& ppx,
              VectorDoubleParticleData& ppy,
              VectorDoubleParticleData& ppz,
              VectorDoubleParticleData& ggx,
              VectorDoubleParticleData& ggy,
              VectorDoubleParticleData& ggz);

} // end namespace


#endif /* MPMUTILS_H_ */
