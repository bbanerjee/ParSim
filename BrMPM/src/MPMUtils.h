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
  void integrate(const VectorIntParticleData& cIdx,
                 const VectorDoubleParticleData& cW,
                 const DoubleParticleData& ppx,
                 const DoubleParticleData& ppy,
                 const DoubleParticleData& ppz,
                 DoubleNodeData& ggx,
                 DoubleNodeData& ggy,
                 DoubleNodeData& ggz);

  void interpolate(const VectorIntParticleData& cIdx,
                   const VectorDoubleParticleData& cW,
                   DoubleParticleData& ppx,
                   DoubleParticleData& ppy,
                   DoubleParticleData& ppz,
                   const DoubleNodeData& ggx,
                   const DoubleNodeData& ggy,
                   const DoubleNodeData& ggz);

  void gradient(const VectorIntParticleData& cIdx,
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
                const DoubleNodeData& ggz);

  void divergence(const VectorIntParticleData& cIdx,
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
                  DoubleNodeData& ggz);

  void gradscalar(const VectorIntParticleData& cIdx,
                  const VectorDoubleParticleData& cGradx,
                  const VectorDoubleParticleData& cGrady,
                  const VectorDoubleParticleData& cGradz,
                  const DoubleParticleData& pp,
                  DoubleNodeData& ggx,
                  DoubleNodeData& ggy,
                  DoubleNodeData& ggz);

  void dotAdd(Matrix3DParticleData& pp,
              const Matrix3DParticleData& qq);

} // end namespace


#endif /* MPMUTILS_H_ */
