/*
 * MPMUtils.h
 *
 *  Created on: 15/10/2013
 *      Author: banerjee
 */

#ifndef MPMUTILS_H_
#define MPMUTILS_H_

#include <MPMDataTypes.h>

using namespace BrMPM;

namespace MPMUtils
{
  // Integrate particle values and move to grid
  template<typename T1, typename T2>
  void integrate(const VectorIntParticleData& cIdx,
                 const VectorDoubleParticleData& cW,
                 const T1& pp,
                 T2& gg);

  void interpolate(const VectorIntParticleData& cIdx,
                   const VectorDoubleParticleData& cW,
                   Vector3DParticleData& pp,
                   const Vector3DNodeData& gg);

  void gradient(const VectorIntParticleData& cIdx,
                const VectorDoubleParticleData& cGradx,
                const VectorDoubleParticleData& cGrady,
                const VectorDoubleParticleData& cGradz,
                Matrix3DParticleData& pp,
                const Vector3DNodeData& gg);

  void divergence(const VectorIntParticleData& cIdx,
                  const VectorDoubleParticleData& cGradx,
                  const VectorDoubleParticleData& cGrady,
                  const VectorDoubleParticleData& cGradz,
                  const Matrix3DParticleData& pp,
                  Vector3DNodeData& gg);

  void gradscalar(const VectorIntParticleData& cIdx,
                  const VectorDoubleParticleData& cGradx,
                  const VectorDoubleParticleData& cGrady,
                  const VectorDoubleParticleData& cGradz,
                  const DoubleParticleData& pp,
                  Vector3DNodeData& gg);

  void dotAdd(Matrix3DParticleData& pp,
              const Matrix3DParticleData& qq);


} // end namespace

// template instantiations
namespace MPMUtils
{
  template<>
  void integrate<DoubleParticleData, DoubleNodeData>(const VectorIntParticleData& cIdx,
                 const VectorDoubleParticleData& cW,
                 const DoubleParticleData& pp,
                 DoubleNodeData& gg);
  template<>
  void integrate<Vector3DParticleData, Vector3DNodeData>(const VectorIntParticleData& cIdx,
                 const VectorDoubleParticleData& cW,
                 const Vector3DParticleData& pp,
                 Vector3DNodeData& gg);

} // end MPMUtils namespace


#endif /* MPMUTILS_H_ */
