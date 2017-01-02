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
