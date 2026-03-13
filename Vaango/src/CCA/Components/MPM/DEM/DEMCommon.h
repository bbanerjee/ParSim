/*
 * The MIT License
 *
 * Copyright (c) 2013-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#ifndef __VAANGO_CCA_COMPONENTS_MPM_DEM_COMMON_H__
#define __VAANGO_CCA_COMPONENTS_MPM_DEM_COMMON_H__

#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Math/Matrix3.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarTypes.h>

namespace Vaango {

using constParticleVarPoint = Uintah::constParticleVariable<Uintah::Point>;
using constParticleVarVector = Uintah::constParticleVariable<Uintah::Vector>;
using constParticleVarMatrix3 = Uintah::constParticleVariable<Uintah::Matrix3>;
using constParticleVarLong64 = Uintah::constParticleVariable<Uintah::long64>;
using constParticleVarDouble = Uintah::constParticleVariable<double>;
using ParticleVarVector = Uintah::ParticleVariable<Uintah::Vector>;
using ParticleVarLong64 = Uintah::ParticleVariable<Uintah::long64>;
using ParticleIDToCurrentIdxMap = std::map<Uintah::long64, int>;

// DEM force calculation data structures
struct DEMParticleInputData {
  explicit DEMParticleInputData(int num_mats)
    : pX_old(num_mats), pX0_old(num_mats), 
      pMass_old(num_mats), pRadius_old(num_mats), pSize_old(num_mats),
      pOrientation_old(num_mats), pVelocity_old(num_mats), pAngVel_old(num_mats),
      pRigidBodyID_old(num_mats) 
  {}

  std::vector<constParticleVarPoint>   pX_old, pX0_old;
  std::vector<constParticleVarDouble>  pMass_old, pRadius_old;
  std::vector<constParticleVarMatrix3> pSize_old, pOrientation_old;
  std::vector<constParticleVarVector>  pVelocity_old, pAngVel_old;
  std::vector<constParticleVarLong64>  pRigidBodyID_old;
};

struct DEMParticleSets {
  explicit DEMParticleSets(int num_mats)
    : psets_all(num_mats, nullptr), psets_real(num_mats, nullptr)
  {}

  std::vector<Uintah::ParticleSubset*>  psets_all;
  std::vector<Uintah::ParticleSubset*>  psets_real;
};

struct DEMParticleOutputData {
  explicit DEMParticleOutputData(int num_mats)
    : pExtForce_new(num_mats), pTorque_new(num_mats), pRigidBodyID_new(num_mats)
  {}

  std::vector<ParticleVarVector> pExtForce_new;
  std::vector<ParticleVarVector> pTorque_new;
  std::vector<ParticleVarLong64> pRigidBodyID_new;
};

struct DEMContactResult {
  Vector totalForce   { 0, 0, 0 };
  Vector arm_i        { 0, 0, 0 };
  Vector arm_j_center { 0, 0, 0 };
  bool   collision    { false };
};

// Contact property helpers 
struct DEMContactProps {
  double kn, kt, gamma, mu;
};

// Hash for IntVector so it can be used as an unordered_map key
struct IntVectorHash {
  size_t operator()(const IntVector& k) const {
    size_t h = std::hash<int>{}(k.x());
    h ^= std::hash<int>{}(k.y()) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= std::hash<int>{}(k.z()) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
  }
};

} // namespace Vaango

#endif //__VAANGO_CCA_COMPONENTS_MPM_DEM_COMMON_H__