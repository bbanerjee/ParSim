/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#ifndef __VAANGO_CCA_COMPONENTS_PERIDYNAMICS_PARTICLE_CREATOR_H__
#define __VAANGO_CCA_COMPONENTS_PERIDYNAMICS_PARTICLE_CREATOR_H__

#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>

#include <Core/Parallel/CrowdMonitor.h>

#include <map>
#include <vector>

namespace Uintah {
class GeometryObject;
class Patch;
class DataWarehouse;
class ParticleSubset;
class VarLabel;
} // namespace Uintah

namespace Vaango {

using particleIndex = int;
using particleId    = int;

class PeridynamicsFlags;
class PeridynamicsMaterial;
class PeridynamicsLabel;

class ParticleCreator
{

public:
  ParticleCreator(PeridynamicsMaterial* matl, const PeridynamicsFlags* flags);
  virtual ~ParticleCreator();

  using VecGeometryObjectSP =
    std::vector<std::shared_ptr<Uintah::GeometryObject>>;

  virtual particleIndex
  createParticles(PeridynamicsMaterial* matl,
                  Uintah::CCVariable<short int>& cellNAPID,
                  const Uintah::Patch*,
                  Uintah::DataWarehouse* new_dw,
                  const VecGeometryObjectSP& objects);

  virtual void
  registerPermanentParticleState(PeridynamicsMaterial* matl);

  std::vector<const Uintah::VarLabel*>
  returnParticleState();

  std::vector<const Uintah::VarLabel*>
  returnParticleStatePreReloc();

protected:
  using GeomName = std::pair<std::string, Uintah::GeometryObject*>;
  using GeomPoint =
    std::map<Uintah::GeometryObject*, std::vector<Uintah::Point>>;
  using GeomScalar = std::map<GeomName, std::vector<double>>;
  using GeomVector = std::map<GeomName, std::vector<Uintah::Vector>>;
  using GeomTensor = std::map<GeomName, std::vector<Uintah::Matrix3>>;

  struct ObjectVars
  {
    GeomPoint points;
    GeomScalar scalars;
    GeomVector vectors;
    GeomTensor tensors;
  };

  struct ParticleVars
  {
    Uintah::ParticleVariable<Uintah::Point> position;
    Uintah::ParticleVariable<Uintah::Vector> pVelocity, pExternalForce;
    Uintah::ParticleVariable<Uintah::Matrix3> pSize;
    Uintah::ParticleVariable<double> pMass, pVolume;
    Uintah::ParticleVariable<Uintah::long64> pParticleID;
    Uintah::ParticleVariable<Uintah::Vector> pDisplacement;
    Uintah::ParticleVariable<double> pHorizon;
    Uintah::ParticleVariable<int> pLoadCurveID;
    Uintah::ParticleVariable<IntVector> pLoadCurveIDVec;
  };

protected:
  virtual Uintah::ParticleSubset*
  allocateVariables(particleIndex numParticles,
                    int dwi,
                    const Uintah::Patch* patch,
                    Uintah::DataWarehouse* new_dw,
                    ParticleVars& pvars);

  virtual particleIndex
  countAndCreateParticles(const Uintah::Patch*,
                          Uintah::GeometryObject* obj,
                          ObjectVars& vars);

  void
  createPoints(const Uintah::Patch* patch,
               Uintah::GeometryObject* obj,
               ObjectVars& vars);

  virtual void
  initializeParticle(const Uintah::Patch* patch,
                     Uintah::GeometryObject* obj,
                     PeridynamicsMaterial* matl,
                     Uintah::Point p,
                     Uintah::IntVector cell_idx,
                     particleIndex i,
                     Uintah::CCVariable<short int>& cellNAPI,
                     ParticleVars& pvars);

  int
  checkForSurface(const Uintah::GeometryPieceP piece,
                  const Uintah::Point p,
                  const Uintah::Vector dxpp);

  int
  getLoadCurveID(const Uintah::Point& pp, const Uintah::Vector& dxpp);

protected:
  std::unique_ptr<PeridynamicsLabel> d_varLabel;
  const PeridynamicsFlags* d_flags;

  std::vector<const Uintah::VarLabel*> particle_state, particle_state_preReloc;
};

} // End of namespace Vaango

#endif // __VAANGO_CCA_COMPONENTS_PERIDYNAMICS_PARTICLE_CREATOR_H__
