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

#ifndef __VAANGO_PARTICLE_CREATOR_H__
#define __VAANGO_PARTICLE_CREATOR_H__

#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Task.h>

#include <Core/Thread/CrowdMonitor.h>

#include <vector>
#include <map>

namespace Uintah {
  class GeometryObject;
  class Patch;
  class DataWarehouse;
  class ParticleSubset;
  class VarLabel;
}

namespace Vaango {

  typedef int particleIndex;
  typedef int particleId;

  class PeridynamicsFlags;
  class PeridynamicsMaterial;
  class PeridynamicsLabel;

  class ParticleCreator {

  public:
    
    ParticleCreator(PeridynamicsMaterial* matl, PeridynamicsFlags* flags);
    virtual ~ParticleCreator();

    virtual particleIndex countParticles(const Uintah::Patch*,
                                         std::vector<Uintah::GeometryObject*>&);

    virtual Uintah::ParticleSubset* createParticles(PeridynamicsMaterial* matl,
                                                    particleIndex numParticles,
                                                    Uintah::CCVariable<short int>& cellNAPID,
                                                    const Uintah::Patch*,
                                                    Uintah::DataWarehouse* new_dw,
                                                    std::vector<Uintah::GeometryObject*>&);

    virtual Uintah::ParticleSubset* allocateVariables(particleIndex numParticles,
                                                      int dwi, 
                                                      const Uintah::Patch* patch,
                                                      Uintah::DataWarehouse* new_dw);

    virtual void allocateVariablesAddRequires(Uintah::Task* task, 
                                              const PeridynamicsMaterial* matl,
                                              const Uintah::PatchSet* patch) const;

    virtual void registerPermanentParticleState(PeridynamicsMaterial* matl);

    std::vector<const Uintah::VarLabel* > returnParticleState();
    std::vector<const Uintah::VarLabel* > returnParticleStatePreReloc();

  protected:

    particleIndex countAndCreateParticles(const Uintah::Patch*,
                                          Uintah::GeometryObject* obj);

    void createPoints(const Uintah::Patch* patch, Uintah::GeometryObject* obj);

    virtual void initializeParticle(const Uintah::Patch* patch,
                                    std::vector<Uintah::GeometryObject*>::const_iterator obj,
                                    PeridynamicsMaterial* matl,
                                    Uintah::Point p, Uintah::IntVector cell_idx,
                                    particleIndex i,
                                    Uintah::CCVariable<short int>& cellNAPI);
    
   int checkForSurface(const Uintah::GeometryPieceP piece, const Uintah::Point p, const Uintah::Vector dxpp);

   int getLoadCurveID(const Uintah::Point& pp, const Uintah::Vector& dxpp);

  protected:

    Uintah::ParticleVariable<Uintah::Point> d_position;
    Uintah::ParticleVariable<Uintah::Vector> d_pvelocity, d_pexternalforce;
    Uintah::ParticleVariable<Uintah::Matrix3> d_psize;
    Uintah::ParticleVariable<double> d_pmass, d_pvolume;
    Uintah::ParticleVariable<Uintah::long64> d_pparticleID;
    Uintah::ParticleVariable<Uintah::Vector> d_pdisp;

    Uintah::ParticleVariable<double> d_pHorizon;
    Uintah::ParticleVariable<int> d_pLoadCurveID;

    PeridynamicsLabel* d_varLabel;
    PeridynamicsFlags* d_flags;

    std::vector<const Uintah::VarLabel* > particle_state, particle_state_preReloc;

    typedef std::vector<Uintah::Point> PointArray;
    typedef std::vector<double> DoubleArray;
    typedef std::vector<Uintah::Vector> VectorArray;

    typedef std::pair<const Uintah::Patch*, Uintah::GeometryObject*> PatchGeometryObjectPair;

    typedef std::map<PatchGeometryObjectPair, PointArray > GeometryPoints;
    typedef std::map<PatchGeometryObjectPair, DoubleArray > GeometryScalars;
    typedef std::map<PatchGeometryObjectPair, VectorArray > GeometryVectors;

    GeometryPoints d_object_points;
    GeometryScalars d_object_vols;
    GeometryVectors d_object_velocity; 
    GeometryVectors d_object_forces; 
    
    mutable Uintah::CrowdMonitor   d_lock;
  };

} // End of namespace Vaango

#endif // __VAANGO_PARTICLE_CREATOR_H__
