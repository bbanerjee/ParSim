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

    virtual particleIndex countParticles(const Uintah::Patch*,
                                         std::vector<Uintah::GeometryObject*>&);

    virtual particleIndex countAndCreateParticles(const Uintah::Patch*,
                                                  Uintah::GeometryObject* obj);

    std::vector<const Uintah::VarLabel* > returnParticleState();
    std::vector<const Uintah::VarLabel* > returnParticleStatePreReloc();

  protected:

    void createPoints(const Uintah::Patch* patch, Uintah::GeometryObject* obj);

    virtual void initializeParticle(const Uintah::Patch* patch,
                                    std::vector<Uintah::GeometryObject*>::const_iterator obj,
                                    PeridynamicsMaterial* matl,
                                    Uintah::Point p, Uintah::IntVector cell_idx,
                                    particleIndex i,
                                    Uintah::CCVariable<short int>& cellNAPI);
    
  protected:

    Uintah::ParticleVariable<Uintah::Point> d_position;
    Uintah::ParticleVariable<Uintah::Vector> d_pvelocity, d_pexternalforce;
    Uintah::ParticleVariable<Uintah::Matrix3> d_psize;
    Uintah::ParticleVariable<double> d_pmass, d_pvolume;
    Uintah::ParticleVariable<Uintah::long64> d_pparticleID;
    Uintah::ParticleVariable<Uintah::Vector> d_pdisp;

    Uintah::ParticleVariable<double> d_pHorizon;

    PeridynamicsLabel* d_varLabel;
    PeridynamicsFlags* d_flags;

    std::vector<const Uintah::VarLabel* > particle_state, particle_state_preReloc;

    typedef std::map<std::pair<const Uintah::Patch*,Uintah::GeometryObject*>,std::vector<Uintah::Point> > GeometryPoints;
    typedef std::map<std::pair<const Uintah::Patch*,Uintah::GeometryObject*>,std::vector<double> > GeometryVolumes;
    typedef std::map<std::pair<const Uintah::Patch*,Uintah::GeometryObject*>,std::vector<Uintah::Vector> > GeometryVectors;

    GeometryPoints d_object_points;
    GeometryVolumes d_object_vols;
    GeometryVectors d_object_velocity; // gcd add
    
    mutable Uintah::CrowdMonitor   d_lock;
  };

} // End of namespace Vaango

#endif // __VAANGO_PARTICLE_CREATOR_H__
