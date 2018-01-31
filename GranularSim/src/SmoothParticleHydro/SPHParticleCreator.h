#ifndef SPH_PARTICLE_CREATOR_H
#define SPH_PARTICLE_CREATOR_H

#include <SmoothParticleHydro/SPHContainers.h>
#include <SmoothParticleHydro/SPHParticle.h>
#include <Core/Geometry/Box.h>

namespace sph {

class SPHParticleCreator
{
public:
  SPHParticleCreator() = default;
  ~SPHParticleCreator() = default;

  template <int dim>
  SPHParticlePArray generateSPHParticleDam(const dem::Box& spatialDomain,
                           const dem::DEMParticlePArray& allDEMParticles);

  template <int dim>
  SPHParticlePArray generateSPHParticleNoBottom(const dem::Box& spatialDomain,
                                   const dem::DEMParticlePArray& allDEMParticles);

  template <int dim>
  SPHParticlePArray generateSPHParticleMiddleLayers(const dem::Box& spatialDomain,
                                   const dem::DEMParticlePArray& allDEMParticles);

  template <int dim>
  REAL computeMass(const double& density, const double& length,
                   const std::size_t& numPts) const;

  template <int dim>
  void createCoords(const dem::Vec& vmin, const dem::Vec& vmax,
                    const REAL& spaceInterval, const int& numLayers,
                    std::vector<REAL>& xCoords, std::vector<REAL>& yCoords,
                    std::vector<REAL>& zCoords) const;

  template <int dim>
  SPHParticlePArray createParticleArray(const REAL& mass, const REAL& density,
                           const std::vector<REAL>& xCoords,
                           const std::vector<REAL>& yCoords,
                           const std::vector<REAL>& zCoords);

  void removeRedundantSPHParticles(SPHParticlePArray& allSPHParticles);
};

} // End namespace sph

#endif