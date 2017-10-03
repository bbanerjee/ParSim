#ifndef DEM_PARTICLE_CREATOR_H
#define DEM_PARTICLE_CREATOR_H

#include <DiscreteElements/DEMContainers.h>
#include <DiscreteElements/DEMParticle.h>
#include <Core/Geometry/Box.h>

namespace dem {

class DEMParticleCreator
{
public:

  DEMParticleCreator() = default;
  ~DEMParticleCreator() = default;

  template <std::size_t layerFlag>
  DEMParticlePArray generateDEMParticles(const dem::Box& allContainer,
                                         dem::Gradation& gradation);

private:

  struct ParticleParameters {
    REAL youngModulus;
    REAL poissonRatio;
    REAL maxDiameter;
    REAL edge;
    REAL offset;
  };

  ParticleParameters getParticleParameters(const dem::Gradation& gradation);

  template <std::size_t layerFlag>
  DEMParticlePArray createParticles(const dem::Box& allContainer, 
    dem::Gradation& gradation, const ParticleParameters& params);

};

} // End namespace dem

#endif