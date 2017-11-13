#ifndef DEM_PARTICLE_CREATOR_H
#define DEM_PARTICLE_CREATOR_H

#include <DiscreteElements/DEMContainers.h>
#include <DiscreteElements/DEMParticle.h>
#include <Core/Geometry/Box.h>
#include <Core/Geometry/FaceEdge.h>

namespace dem {

class DEMParticleCreator
{
public:

  struct ParticleParameters {
    REAL youngModulus;
    REAL poissonRatio;
    REAL maxDiameter;
    REAL edge;
    REAL offset;
  };

  DEMParticleCreator() = default;
  ~DEMParticleCreator() = default;

  DEMParticlePArray generateDEMParticles(std::size_t layerFlag,
                                         DEMParticle::DEMParticleShape,
                                         const dem::Box& spatialDomain,
                                         dem::Gradation& gradation);

  DEMParticlePArray generatePeriodicDEMParticles(DEMParticlePArray& parts,
                                                 dem::Box& spatialDomain,
                                                 REAL marginFactor = 2,
                                                 REAL faceShiftFactor = 0);

private:

  ParticleParameters getParticleParameters(const dem::Gradation& gradation);

  template <std::size_t layerFlag>
  DEMParticlePArray generateDEMParticles(DEMParticle::DEMParticleShape type,
                                         const dem::Box& spatialDomain,
                                         dem::Gradation& gradation);

  template <std::size_t layerFlag>
  DEMParticlePArray createParticles(DEMParticle::DEMParticleShape type,
                                    const dem::Box& spatialDomain, 
                                    dem::Gradation& gradation, 
                                    const ParticleParameters& params);

  using IntersectionStatus = std::pair<bool, std::pair<Face::Location, int>>;
  void addExtraTranslations(const Vec& shift, 
                            REAL boundaryMargin,
                            Vec inPlaneDiag,
                            REAL widthX, REAL widthY, REAL widthZ,
                            const Face& face,
                            const IntersectionStatus status, 
                            std::vector<Vec>& translations);

  void removeDuplicates(DEMParticlePArray& input);
};

} // End namespace dem

#endif
