#ifndef DEM_PARTICLE_CREATOR_H
#define DEM_PARTICLE_CREATOR_H

#include <DiscreteElements/DEMContainers.h>
#include <DiscreteElements/DEMParticle.h>
#include <Core/Geometry/Box.h>
#include <Core/Geometry/OrientedBox.h>
#include <Core/Geometry/FaceEdge.h>

namespace dem {

class DEMParticleCreator
{
public:

  struct ParticleParameters {
    int numLayers;
    REAL youngModulus;
    REAL poissonRatio;
    REAL maxDiameter;
    REAL edge;
    REAL offset;
    bool randomOrientation;
    bool randomRadiusRatio;
    bool randomVelocity;

    friend
    std::ostream& operator<<(std::ostream& os, const ParticleParameters& p) {
      os << "layers = " << p.numLayers
         << " ym = " << p.youngModulus
         << " pr = " << p.poissonRatio
         << " dia = " << p.maxDiameter
         << " edge = " << p.edge
         << " offset = " << p.offset
         << " orient = " << std::boolalpha << p.randomOrientation
         << " ratio = " << std::boolalpha << p.randomRadiusRatio
         << " vel = " << std::boolalpha << p.randomVelocity << "\n";
      return os;
    }
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

  DEMParticlePArray updatePeriodicDEMParticles(const OrientedBox& periodicDomain, 
                                               DEMParticlePArray& particles);

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
  void addExtraTranslations(const Vec& shift, 
                            const Face& face,
                            const IntersectionStatus status, 
                            std::vector<Vec>& translations);

  void removeDuplicates(DEMParticlePArray& input);
};

} // End namespace dem

#endif
