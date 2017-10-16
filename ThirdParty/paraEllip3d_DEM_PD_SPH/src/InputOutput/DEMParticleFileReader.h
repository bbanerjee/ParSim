#ifndef ELLIP3D_PARTICLE_FILE_READER_H
#define ELLIP3D_PARTICLE_FILE_READER_H

#include <Core/Types/RealTypes.h>
#include <DiscreteElements/DEMContainers.h>
#include <DiscreteElements/Gradation.h>
#include <InputOutput/zenxml/xml.h>
#include <vector>

namespace dem {

class DEMParticleFileReader
{

public:
  DEMParticleFileReader() = default;
  ~DEMParticleFileReader() = default;

  void read(const std::string& fileName, const REAL& youngModulus,
            const REAL& poissonRatio, bool doInitialize,
            DEMParticlePArray& particles, Gradation& gradation);

private:
  REAL d_youngModulus;
  REAL d_poissonRatio;
  bool d_doInitialize;

  void readParticlesText(const std::string& inputParticle,
                         DEMParticlePArray& particles, Gradation& gradation) const;

  bool readParticlesXML(const std::string& inputFileName,
                        DEMParticlePArray& particles, Gradation& gradation) const;

  template <typename T>
  bool readParticleValues(zen::XmlIn& ps, const std::string& name,
                          const std::string& particleType,
                          std::vector<T>& output) const;

  void updateTotalMass(DEMParticlePArray& particles);

  DEMParticleFileReader(DEMParticleFileReader const&) = delete; // don't implement
  void operator=(DEMParticleFileReader const&) = delete;     // don't implement
};
}
#endif
