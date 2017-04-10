#ifndef ELLIP3D_PARTICLE_FILE_READER_H
#define ELLIP3D_PARTICLE_FILE_READER_H

#include <Core/Types/realtypes.h>
#include <DiscreteElements/Containers.h>
#include <DiscreteElements/Gradation.h>
#include <InputOutput/zenxml/xml.h>
#include <vector>

namespace dem {

class ParticleFileReader
{

public:
  ParticleFileReader() = default;
  ~ParticleFileReader() = default;

  void read(const std::string& fileName, const REAL& youngModulus,
            const REAL& poissonRatio, bool doInitialize,
            ParticlePArray& particles, Gradation& gradation);

private:
  REAL d_youngModulus;
  REAL d_poissonRatio;
  bool d_doInitialize;

  void readParticlesText(const std::string& inputParticle,
                         ParticlePArray& particles, Gradation& gradation) const;

  bool readParticlesXML(const std::string& inputFileName,
                        ParticlePArray& particles, Gradation& gradation) const;

  template <typename T>
  bool readParticleValues(zen::XmlIn& ps, const std::string& name,
                          const std::string& particleType,
                          std::vector<T>& output) const;

  template <typename T>
  bool decodeAndUncompress(const std::string& inputStr,
                           const int& numComponents,
                           std::vector<T>& output) const;

  template <typename T>
  T convert(const std::string& str) const;

  template <typename T>
  std::vector<T> convertStrArray(const std::string& str) const;

  ParticleFileReader(ParticleFileReader const&) = delete; // don't implement
  void operator=(ParticleFileReader const&) = delete;     // don't implement
};
}
#endif
