#ifndef ELLIP3D_PERI_PARTICLE_FILE_READER_H
#define ELLIP3D_PERI_PARTICLE_FILE_READER_H

#include <Core/Types/realtypes.h>
#include <Peridynamics/PeriContainers.h>
#include <InputOutput/zenxml/xml.h>
#include <vector>

namespace pd {

class PeriParticleFileReader
{

public:
  PeriParticleFileReader() = default;
  ~PeriParticleFileReader() = default;

  void read(const std::string& fileName,
            PeriParticlePArray& particles,
            PeriElements& connectivity) const;

private:

  void readPeriParticlesText(const std::string& inputFileName,
                             PeriParticlePArray& particles,
                             PeriElements& connectivity) const;

  bool readPeriParticlesXML(const std::string& inputFileName,
                            PeriParticlePArray& particles,
                            PeriElements& connectivity) const;

  template <typename T>
  bool readPeriParticleValues(zen::XmlIn& ps, const std::string& name,
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

  PeriParticleFileReader(PeriParticleFileReader const&) = delete; // don't implement
  void operator=(PeriParticleFileReader const&) = delete;     // don't implement
};
}
#endif
