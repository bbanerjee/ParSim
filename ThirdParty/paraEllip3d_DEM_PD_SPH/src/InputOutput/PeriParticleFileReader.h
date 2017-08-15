#ifndef ELLIP3D_PERI_PARTICLE_FILE_READER_H
#define ELLIP3D_PERI_PARTICLE_FILE_READER_H

#include <Core/Types/realtypes.h>
#include <Core/Types/integertypes.h>
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
            PeriElementArray& connectivity) const;

private:

  bool checkAbaqusFileFormat(const std::string& fileName) const;

  void readPeriParticlesText(const std::string& inputFileName,
                             PeriParticlePArray& particles,
                             PeriElementArray& connectivity) const;

  bool readPeriParticlesAbaqus(const std::string& inputFileName,
                               PeriParticlePArray& particles,
                               PeriElementArray& connectivity) const;

  struct VolumeElement
  {
    VolumeElement(ElementID id, std::vector<ParticleID> nodes) {
      id_ = id;
      node1_ = nodes[0]; node2_ = nodes[1];
      node3_ = nodes[2]; node4_ = nodes[3];
    }

    ElementID id_;
    ParticleID node1_;
    ParticleID node2_;
    ParticleID node3_;
    ParticleID node4_;
  };

  struct MeshNode
  {
    MeshNode(ParticleID id, double x, double y, double z) {
      id_ = id; x_ = x; y_ = y; z_ = z; 
    }

    ParticleID id_;
    double x_;
    double y_;
    double z_;
  };

  void readAbaqusMeshNode(const std::string& inputLine,
                          std::vector<MeshNode>& nodes) const;
  void readAbaqusMeshVolumeElement(const std::string& inputLine,
                                   std::vector<VolumeElement>& elements) const;

  PeriParticleFileReader(PeriParticleFileReader const&) = delete; // don't implement
  void operator=(PeriParticleFileReader const&) = delete;     // don't implement

};
}
#endif
