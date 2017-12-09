#include <InputOutput/DEMParticleFileReader.h>
#include <DiscreteElements/DEMParticle.h>
#include <Core/Util/Utility.h>
#include <gtest/gtest.h>

using namespace dem;

TEST(DEMParticleFileReaderTest, readVTK) {

  // Some parameters needed by DEMParticle
  dem::InputParameter::get().addParameter("young", 1.0e9);
  dem::InputParameter::get().addParameter("poisson", 0.3);
  dem::InputParameter::get().addParameter("specificG", 1.5);
  dem::InputParameter::get().addParameter("gravAccel", 9.81);
  dem::InputParameter::get().addParameter("gravScale", 1);

  REAL youngModulus = util::getParam<REAL>("young");
  REAL poissonRatio = util::getParam<REAL>("poisson");

  // Read particles
  DEMParticleFileReader reader;
  DEMParticlePArray readParticles;

  reader.readVTK("particle_00000.vtu", youngModulus, poissonRatio, readParticles);
  for (auto& particle : readParticles) {
    std::cout << *particle <<"\n";
  }
}

