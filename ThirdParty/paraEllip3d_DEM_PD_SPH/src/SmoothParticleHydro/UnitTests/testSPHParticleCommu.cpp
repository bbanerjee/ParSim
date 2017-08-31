#include <SmoothParticleHydro/SmoothParticleHydro.h>
#include <SmoothParticleHydro/SPHParticleCreator.h>
#include <Core/Util/Utility.h>

#include <gtest/gtest.h>
#include <boost/mpi.hpp>

using namespace sph;

class MPIEnvironment : public ::testing::Environment
{
public:
  // Override this to define how to set up the environment.
  virtual void SetUp() {
    char** argv;
    int argc = 0;
    int mpiError = MPI_Init(&argc, &argv);
    ASSERT_FALSE(mpiError);
  }
  // Override this to define how to tear down the environment.
  virtual void TearDown() {
    int mpiError = MPI_Finalize();
    ASSERT_FALSE(mpiError);
  }

  virtual ~MPIEnvironment() {}

};

int main(int argc, char* argv[]) 
{
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new MPIEnvironment);
  return RUN_ALL_TESTS();
}

TEST(SPHParticleCommuTest, communicate)
{
  // Set up communicator
  boost::mpi::communicator boostWorld;
  dem::InputParameter::get().addParameter("mpiProcX", 2);
  dem::InputParameter::get().addParameter("mpiProcY", 1);
  dem::InputParameter::get().addParameter("mpiProcZ", 1);

  // Setup the parameters
  dem::InputParameter::get().addParameter("specificG", 1.1);
  dem::InputParameter::get().addParameter("gravAccel", 10);
  dem::InputParameter::get().addParameter("gravScale", 1);
  dem::InputParameter::get().addParameter("spaceInterval", 0.05);
  dem::InputParameter::get().addParameter("numLayers", 1);
  dem::InputParameter::get().addParameter("gamma", 1.4);
  dem::InputParameter::get().addParameter("P0", 1.0e5);
  dem::InputParameter::get().addParameter("SPHInitialDensity", 1000);
  dem::InputParameter::get().addParameter("nu", 1);


  SPHParticlePArray sph_particles;

  if (boostWorld.rank() == 0) {
    dem::DEMParticlePArray allDEMParticles;

    SPHParticleCreator creator;
    dem::Box container(dem::Vec(0, 0, 0), dem::Vec(0.1, 0.1, 0.1));
    sph_particles = creator.generateSPHParticleNoBottom<2>(container, allDEMParticles);

    EXPECT_EQ(sph_particles.size(), 25);
  }


  // Set up SPH domain
  dem::Box allContainer(dem::Vec(-0.01, -0.01, -0.01), dem::Vec(0.11, 0.11, 0.11));
  REAL spaceInterval = util::getParam<REAL>("spaceInterval");
  int numLayers = util::getParam<int>("numLayers");
  REAL ghostWidth = spaceInterval;
  REAL bufferLength = spaceInterval*numLayers;

  // Scatter the set of particles
  SmoothParticleHydro sph;
  sph.setCommunicator(boostWorld);
  if (boostWorld.rank() == 0) {
    sph.setAllSPHParticleVec(sph_particles);
    EXPECT_EQ(sph.getAllSPHParticleVec().size(), 25);
  }
  sph.scatterSPHParticle(allContainer, ghostWidth, bufferLength);

  // Communicate ghost regions
  int iteration = 0;
  sph.commuSPHParticle(iteration, ghostWidth);
  if (boostWorld.rank() == 0) {
    EXPECT_EQ(sph.getMergedSPHParticleVec().size(), 15);
  } else {
    EXPECT_EQ(sph.getMergedSPHParticleVec().size(), 20);
  }

}