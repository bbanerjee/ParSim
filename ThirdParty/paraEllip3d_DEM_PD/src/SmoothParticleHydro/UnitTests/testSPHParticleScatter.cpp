#include <SmoothParticleHydro/SmoothParticleHydro.h>

#include <gtest/gtest.h>
#include <boost/mpi.hpp>

using namespace sph;

int main(int argc, char* argv[]) 
{
  ::testing::InitGoogleTest(&argc, argv);
  boost::mpi::environment(argc, argv);
  return RUN_ALL_TESTS();
}

TEST(SPHParticleScatterTest, scatter)
{
  boost::mpi::communicator boostworld;
  SmoothParticleHydro sph;
}
