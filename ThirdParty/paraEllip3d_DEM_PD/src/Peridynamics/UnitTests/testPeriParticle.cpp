#include <Peridynamics/PeriParticle.h>
#include <gtest/gtest.h>

using namespace pd;

TEST(PeriParticleTest, accFunctions) {

  PeriParticle particle;
  EXPECT_EQ(particle.getId(), 0);
}