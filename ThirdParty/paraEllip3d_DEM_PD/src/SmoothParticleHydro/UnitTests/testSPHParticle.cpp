#include <SmoothParticleHydro/SPHParticle.h>
#include <Core/Util/Utility.h>
#include <InputOutput/InputParameter.h>
#include <gtest/gtest.h>

using namespace sph;

TEST(SPHParticleTest, construction0) {

  SPHParticle particle;
  EXPECT_EQ(particle.getId(), 0);
  EXPECT_DOUBLE_EQ(particle.getMass(), 0);
  EXPECT_DOUBLE_EQ(particle.getDensity(), 0);
  EXPECT_DOUBLE_EQ(particle.getVolume(), 0);
  EXPECT_DOUBLE_EQ(particle.getPressure(), 0);
  EXPECT_DOUBLE_EQ(particle.getViscosity(), 0);
  EXPECT_DOUBLE_EQ(particle.getVelocity().x(), 0);
  EXPECT_DOUBLE_EQ(particle.getVelocity().y(), 0);
  EXPECT_DOUBLE_EQ(particle.getVelocity().z(), 0);
  EXPECT_EQ(particle.getType(), SPHParticleType::NONE);

}

TEST(SPHParticleTest, construction1) {

  // Setup the parameters that are used by the constructor
  dem::InputParameter::get().addParameter("P0", 1600.0);
  dem::InputParameter::get().addParameter("SPHInitialDensity", 1100.0);
  dem::InputParameter::get().addParameter("gamma", 1.4);
  dem::InputParameter::get().addParameter("nu", 100.0);
  ParticleID id = 100;
  REAL mass = 1.5;
  REAL density = 1200.0;
  REAL x = 1.0;
  REAL y = -1.2;
  REAL z = 5.2;
  dem::Vec local = 0;
  SPHParticleType type = SPHParticleType::FREE;
  SPHParticle particle(id, mass, density, x, y, z, local, type);
  EXPECT_EQ(particle.getId(), 100);
  EXPECT_DOUBLE_EQ(particle.getMass(), 1.5);
  EXPECT_DOUBLE_EQ(particle.getDensity(), 1200.0);
  EXPECT_DOUBLE_EQ(particle.getVolume(), 0.00125);
  EXPECT_NEAR(particle.getPressure(), 207.273863, 1.0e-6);
  EXPECT_DOUBLE_EQ(particle.getViscosity(), 120000);
  EXPECT_DOUBLE_EQ(particle.getVelocity().x(), 0);
  EXPECT_DOUBLE_EQ(particle.getVelocity().y(), 0);
  EXPECT_DOUBLE_EQ(particle.getVelocity().z(), 0);
  EXPECT_EQ(particle.getDEMParticle(), nullptr);
  EXPECT_EQ(particle.getType(), SPHParticleType::FREE);

}