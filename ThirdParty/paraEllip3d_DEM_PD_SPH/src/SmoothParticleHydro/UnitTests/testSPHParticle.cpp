#include <Core/Util/Utility.h>
#include <DiscreteElements/DEMParticle.h>
#include <InputOutput/InputParameter.h>
#include <SmoothParticleHydro/SPHParticle.h>
#include <gtest/gtest.h>

using namespace sph;

TEST(SPHParticleTest, construction0)
{

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

TEST(SPHParticleTest, construction1)
{

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

TEST(SPHParticleTest, insideDEMParticle)
{

  SPHParticle particle;
  particle.setInitialPos(dem::Vec(0, 0, 0));

  dem::DEMParticleP dem_particle = std::make_shared<dem::DEMParticle>();

  REAL a = 1;
  REAL b = 1;
  REAL c = 1;

  dem::Vec adir(0, dem::Pi / 2, dem::Pi / 2);
  dem::Vec bdir(dem::Pi / 2, 0, dem::Pi / 2);
  dem::Vec cdir(dem::Pi / 2, dem::Pi / 2, 0);

  dem::Vec pos(0, 0, 0);

  dem_particle->setA(a);
  dem_particle->setB(b);
  dem_particle->setC(c);

  dem_particle->setCurrPos(pos);

  dem_particle->setCurrDirecA(adir);
  dem_particle->setCurrDirecB(bdir);
  dem_particle->setCurrDirecC(cdir);

  REAL buffer = 0;
  dem::Vec localCoord = 0;
  bool inside = false;
  bool insideGhostLayer = false;

  buffer = -a;
  particle.setInitialPos(dem::Vec(0.5 * a, 0, 0));
  inside = particle.isInsideDEMParticle<2>(buffer, dem_particle, localCoord,
                                           insideGhostLayer);
  EXPECT_EQ(inside, true);
  EXPECT_EQ(insideGhostLayer, false);

  buffer = 0;
  particle.setInitialPos(dem::Vec(0, a, 0));
  inside = particle.isInsideDEMParticle<3>(buffer, dem_particle, localCoord,
                                           insideGhostLayer);
  EXPECT_EQ(inside, true);
  EXPECT_EQ(insideGhostLayer, false);

  buffer = a;
  particle.setInitialPos(dem::Vec(0, 0, 2.0 * a));
  inside = particle.isInsideDEMParticle<2>(buffer, dem_particle, localCoord,
                                           insideGhostLayer);
  EXPECT_EQ(inside, false);
  EXPECT_EQ(insideGhostLayer, false);

  buffer = 2 * a;
  particle.setInitialPos(dem::Vec(0, 0, 2.0 * a));
  inside = particle.isInsideDEMParticle<3>(buffer, dem_particle, localCoord,
                                           insideGhostLayer);
  EXPECT_EQ(inside, false);
  EXPECT_EQ(insideGhostLayer, false);

  buffer = 0.2 * a;
  particle.setInitialPos(dem::Vec(0.9 * a, 0, 0));
  inside = particle.isInsideDEMParticle<2>(buffer, dem_particle, localCoord,
                                           insideGhostLayer);
  EXPECT_EQ(inside, true);
  EXPECT_EQ(insideGhostLayer, true);

  particle.setInitialPos(dem::Vec(0, 0.9 * a, 0));
  inside = particle.isInsideDEMParticle<3>(buffer, dem_particle, localCoord,
                                           insideGhostLayer);
  EXPECT_EQ(inside, true);
  EXPECT_EQ(insideGhostLayer, true);

  particle.setInitialPos(dem::Vec(0, 0, 0.9 * a));
  inside = particle.isInsideDEMParticle<2>(buffer, dem_particle, localCoord,
                                           insideGhostLayer);
  EXPECT_EQ(inside, true);
  EXPECT_EQ(insideGhostLayer, true);

  buffer = 0.3 * a;
  particle.setInitialPos(dem::Vec(0.81 * a * cos(dem::Pi / 3),
                                  0.81 * a * cos(dem::Pi / 3),
                                  0.9 * a * cos(dem::Pi / 3)));
  inside = particle.isInsideDEMParticle<3>(buffer, dem_particle, localCoord,
                                           insideGhostLayer);
  EXPECT_EQ(inside, true);
  EXPECT_EQ(insideGhostLayer, true);

  buffer = 0.3001;
  particle.setInitialPos(dem::Vec(0.7, 0.7, 0.7));
  inside = particle.isInsideDEMParticle<3>(buffer, dem_particle, localCoord,
                                           insideGhostLayer);
  EXPECT_EQ(inside, false);
  EXPECT_EQ(insideGhostLayer, false);

  /*
  std::cout << "Pt: " << point
            << "Particle center: " << pos
            << "LocalCoord = " << localCoord
            << "inside = " << std::boolalpha << inside << std::endl;
  */
}
