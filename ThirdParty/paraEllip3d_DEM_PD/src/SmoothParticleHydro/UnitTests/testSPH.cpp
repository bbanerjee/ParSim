#include <SmoothParticleHydro/SmoothParticleHydro.h>
#include <Core/Util/Utility.h>
#include <InputOutput/InputParameter.h>
#include <gtest/gtest.h>

using namespace sph;

TEST(SPHTest, computeMass) {

  SmoothParticleHydro sph;

  REAL density = 1000;
  REAL length = 2.5;
  std::size_t numPts = 2500;

  REAL mass = sph.computeMass<2>(density, length, numPts);
  EXPECT_DOUBLE_EQ(mass, 0.001);

  mass = sph.computeMass<3>(density, length, numPts);
  EXPECT_DOUBLE_EQ(mass, 1.0e-6);
}

TEST(SPHTest, createCoords) {

  SmoothParticleHydro sph;

  dem::Vec vmin(0.4, 0.4, 0.4);
  dem::Vec vmax(1, 2, 3);
  REAL spaceInterval = 0.2;
  int numLayers = 2;

  // Create an linearly spaced array of x/y/zcoords from x/y/zminBuffered to x/y/zmaxBuffered
  std::vector<REAL> xCoords, yCoords, zCoords;
  sph.createCoords<2>(vmin, vmax, spaceInterval, numLayers, xCoords, yCoords, zCoords);

  EXPECT_EQ(xCoords.size(), 8);
  EXPECT_EQ(yCoords.size(), 0);
  EXPECT_EQ(zCoords.size(), 18);

  EXPECT_DOUBLE_EQ(xCoords[0], 0);
  EXPECT_DOUBLE_EQ(xCoords[7], 1.4);
  EXPECT_DOUBLE_EQ(zCoords[0], 0);
  EXPECT_DOUBLE_EQ(zCoords[17], 3.4);

  spaceInterval = 0.23;
  numLayers = 2;
  xCoords.clear();
  yCoords.clear();
  zCoords.clear();
  sph.createCoords<3>(vmin, vmax, spaceInterval, numLayers, xCoords, yCoords, zCoords);

  /*
  std::cout << "z: ";
  for (auto z: zCoords) {
    std::cout << z << " ";
  }
  std::cout << std::endl;
  */

  EXPECT_EQ(xCoords.size(), 8);
  EXPECT_EQ(yCoords.size(), 12);
  EXPECT_EQ(zCoords.size(), 16);

  EXPECT_DOUBLE_EQ(xCoords[0], -0.06);
  EXPECT_NEAR(xCoords[6], 1.32, 1.0e-6);
  EXPECT_DOUBLE_EQ(yCoords[0], -0.06);
  EXPECT_NEAR(yCoords[10], 2.24, 1.0e-6);
  EXPECT_DOUBLE_EQ(zCoords[0], -0.06);
  EXPECT_NEAR(zCoords[15], 3.39, 1.0e-6);
}

TEST(SPHTest, createParticleArray) {

  dem::InputParameter::get().addParameter("P0", 1600.0);
  dem::InputParameter::get().addParameter("SPHInitialDensity", 1100.0);
  dem::InputParameter::get().addParameter("gamma", 1.4);
  dem::InputParameter::get().addParameter("nu", 100.0);

  SmoothParticleHydro sph;

  dem::Vec vmin(0, 0, 0);
  dem::Vec vmax(1, 0, 1);
  REAL spaceInterval = 0.2;
  int numLayers = 2;

  std::vector<REAL> xCoords, yCoords, zCoords;
  sph.createCoords<2>(vmin, vmax, spaceInterval, numLayers, xCoords, yCoords, zCoords);

  REAL mass = 0.001;
  REAL density = 1000.0;
  sph.createParticleArray<2>(mass, density, xCoords, yCoords, zCoords);

  SPHParticlePArray particles = sph.getAllSPHParticleVec();
  EXPECT_EQ(particles.size(), 100);
  EXPECT_NEAR(particles[0]->getInitPosition().x(), -0.4, 1.0e-16);
  EXPECT_DOUBLE_EQ(particles[0]->getInitPosition().y(), 0);
  EXPECT_NEAR(particles[0]->getInitPosition().z(), -0.4, 1.0e-16);
  EXPECT_NEAR(particles[99]->getInitPosition().x(), 1.4, 1.0e-16);
  EXPECT_DOUBLE_EQ(particles[99]->getInitPosition().y(), 0);
  EXPECT_NEAR(particles[99]->getInitPosition().z(), 1.4, 1.0e-16);

  vmax = dem::Vec(1, 0.1, 1);
  xCoords.clear();
  yCoords.clear();
  zCoords.clear();
  sph.createCoords<3>(vmin, vmax, spaceInterval, numLayers, xCoords, yCoords, zCoords);

  /*
  std::cout << "x, y, z = " << xCoords.size() << "," << yCoords.size()
            << "," << zCoords.size() << std::endl;
  std::cout << "y: ";
  for (auto y: yCoords) {
    std::cout << y << " ";
  }
  std::cout << std::endl;
  */

  sph.createParticleArray<3>(mass, density, xCoords, yCoords, zCoords);
  particles = sph.getAllSPHParticleVec();
  std::cout << "size = " << particles.size() << std::endl;
  EXPECT_EQ(particles.size(), 600);
  EXPECT_NEAR(particles[0]->getInitPosition().x(), -0.4, 1.0e-16);
  EXPECT_NEAR(particles[0]->getInitPosition().y(), -0.4, 1.0e-16);
  EXPECT_NEAR(particles[0]->getInitPosition().z(), -0.4, 1.0e-16);
  EXPECT_NEAR(particles[599]->getInitPosition().x(), 1.4, 1.0e-16);
  EXPECT_NEAR(particles[599]->getInitPosition().y(), 0.6, 1.0e-16);
  EXPECT_NEAR(particles[599]->getInitPosition().z(), 1.4, 1.0e-16);
}
