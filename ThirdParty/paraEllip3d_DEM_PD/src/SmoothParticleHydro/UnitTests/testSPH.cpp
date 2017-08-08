#include <Core/Util/Utility.h>
#include <DiscreteElements/DEMParticle.h>
#include <InputOutput/InputParameter.h>
#include <SmoothParticleHydro/SmoothParticleHydro.h>
#include <gtest/gtest.h>

using namespace sph;

TEST(SPHTest, computeMass)
{

  SmoothParticleHydro sph;

  REAL density = 1000;
  REAL length = 2.5;
  std::size_t numPts = 2500;

  REAL mass = sph.computeMass<2>(density, length, numPts);
  EXPECT_DOUBLE_EQ(mass, 0.001);

  mass = sph.computeMass<3>(density, length, numPts);
  EXPECT_DOUBLE_EQ(mass, 1.0e-6);
}

TEST(SPHTest, createCoords)
{

  SmoothParticleHydro sph;

  dem::Vec vmin(0.4, 0.4, 0.4);
  dem::Vec vmax(1, 2, 3);
  REAL spaceInterval = 0.2;
  int numLayers = 2;

  // Create an linearly spaced array of x/y/zcoords from x/y/zminBuffered to
  // x/y/zmaxBuffered
  std::vector<REAL> xCoords, yCoords, zCoords;
  sph.createCoords<2>(vmin, vmax, spaceInterval, numLayers, xCoords, yCoords,
                      zCoords);

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
  sph.createCoords<3>(vmin, vmax, spaceInterval, numLayers, xCoords, yCoords,
                      zCoords);

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

TEST(SPHTest, createParticleArray)
{

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
  sph.createCoords<2>(vmin, vmax, spaceInterval, numLayers, xCoords, yCoords,
                      zCoords);

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
  sph.createCoords<3>(vmin, vmax, spaceInterval, numLayers, xCoords, yCoords,
                      zCoords);

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
  EXPECT_EQ(particles.size(), 600);
  EXPECT_NEAR(particles[0]->getInitPosition().x(), -0.4, 1.0e-16);
  EXPECT_NEAR(particles[0]->getInitPosition().y(), -0.4, 1.0e-16);
  EXPECT_NEAR(particles[0]->getInitPosition().z(), -0.4, 1.0e-16);
  EXPECT_NEAR(particles[599]->getInitPosition().x(), 1.4, 1.0e-16);
  EXPECT_NEAR(particles[599]->getInitPosition().y(), 0.6, 1.0e-16);
  EXPECT_NEAR(particles[599]->getInitPosition().z(), 1.4, 1.0e-16);
}

TEST(SPHTest, generateSPHParticle)
{

  // Setup the parameters that are used by the constructor
  dem::InputParameter::get().addParameter("specificG", 1.1);
  dem::InputParameter::get().addParameter("gravAccel", 10);
  dem::InputParameter::get().addParameter("gravScale", 1);
  dem::InputParameter::get().addParameter("waterLength", 1);
  dem::InputParameter::get().addParameter("nSPHPoint", 20);
  dem::InputParameter::get().addParameter("numLayers", 1);
  dem::InputParameter::get().addParameter("gamma", 1.4);
  dem::InputParameter::get().addParameter("P0", 1.0e5);
  dem::InputParameter::get().addParameter("SPHInitialDensity", 1000);

  // Create three DEM particles
  std::vector<ParticleID> dem_id = { 0, 1, 2 };
  std::vector<REAL> dem_radii = { 1.15, 0.5,  0.69, 1.15, 0.5,
                                  0.69, 1.15, 0.5,  0.69 };
  std::vector<REAL> dem_axle_a = { 2.250054,  1.570796, 0.6792575,
                                   2.536687,  1.570796, 0.9658907,
                                   0.8923825, 1.570796, 0.6784139 };
  std::vector<REAL> dem_axle_b = { 1.570796, 4.61936e-07,  1.570796,
                                   1.570797, 1.803533e-06, 1.570798,
                                   1.570799, 2.865908e-06, 1.570795 };
  std::vector<REAL> dem_axle_c = { 2.462335, 1.570796, 2.250054,
                                   2.175704, 1.570796, 2.536685,
                                   2.463179, 1.570796, 0.8923828 };
  std::vector<REAL> dem_pos = { 49.73174, 0,        0.9936725, 50.90096, 0,
                                5.320529, 30.81142, 0,         3.898033 };
  std::vector<REAL> dem_vel = { 0.02008031, 0, 0.05505538,
                                -0.3014865, 0, -0.1031938,
                                -0.1538927, 0, 0.04640176 };
  std::vector<REAL> dem_omega = { 0, 0.02781917, 0,          0, -0.07066198,
                                  0, 0,          -0.2856931, 0 };
  std::vector<REAL> dem_force = { -491684,   0,       -126938, -4413.702, 0,
                                  -10805.48, 22966.5, 0,       6613.044 };
  std::vector<REAL> dem_moment = { 0.1688596,  -118712,   -0.3092984,
                                   -0.1174683, -36733.58, -0.05749536,
                                   0.0414009,  -5030.464, -0.0190764 };

  dem::DEMParticlePArray allDEMParticles;
  for (auto id : dem_id) {
    std::size_t startIndex = static_cast<std::size_t>(3 * id);
    std::size_t endIndex = static_cast<std::size_t>(3 * id + 2);
    dem::Vec radii(dem_radii.cbegin() + startIndex,
                   dem_radii.cbegin() + endIndex);
    dem::Vec axle_a(dem_axle_a.cbegin() + startIndex,
                    dem_axle_a.cbegin() + endIndex);
    dem::Vec axle_b(dem_axle_b.cbegin() + startIndex,
                    dem_axle_b.cbegin() + endIndex);
    dem::Vec axle_c(dem_axle_c.cbegin() + startIndex,
                    dem_axle_c.cbegin() + endIndex);
    dem::Vec pos(dem_pos.cbegin() + startIndex, dem_pos.cbegin() + endIndex);
    dem::Vec vel(dem_vel.cbegin() + startIndex, dem_vel.cbegin() + endIndex);
    dem::Vec omega(dem_omega.cbegin() + startIndex,
                   dem_omega.cbegin() + endIndex);
    dem::Vec force(dem_force.cbegin() + startIndex,
                   dem_force.cbegin() + endIndex);
    dem::Vec moment(dem_moment.cbegin() + startIndex,
                    dem_moment.cbegin() + endIndex);
    // std::cout << "id = " << id << " r = " << radii << " a = " << axle_a
    //          << " b = " << axle_b << " c = " << axle_c << " pos = " << pos
    //          << " vel = " << vel << " omega = " << omega << " force = " <<
    //          force
    //          << " moment = " << moment << std::endl;
    dem::DEMParticleP pt = std::make_shared<dem::DEMParticle>(
      id, 0, radii, pos, axle_a, axle_b, axle_c, 1.0e9, 0.3);
    pt->setPrevVeloc(vel);
    pt->setCurrVeloc(vel);
    pt->setPrevOmega(omega);
    pt->setCurrOmega(omega);
    pt->setForce(force);
    pt->setMoment(moment);
    allDEMParticles.push_back(pt);
  }

  SmoothParticleHydro sph;
  dem::Box container(dem::Vec(30.0, -2.0, 0.95), dem::Vec(52.0, 2.0, 5.5));
  sph.generateSPHParticle<2>(container, allDEMParticles);

  SPHParticlePArray sph_particles = sph.getAllSPHParticleVec();
  // std::cout << "Size:" << sph_particles.size() << "\n";
  std::size_t noneP = 0;
  std::size_t ghostP = 0;
  std::size_t freeP = 0;
  std::size_t boundaryP = 0;
  for (auto pt : sph_particles) {
    switch (pt->getType()) {
      case SPHParticleType::GHOST:
        ghostP++;
        // std::cout << "pt " << pt->getId() << " has type GHOST" << std::endl;
        break;
      case SPHParticleType::FREE:
        freeP++;
        // std::cout << "pt " << pt->getId() << " has type FREE" << std::endl;
        break;
      case SPHParticleType::BOUNDARY:
        boundaryP++;
        // std::cout << "pt " << pt->getId() << " has type BOUNDARY" <<
        // std::endl;
        break;
      default:
        noneP++;
        // std::cout << "pt " << pt->getId() << " has type NONE" << std::endl;
        break;
    }
  }

  std::size_t testIndex = round(sph_particles.size() / 2);
  SPHParticleP testParticle = sph_particles[testIndex];
  EXPECT_EQ(sph_particles.size(), 36376);
  EXPECT_EQ(noneP, 0);
  EXPECT_EQ(ghostP, 929);
  EXPECT_EQ(freeP, 34890);
  EXPECT_EQ(boundaryP, 557);
  EXPECT_EQ(testParticle->getId(), 18475);
  EXPECT_EQ(testParticle->getType(), SPHParticleType::FREE);
  EXPECT_NEAR(testParticle->getInitPosition().x(), 49.526315789, 1.0e-6);
  EXPECT_NEAR(testParticle->getInitPosition().y(), 0, 1.0e-6);
  EXPECT_NEAR(testParticle->getInitPosition().z(), 3.1605263157, 1.0e-6);

  /*
  std::cout << "none = " << noneP << " ghost = " << ghostP
            << " free = " << freeP << " boundary = " << boundaryP << std::endl;
  std::cout << " Test: " << testIndex
            << " Id = " << testParticle->getId()
            << " Type = " << static_cast<int>(testParticle->getType())
            << " Pos = " << testParticle->getInitPosition() << std::endl;
  */

  sph_particles.clear();
  sph.clearAllSPHParticleVec();
  noneP = 0;
  ghostP = 0;
  freeP = 0;
  boundaryP = 0;

  dem::Box container3D(dem::Vec(30.0, -1.0, 3.5), dem::Vec(35.0, 1.0, 4.0));
  sph.generateSPHParticle<3>(container3D, allDEMParticles);
  sph_particles = sph.getAllSPHParticleVec();
  // std::cout << "Size:" << sph_particles.size() << "\n";
  for (auto pt : sph_particles) {
    switch (pt->getType()) {
      case SPHParticleType::GHOST:
        ghostP++;
        // std::cout << "pt " << pt->getId() << " has type GHOST" << std::endl;
        break;
      case SPHParticleType::FREE:
        freeP++;
        // std::cout << "pt " << pt->getId() << " has type FREE" << std::endl;
        break;
      case SPHParticleType::BOUNDARY:
        boundaryP++;
        // std::cout << "pt " << pt->getId() << " has type BOUNDARY" <<
        // std::endl;
        break;
      default:
        noneP++;
        // std::cout << "pt " << pt->getId() << " has type NONE" << std::endl;
        break;
    }
  }
  testIndex = round(sph_particles.size() / 2);
  testParticle = sph_particles[testIndex];

  EXPECT_EQ(sph_particles.size(), 46502);
  EXPECT_EQ(noneP, 0);
  EXPECT_EQ(ghostP, 3408);
  EXPECT_EQ(freeP, 36441);
  EXPECT_EQ(boundaryP, 6653);
  EXPECT_EQ(testParticle->getId(), 24042);
  EXPECT_EQ(testParticle->getType(), SPHParticleType::BOUNDARY);
  EXPECT_NEAR(testParticle->getInitPosition().x(), 31.631578947, 1.0e-6);
  EXPECT_NEAR(testParticle->getInitPosition().y(), 1.0526315789, 1.0e-6);
  EXPECT_NEAR(testParticle->getInitPosition().z(), 3.7105263157, 1.0e-6);

  /*
  std::cout << "none = " << noneP << " ghost = " << ghostP
            << " free = " << freeP << " boundary = " << boundaryP << std::endl;
  std::cout << " Test: " << testIndex
            << " Id = " << testParticle->getId()
            << " Type = " << static_cast<int>(testParticle->getType())
            << " Pos = " << testParticle->getInitPosition() << std::endl;
  */
}

TEST(SPHTest, generateSPHParticleNoBottom)
{

  // Setup the parameters that are used by the constructor
  dem::InputParameter::get().addParameter("specificG", 1.1);
  dem::InputParameter::get().addParameter("gravAccel", 10);
  dem::InputParameter::get().addParameter("gravScale", 1);
  dem::InputParameter::get().addParameter("spaceInterval", 0.0526316);
  dem::InputParameter::get().addParameter("numLayers", 1);
  dem::InputParameter::get().addParameter("gamma", 1.4);
  dem::InputParameter::get().addParameter("P0", 1.0e5);
  dem::InputParameter::get().addParameter("SPHInitialDensity", 1000);

  // Create three DEM particles
  std::vector<ParticleID> dem_id = { 0, 1, 2 };
  std::vector<REAL> dem_radii = { 1.15, 0.5,  0.69, 1.15, 0.5,
                                  0.69, 1.15, 0.5,  0.69 };
  std::vector<REAL> dem_axle_a = { 2.250054,  1.570796, 0.6792575,
                                   2.536687,  1.570796, 0.9658907,
                                   0.8923825, 1.570796, 0.6784139 };
  std::vector<REAL> dem_axle_b = { 1.570796, 4.61936e-07,  1.570796,
                                   1.570797, 1.803533e-06, 1.570798,
                                   1.570799, 2.865908e-06, 1.570795 };
  std::vector<REAL> dem_axle_c = { 2.462335, 1.570796, 2.250054,
                                   2.175704, 1.570796, 2.536685,
                                   2.463179, 1.570796, 0.8923828 };
  std::vector<REAL> dem_pos = { 49.73174, 0,        0.9936725, 50.90096, 0,
                                5.320529, 30.81142, 0,         3.898033 };
  std::vector<REAL> dem_vel = { 0.02008031, 0, 0.05505538,
                                -0.3014865, 0, -0.1031938,
                                -0.1538927, 0, 0.04640176 };
  std::vector<REAL> dem_omega = { 0, 0.02781917, 0,          0, -0.07066198,
                                  0, 0,          -0.2856931, 0 };
  std::vector<REAL> dem_force = { -491684,   0,       -126938, -4413.702, 0,
                                  -10805.48, 22966.5, 0,       6613.044 };
  std::vector<REAL> dem_moment = { 0.1688596,  -118712,   -0.3092984,
                                   -0.1174683, -36733.58, -0.05749536,
                                   0.0414009,  -5030.464, -0.0190764 };

  dem::DEMParticlePArray allDEMParticles;
  for (auto id : dem_id) {
    std::size_t startIndex = static_cast<std::size_t>(3 * id);
    std::size_t endIndex = static_cast<std::size_t>(3 * id + 2);
    dem::Vec radii(dem_radii.cbegin() + startIndex,
                   dem_radii.cbegin() + endIndex);
    dem::Vec axle_a(dem_axle_a.cbegin() + startIndex,
                    dem_axle_a.cbegin() + endIndex);
    dem::Vec axle_b(dem_axle_b.cbegin() + startIndex,
                    dem_axle_b.cbegin() + endIndex);
    dem::Vec axle_c(dem_axle_c.cbegin() + startIndex,
                    dem_axle_c.cbegin() + endIndex);
    dem::Vec pos(dem_pos.cbegin() + startIndex, dem_pos.cbegin() + endIndex);
    dem::Vec vel(dem_vel.cbegin() + startIndex, dem_vel.cbegin() + endIndex);
    dem::Vec omega(dem_omega.cbegin() + startIndex,
                   dem_omega.cbegin() + endIndex);
    dem::Vec force(dem_force.cbegin() + startIndex,
                   dem_force.cbegin() + endIndex);
    dem::Vec moment(dem_moment.cbegin() + startIndex,
                    dem_moment.cbegin() + endIndex);
    // std::cout << "id = " << id << " r = " << radii << " a = " << axle_a
    //          << " b = " << axle_b << " c = " << axle_c << " pos = " << pos
    //          << " vel = " << vel << " omega = " << omega 
    //          << " force = " << force
    //          << " moment = " << moment << std::endl;
    dem::DEMParticleP pt = std::make_shared<dem::DEMParticle>(
      id, 0, radii, pos, axle_a, axle_b, axle_c, 1.0e9, 0.3);
    pt->setPrevVeloc(vel);
    pt->setCurrVeloc(vel);
    pt->setPrevOmega(omega);
    pt->setCurrOmega(omega);
    pt->setForce(force);
    pt->setMoment(moment);
    allDEMParticles.push_back(pt);
  }

  SmoothParticleHydro sph;
  dem::Box container(dem::Vec(30.0, -2.0, 0.95), dem::Vec(52.0, 2.0, 5.5));
  sph.generateSPHParticleNoBottom<2>(container, allDEMParticles);

  SPHParticlePArray sph_particles = sph.getAllSPHParticleVec();
  //std::cout << "Size:" << sph_particles.size() << "\n";
  std::size_t noneP = 0;
  std::size_t ghostP = 0;
  std::size_t freeP = 0;
  std::size_t boundaryP = 0;
  for (auto pt : sph_particles) {
    switch (pt->getType()) {
      case SPHParticleType::GHOST:
        ghostP++;
        // std::cout << "pt " << pt->getId() << " has type GHOST" << std::endl;
        break;
      case SPHParticleType::FREE:
        freeP++;
        // std::cout << "pt " << pt->getId() << " has type FREE" << std::endl;
        break;
      case SPHParticleType::BOUNDARY:
        boundaryP++;
        // std::cout << "pt " << pt->getId() << " has type BOUNDARY" <<
        // std::endl;
        break;
      default:
        noneP++;
        // std::cout << "pt " << pt->getId() << " has type NONE" << std::endl;
        break;
    }
  }

  std::size_t testIndex = round(sph_particles.size() / 2);
  SPHParticleP testParticle = sph_particles[testIndex];
  EXPECT_EQ(sph_particles.size(), 37469);
  EXPECT_EQ(noneP, 0);
  EXPECT_EQ(ghostP, 2022);
  EXPECT_EQ(freeP, 35279);
  EXPECT_EQ(boundaryP, 168);
  EXPECT_EQ(testParticle->getId(), 18734);
  EXPECT_EQ(testParticle->getType(), SPHParticleType::FREE);
  EXPECT_NEAR(testParticle->getInitPosition().x(), 41.000004399, 1.0e-6);
  EXPECT_NEAR(testParticle->getInitPosition().y(), 0, 1.0e-6);
  EXPECT_NEAR(testParticle->getInitPosition().z(), 3.2131588, 1.0e-6);

  /*
  std::cout << "none = " << noneP << " ghost = " << ghostP
            << " free = " << freeP << " boundary = " << boundaryP << std::endl;
  std::cout << " Test: " << testIndex
            << " Id = " << testParticle->getId()
            << " Type = " << static_cast<int>(testParticle->getType())
            << " Pos = " << testParticle->getInitPosition() << std::endl;
  */

  sph_particles.clear();
  sph.clearAllSPHParticleVec();
  noneP = 0;
  ghostP = 0;
  freeP = 0;
  boundaryP = 0;

  dem::Box container3D(dem::Vec(30.0, -1.0, 3.5), dem::Vec(35.0, 1.0, 4.0));
  sph.generateSPHParticleNoBottom<3>(container3D, allDEMParticles);
  sph_particles = sph.getAllSPHParticleVec();
  // std::cout << "Size:" << sph_particles.size() << "\n";
  for (auto pt : sph_particles) {
    switch (pt->getType()) {
      case SPHParticleType::GHOST:
        ghostP++;
        // std::cout << "pt " << pt->getId() << " has type GHOST" << std::endl;
        break;
      case SPHParticleType::FREE:
        freeP++;
        // std::cout << "pt " << pt->getId() << " has type FREE" << std::endl;
        break;
      case SPHParticleType::BOUNDARY:
        boundaryP++;
        // std::cout << "pt " << pt->getId() << " has type BOUNDARY" <<
        // std::endl;
        break;
      default:
        noneP++;
        // std::cout << "pt " << pt->getId() << " has type NONE" << std::endl;
        break;
    }
  }

  testIndex = round(sph_particles.size() / 2);
  testParticle = sph_particles[testIndex];

  EXPECT_EQ(sph_particles.size(), 48216);
  EXPECT_EQ(noneP, 0);
  EXPECT_EQ(ghostP, 5122);
  EXPECT_EQ(freeP, 39829);
  EXPECT_EQ(boundaryP, 3265);
  EXPECT_EQ(testParticle->getId(), 24108);
  EXPECT_EQ(testParticle->getType(), SPHParticleType::BOUNDARY);
  EXPECT_NEAR(testParticle->getInitPosition().x(), 29.947368399, 1.0e-6);
  EXPECT_NEAR(testParticle->getInitPosition().y(), -1.0526316, 1.0e-6);
  EXPECT_NEAR(testParticle->getInitPosition().z(), 3.763158, 1.0e-6);

  /*
  std::cout << "none = " << noneP << " ghost = " << ghostP
            << " free = " << freeP << " boundary = " << boundaryP << std::endl;
  std::cout << " Test: " << testIndex
            << " Id = " << testParticle->getId()
            << " Type = " << static_cast<int>(testParticle->getType())
            << " Pos = " << testParticle->getInitPosition() << std::endl;
  */
}