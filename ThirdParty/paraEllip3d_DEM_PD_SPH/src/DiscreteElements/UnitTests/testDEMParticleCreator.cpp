#include <DiscreteElements/DEMParticleCreator.h>
#include <Core/Util/Utility.h>
#include <InputOutput/InputParameter.h>
#include <InputOutput/OutputTecplot.h>
#include <gtest/gtest.h>

using namespace dem;

TEST(DEMParticleCreatorTest, periodic) {

  Box boundingBox(Vec(-10, -10, -10), Vec(10, 10, 10));

  DEMParticlePArray particles;

  REAL a = 1;
  REAL b = 2;
  REAL c = 3;
  dem::Vec adir(0.00000,  1.00000,  0.00000);
  dem::Vec bdir(-1.00000,  0.00000,  0.00000);
  dem::Vec cdir(0.00000,  0.00000,  1.00000);
  dem::Vec pos(2, 1.5, 1.0);

  DEMParticle particle1;
  particle1.setRadiusA(a);
  particle1.setRadiusB(b);
  particle1.setRadiusC(c);
  particle1.setCurrentPosition(pos);
  particle1.setCurrentAnglesAxisA(vacos(adir));
  particle1.setCurrentAnglesAxisB(vacos(bdir));
  particle1.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle1, 
                                                   [](DEMParticle*){}));

  adir = Vec(0.00000,  1.00000,  0.00000);
  bdir = Vec(-1.00000,  0.00000,  0.00000);
  cdir = Vec(0.00000,  0.00000,  1.00000);
  pos = Vec(10, 1.5, 1.0);
  DEMParticle particle2;
  particle2.setRadiusA(1.9*a);
  particle2.setRadiusB(0.8*b);
  particle2.setRadiusC(0.3*c);
  particle2.setCurrentPosition(pos);
  particle2.setCurrentAnglesAxisA(vacos(adir));
  particle2.setCurrentAnglesAxisB(vacos(bdir));
  particle2.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle2, 
                                                   [](DEMParticle*){}));

  adir = Vec(1,  0,  0);
  bdir = Vec(0,  1,  0);
  cdir = Vec(0,  0,  1);
  pos = Vec(-10, -1.5, 1);
  DEMParticle particle3;
  particle3.setRadiusA(1.1*a);
  particle3.setRadiusB(0.6*b);
  particle3.setRadiusC(0.5*c);
  particle3.setCurrentPosition(pos);
  particle3.setCurrentAnglesAxisA(vacos(adir));
  particle3.setCurrentAnglesAxisB(vacos(bdir));
  particle3.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle3, 
                                                   [](DEMParticle*){}));

  adir = Vec(0.00000,  1.00000,  0.00000);
  bdir = Vec(-1.00000,  0.00000,  0.00000);
  cdir = Vec(0.00000,  0.00000,  1.00000);
  pos = Vec(-0.5, -10.5, 1);
  DEMParticle particle4;
  particle4.setRadiusA(1.5*a);
  particle4.setRadiusB(0.2*b);
  particle4.setRadiusC(0.3*c);
  particle4.setCurrentPosition(pos);
  particle4.setCurrentAnglesAxisA(vacos(adir));
  particle4.setCurrentAnglesAxisB(vacos(bdir));
  particle4.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle4, 
                                                   [](DEMParticle*){}));


  adir = Vec(0.76604, -0.64279,  0.00000);
  bdir = Vec(0.64279,  0.76604,  0.00000);
  cdir = Vec(0.00000,  0.00000,  1.00000);
  pos = Vec(10.6, -10, 1);
  DEMParticle particle5;
  particle5.setRadiusA(0.9*a);
  particle5.setRadiusB(0.8*b);
  particle5.setRadiusC(0.6*c);
  particle5.setCurrentPosition(pos);
  particle5.setCurrentAnglesAxisA(vacos(adir));
  particle5.setCurrentAnglesAxisB(vacos(bdir));
  particle5.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle5, 
                                                   [](DEMParticle*){}));

  adir = Vec(0.76604,  0.64279,  0.00000);
  bdir = Vec(-0.64279,  0.76604,  0.00000);
  cdir = Vec(0.00000,  0.00000,  1.00000);
  pos = Vec(10, 10, 1);
  DEMParticle particle6;
  particle6.setRadiusA(1.7*a);
  particle6.setRadiusB(0.5*b);
  particle6.setRadiusC(0.1*c);
  particle6.setCurrentPosition(pos);
  particle6.setCurrentAnglesAxisA(vacos(adir));
  particle6.setCurrentAnglesAxisB(vacos(bdir));
  particle6.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle6, 
                                                   [](DEMParticle*){}));

  adir = Vec(0.93969, -0.34202,  0.00000);
  bdir = Vec(0.34202,  0.93969,  0.00000);
  cdir = Vec(0.00000,  0.00000,  1.00000);
  pos = Vec(-10, 10, 1);
  DEMParticle particle7;
  particle7.setRadiusA(1.2*a);
  particle7.setRadiusB(0.8*b);
  particle7.setRadiusC(0.8*c);
  particle7.setCurrentPosition(pos);
  particle7.setCurrentAnglesAxisA(vacos(adir));
  particle7.setCurrentAnglesAxisB(vacos(bdir));
  particle7.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle7, 
                                                   [](DEMParticle*){}));

  adir = Vec(0.93969,  0.34202,  0.00000);
  bdir = Vec(-0.34202,  0.93969,  0.00000);
  cdir = Vec(0.00000,  0.00000,  1.00000);
  pos = Vec(-10, -10, 1);
  DEMParticle particle8;
  particle8.setRadiusA(1.4*a);
  particle8.setRadiusB(0.1*b);
  particle8.setRadiusC(0.5*c);
  particle8.setCurrentPosition(pos);
  particle8.setCurrentAnglesAxisA(vacos(adir));
  particle8.setCurrentAnglesAxisB(vacos(bdir));
  particle8.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle8, 
                                                   [](DEMParticle*){}));

  DEMParticle particle9;
  particle9.setRadiusA(0.3);
  particle9.setRadiusB(0.3);
  particle9.setRadiusC(0.3);
  pos = Vec(10, 10, 10);
  particle9.setCurrentPosition(pos);
  particle9.setCurrentAnglesAxisA(vacos(adir));
  particle9.setCurrentAnglesAxisB(vacos(bdir));
  particle9.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle9, 
                                                   [](DEMParticle*){}));

  DEMParticle particle90;
  particle90.setRadiusA(0.5);
  particle90.setRadiusB(0.5);
  particle90.setRadiusC(0.5);
  pos = Vec(-10, 10, 10);
  particle90.setCurrentPosition(pos);
  particle90.setCurrentAnglesAxisA(vacos(adir));
  particle90.setCurrentAnglesAxisB(vacos(bdir));
  particle90.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle90, 
                                                   [](DEMParticle*){}));

  DEMParticle particle91;
  particle91.setRadiusA(0.7);
  particle91.setRadiusB(0.7);
  particle91.setRadiusC(0.7);
  pos = Vec(10, -10, 10);
  particle91.setCurrentPosition(pos);
  particle91.setCurrentAnglesAxisA(vacos(adir));
  particle91.setCurrentAnglesAxisB(vacos(bdir));
  particle91.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle91, 
                                                   [](DEMParticle*){}));

  DEMParticle particle92;
  particle92.setRadiusA(0.9);
  particle92.setRadiusB(0.9);
  particle92.setRadiusC(0.9);
  pos = Vec(-10, -10, 10);
  particle92.setCurrentPosition(pos);
  particle92.setCurrentAnglesAxisA(vacos(adir));
  particle92.setCurrentAnglesAxisB(vacos(bdir));
  particle92.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle92, 
                                                   [](DEMParticle*){}));

  DEMParticle particle93;
  particle93.setRadiusA(1.1);
  particle93.setRadiusB(1.1);
  particle93.setRadiusC(1.1);
  pos = Vec(-10, 10, -10);
  particle93.setCurrentPosition(pos);
  particle93.setCurrentAnglesAxisA(vacos(adir));
  particle93.setCurrentAnglesAxisB(vacos(bdir));
  particle93.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle93, 
                                                   [](DEMParticle*){}));

  DEMParticle particle94;
  particle94.setRadiusA(1.3);
  particle94.setRadiusB(1.3);
  particle94.setRadiusC(1.3);
  pos = Vec(10, -10, -10);
  particle94.setCurrentPosition(pos);
  particle94.setCurrentAnglesAxisA(vacos(adir));
  particle94.setCurrentAnglesAxisB(vacos(bdir));
  particle94.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle94, 
                                                   [](DEMParticle*){}));

  DEMParticle particle95;
  particle95.setRadiusA(1.5);
  particle95.setRadiusB(1.5);
  particle95.setRadiusC(1.5);
  pos = Vec(10, 10, -10);
  particle95.setCurrentPosition(pos);
  particle95.setCurrentAnglesAxisA(vacos(adir));
  particle95.setCurrentAnglesAxisB(vacos(bdir));
  particle95.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle95, 
                                                   [](DEMParticle*){}));

  adir = Vec(0.76604,  0.64279,  0.00000);
  bdir = Vec(-0.64279,  0.76604,  0.00000);
  cdir = Vec(0.00000,  0.00000,  1.00000);
  pos = Vec(-10, -10, -10);
  DEMParticle particle10;
  particle10.setRadiusA(a);
  particle10.setRadiusB(b);
  particle10.setRadiusC(1.2*c);
  particle10.setCurrentPosition(pos);
  particle10.setCurrentAnglesAxisA(vacos(adir));
  particle10.setCurrentAnglesAxisB(vacos(bdir));
  particle10.setCurrentAnglesAxisC(vacos(cdir));

  particles.push_back(std::shared_ptr<DEMParticle>(&particle10, 
                                                   [](DEMParticle*){}));

  /*
  for (const auto particle : particles) {
    std::cout << *particle;
  }
  */

  /*
  DEMParticleCreator creator;
  DEMParticlePArray periodic = 
    creator.generatePeriodicDEMParticles(particles, boundingBox);

  particles.insert(particles.end(), periodic.begin(), periodic.end());
  
  auto outputFolder = util::createOutputFolder("test_particles");
  OutputTecplot<DEMParticlePArray> writer(outputFolder, 0);
  writer.writeParticles(&particles, 0);

  outputFolder = util::createOutputFolder("periodic_particles");
  OutputTecplot<DEMParticlePArray> writer1(outputFolder, 0);
  writer1.writeParticles(&periodic, 0);

  //EXPECT_EQ(periodic.size(), 11);
  for (const auto particle : periodic) {
    std::cout << *particle;
  }
  */
}

TEST(DEMParticleCreatorTest, periodicVertex) {

  Box boundingBox(Vec(-10, -10, -10), Vec(10, 10, 10));

  dem::Vec adir(0.00000,  1.00000,  0.00000);
  dem::Vec bdir(-1.00000,  0.00000,  0.00000);
  dem::Vec cdir(0.00000,  0.00000,  1.00000);
  dem::Vec pos(2, 1.5, 1.0);

  DEMParticlePArray particles;

  DEMParticle particle9;
  particle9.setRadiusA(0.3);
  particle9.setRadiusB(0.3);
  particle9.setRadiusC(0.3);
  pos = Vec(10, 10, 10);
  particle9.setCurrentPosition(pos);
  particle9.setCurrentAnglesAxisA(vacos(adir));
  particle9.setCurrentAnglesAxisB(vacos(bdir));
  particle9.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle9, 
                                                   [](DEMParticle*){}));

  DEMParticle particle90;
  particle90.setRadiusA(0.5);
  particle90.setRadiusB(0.5);
  particle90.setRadiusC(0.5);
  pos = Vec(-10, 10, 10);
  particle90.setCurrentPosition(pos);
  particle90.setCurrentAnglesAxisA(vacos(adir));
  particle90.setCurrentAnglesAxisB(vacos(bdir));
  particle90.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle90, 
                                                   [](DEMParticle*){}));

  DEMParticle particle91;
  particle91.setRadiusA(0.7);
  particle91.setRadiusB(0.7);
  particle91.setRadiusC(0.7);
  pos = Vec(10, -10, 10);
  particle91.setCurrentPosition(pos);
  particle91.setCurrentAnglesAxisA(vacos(adir));
  particle91.setCurrentAnglesAxisB(vacos(bdir));
  particle91.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle91, 
                                                   [](DEMParticle*){}));

  DEMParticle particle92;
  particle92.setRadiusA(0.9);
  particle92.setRadiusB(0.9);
  particle92.setRadiusC(0.9);
  pos = Vec(-10, -10, 10);
  particle92.setCurrentPosition(pos);
  particle92.setCurrentAnglesAxisA(vacos(adir));
  particle92.setCurrentAnglesAxisB(vacos(bdir));
  particle92.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle92, 
                                                   [](DEMParticle*){}));

  DEMParticle particle93;
  particle93.setRadiusA(1.1);
  particle93.setRadiusB(1.1);
  particle93.setRadiusC(1.1);
  pos = Vec(-10, 10, -10);
  particle93.setCurrentPosition(pos);
  particle93.setCurrentAnglesAxisA(vacos(adir));
  particle93.setCurrentAnglesAxisB(vacos(bdir));
  particle93.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle93, 
                                                   [](DEMParticle*){}));

  DEMParticle particle94;
  particle94.setRadiusA(1.3);
  particle94.setRadiusB(1.3);
  particle94.setRadiusC(1.3);
  pos = Vec(10, -10, -10);
  particle94.setCurrentPosition(pos);
  particle94.setCurrentAnglesAxisA(vacos(adir));
  particle94.setCurrentAnglesAxisB(vacos(bdir));
  particle94.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle94, 
                                                   [](DEMParticle*){}));

  DEMParticle particle95;
  particle95.setRadiusA(1.5);
  particle95.setRadiusB(1.5);
  particle95.setRadiusC(1.5);
  pos = Vec(10, 10, -10);
  particle95.setCurrentPosition(pos);
  particle95.setCurrentAnglesAxisA(vacos(adir));
  particle95.setCurrentAnglesAxisB(vacos(bdir));
  particle95.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle95, 
                                                   [](DEMParticle*){}));

  REAL a = 1;
  REAL b = 2;
  REAL c = 3;
  adir = Vec(0.76604,  0.64279,  0.00000);
  bdir = Vec(-0.64279,  0.76604,  0.00000);
  cdir = Vec(0.00000,  0.00000,  1.00000);
  pos = Vec(-10, -10, -10);
  DEMParticle particle10;
  particle10.setRadiusA(a);
  particle10.setRadiusB(b);
  particle10.setRadiusC(1.2*c);
  particle10.setCurrentPosition(pos);
  particle10.setCurrentAnglesAxisA(vacos(adir));
  particle10.setCurrentAnglesAxisB(vacos(bdir));
  particle10.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle10, 
                                                   [](DEMParticle*){}));

  DEMParticleCreator creator;
  DEMParticlePArray periodic = 
    creator.generatePeriodicDEMParticles(particles, boundingBox);

  particles.insert(particles.end(), periodic.begin(), periodic.end());
  
  auto outputFolder = util::createOutputFolder("test_particles_vertex");
  OutputTecplot<DEMParticlePArray> writer(outputFolder, 0);
  writer.writeParticles(&particles, 0);

  //EXPECT_EQ(periodic.size(), 11);
  /*
  for (const auto particle : periodic) {
    std::cout << *particle;
  }
  */
}

