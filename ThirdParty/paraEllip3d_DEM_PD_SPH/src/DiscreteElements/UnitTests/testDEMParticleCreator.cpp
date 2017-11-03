#include <DiscreteElements/DEMParticleCreator.h>
#include <Core/Util/Utility.h>
#include <InputOutput/InputParameter.h>
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
  particle2.setRadiusA(a);
  particle2.setRadiusB(b);
  particle2.setRadiusC(c);
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
  particle3.setRadiusA(a);
  particle3.setRadiusB(b);
  particle3.setRadiusC(c);
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
  particle4.setRadiusA(a);
  particle4.setRadiusB(b);
  particle4.setRadiusC(c);
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
  particle5.setRadiusA(a);
  particle5.setRadiusB(b);
  particle5.setRadiusC(c);
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
  particle6.setRadiusA(a);
  particle6.setRadiusB(b);
  particle6.setRadiusC(c);
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
  particle7.setRadiusA(a);
  particle7.setRadiusB(b);
  particle7.setRadiusC(c);
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
  particle8.setRadiusA(a);
  particle8.setRadiusB(b);
  particle8.setRadiusC(c);
  particle8.setCurrentPosition(pos);
  particle8.setCurrentAnglesAxisA(vacos(adir));
  particle8.setCurrentAnglesAxisB(vacos(bdir));
  particle8.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle8, 
                                                   [](DEMParticle*){}));

  DEMParticle particle9;
  particle9.setRadiusA(0.5*a);
  particle9.setRadiusB(0.5*b);
  particle9.setRadiusC(0.5*c);
  particle9.setCurrentPosition(pos);
  particle9.setCurrentAnglesAxisA(vacos(adir));
  particle9.setCurrentAnglesAxisB(vacos(bdir));
  particle9.setCurrentAnglesAxisC(vacos(cdir));
  particles.push_back(std::shared_ptr<DEMParticle>(&particle9, 
                                                   [](DEMParticle*){}));

  DEMParticleCreator creator;
  DEMParticlePArray periodic = 
    creator.generatePeriodicDEMParticles(particles, boundingBox);
  
  for (const auto particle : periodic) {
    std::cout << *particle;
  }
}

