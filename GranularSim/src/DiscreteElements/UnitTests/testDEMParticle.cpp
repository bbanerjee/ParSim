#include <DiscreteElements/DEMParticle.h>
#include <Core/Util/Utility.h>
#include <InputOutput/InputParameter.h>
#include <gtest/gtest.h>

using namespace dem;

// Spherical particle - axis aligned
TEST(DEMParticleTest, containsPoint1) {

  DEMParticle particle;

  REAL a = 1;
  REAL b = 1;
  REAL c = 1;
  
  dem::Vec adir(0, dem::Pi/2, dem::Pi/2);
  dem::Vec bdir(dem::Pi/2, 0, dem::Pi/2);
  dem::Vec cdir(dem::Pi/2, dem::Pi/2, 0);

  dem::Vec pos(0,0,0);

  particle.setRadiusA(a);
  particle.setRadiusB(b);
  particle.setRadiusC(c);

  particle.setCurrentPosition(pos);

  particle.setCurrentAnglesAxisA(adir);
  particle.setCurrentAnglesAxisB(bdir);
  particle.setCurrentAnglesAxisC(cdir);

  REAL buffer = 0;
  dem::Vec point = 0;
  dem::Vec localCoord = 0;
  bool inside = false;
  bool insideGhostLayer = false;

  buffer = -a;
  point = dem::Vec(0.5*a, 0, 0);
  inside = particle.containsPoint(point, pos, buffer, localCoord, insideGhostLayer);
  EXPECT_EQ(inside, true);
  EXPECT_EQ(insideGhostLayer, false);

  buffer = 0;
  point = dem::Vec(0, a, 0);
  inside = particle.containsPoint(point, pos, buffer, localCoord, insideGhostLayer);
  EXPECT_EQ(inside, true);
  EXPECT_EQ(insideGhostLayer, false);

  buffer = a;
  point = dem::Vec(0, 0, 2.0*a);
  inside = particle.containsPoint(point, pos, buffer, localCoord, insideGhostLayer);
  EXPECT_EQ(inside, false);
  EXPECT_EQ(insideGhostLayer, false);

  buffer = 2*a;
  point = dem::Vec(0, 0, 2.0*a);
  inside = particle.containsPoint(point, pos, buffer, localCoord, insideGhostLayer);
  EXPECT_EQ(inside, false);
  EXPECT_EQ(insideGhostLayer, false);

  buffer = 0.2*a;
  point = dem::Vec(0.9*a, 0, 0);
  inside = particle.containsPoint(point, pos, buffer, localCoord, insideGhostLayer);
  EXPECT_EQ(inside, true);
  EXPECT_EQ(insideGhostLayer, true);

  point = dem::Vec(0, 0.9*a, 0);
  inside = particle.containsPoint(point, pos, buffer, localCoord, insideGhostLayer);
  EXPECT_EQ(inside, true);
  EXPECT_EQ(insideGhostLayer, true);

  point = dem::Vec(0, 0, 0.9*a);
  inside = particle.containsPoint(point, pos, buffer, localCoord, insideGhostLayer);
  EXPECT_EQ(inside, true);
  EXPECT_EQ(insideGhostLayer, true);

  buffer = 0.3*a;
  point = dem::Vec(0.81*a*cos(dem::Pi/3), 0.81*a*cos(dem::Pi/3), 0.9*a*cos(dem::Pi/3));
  inside = particle.containsPoint(point, pos, buffer, localCoord, insideGhostLayer);
  EXPECT_EQ(inside, true);
  EXPECT_EQ(insideGhostLayer, true);

  buffer = 0.3001;
  point = dem::Vec(0.7, 0.7, 0.7);
  inside = particle.containsPoint(point, pos, buffer, localCoord, insideGhostLayer);
  EXPECT_EQ(inside, false);
  EXPECT_EQ(insideGhostLayer, false);

  /*
  std::cout << "Pt: " << point
            << "Particle center: " << pos
            << "LocalCoord = " << localCoord
            << "inside = " << std::boolalpha << inside << std::endl;
  */

}

