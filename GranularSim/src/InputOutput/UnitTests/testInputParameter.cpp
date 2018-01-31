#include <InputOutput/InputParameter.h>
#include <Core/Util/Utility.h>
#include <DiscreteElements/DEMBoundaryConditionUtils.h>
#include <gtest/gtest.h>

using namespace dem;

TEST(InputParameterTest, testInputXMLFileReader) {

  std::string filename = "input001_p1.xml";
  InputParameter::get().readInXML(filename);
  //InputParameter::get().writeOut();

  auto param = dem::InputParameter::get().param;
  //for (const auto& p : param) {
  //  std::cout << p.first << "  " << p.second << std::endl;
  //}
  //std::cout << "size = " << param.size() << "\n";
  EXPECT_EQ(param.size(), 88);

  EXPECT_NEAR(util::getParam<REAL>("Aphi"), 1.55678, 1.0e-5);
  EXPECT_NEAR(util::getParam<REAL>("Apsi"), 1.55678, 1.0e-5);
  EXPECT_NEAR(util::getParam<REAL>("Bphi"), 1.41767, 1.0e-5);
  EXPECT_NEAR(util::getParam<REAL>("Bpsi"), 1.41767, 1.0e-5);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("Xmax"), 100);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("Xmin"), -20);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("Ymax"), 20);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("Ymin"), -20);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("Zmax"), 50);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("Zmin"), 0);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("beta"), -1);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("bodyDensity"), 5095.5);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("bondStretchLimit"), 100);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("boundaryFric"), 0.5);
  EXPECT_NEAR(util::getParam<REAL>("c"), 208848, 1.0);
  EXPECT_NEAR(util::getParam<REAL>("chi"), 208848, 1.0);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("contactCohesion"), 0);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("contactDamp"), 0.55);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("contactFric"), 0.5);
  EXPECT_EQ(util::getParam<int>("demToInitParticle"), 1);
  EXPECT_EQ(util::getParam<int>("endSnap"), 5);
  EXPECT_EQ(util::getParam<int>("endStep"), 5);
  EXPECT_EQ(util::getParam<int>("fixRadius"), 0);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("forceDamp"), 0);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("gravAccel"), 9.8);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("gravScale"), 0);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("hchi"), 0);
  EXPECT_NEAR(util::getParam<REAL>("kBulk"), 3.33333e+07, 1.0e2);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("kappa"), -50000);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("lambda"), 2e+07);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("massScale"), 1);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("maxAllowableRelativeOverlap"), 0.01);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("minMeasurableOverlap"), 1e-08);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("memYoung"), 1.4e+06);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("minAllowableRelativeOverlap"), 1e-06);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("mntScale"), 1);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("momentDamp"), 0);
  EXPECT_EQ(util::getParam<int>("mpiProcX"), 1);
  EXPECT_EQ(util::getParam<int>("mpiProcY"), 1);
  EXPECT_EQ(util::getParam<int>("mpiProcZ"), 1);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("mu"), 2e+07);
  EXPECT_EQ(util::getParam<int>("ompThreads"), 1);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("periDensity"), 1250);
  EXPECT_EQ(util::getParam<int>("periFixCentroidX"), 0);
  EXPECT_EQ(util::getParam<int>("periFixCentroidY"), 0);
  EXPECT_EQ(util::getParam<int>("periFixCentroidZ"), 0);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("periForce"), 0);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("periGeomScaleFac"), 1);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("periPoisson"), 0.25);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("periReflVecX"), 1);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("periReflVecY"), 1);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("periReflVecZ"), 1);
  EXPECT_EQ(util::getParam<int>("periToInitParticle"), 1);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("periTransVecX"), 50);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("periTransVecY"), 0);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("periTransVecZ"), 0);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("periYoung") ,5e+07);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("periodicBoundaryFaceShiftFactor"), 0);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("periodicBoundaryMarginFactor"), 2);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("phi"), 0.738663);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("pileRate"), 0);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("poisson"), 0.25);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("psi"), 0.738663);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("rEllip"), 1);
  EXPECT_EQ(util::getParam<int>("rampStep"), 1);
  EXPECT_EQ(util::getParam<int>("removeDEMParticles"), 1);
  EXPECT_EQ(util::getParam<int>("removePeriParticles"), 0);
  EXPECT_EQ(util::getParam<int>("simuType"), 1);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("specificG"), 2.65);
  EXPECT_EQ(util::getParam<int>("startSnap"), 1);
  EXPECT_EQ(util::getParam<int>("startStep"), 1);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("tangentModulus11"), 6e+07);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("tangentModulus12"), 2e+07);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("tangentModulus13"), 2e+07);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("tangentModulus21"), 2e+07);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("tangentModulus22"), 6e+07);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("tangentModulus23"), 2e+07);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("tangentModulus31"), 2e+07);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("tangentModulus32"), 2e+07);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("tangentModulus33"), 6e+07);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("tangentModulus44"), 2e+07);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("tangentModulus55"), 2e+07);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("tangentModulus66"), 2e+07);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("timeAccrued"), 0);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("timeStep"), 1e-06);
  EXPECT_EQ(util::getParam<int>("typeConstitutive"), 1);
  EXPECT_DOUBLE_EQ(util::getParam<REAL>("young"), 4.5e+10);
  EXPECT_EQ(util::getFilename("demBoundaryConditionFilename"), "none");
  EXPECT_EQ(util::getFilename("boundaryFilename"), "input001_boundary.xml");
  EXPECT_EQ(util::getFilename("outputFolder"), "deposit.pe3d");
  EXPECT_EQ(util::getFilename("particleFilename"), "input001_particles.xml");
  EXPECT_EQ(util::getFilename("periFilename"), "peri_part.inp");
  EXPECT_EQ(util::getFilename("generatedBoundaryOutputFilename"), "generated_periodic_boundary");
  EXPECT_EQ(util::getFilename("generatedParticleOutputFilename"), "generated_periodic_particles");
}

