#include <DiscreteElements/DEMPeriodicParticleBC.h>
#include <Core/MechanicsConcepts/Deformations.h>
#include <Core/MechanicsConcepts/StrainTensors.h>
#include <gtest/gtest.h>

using namespace dem;

TEST(DEMPeriodicParticleBCTest, displacementBC) {

  int num_points = 3;
  //int num_components = 3;
  std::vector<double> times(num_points);
  std::vector<Vec> values(num_points);
  for (int ii = 0; ii < num_points; ii++) {
    times[ii] = ii+1;
    values[ii] = Vec(5*(ii+1), 6*(ii+1), 7*(ii+1));
  }

  /* Read xml data */
  zen::XmlDoc doc;
  zen::XmlOut xml(doc);
  std::string timeStr = "";
  std::string valueStr = "";
  for (int ii = 0; ii < num_points; ii++) {
    timeStr += std::to_string(times[ii]) + " ";
    valueStr += std::to_string(values[ii].x()) + " " +
                std::to_string(values[ii].y()) + " " +
                std::to_string(values[ii].z()) + " ";
  }
  xml["time"](timeStr);
  xml["value"](valueStr);
  zen::XmlIn in_xml(doc);

  //std::cout << zen::serialize(doc) << std::endl;

  DEMPeriodicParticleBC<Displacement, 3> disp_bc_xml(in_xml);
  Displacement disp = disp_bc_xml.getBCValue(1.9);
  //std::cout << disp_bc_xml;
  //std::cout << disp << "\n";
  EXPECT_DOUBLE_EQ(disp[0], 9.5);
  EXPECT_DOUBLE_EQ(disp[1], 11.4);
  EXPECT_DOUBLE_EQ(disp[2], 13.3);

  DisplacementRate vel = disp_bc_xml.getBCRate(-1.0);
  EXPECT_DOUBLE_EQ(vel[0], 0);
  EXPECT_DOUBLE_EQ(vel[1], 0);
  EXPECT_DOUBLE_EQ(vel[2], 0);

  vel = disp_bc_xml.getBCRate(5.0);
  EXPECT_DOUBLE_EQ(vel[0], 0);
  EXPECT_DOUBLE_EQ(vel[1], 0);
  EXPECT_DOUBLE_EQ(vel[2], 0);

  vel = disp_bc_xml.getBCRate(1.9);
  EXPECT_NEAR(vel[0], 19.99999999, 1.0e-6);
  EXPECT_NEAR(vel[1], 23.99999999, 1.0e-6);
  EXPECT_NEAR(vel[2], 28.00000000, 1.0e-6);
  //std::cout << vel << "\n";
}

TEST(DEMPeriodicParticleBCTest, deformationGradientBC) {

  int num_points = 3;
  constexpr int num_components = 9;
  std::vector<double> times(num_points);
  std::vector<std::array<double, num_components>> values(num_points);
  for (int ii = 0; ii < num_points; ii++) {
    times[ii] = ii+1;
    values[ii] = {{5.0*(ii+1), 6.0*(ii+1), 7.0*(ii+1),
                   0.1*(ii+1), 0.2*(ii+1), 0.3*(ii+1),
                   1.1*(ii+1), 1.2*(ii+1), 1.3*(ii+1)}};
  }

  /* Read xml data */
  zen::XmlDoc doc;
  zen::XmlOut xml(doc);
  std::string timeStr = "";
  std::string valueStr = "";
  for (int ii = 0; ii < num_points; ii++) {
    timeStr += std::to_string(times[ii]) + " ";
    for (int jj = 0; jj < num_components; jj++) {
      valueStr += std::to_string(values[ii][jj]) + " ";
    }
    valueStr += "\n";
  }
  xml["time"](timeStr);
  xml["value"](valueStr);
  zen::XmlIn in_xml(doc);

  //std::cout << zen::serialize(doc) << std::endl;

  DEMPeriodicParticleBC<DeformationGradient, 9> defgrad_bc_xml(in_xml);
  DeformationGradient defgrad = defgrad_bc_xml.getBCValue(1.9);
  //std::cout << defgrad_bc_xml;
  //std::cout << defgrad << "\n";

  DEMPeriodicParticleBC<DeformationGradient, 9> defgrad_bc_xml1;
  defgrad_bc_xml1.read(in_xml);
  defgrad = defgrad_bc_xml1.getBCValue(1.9);
  //std::cout << defgrad_bc_xml1;
  EXPECT_DOUBLE_EQ(defgrad(0,0), 9.5);
  EXPECT_DOUBLE_EQ(defgrad(0,1), 11.4);
  EXPECT_DOUBLE_EQ(defgrad(0,2), 13.3);
  EXPECT_DOUBLE_EQ(defgrad(1,0), 0.19);
  EXPECT_DOUBLE_EQ(defgrad(1,1), 0.38);
  EXPECT_DOUBLE_EQ(defgrad(1,2), 0.57);
  EXPECT_DOUBLE_EQ(defgrad(2,0), 2.09);
  EXPECT_DOUBLE_EQ(defgrad(2,1), 2.28);
  EXPECT_DOUBLE_EQ(defgrad(2,2), 2.47);

  DeformationGradientRate Fdot = defgrad_bc_xml.getBCRate(-1.0);
  EXPECT_DOUBLE_EQ(Fdot(0,0), 0);
  EXPECT_DOUBLE_EQ(Fdot(0,1), 0);
  EXPECT_DOUBLE_EQ(Fdot(0,2), 0);

  Fdot = defgrad_bc_xml.getBCRate(5.0);
  EXPECT_DOUBLE_EQ(Fdot(0,0), 0);
  EXPECT_DOUBLE_EQ(Fdot(0,1), 0);
  EXPECT_DOUBLE_EQ(Fdot(0,2), 0);

  Fdot = defgrad_bc_xml.getBCRate(1.9);
  EXPECT_NEAR(Fdot(0,0), 19.99999999, 1.0e-6);
  EXPECT_NEAR(Fdot(0,1), 23.99999999, 1.0e-6);
  EXPECT_NEAR(Fdot(0,2), 28.00000000, 1.0e-6);
  //std::cout << Fdot << "\n";
}

TEST(DEMPeriodicParticleBCTest, axisymmetricStrainBC) {

  int num_points = 3;
  constexpr int num_components = 4;
  std::vector<double> times(num_points);
  std::vector<std::array<double, num_components>> values(num_points);
  for (int ii = 0; ii < num_points; ii++) {
    times[ii] = ii+1;
    values[ii] = {{5.0*(ii+1), 6.0*(ii+1), 7.0*(ii+1),
                   0.1*(ii+1)}};
  }

  /* Read xml data */
  zen::XmlDoc doc;
  zen::XmlOut xml(doc);
  std::string timeStr = "";
  std::string valueStr = "";
  for (int ii = 0; ii < num_points; ii++) {
    timeStr += std::to_string(times[ii]) + " ";
    for (int jj = 0; jj < num_components; jj++) {
      valueStr += std::to_string(values[ii][jj]) + " ";
    }
    valueStr += "\n";
  }
  xml["time"](timeStr);
  xml["value"](valueStr);
  zen::XmlIn in_xml(doc);

  //std::cout << zen::serialize(doc) << std::endl;

  DEMPeriodicParticleBC<AxisymmetricStrain, 4> axistrain_bc_xml(in_xml);
  AxisymmetricStrain axistrain = axistrain_bc_xml.getBCValue(1.9);
  //std::cout << axistrain_bc_xml;
  //std::cout << axistrain << "\n";
  EXPECT_DOUBLE_EQ(axistrain[0], 9.5);
  EXPECT_DOUBLE_EQ(axistrain[1], 11.4);
  EXPECT_DOUBLE_EQ(axistrain[2], 13.3);
  EXPECT_DOUBLE_EQ(axistrain[3], 0.19);

  AxisymmetricStrainRate epsdot = axistrain_bc_xml.getBCRate(-1.0);
  EXPECT_DOUBLE_EQ(epsdot[0], 0);
  EXPECT_DOUBLE_EQ(epsdot[1], 0);
  EXPECT_DOUBLE_EQ(epsdot[2], 0);

  epsdot = axistrain_bc_xml.getBCRate(5.0);
  EXPECT_DOUBLE_EQ(epsdot[0], 0);
  EXPECT_DOUBLE_EQ(epsdot[1], 0);
  EXPECT_DOUBLE_EQ(epsdot[2], 0);

  epsdot = axistrain_bc_xml.getBCRate(1.9);
  EXPECT_NEAR(epsdot[0], 19.99999999, 1.0e-6);
  EXPECT_NEAR(epsdot[1], 23.99999999, 1.0e-6);
  EXPECT_NEAR(epsdot[2], 28.00000000, 1.0e-6);
  //std::cout << epsdot << "\n";
}