#include <Boundary/PeriodicParticleBC.h>
#include <Core/MechanicsConcepts/Deformations.h>
#include <Core/MechanicsConcepts/StrainTensors.h>
#include <gtest/gtest.h>

using namespace dem;

TEST(PeriodicParticleBCTest, displacementBC) {

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

  std::cout << zen::serialize(doc) << std::endl;

  PeriodicParticleBC<Displacement, 3> disp_bc_xml(in_xml);
  Displacement disp = disp_bc_xml.getBCValue(1.9);
  std::cout << disp_bc_xml;
  std::cout << disp << "\n";
}

TEST(PeriodicParticleBCTest, deformationGradientBC) {

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

  std::cout << zen::serialize(doc) << std::endl;

  PeriodicParticleBC<DeformationGradient, 9> defgrad_bc_xml(in_xml);
  DeformationGradient defgrad = defgrad_bc_xml.getBCValue(1.9);
  std::cout << defgrad_bc_xml;
  std::cout << defgrad << "\n";
}

TEST(PeriodicParticleBCTest, axisymmetricStrainBC) {

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

  std::cout << zen::serialize(doc) << std::endl;

  PeriodicParticleBC<AxisymmetricStrain, 4> axistrain_bc_xml(in_xml);
  AxisymmetricStrain axistrain = axistrain_bc_xml.getBCValue(1.9);
  std::cout << axistrain_bc_xml;
  std::cout << axistrain << "\n";
}