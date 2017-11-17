#include <DiscreteElements/DEMBoundaryConditions.h>
#include <Core/MechanicsConcepts/Deformations.h>
#include <Core/MechanicsConcepts/StrainTensors.h>
#include <gtest/gtest.h>

using namespace dem;

TEST(DEMBoundaryConditionsTest, displacementBC) {

  std::string filename = "input_dem_particle_disp_bc.xml";

  DEMBoundaryConditions disp_bc;
  disp_bc.read(filename);
  Displacement disp = disp_bc.getDisplacement(1.9);
  //std::cout << disp << "\n";
  EXPECT_DOUBLE_EQ(disp.data()[0], 9.5);
  EXPECT_DOUBLE_EQ(disp.data()[1], 11.4);
  EXPECT_DOUBLE_EQ(disp.data()[2], 13.3);
}

TEST(DEMBoundaryConditionsTest, deformationGradientBC) {

  std::string filename = "input_dem_particle_defgrad_bc.xml";

  DEMBoundaryConditions defgrad_bc;
  defgrad_bc.read(filename);
  DeformationGradient defgrad = defgrad_bc.getDeformationGradient(1.9);
  //std::cout << defgrad << "\n";
  EXPECT_DOUBLE_EQ(defgrad.data()[0], 9.5);
  EXPECT_DOUBLE_EQ(defgrad.data()[1], 11.4);
  EXPECT_DOUBLE_EQ(defgrad.data()[2], 13.3);
  EXPECT_DOUBLE_EQ(defgrad.data()[3], 0.19);
  EXPECT_DOUBLE_EQ(defgrad.data()[4], 0.38);
  EXPECT_DOUBLE_EQ(defgrad.data()[5], 0.57);
  EXPECT_DOUBLE_EQ(defgrad.data()[6], 2.09);
  EXPECT_DOUBLE_EQ(defgrad.data()[7], 2.28);
  EXPECT_DOUBLE_EQ(defgrad.data()[8], 2.47);
}

TEST(DEMBoundaryConditionsTest, axisymmetricStrainBC) {

  std::string filename = "input_dem_particle_axisymm_strain_bc.xml";

  DEMBoundaryConditions axisymm_bc;
  axisymm_bc.read(filename);
  AxisymmetricStrain axistrain = axisymm_bc.getAxisymmetricStrain(1.9);
  //std::cout << axistrain << "\n";
  EXPECT_DOUBLE_EQ(axistrain.data()[0], 9.5);
  EXPECT_DOUBLE_EQ(axistrain.data()[1], 11.4);
  EXPECT_DOUBLE_EQ(axistrain.data()[2], 13.3);
  EXPECT_DOUBLE_EQ(axistrain.data()[3], 0.19);
}