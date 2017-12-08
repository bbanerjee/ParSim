#include <Simulations/DEM/PeriodicBCComputeStressStrain.h>
#include <Core/Util/Utility.h>
#include <InputOutput/InputParameter.h>
#include <InputOutput/OutputTecplot.h>
#include <gtest/gtest.h>

using namespace dem;

TEST(DEM4002Test, computeStressStrain) {
  DiscreteElements dem;
  dem::InputParameter::get().addFilename("inputDataDirectory", "axisymmetric_strain_pb_two.000");

  PeriodicBCComputeStressStrain command;
  command.execute(&dem);
}

