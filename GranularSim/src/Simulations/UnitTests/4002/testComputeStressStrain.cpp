#include <Simulations/DEM/PeriodicBCComputeStressStrain.h>
#include <Core/Util/Utility.h>
#include <InputOutput/InputParameter.h>
#include <InputOutput/OutputTecplot.h>
#include <gtest/gtest.h>

using namespace dem;

TEST(DEM4002Test, computeStressStrain) {
  DiscreteElements dem;
  dem::InputParameter::get().addFilename("inputDataDirectory", 
                                         "axisymmetric_strain_pb_seven.000");

  dem::InputParameter::get().addParameter("young", 1.0e9);
  dem::InputParameter::get().addParameter("poisson", 0.3);
  dem::InputParameter::get().addParameter("specificG", 1.5);
  dem::InputParameter::get().addParameter("gravAccel", 9.81);
  dem::InputParameter::get().addParameter("gravScale", 1);

  PeriodicBCComputeStressStrain command;
  command.execute(&dem);
}

