#include <Simulations/DEM/ProceedFromPreset.h>
#include <Boundary/BoundaryReader.h>
#include <InputOutput/InputParameter.h>

#include <iostream>
#include <string>

using namespace dem;
void
ProceedFromPreset::execute(DiscreteElements* dem)
{
  std::string boundaryFilename =
    InputParameter::get().datafile["boundaryFilename"];
  std::string particleFilename =
    InputParameter::get().datafile["particleFilename"];
  //std::cout << "Called ProceedFromPresetCommand with: "
  //          << " BoundaryFile = " << boundaryFilename << " and "
  //          << " ParticleFile = " << particleFilename << std::endl;
  dem->deposit(boundaryFilename, particleFilename);
}
