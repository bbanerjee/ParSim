#include <Simulations/DEM/ProceedFromPreset.h>
#include <Boundary/BoundaryReader.h>
#include <InputOutput/InputParameter.h>

#include <iostream>
#include <string>

using namespace dem;
void
ProceedFromPreset::execute(DiscreteElements* dem)
{
  std::string boundFile =
    InputParameter::get().datafile["boundaryFile"];
  std::string partFile =
    InputParameter::get().datafile["particleFile"];
  //std::cout << "Called ProceedFromPresetCommand with: "
  //          << " BoundaryFile = " << boundFile << " and "
  //          << " ParticleFile = " << partFile << std::endl;
  // dem->deposit(boundFile.c_str(), partFile.c_str());
  dem->deposit(boundFile, partFile.c_str());
}
