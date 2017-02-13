#include <Simulations/ProceedFromPreset.h>
#include <InputOutput/Parameter.h>
#include <Boundary/BoundaryReader.h>

#include <string>
#include <iostream>

using namespace dem;
void
ProceedFromPreset::execute(Assembly* assembly)
{

  std::string boundFile =
    dem::Parameter::getSingleton().datafile["boundaryFile"];
  std::string partFile =
    dem::Parameter::getSingleton().datafile["particleFile"];
  std::cout << "Called ProceedFromPresetCommand with: "
            << " BoundaryFile = " << boundFile << " and "
            << " ParticleFile = " << partFile << std::endl;
  //assembly->deposit(boundFile.c_str(), partFile.c_str());
  assembly->deposit(boundFile, partFile.c_str());

}
