#include <Commands/ProceedFromPresetCommand.h>
#include <InputOutput/Parameter.h>

#include <string>
#include <iostream>

using namespace dem;
void
ProceedFromPresetCommand::execute(Assembly* assembly)
{

  std::string boundFile =
    dem::Parameter::getSingleton().datafile["boundaryFile"];
  std::string partFile =
    dem::Parameter::getSingleton().datafile["particleFile"];
  std::cout << "Called ProceedFromPresetCommand with: "
            << " BoundaryFile = " << boundFile << " and "
            << " ParticleFile = " << partFile << std::endl;
  assembly->deposit(boundFile.c_str(), partFile.c_str());

}
