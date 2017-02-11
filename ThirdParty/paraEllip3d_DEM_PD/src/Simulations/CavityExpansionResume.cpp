#include <Commands/CavityExpansionResumeCommand.h>

using namespace dem;
void
CavityExpansionResumeCommand::execute(Assembly* assembly)
{
  assembly->deposit(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str(),
          dem::Parameter::getSingleton().datafile["particleFile"].c_str());
}
