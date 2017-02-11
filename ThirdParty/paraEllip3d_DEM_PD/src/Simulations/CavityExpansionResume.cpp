#include <Simulations/CavityExpansionResume.h>

using namespace dem;
void
CavityExpansionResume::execute(Assembly* assembly)
{
  assembly->deposit(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str(),
          dem::Parameter::getSingleton().datafile["particleFile"].c_str());
}
