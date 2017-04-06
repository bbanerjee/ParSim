#include <Simulations/CavityExpansionResume.h>

using namespace dem;
void
CavityExpansionResume::execute(Assembly* assembly)
{
  assembly->deposit(
    Parameter::get().datafile["boundaryFile"].c_str(),
    Parameter::get().datafile["particleFile"].c_str());
}
