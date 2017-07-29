#include <Simulations/CavityExpansionResume.h>

using namespace dem;
void
CavityExpansionResume::execute(DiscreteElements* dem)
{
  dem->deposit(
    Parameter::get().datafile["boundaryFile"].c_str(),
    Parameter::get().datafile["particleFile"].c_str());
}
