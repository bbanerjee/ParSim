#include <Simulations/CavityExpansionResume.h>

using namespace dem;
void
CavityExpansionResume::execute(DiscreteElements* dem)
{
  dem->deposit(
    InputParameter::get().datafile["boundaryFile"].c_str(),
    InputParameter::get().datafile["particleFile"].c_str());
}
