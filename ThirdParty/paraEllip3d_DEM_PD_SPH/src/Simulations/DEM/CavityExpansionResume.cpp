#include <Simulations/DEM/CavityExpansionResume.h>

using namespace dem;
void
CavityExpansionResume::execute(DiscreteElements* dem)
{
  dem->deposit(
    InputParameter::get().datafile["boundaryFilename"].c_str(),
    InputParameter::get().datafile["particleFilename"].c_str());
}
