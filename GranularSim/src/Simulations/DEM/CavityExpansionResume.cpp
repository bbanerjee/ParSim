#include <Simulations/DEM/CavityExpansionResume.h>

using namespace dem;
void
CavityExpansionResume::execute(DiscreteElements* dem)
{
  dem->allowPatchDomainResize(Boundary::BoundaryID::ZPLUS);

  dem->deposit(
    util::getFilename("boundaryFilename").c_str(),
    util::getFilename("particleFilename").c_str());
}
