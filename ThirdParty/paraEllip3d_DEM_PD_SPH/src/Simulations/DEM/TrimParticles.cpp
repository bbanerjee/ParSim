#include <Simulations/DEM/TrimParticles.h>

using namespace dem;
void
TrimParticles::execute(DiscreteElements* dem)
{
  if (dem->getMPIRank() == 0) {
    dem->allowPatchDomainResize(Boundary::BoundaryID::ZPLUS);
    dem->readBoundary( util::getFilename("boundaryFilename").c_str());
    dem->trim( true, util::getFilename("particleFilename").c_str(),
      "trim_particle_end");
  }
}
