#include <Simulations/DEM/TrimParticles.h>

using namespace dem;
void
TrimParticles::execute(DiscreteElements* dem)
{
  if (dem->getMPIRank() == 0) {
    dem->readBoundary(
      InputParameter::get().datafile["boundaryFilename"].c_str());
    dem->trim(
      true, InputParameter::get().datafile["particleFilename"].c_str(),
      "trim_particle_end");
  }
}
