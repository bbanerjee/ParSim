#include <Simulations/TrimParticles.h>

using namespace dem;
void
TrimParticles::execute(Assembly* assembly)
{
  if (assembly->getMPIRank() == 0) {
    assembly->readBoundary(
      Parameter::get().datafile["boundaryFile"].c_str());
    assembly->trim(
      true, Parameter::get().datafile["particleFile"].c_str(),
      "trim_particle_end");
  }
}
