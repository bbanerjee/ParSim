#include <Simulations/TrimParticles.h>

using namespace dem;
void
TrimParticles::execute(Assembly* assembly)
{
  if (assembly->getMPIRank() == 0) {
    assembly->readBoundary(
      dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
    assembly->trim(
      true, dem::Parameter::getSingleton().datafile["particleFile"].c_str(),
      "trim_particle_end");
  }
}
