#include <Commands/CavityExpansionCommand.h>

using namespace dem;
void
CavityExpansionCommand::execute(Assembly* assembly)
{
  if (assembly->getMPIRank() == 0) {
    const char* inputParticle =
      dem::Parameter::getSingleton().datafile["particleFile"].c_str();
    REAL percent = dem::Parameter::getSingleton().parameter["expandPercent"];
    assembly->readParticle(inputParticle);

    REAL x1 = dem::Parameter::getSingleton().parameter["cavityMinX"];
    REAL y1 = dem::Parameter::getSingleton().parameter["cavityMinY"];
    REAL z1 = dem::Parameter::getSingleton().parameter["cavityMinZ"];
    REAL x2 = dem::Parameter::getSingleton().parameter["cavityMaxX"];
    REAL y2 = dem::Parameter::getSingleton().parameter["cavityMaxY"];
    REAL z2 = dem::Parameter::getSingleton().parameter["cavityMaxZ"];

    ParticlePArray cavityParticleVec;
    Vec center;

    for (const auto & it : assembly->getAllParticleVec()) {
      center = it->getCurrPos();
      if (center.getX() > x1 && center.getX() < x2 && center.getY() > y1 &&
          center.getY() < y2 && center.getZ() > z1 && center.getZ() < z2) {
        it->expand(percent);
        cavityParticleVec.push_back(it);
      }
    }

    assembly->printParticle("cavity_particle_ini", cavityParticleVec);
    assembly->printParticle("expand_particle_ini");
  }

  assembly->deposit(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str(),
          "expand_particle_ini");
}
