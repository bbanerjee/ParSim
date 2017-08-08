#include <Simulations/DEM/CavityExpansion.h>
#include <Core/Util/Utility.h>

using namespace dem;

void
CavityExpansion::execute(DiscreteElements* dem)
{
  if (dem->getMPIRank() == 0) {
    auto inputParticle = InputParameter::get().datafile["particleFile"];
    REAL percent = util::getParam<REAL>("expandPercent");
    dem->readParticles(inputParticle);

    REAL x1 = util::getParam<REAL>("cavityMinX");
    REAL y1 = util::getParam<REAL>("cavityMinY");
    REAL z1 = util::getParam<REAL>("cavityMinZ");
    REAL x2 = util::getParam<REAL>("cavityMaxX");
    REAL y2 = util::getParam<REAL>("cavityMaxY");
    REAL z2 = util::getParam<REAL>("cavityMaxZ");

    DEMParticlePArray cavityParticleVec;
    Vec center;

    for (const auto& it : dem->getAllDEMParticleVec()) {
      center = it->currentPosition();
      if (center.x() > x1 && center.x() < x2 && center.y() > y1 &&
          center.y() < y2 && center.z() > z1 && center.z() < z2) {
        it->expand(percent);
        cavityParticleVec.push_back(it);
      }
    }

    dem->printParticle("cavity_particle_ini", cavityParticleVec, 0);
    dem->printParticle("expand_particle_ini", 0);
  }

  dem->deposit(
    InputParameter::get().datafile["boundaryFile"],
    "expand_particle_ini");
}
