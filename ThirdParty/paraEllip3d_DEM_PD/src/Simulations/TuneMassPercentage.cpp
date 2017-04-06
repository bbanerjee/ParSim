#include <Simulations/TuneMassPercentage.h>
#include <Core/Util/Utility.h>

using namespace dem;

// input:   number percentage smaller from data file
// output:  mass percentage smaller to disk file debugInf
// purpose: let mass percentage smaller satisfy particle size distribution curve
// method:  use trial and error method on number percentage until mass
// percentage is satisfied
void
TuneMassPercentage::execute(Assembly* assembly)
{
  if (assembly->getMPIRank() == 0) {
    REAL minX = util::getParam<REAL>("minX");
    REAL minY = util::getParam<REAL>("minY");
    REAL minZ = util::getParam<REAL>("minZ");
    REAL maxX = util::getParam<REAL>("maxX");
    REAL maxY = util::getParam<REAL>("maxY");
    REAL maxZ = util::getParam<REAL>("maxZ");
    auto particleLayers = util::getParam<std::size_t>("particleLayers");

    assembly->setContainer(Box(minX, minY, minZ, maxX, maxY, maxZ));

    assembly->buildBoundary(5, "deposit_boundary_ini");

    auto sieveNum = util::getParam<std::size_t>("sieveNum");
    std::vector<REAL> percent(sieveNum), size(sieveNum);
    std::vector<std::pair<REAL, REAL>>& grada =
      Parameter::get().gradation;
    assert(grada.size() == sieveNum);
    for (std::size_t i = 0; i < sieveNum; ++i) {
      percent[i] = grada[i].first;
      size[i] = grada[i].second;
    }
    REAL ratioBA = util::getParam<REAL>("ratioBA");
    REAL ratioCA = util::getParam<REAL>("ratioCA");
    assembly->setGradation(
      Gradation(sieveNum, percent, size, ratioBA, ratioCA));

    assembly->generateParticle(particleLayers, "float_particle_ini");

    // statistics of mass distribution
    Gradation massGrad = assembly->getGradation();
    std::vector<REAL>& massPercent = massGrad.getPercent();
    std::vector<REAL>& massSize = massGrad.getSize();
    for (double& i : massPercent)
      i = 0;

    for (const auto& itr : assembly->getAllParticleVec())
      for (int i = massPercent.size() - 1; i >= 0;
           --i) { // do not use size_t for descending series
        if (itr->getA() <= massSize[i])
          massPercent[i] += itr->getMass();
      }
    REAL totalMass = massPercent[0];
    for (double& i : massPercent)
      i /= totalMass;
    debugInf << std::endl
             << "mass percentage of particles:" << std::endl
             << std::setw(OWID) << massPercent.size() << std::endl;
    for (std::size_t i = 0; i < massPercent.size(); ++i)
      debugInf << std::setw(OWID) << massPercent[i] << std::setw(OWID)
               << massSize[i] << std::endl;
  }
}
