#include <Commands/TuneMassPercentageCommand.h>

using namespace dem;


// input:   number percentage smaller from data file
// output:  mass percentage smaller to disk file debugInf
// purpose: let mass percentage smaller satisfy particle size distribution curve
// method:  use trial and error method on number percentage until mass
// percentage is satisfied
void
TuneMassPercentageCommand::execute(Assembly* assembly)
{
  if (assembly->getMPIRank() == 0) {
    REAL minX = dem::Parameter::getSingleton().parameter["minX"];
    REAL minY = dem::Parameter::getSingleton().parameter["minY"];
    REAL minZ = dem::Parameter::getSingleton().parameter["minZ"];
    REAL maxX = dem::Parameter::getSingleton().parameter["maxX"];
    REAL maxY = dem::Parameter::getSingleton().parameter["maxY"];
    REAL maxZ = dem::Parameter::getSingleton().parameter["maxZ"];
    std::size_t particleLayers =
      dem::Parameter::getSingleton().parameter["particleLayers"];

    assembly->setContainer(Box(minX, minY, minZ, maxX, maxY, maxZ));

    assembly->buildBoundary(5, "deposit_boundary_ini");

    std::size_t sieveNum = static_cast<std::size_t>(
      dem::Parameter::getSingleton().parameter["sieveNum"]);
    std::vector<REAL> percent(sieveNum), size(sieveNum);
    std::vector<std::pair<REAL, REAL>>& grada =
      dem::Parameter::getSingleton().gradation;
    assert(grada.size() == sieveNum);
    for (std::size_t i = 0; i < sieveNum; ++i) {
      percent[i] = grada[i].first;
      size[i] = grada[i].second;
    }
    REAL ratioBA = dem::Parameter::getSingleton().parameter["ratioBA"];
    REAL ratioCA = dem::Parameter::getSingleton().parameter["ratioCA"];
    assembly->setGradation(Gradation(sieveNum, percent, size, ratioBA, ratioCA));

    assembly->generateParticle(particleLayers, "float_particle_ini");

    // statistics of mass distribution
    Gradation massGrad = assembly->getGradation();
    std::vector<REAL>& massPercent = massGrad.getPercent();
    std::vector<REAL>& massSize = massGrad.getSize();
    for (double & i : massPercent)
      i = 0;

    for (const auto & itr : assembly->getAllParticleVec())
      for (int i = massPercent.size() - 1; i >= 0;
           --i) { // do not use size_t for descending series
        if (itr->getA() <= massSize[i])
          massPercent[i] += itr->getMass();
      }
    REAL totalMass = massPercent[0];
    for (double & i : massPercent)
      i /= totalMass;
    debugInf << std::endl
             << "mass percentage of particles:" << std::endl
             << std::setw(OWID) << massPercent.size() << std::endl;
    for (std::size_t i = 0; i < massPercent.size(); ++i)
      debugInf << std::setw(OWID) << massPercent[i] << std::setw(OWID)
               << massSize[i] << std::endl;
  }
}
