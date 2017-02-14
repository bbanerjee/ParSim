#include <Simulations/DepositIntoContainer.h>

using namespace dem;

void
DepositIntoContainer::execute(Assembly* assembly)
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
    assembly->setGradation(
      Gradation(sieveNum, percent, size, ratioBA, ratioCA));

    assembly->generateParticle(particleLayers, "float_particle_ini");
  }

  assembly->deposit("deposit_boundary_ini", "float_particle_ini");

  if (assembly->getMPIRank() == 0) {
    const Box allContainer = assembly->getAllContainer();
    assembly->setContainer(Box(
      allContainer.getMinCorner().getX(), allContainer.getMinCorner().getY(),
      allContainer.getMinCorner().getZ(), allContainer.getMaxCorner().getX(),
      allContainer.getMaxCorner().getY(),
      dem::Parameter::getSingleton().parameter["trimHeight"]));
    assembly->buildBoundary(6, "trim_boundary_ini");
    char cstr[50];
    std::size_t endSnap = static_cast<std::size_t>(
      dem::Parameter::getSingleton().parameter["endSnap"]);
    assembly->trim(
      false, Assembly::combineString(cstr, "deposit_particle_", endSnap, 3),
      "trim_particle_ini");
  }
}
