#include <Simulations/DepositIntoContainer.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;

void
DepositIntoContainer::execute(Assembly* assembly)
{
  if (assembly->getMPIRank() == 0) {
    REAL minX = util::getParam<REAL>("minX");
    REAL minY = util::getParam<REAL>("minY");
    REAL minZ = util::getParam<REAL>("minZ");
    REAL maxX = util::getParam<REAL>("maxX");
    REAL maxY = util::getParam<REAL>("maxY");
    REAL maxZ = util::getParam<REAL>("maxZ");
    auto particleLayers = util::getParam<size_t>("particleLayers");

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
  }

  assembly->deposit("deposit_boundary_ini", "float_particle_ini");

  if (assembly->getMPIRank() == 0) {
    const Box allContainer = assembly->getAllContainer();
    assembly->setContainer(Box(
      allContainer.getMinCorner().x(), allContainer.getMinCorner().y(),
      allContainer.getMinCorner().z(), allContainer.getMaxCorner().x(),
      allContainer.getMaxCorner().y(),
      util::getParam<REAL>("trimHeight")));
    assembly->buildBoundary(6, "trim_boundary_ini");
    auto endSnap = util::getParam<std::size_t>("endSnap");
    assembly->trim(
      false, combine(".", "deposit_particle_", endSnap, 3),
      "trim_particle_ini");
  }
}
