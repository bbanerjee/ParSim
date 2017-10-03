#include <Simulations/DEM/DepositIntoContainer.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;

void
DepositIntoContainer::execute(DiscreteElements* dem)
{
  if (dem->getMPIRank() == 0) {
    REAL minX = util::getParam<REAL>("minX");
    REAL minY = util::getParam<REAL>("minY");
    REAL minZ = util::getParam<REAL>("minZ");
    REAL maxX = util::getParam<REAL>("maxX");
    REAL maxY = util::getParam<REAL>("maxY");
    REAL maxZ = util::getParam<REAL>("maxZ");
    auto particleLayers = util::getParam<std::size_t>("particleLayers");

    dem->setContainer(Box(minX, minY, minZ, maxX, maxY, maxZ));

    dem->buildBoundary(5, "deposit_boundary_ini");

    auto sieveNum = util::getParam<std::size_t>("sieveNum");
    std::vector<REAL> percent(sieveNum), size(sieveNum);
    std::vector<std::pair<REAL, REAL>>& grada =
      InputParameter::get().gradation;
    assert(grada.size() == sieveNum);
    for (std::size_t i = 0; i < sieveNum; ++i) {
      percent[i] = grada[i].first;
      size[i] = grada[i].second;
    }
    REAL ratioBA = util::getParam<REAL>("ratioBA");
    REAL ratioCA = util::getParam<REAL>("ratioCA");
    dem->setGradation(
      Gradation(sieveNum, percent, size, ratioBA, ratioCA));

    dem->generateParticle(particleLayers, "float_particle_ini");
  }

  dem->deposit("deposit_boundary_ini", "float_particle_ini");

  if (dem->getMPIRank() == 0) {
    const Box allContainer = dem->getAllContainer();
    dem->setContainer(Box(
      allContainer.getMinCorner().x(), allContainer.getMinCorner().y(),
      allContainer.getMinCorner().z(), allContainer.getMaxCorner().x(),
      allContainer.getMaxCorner().y(),
      util::getParam<REAL>("trimHeight")));
    dem->buildBoundary(6, "trim_boundary_ini");
    auto endSnap = util::getParam<std::size_t>("endSnap");
    dem->trim(
      false, combine(".", "deposit_particle_", endSnap, 3),
      "trim_particle_ini");
  }
}
