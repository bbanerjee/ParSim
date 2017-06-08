#include <Simulations/DepositIntoContainerResume.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;

void
DepositIntoContainerResume::execute(Assembly* assembly)
{
  assembly->deposit(
    Parameter::get().datafile["boundaryFile"],
    Parameter::get().datafile["particleFile"]);

  if (assembly->getMPIRank() == 0) {
    const Box& allContainer = assembly->getAllContainer();
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
