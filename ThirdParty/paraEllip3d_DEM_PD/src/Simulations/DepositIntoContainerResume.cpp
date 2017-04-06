#include <Simulations/DepositIntoContainerResume.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;

void
DepositIntoContainerResume::execute(Assembly* assembly)
{
  assembly->deposit(
    dem::Parameter::getSingleton().datafile["boundaryFile"].c_str(),
    dem::Parameter::getSingleton().datafile["particleFile"].c_str());

  if (assembly->getMPIRank() == 0) {
    const Box& allContainer = assembly->getAllContainer();
    assembly->setContainer(Box(
      allContainer.getMinCorner().x(), allContainer.getMinCorner().y(),
      allContainer.getMinCorner().z(), allContainer.getMaxCorner().x(),
      allContainer.getMaxCorner().y(),
      dem::Parameter::getSingleton().parameter["trimHeight"]));
    assembly->buildBoundary(6, "trim_boundary_ini");
    std::size_t endSnap = static_cast<std::size_t>(
      dem::Parameter::getSingleton().parameter["endSnap"]);
    assembly->trim(
      false, combine( "deposit_particle_", endSnap, 3),
      "trim_particle_ini");
  }
}
