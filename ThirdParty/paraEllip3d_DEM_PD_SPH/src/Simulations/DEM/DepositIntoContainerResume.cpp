#include <Simulations/DEM/DepositIntoContainerResume.h>
#include <Boundary/BoundaryFileWriter.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;

void
DepositIntoContainerResume::execute(DiscreteElements* dem)
{
  dem->deposit(
    InputParameter::get().datafile["boundaryFilename"],
    InputParameter::get().datafile["particleFilename"]);

  if (dem->getMPIRank() == 0) {
    const Box& allContainer = dem->getAllContainer();
    dem->setContainer(Box(
      allContainer.getMinCorner().x(), allContainer.getMinCorner().y(),
      allContainer.getMinCorner().z(), allContainer.getMaxCorner().x(),
      allContainer.getMaxCorner().y(),
      util::getParam<REAL>("trimHeight")));

    BoundaryFileWriter boundaryWriter;
    boundaryWriter.writeXML(6, "trim_boundary_ini.xml", dem->getAllContainer());
    boundaryWriter.writeCSV(6, "trim_boundary_ini.csv", dem->getAllContainer());

    auto endSnap = util::getParam<std::size_t>("endSnap");
    dem->trim(
      false, combine(".", "deposit_particle_", endSnap, 3),
      "trim_particle_ini");
  }
}
