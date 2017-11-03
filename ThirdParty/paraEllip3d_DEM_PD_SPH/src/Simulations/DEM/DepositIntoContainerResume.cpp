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
    const Box& spatialDomain = dem->getSpatialDomain();
    dem->setSpatialDomain(Box(
      spatialDomain.minCorner().x(), spatialDomain.minCorner().y(),
      spatialDomain.minCorner().z(), spatialDomain.maxCorner().x(),
      spatialDomain.maxCorner().y(),
      util::getParam<REAL>("trimHeight")));

    BoundaryFileWriter boundaryWriter;
    boundaryWriter.writeXML(6, "trim_boundary_ini.xml", dem->getSpatialDomain());
    boundaryWriter.writeCSV(6, "trim_boundary_ini.csv", dem->getSpatialDomain());

    dem->trim(false, "dummy", "trim_particle_ini");
  }
}
