#include <Simulations/DEM/DepositIntoContainer.h>
#include <DiscreteElements/DEMParticleCreator.h>
#include <InputOutput/DEMParticleFileWriter.h>
#include <Boundary/BoundaryFileWriter.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;


void
DepositIntoContainer::execute(DiscreteElements* dem)
{
  std::string particleFile = "float_particle_ini.xml";
  std::string boundaryFile = "deposit_boundary_ini.xml";

  if (dem->getMPIRank() == 0) {

    // Set up gradation
    Gradation gradation;
    gradation.initializeFromInputParameters();
    dem->setGradation(gradation);

    // Create and save particle and boundaries
    dem->createAndSaveParticlesAndBoundaries(boundaryFile, particleFile);
  }

  // Do gravity deposition
  dem->deposit(boundaryFile, particleFile);

  // Reduce the size of the domain for further calculations
  if (dem->getMPIRank() == 0) {
    const Box allContainer = dem->getAllContainer();
    Box trimmedContainer(allContainer.getMinCorner().x(), 
      allContainer.getMinCorner().y(), allContainer.getMinCorner().z(), 
      allContainer.getMaxCorner().x(), allContainer.getMaxCorner().y(),
      util::getParam<REAL>("trimHeight"));

    dem->setContainer(trimmedContainer);

    BoundaryFileWriter boundaryWriter;
    boundaryWriter.writeXML(6, "trim_boundary_ini.xml", trimmedContainer);
    boundaryWriter.writeCSV(6, "trim_boundary_ini.csv", trimmedContainer);

    dem->trim(false, "dummy", "trim_particle_ini");
  }
}
