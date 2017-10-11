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
  if (dem->getMPIRank() == 0) {
    std::string boundaryFilename =
      InputParameter::get().datafile["boundaryFilename"];
    dem->readBoundary(boundaryFilename);

    //REAL minX = util::getParam<REAL>("minX");
    //REAL minY = util::getParam<REAL>("minY");
    //REAL minZ = util::getParam<REAL>("minZ");
    //REAL maxX = util::getParam<REAL>("maxX");
    //REAL maxY = util::getParam<REAL>("maxY");
    //REAL maxZ = util::getParam<REAL>("maxZ");
    //dem->setContainer(Box(minX, minY, minZ, maxX, maxY, maxZ));
    Box container = dem->getAllContainer();

    BoundaryFileWriter boundaryWriter;
    boundaryWriter.writeXML(5, "deposit_boundary_ini.xml", container);
    boundaryWriter.writeCSV(5, "deposit_boundary_ini.txt", container);

    Gradation gradation;
    gradation.initializeFromInputParameters();
    dem->setGradation(gradation);

    auto layerFlag = util::getParam<std::size_t>("particleLayers");
    DEMParticle::DEMParticleShape particleType = 
      DEMParticle::DEMParticleShape::ELLIPSOID;

    DEMParticleCreator creator;
    DEMParticlePArray particles = 
      creator.generateDEMParticles(layerFlag, particleType, container, gradation);
    dem->setAllDEMParticleVec(particles);

    DEMParticleFileWriter writer;
    writer.writeCSV(particles, gradation, "float_particle_ini.csv");
    writer.writeXML(particles, gradation, "float_particle_ini.xml");
  }

  dem->deposit("deposit_boundary_ini.xml", "float_particle_ini.xml");

  if (dem->getMPIRank() == 0) {
    const Box allContainer = dem->getAllContainer();
    Box trimmedContainer(allContainer.getMinCorner().x(), 
      allContainer.getMinCorner().y(), allContainer.getMinCorner().z(), 
      allContainer.getMaxCorner().x(), allContainer.getMaxCorner().y(),
      util::getParam<REAL>("trimHeight"));

    dem->setContainer(trimmedContainer);

    BoundaryFileWriter boundaryWriter;
    boundaryWriter.writeXML(6, "trim_boundary_ini.xml", trimmedContainer);
    boundaryWriter.writeCSV(6, "trim_boundary_ini.txt", trimmedContainer);

    dem->trim(false, "dummy", "trim_particle_ini");
  }
}
