#include <Simulations/DEM/DepositIntoContainer.h>
#include <DiscreteElements/DEMParticleCreator.h>
#include <Boundary/BoundaryFileWriter.h>
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
    dem->setContainer(Box(minX, minY, minZ, maxX, maxY, maxZ));
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

    dem->printParticle("float_particle_ini.xml", 0);
  }

  dem->deposit("deposit_boundary_ini", "float_particle_ini");

  if (dem->getMPIRank() == 0) {
    const Box allContainer = dem->getAllContainer();
    dem->setContainer(Box(
      allContainer.getMinCorner().x(), allContainer.getMinCorner().y(),
      allContainer.getMinCorner().z(), allContainer.getMaxCorner().x(),
      allContainer.getMaxCorner().y(),
      util::getParam<REAL>("trimHeight")));

    BoundaryFileWriter boundaryWriter;
    boundaryWriter.writeXML(6, "trim_boundary_ini.xml", dem->getAllContainer());
    boundaryWriter.writeCSV(6, "trim_boundary_ini.txt", dem->getAllContainer());

    auto endSnap = util::getParam<std::size_t>("endSnap");
    dem->trim(
      false, combine(".", "deposit_particle_", endSnap, 3),
      "trim_particle_ini");
  }
}
