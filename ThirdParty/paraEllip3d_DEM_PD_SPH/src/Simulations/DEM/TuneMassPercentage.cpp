#include <Simulations/DEM/TuneMassPercentage.h>
#include <DiscreteElements/DEMParticleCreator.h>
#include <InputOutput/DEMParticleFileWriter.h>
#include <Boundary/BoundaryFileWriter.h>
#include <Core/Util/Utility.h>

using namespace dem;

// input:   number percentage smaller from data file
// output:  mass percentage smaller to disk file debugInf
// purpose: let mass percentage smaller satisfy particle size distribution curve
// method:  use trial and error method on number percentage until mass
// percentage is satisfied
void
TuneMassPercentage::execute(DiscreteElements* dem)
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
    boundaryWriter.writeXML(5, "deposit_boundary_ini.xml", dem->getAllContainer());
    boundaryWriter.writeCSV(5, "deposit_boundary_ini.txt", dem->getAllContainer());

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

    // statistics of mass distribution
    Gradation massGrad = dem->getGradation();
    std::vector<REAL>& massPercent = massGrad.getPercent();
    std::vector<REAL>& massSize = massGrad.getSize();
    for (double& i : massPercent)
      i = 0;

    for (const auto& itr : dem->getAllDEMParticleVec())
      for (int i = massPercent.size() - 1; i >= 0;
           --i) { // do not use size_t for descending series
        if (itr->getA() <= massSize[i])
          massPercent[i] += itr->getMass();
      }
    REAL totalMass = massPercent[0];
    for (double& i : massPercent)
      i /= totalMass;
    debugInf << std::endl
             << "mass percentage of particles:" << std::endl
             << std::setw(OWID) << massPercent.size() << std::endl;
    for (std::size_t i = 0; i < massPercent.size(); ++i)
      debugInf << std::setw(OWID) << massPercent[i] << std::setw(OWID)
               << massSize[i] << std::endl;
  }
}
