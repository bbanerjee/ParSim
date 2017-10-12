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

    std::string boundaryFile = "deposit_boundary_ini.xml";
    std::string particleFile = "float_particle_ini.xml";

    // Set up gradation
    Gradation gradation;
    gradation.initializeFromInputParameters();
    dem->setGradation(gradation);

    // Create and save particle and boundaries
    dem->createAndSaveParticlesAndBoundaries(boundaryFile, particleFile);

    // statistics of mass distribution
    Gradation massGrad = dem->getGradation();
    std::vector<REAL>& massPercent = massGrad.getPercent();
    std::vector<REAL>& massSize = massGrad.getSize();
    for (double& i : massPercent)
      i = 0;

    for (const auto& particle : dem->getAllDEMParticleVec()) {
      for (int i = massPercent.size() - 1; i >= 0; --i) { 
        if (particle->getA() <= massSize[i])
          massPercent[i] += particle->getMass();
      }
    }

    REAL totalMass = massPercent[0];
    for (double& i : massPercent) {
      i /= totalMass;
    }
    debugInf << std::endl
             << "mass percentage of particles:" << std::endl
             << std::setw(OWID) << massPercent.size() << std::endl;
    for (std::size_t i = 0; i < massPercent.size(); ++i) {
      debugInf << std::setw(OWID) << massPercent[i] << std::setw(OWID)
               << massSize[i] << std::endl;
    }
  }
}
