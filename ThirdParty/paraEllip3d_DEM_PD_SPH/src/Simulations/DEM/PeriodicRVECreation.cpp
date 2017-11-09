#include <Simulations/DEM/PeriodicRVECreation.h>
#include <DiscreteElements/DEMParticleCreator.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;

void
PeriodicRVECreation::execute(DiscreteElements* dem)
{
  if (dem->getMPIRank() == 0) {

    // Read the particles and boundary from file
    std::string boundaryFile = InputParameter::get().datafile["boundaryFilename"];
    std::string particleFile = InputParameter::get().datafile["particleFilename"];
    dem->readBoundary(boundaryFile);
    dem->readParticles(particleFile);

    // Create the output writer in the master process
    // <outputFolder> foldername </outputFolder>
    std::string outputFolder(".");
    auto folderName =  dem::InputParameter::get().datafile["outputFolder"];
    outputFolder = util::createOutputFolder(folderName);
    dem->createOutputWriter(outputFolder, 0);

    auto& particles = dem->getModifiableAllParticleVec(); 
    DEMParticleCreator creator;
    DEMParticlePArray periodic = 
      creator.generatePeriodicDEMParticles(particles, dem->getSpatialDomain(), 1, 5);
  
    particles.insert(particles.end(), periodic.begin(), periodic.end());
    std::cout << "periodic = " << periodic.size() << "\n";
    std::cout << "total particles = " << particles.size() << "\n";

    dem->writeBoundaryToFile();
    dem->writeParticlesToFile(0);
    dem->printParticle("generated_particles.csv", 0);
  }
}
