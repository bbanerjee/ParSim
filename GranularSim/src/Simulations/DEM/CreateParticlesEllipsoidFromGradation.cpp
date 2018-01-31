#include <Simulations/DEM/CreateParticlesEllipsoidFromGradation.h>
#include <DiscreteElements/DEMParticleCreator.h>
#include <InputOutput/DEMParticleFileWriter.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;

void
CreateParticlesEllipsoidFromGradation::execute(DiscreteElements* dem)
{
  //std::string particleFile = "generated_spherical_particles.xml";
  //std::string boundaryFile = "generated_boundary_file.xml";

  auto boundaryFile = util::getFilename("generatedBoundaryOutputFilename");
  auto particleFile = util::getFilename("generatedParticleOutputFilename");

  if (dem->getMPIRank() == 0) {

    // Set up gradation
    Gradation gradation;
    gradation.initializeFromInputParameters();
    dem->setGradation(gradation);

    // Create and save particle and boundaries
    dem->createAndSaveParticlesAndBoundaries(boundaryFile, particleFile);
  }

}
