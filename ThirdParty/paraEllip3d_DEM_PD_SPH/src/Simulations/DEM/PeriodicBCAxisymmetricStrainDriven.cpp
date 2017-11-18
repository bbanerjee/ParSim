#include <Simulations/DEM/PeriodicBCAxisymmetricStrainDriven.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;

void
PeriodicBCAxisymmetricStrainDriven::execute(DiscreteElements* dem)
{
  // Read the boundary conditions in all processes
  // The file should be small enough
  std::string bcFile = util::getFilename("demBoundaryConditionFilename");
  dem->readBoundaryConditions(bcFile);

  std::ofstream progressInf;
  std::ofstream balancedInf;

  std::string outputFolder(".");
  REAL distX = 0, distY = 0, distZ = 0;

  if (dem->getMPIRank() == 0) {
    std::string boundaryFile = util::getFilename("boundaryFilename");
    std::string particleFile = util::getFilename("particleFilename");
    dem->readBoundary(boundaryFile);
    dem->readParticles(particleFile);

    // Create the output writer in the master process
    // <outputFolder> filename.pe3d </outputFolder>
    auto folderName =  dem::InputParameter::get().datafile["outputFolder"];
    outputFolder = util::createOutputFolder(folderName);
    dem->createOutputWriter(outputFolder, 0);

    dem->openProgressOutputFile(progressInf, folderName + ".progress");
    dem->openProgressOutputFile(balancedInf, folderName + ".balanced");

    dem->writeBoundaryToFile();
    dem->writePatchGridToFile();
    dem->writeParticlesToFile(0);
    dem->printBoundary();
    dem->printBoundaryContacts();
    dem->getStartDimension(distX, distY, distZ);
  }

  // Broadcast the output folder to all processes
  broadcast(dem->getMPIWorld(), outputFolder, 0);

  // Scatter the particles
  dem->scatterParticles();

  auto startStep = util::getParam<std::size_t>("startStep");
  auto endStep = util::getParam<std::size_t>("endStep");
  auto startSnap = util::getParam<std::size_t>("startSnap");
  auto endSnap = util::getParam<std::size_t>("endSnap");

  auto netStep = endStep - startStep + 1;
  auto netSnap = endSnap - startSnap + 1;

  auto iterSnap = startSnap;
  auto iteration = startStep;

  auto startTime = util::getParam<REAL>("startTime");
  auto endTime = util::getParam<REAL>("endTime");
  auto timeStep = util::getParam<REAL>("timeStep");

  auto curTime = startTime;

  while (iteration <= endStep && curTime < endTime) {

    // Communicate ghost particles to patches
    dem->communicateGhostParticles(iteration);

    // Apply the particle boundary conditions
    dem->applyPatchParticleBC(curTime);

    timeStep = dem->calcTimeStep(); 
    dem->findContact();
    if (dem->isBoundaryProcess()) {
      dem->findBoundaryContacts();
    }
    dem->clearContactForce();
    dem->internalForce();
    if (dem->isBoundaryProcess()) {
      dem->boundaryForce();
    }
    dem->updateParticles();
    dem->gatherBoundaryContacts(); 
    dem->updatePatchBox();

    if (iteration % (netStep / netSnap) == 0) {

      dem->gatherParticles();
      dem->gatherEnergy();

      if (dem->getMPIRank() == 0) {
        dem->updateFileNames(iterSnap);
        dem->writeBoundaryToFile();
        dem->printBoundary();
        dem->writePatchGridToFile();
        dem->writeParticlesToFile(iterSnap);
        dem->printBoundaryContacts();
        dem->appendToProgressOutputFile(progressInf, distX, distY, distZ);
      }

      dem->printContact(combine(outputFolder, "contact_", iterSnap, 3));
      ++iterSnap;
    }

    dem ->releaseReceivedParticles();
    dem->migrateParticles();

    ++iteration;
    curTime += timeStep;
  }

  if (dem->getMPIRank() == 0) {
    dem->closeProgressOutputFile(progressInf);
    dem->closeProgressOutputFile(balancedInf);
  }
}
