#include <Simulations/DEM/TriaxialLoading.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;
void
TriaxialLoading::execute(DiscreteElements* dem)
{
  std::ofstream progressInf;

  if (dem->getMPIRank() == 0) {
    dem->readBoundary(
      util::getFilename("boundaryFilename"));
    dem->readParticles(
      util::getFilename("particleFilename"));
    dem->openProgressOutputFile(progressInf, "triaxial_progress");
  }
  dem->scatterParticles();

  auto startStep = util::getParam<std::size_t>("startStep");
  auto endStep = util::getParam<std::size_t>("endStep");
  auto startSnap = util::getParam<std::size_t>("startSnap");
  auto endSnap = util::getParam<std::size_t>("endSnap");
  std::size_t netStep = endStep - startStep + 1;
  std::size_t netSnap = endSnap - startSnap + 1;
  REAL sigmaConf = util::getParam<REAL>("sigmaConf");
  auto timeStep = util::getParam<REAL>("timeStep");

  REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
  auto iteration = startStep;
  std::size_t iterSnap = startSnap;
  REAL distX, distY, distZ;
  std::string outputFolder(".");
  if (dem->getMPIRank() == 0) {

    // Create the output writer in the master process
    // <outputFolder> triaxial.pe3d </outputFolder>
    auto folderName =  dem::InputParameter::get().datafile["outputFolder"];
    outputFolder = util::createOutputFolder(folderName);
    //std::cout << "Output folder = " << outputFolder << "\n";
    dem->createOutputWriter(outputFolder, iterSnap-1);

    dem->writeBoundaryToFile();
    dem->writePatchGridToFile();
    dem->writeParticlesToFile(iterSnap);
    dem->printBoundaryContacts();
    dem->printBoundary();
    dem->getStartDimension(distX, distY, distZ);
  }
  if (dem->getMPIRank() == 0)
    debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT"
             << std::setw(OWID) << "migraT" << std::setw(OWID) << "totalT"
             << std::setw(OWID) << "overhead%" << std::endl;

  // Broadcast the output folder to all processes
  broadcast(dem->getMPIWorld(), outputFolder, 0);

  while (iteration <= endStep) {
    commuT = migraT = gatherT = totalT = 0;
    time0 = MPI_Wtime();
    dem->communicateGhostParticles();
    time2 = MPI_Wtime();
    commuT = time2 - time0;

    // displacement control relies on constant time step, so do not call
    // calcTimeStep().
    // calcTimeStep(); // use values from last step, must call before findConact
    dem->findContact(iteration);
    if (dem->isBoundaryProcess())
      dem->findBoundaryContacts(iteration);

    dem->initializeForces();
    dem->internalForce(timeStep, iteration);
    if (dem->isBoundaryProcess())
      dem->boundaryForce(timeStep, iteration);

    dem->updateParticles(timeStep, iteration);
    dem->gatherBoundaryContacts(); // must call before updateBoundary
    dem->updateBoundary(sigmaConf, "triaxial");
    dem->updatePatchBox();

    if (iteration % (netStep / netSnap) == 0) {
      time1 = MPI_Wtime();
      dem->gatherParticles();
      dem->gatherEnergy();
      time2 = MPI_Wtime();
      gatherT = time2 - time1;

      if (dem->getMPIRank() == 0) {
        dem->updateFileNames(iterSnap);
        dem->writeBoundaryToFile();
        dem->writePatchGridToFile();
        dem->writeParticlesToFile(iterSnap);
        dem->printBoundaryContacts();
        dem->printBoundary();
        // dem->appendToProgressOutputFile(progressInf, distX, distY, distZ); //
        // redundant
      }
      dem->printContact(combine(outputFolder, "triaxial_contact_", iterSnap, 3));
      ++iterSnap;
    }

    dem->releaseReceivedParticles(); // late release because
                                     // dem->printContact refers to
                                     // received particles
    time1 = MPI_Wtime();
    dem->migrateParticles(iteration);
    time2 = MPI_Wtime();
    migraT = time2 - time1;
    totalT = time2 - time0;
    if (dem->getMPIRank() == 0 &&
        (iteration + 1) % (netStep / netSnap) ==
          0) // ignore gather and dem->print time at this step
      debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT
               << std::setw(OWID) << migraT << std::setw(OWID) << totalT
               << std::setw(OWID) << (commuT + migraT) / totalT * 100
               << std::endl;

    if (dem->getMPIRank() == 0 && iteration % 10 == 0)
      dem->appendToProgressOutputFile(progressInf, timeStep, distX, distY, distZ);

    // no break condition, just through top/bottom displacement control
    ++iteration;
  }

  if (dem->getMPIRank() == 0) {
    dem->updateFileNames(iterSnap, ".end");
    dem->writeParticlesToFile(iterSnap);
    dem->printBoundaryContacts();
    dem->printBoundary();
    dem->appendToProgressOutputFile(progressInf, timeStep, distX, distY, distZ);
  }

  if (dem->getMPIRank() == 0)
    dem->closeProgressOutputFile(progressInf);
}
