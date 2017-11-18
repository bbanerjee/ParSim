#include <Simulations/DEM/TractionLoading.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;
void
TractionLoading::execute(DiscreteElements* dem)
{
  std::ofstream progressInf;
  std::ofstream balancedInf;

  if (dem->getMPIRank() == 0) {
    std::string boundaryFile = InputParameter::get().datafile["boundaryFilename"];
    std::string particleFile = InputParameter::get().datafile["particleFilename"];
    dem->readBoundary(boundaryFile);
    dem->readParticles(particleFile);
    dem->openProgressOutputFile(progressInf, "traction_progress");
    dem->openProgressOutputFile(balancedInf, "traction_balanced");
  }

  dem->scatterParticles();

  auto startStep = util::getParam<std::size_t>("startStep");
  auto endStep = util::getParam<std::size_t>("endStep");
  auto startSnapshot = util::getParam<std::size_t>("startSnap");
  auto endSnapshot = util::getParam<std::size_t>("endSnap");

  auto numSteps = endStep - startStep + 1;
  auto numSnapshots = endSnapshot - startSnapshot + 1;

  auto snapshot = startSnapshot;

  std::string outputFolder(".");
  REAL distX = 0, distY = 0, distZ = 0;
  if (dem->getMPIRank() == 0) {
    // Create the output writer in the master process
    // <outputFolder> isotropic.pe3d </outputFolder>
    auto folderName =  dem::InputParameter::get().datafile["outputFolder"];
    outputFolder = util::createOutputFolder(folderName);
    dem->createOutputWriter(outputFolder, snapshot-1);

    dem->writeBoundaryToFile();
    dem->writePatchGridToFile();
    dem->writeParticlesToFile(snapshot);
    dem->printBoundary();
    dem->printBoundaryContacts();
    dem->getStartDimension(distX, distY, distZ);

    debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT"
             << std::setw(OWID) << "migraT" << std::setw(OWID) << "totalT"
             << std::setw(OWID) << "overhead%" << std::endl;
  }

  // Broadcast the output folder to all processes
  broadcast(dem->getMPIWorld(), outputFolder, 0);

  auto currentTime = util::getParam<REAL>("timeAccrued");
  auto deltaT = util::getParam<REAL>("timeStep");
  auto maxTime = currentTime + deltaT * numSteps;
  auto iteration = startStep;
  while (currentTime < maxTime) {

    auto t0 = MPI_Wtime();
    dem->communicateGhostParticles(iteration);
    auto t1 = MPI_Wtime();
    auto commuT = t1 - t0;

    deltaT = dem->calcTimeStep(); 

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

    dem->gatherBoundaryContacts(); // must call before updateBoundary

    auto mass = dem->getTotalMassFromPatchParticleData();
    dem->updateBoundary(currentTime, deltaT, mass);

    dem->updatePatchBox();

    if (iteration % (numSteps / numSnapshots) == 0) {

      //auto t2 = MPI_Wtime();
      dem->gatherParticles();
      dem->gatherEnergy();
      //auto t3 = MPI_Wtime();
      //auto gatherT = t3 - t2;

      if (dem->getMPIRank() == 0) {
        dem->updateFileNames(snapshot);
        dem->writeBoundaryToFile();
        dem->printBoundary();
        dem->writePatchGridToFile();
        dem->writeParticlesToFile(snapshot);
        dem->printBoundaryContacts();
        dem->appendToProgressOutputFile(progressInf, distX, distY, distZ);
      }
      dem->printContact(combine(outputFolder, "traction", snapshot, 3));
      ++snapshot;
    }

    dem ->releaseReceivedParticles();

    auto t4 = MPI_Wtime();
    dem->migrateParticles();
    auto t5 = MPI_Wtime();
    auto migraT = t5 - t4;
    auto totalT = t5 - t0;

    if (dem->getMPIRank() == 0 && (iteration + 1) % (numSteps / numSnapshots) == 0) {
      debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT
               << std::setw(OWID) << migraT << std::setw(OWID) << totalT
               << std::setw(OWID) << (commuT + migraT) / totalT * 100
               << std::endl;
    }

    if (dem->areBoundaryTractionsEquilibrated(currentTime)) {
      if (dem->getMPIRank() == 0) {
        dem->updateFileNames(snapshot, ".end");
        dem->writeParticlesToFile(snapshot);
        dem->printBoundaryContacts();
        dem->printBoundary();
        dem->appendToProgressOutputFile(balancedInf, distX, distY, distZ);
      }
      break;
    }

    ++iteration;
    currentTime += deltaT;
  }

  if (dem->getMPIRank() == 0) {
    dem->closeProgressOutputFile(progressInf);
    dem->closeProgressOutputFile(balancedInf);
  }
}
