#include <Simulations/DEM/IsotropicLoading.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;
void
IsotropicLoading::execute(DiscreteElements* dem)
{
  std::ofstream progressInf;
  std::ofstream balancedInf;

  if (dem->getMPIRank() == 0) {
    std::string boundaryFile = util::getFilename("boundaryFilename");
    std::string particleFile = util::getFilename("particleFilename");
    dem->readBoundary(boundaryFile);
    dem->readParticles(particleFile);
    dem->openProgressOutputFile(progressInf, "isotropic_progress");
    dem->openProgressOutputFile(balancedInf, "isotropic_balanced");
  }

  dem->scatterParticles();

  auto startStep = util::getParam<std::size_t>("startStep");
  auto endStep = util::getParam<std::size_t>("endStep");
  auto startSnap = util::getParam<std::size_t>("startSnap");
  auto endSnap = util::getParam<std::size_t>("endSnap");

  auto timeStep = util::getParam<REAL>("timeStep");
  auto currentTime = startStep * timeStep;

  auto netStep = endStep - startStep + 1;
  auto netSnap = endSnap - startSnap + 1;

  auto iterSnap = startSnap;

  std::string outputFolder(".");
  REAL distX = 0, distY = 0, distZ = 0;
  if (dem->getMPIRank() == 0) {
    // Create the output writer in the master process
    // <outputFolder> isotropic.pe3d </outputFolder>
    auto folderName =  dem::InputParameter::get().datafile["outputFolder"];
    outputFolder = util::createOutputFolder(folderName);
    dem->createOutputWriter(outputFolder, iterSnap-1);

    dem->writeBoundaryToFile(currentTime);
    dem->writePatchGridToFile(currentTime);
    dem->writeParticlesToFile(iterSnap, currentTime);
    dem->printBoundary();
    dem->printBoundaryContacts();
    dem->getStartDimension(distX, distY, distZ);

    debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT"
             << std::setw(OWID) << "migraT" << std::setw(OWID) << "totalT"
             << std::setw(OWID) << "overhead%" << std::endl;
  }

  // Broadcast the output folder to all processes
  broadcast(dem->getMPIWorld(), outputFolder, 0);
  dem->printContact();

  auto sigmaEnd = util::getParam<REAL>("sigmaEnd");
  auto sigmaDiv = util::getParam<REAL>("sigmaDiv");
  std::vector<REAL>& sigmaPath = InputParameter::get().sigmaPath;

  auto isotropicType = util::getParam<std::size_t>("isotropicType");
  REAL sigmaInc = 0, sigmaVar = 0;
  auto sigma_index = 0u;
  if (isotropicType == 1)
    sigmaVar = sigmaEnd;
  else if (isotropicType == 2) {
    REAL sigmaStart = util::getParam<REAL>("sigmaStart");
    sigmaInc = (sigmaEnd - sigmaStart) / sigmaDiv;
    sigmaVar = sigmaStart;
  } else if (isotropicType == 3) {
    sigmaVar = sigmaPath[sigma_index];
    sigmaInc = (sigmaPath[sigma_index + 1] - sigmaPath[sigma_index]) / sigmaDiv;
    sigmaEnd = sigmaPath[sigmaPath.size() - 1];
  }

  auto iteration = startStep;

  while (iteration <= endStep) {

    auto t0 = MPI_Wtime();
    dem->communicateGhostParticles(iteration);
    auto t1 = MPI_Wtime();
    auto commuT = t1 - t0;

    timeStep = dem->calcTimeStep(timeStep); 
    currentTime += timeStep;

    dem->findContact(iteration);

    if (dem->isBoundaryProcess()) {
      dem->findBoundaryContacts(iteration);
    }

    dem->initializeForces();

    dem->internalForce(timeStep, iteration);

    if (dem->isBoundaryProcess()) {
      dem->boundaryForce(timeStep, iteration);
    }

    dem->updateParticles(timeStep, iteration);

    dem->gatherBoundaryContacts(); // must call before updateBoundary

    dem->updateBoundary(sigmaVar, "isotropic");

    dem->updatePatchBox();

    if (iteration % (netStep / netSnap) == 0) {

      //auto t2 = MPI_Wtime();
      dem->gatherParticles();
      dem->gatherEnergy();
      //auto t3 = MPI_Wtime();
      //auto gatherT = t3 - t2;

      if (dem->getMPIRank() == 0) {
        dem->updateFileNames(iterSnap);
        dem->writeBoundaryToFile(currentTime);
        dem->printBoundary();
        dem->writePatchGridToFile(currentTime);
        dem->writeParticlesToFile(iterSnap, currentTime);
        dem->printBoundaryContacts();
        dem->appendToProgressOutputFile(progressInf, iteration, timeStep, distX, distY, distZ);
      }
      dem->printContact();
      ++iterSnap;
    }

    dem ->releaseReceivedParticles();

    auto t4 = MPI_Wtime();
    dem->migrateParticles(iteration);
    auto t5 = MPI_Wtime();
    auto migraT = t5 - t4;
    auto totalT = t5 - t0;

    if (dem->getMPIRank() == 0 && (iteration + 1) % (netStep / netSnap) == 0) {
      debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT
               << std::setw(OWID) << migraT << std::setw(OWID) << totalT
               << std::setw(OWID) << (commuT + migraT) / totalT * 100
               << std::endl;
    }

    if (isotropicType == 1) {
      if (dem->areBoundaryTractionsEquilibrated(sigmaVar, "isotropic")) {
        if (dem->getMPIRank() == 0) {
          dem->updateFileNames(iterSnap, ".end");
          dem->writeParticlesToFile(iterSnap, currentTime);
          dem->printBoundaryContacts();
          dem->printBoundary();
          dem->appendToProgressOutputFile(balancedInf, iteration, timeStep,  distX, distY, distZ);
        }
        break;
      }
    } else if (isotropicType == 2) {
      if (dem->areBoundaryTractionsEquilibrated(sigmaVar, "isotropic")) {
        if (dem->getMPIRank() == 0)
          dem->appendToProgressOutputFile(balancedInf, iteration, timeStep, distX, distY, distZ);
        sigmaVar += sigmaInc;
      }
      if (dem->areBoundaryTractionsEquilibrated(sigmaEnd, "isotropic")) {
        if (dem->getMPIRank() == 0) {
          dem->updateFileNames(iterSnap, ".end");
          dem->writeParticlesToFile(iterSnap, currentTime);
          dem->printBoundaryContacts();
          dem->printBoundary();
          dem->appendToProgressOutputFile(balancedInf, iteration, timeStep, distX, distY, distZ);
        }
        break;
      }
    } else if (isotropicType == 3) {
      if (dem->areBoundaryTractionsEquilibrated(sigmaVar, "isotropic")) {
        if (dem->getMPIRank() == 0)
          dem->appendToProgressOutputFile(balancedInf, iteration, timeStep, distX, distY, distZ);
        sigmaVar += sigmaInc;
        if (sigmaVar == sigmaPath[sigma_index + 1]) {
          sigmaVar = sigmaPath[++sigma_index];
          sigmaInc = (sigmaPath[sigma_index + 1] - sigmaPath[sigma_index]) / sigmaDiv;
        }
      }
      if (dem->areBoundaryTractionsEquilibrated(sigmaEnd, "isotropic")) {
        if (dem->getMPIRank() == 0) {
          dem->updateFileNames(iterSnap, ".end");
          dem->writeParticlesToFile(iterSnap, currentTime);
          dem->printBoundaryContacts();
          dem->printBoundary();
          dem->appendToProgressOutputFile(balancedInf, iteration, timeStep,  distX, distY, distZ);
        }
        break;
      }
    }

    ++iteration;
  }

  if (dem->getMPIRank() == 0) {
    dem->closeProgressOutputFile(progressInf);
    dem->closeProgressOutputFile(balancedInf);
  }
}
