#include <Simulations/DEM/OedometerLoading.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;

void
OedometerLoading::execute(DiscreteElements* dem)
{
  std::ofstream progressInf;
  std::ofstream balancedInf;

  auto odometerType = util::getParam<std::size_t>("odometerType");
  if (dem->getMPIRank() == 0) {
    dem->readBoundary(
      InputParameter::get().datafile["boundaryFilename"]);
    dem->readParticles(
      InputParameter::get().datafile["particleFilename"]);
    dem->openProgressOutputFile(progressInf, "odometer_progress");
    dem->openProgressOutputFile(balancedInf, "odometer_balanced");
  }
  dem->scatterParticles();

  auto startStep = util::getParam<std::size_t>("startStep");
  auto endStep = util::getParam<std::size_t>("endStep");
  auto startSnap = util::getParam<std::size_t>("startSnap");
  auto endSnap = util::getParam<std::size_t>("endSnap");
  std::size_t netStep = endStep - startStep + 1;
  std::size_t netSnap = endSnap - startSnap + 1;
  timeStep = util::getParam<REAL>("timeStep");

  REAL sigmaEnd, sigmaInc, sigmaVar;
  std::size_t sigmaDiv;

  sigmaEnd = util::getParam<REAL>("sigmaEnd");
  sigmaDiv = util::getParam<REAL>("sigmaDiv");
  std::vector<REAL>& sigmaPath = InputParameter::get().sigmaPath;
  std::size_t sigma_i = 0;

  if (odometerType == 1) {
    REAL sigmaStart = util::getParam<REAL>("sigmaStart");
    sigmaInc = (sigmaEnd - sigmaStart) / sigmaDiv;
    sigmaVar = sigmaStart;
  } else if (odometerType == 2) {
    sigmaVar = sigmaPath[sigma_i];
    sigmaInc = (sigmaPath[sigma_i + 1] - sigmaPath[sigma_i]) / sigmaDiv;
    sigmaEnd = sigmaPath[sigmaPath.size() - 1];
  }

  REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
  iteration = startStep;
  std::size_t iterSnap = startSnap;
  REAL distX, distY, distZ;

  std::string outputFolder(".");
  if (dem->getMPIRank() == 0) {

    // Create the output writer in the master process
    // <outputFolder> oedometer.pe3d </outputFolder>
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

    dem->calcTimeStep(); // use values from last step, must call before
                              // findConact
    dem->findContact();
    if (dem->isBoundaryProcess())
      dem->findBoundaryContacts();

    dem->clearContactForce();
    dem->internalForce();
    if (dem->isBoundaryProcess())
      dem->boundaryForce();

    dem->updateParticles();
    dem->gatherBoundaryContacts(); // must call before updateBoundary
    dem->updateBoundary(sigmaVar, "odometer");
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
        dem->appendToProgressOutputFile(progressInf, distX, distY, distZ);
      }
      dem->printContact(combine(outputFolder, "odometer_contact_", iterSnap, 3));
      ++iterSnap;
    }

    dem
      ->releaseReceivedParticles(); // late release because printContact refers to
                               // received particles
    time1 = MPI_Wtime();
    dem->migrateParticles();
    time2 = MPI_Wtime();
    migraT = time2 - time1;
    totalT = time2 - time0;
    if (dem->getMPIRank() == 0 &&
        (iteration + 1) % (netStep / netSnap) ==
          0) // ignore gather and print time at this step
      debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT
               << std::setw(OWID) << migraT << std::setw(OWID) << totalT
               << std::setw(OWID) << (commuT + migraT) / totalT * 100
               << std::endl;

    if (odometerType == 1) {
      if (dem->areBoundaryTractionsEquilibrated(sigmaVar, "odometer")) {
        if (dem->getMPIRank() == 0)
          dem->appendToProgressOutputFile(balancedInf, distX, distY, distZ);
        sigmaVar += sigmaInc;
      }
      if (dem->areBoundaryTractionsEquilibrated(sigmaEnd, "odometer")) {
        if (dem->getMPIRank() == 0) {
          dem->updateFileNames(iterSnap, ".end");
          dem->writeParticlesToFile(iterSnap);
          dem->printBoundaryContacts();
          dem->printBoundary();
          dem->appendToProgressOutputFile(balancedInf, distX, distY, distZ);
        }
        break;
      }
    } else if (odometerType == 2) {
      if (dem->areBoundaryTractionsEquilibrated(sigmaVar, "odometer")) {
        if (dem->getMPIRank() == 0)
          dem->appendToProgressOutputFile(balancedInf, distX, distY, distZ);
        sigmaVar += sigmaInc;
        if (sigmaVar == sigmaPath[sigma_i + 1]) {
          sigmaVar = sigmaPath[++sigma_i];
          sigmaInc = (sigmaPath[sigma_i + 1] - sigmaPath[sigma_i]) / sigmaDiv;
        }
      }
      if (dem->areBoundaryTractionsEquilibrated(sigmaEnd, "odometer")) {
        if (dem->getMPIRank() == 0) {
          dem->updateFileNames(iterSnap, ".end");
          dem->writeParticlesToFile(iterSnap);
          dem->printBoundaryContacts();
          dem->printBoundary();
          dem->appendToProgressOutputFile(balancedInf, distX, distY, distZ);
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
