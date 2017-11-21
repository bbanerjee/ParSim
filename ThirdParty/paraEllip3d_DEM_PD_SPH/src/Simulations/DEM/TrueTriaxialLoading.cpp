#include <Simulations/DEM/TrueTriaxialLoading.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;
void
TrueTriaxialLoading::execute(DiscreteElements* dem)
{
  std::ofstream progressInf;
  std::ofstream balancedInf;

  auto trueTriaxialType = util::getParam<std::size_t>("trueTriaxialType");
  if (dem->getMPIRank() == 0) {
    dem->readBoundary(
      util::getFilename("boundaryFilename"));
    dem->readParticles(
      util::getFilename("particleFilename"));
    dem->openProgressOutputFile(progressInf, "trueTriaxial_progress");
    dem->openProgressOutputFile(balancedInf, "trueTriaxial_balanced");
  }
  dem->scatterParticles();

  auto startStep = util::getParam<std::size_t>("startStep");
  auto endStep = util::getParam<std::size_t>("endStep");
  auto startSnap = util::getParam<std::size_t>("startSnap");
  auto endSnap = util::getParam<std::size_t>("endSnap");
  std::size_t netStep = endStep - startStep + 1;
  std::size_t netSnap = endSnap - startSnap + 1;
  timeStep = util::getParam<REAL>("timeStep");

  REAL sigmaStart, sigmaEndZ, sigmaEndX, sigmaEndY;
  REAL sigmaDiv, sigmaIncZ, sigmaIncX, sigmaIncY, sigmaVarZ, sigmaVarX,
    sigmaVarY;
  REAL sigmaInit[3], sigmaEnd, sigmaInc, sigmaVar;
  std::size_t changeDirc;
  sigmaDiv = util::getParam<REAL>("sigmaDiv");

  if (trueTriaxialType == 1) {
    sigmaStart = util::getParam<REAL>("sigmaStart");
    sigmaEndZ = util::getParam<REAL>("sigmaEndZ");
    sigmaEndX = util::getParam<REAL>("sigmaEndX");
    sigmaEndY = util::getParam<REAL>("sigmaEndY");
    sigmaIncZ = (sigmaEndZ - sigmaStart) / sigmaDiv;
    sigmaIncX = (sigmaEndX - sigmaStart) / sigmaDiv;
    sigmaIncY = (sigmaEndY - sigmaStart) / sigmaDiv;
    sigmaVarZ = sigmaStart;
    sigmaVarX = sigmaStart;
    sigmaVarY = sigmaStart;
  } else if (trueTriaxialType == 2) {
    sigmaInit[0] = util::getParam<REAL>("sigmaStartX");
    sigmaInit[1] = util::getParam<REAL>("sigmaStartY");
    sigmaInit[2] = util::getParam<REAL>("sigmaStartZ");
    sigmaEnd = util::getParam<REAL>("sigmaEnd");
    changeDirc = util::getParam<size_t>("changeDirc");
    sigmaInc = (sigmaEnd - sigmaInit[changeDirc]) / sigmaDiv;
    sigmaVar = sigmaInit[changeDirc];
  }

  REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
  iteration = startStep;
  std::size_t iterSnap = startSnap;
  REAL distX, distY, distZ;
  std::string outputFolder(".");
  if (dem->getMPIRank() == 0) {

    // Create the output writer in the master process
    // <outputFolder> truetriaxial.pe3d </outputFolder>
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

    if (trueTriaxialType == 1)
      dem->updateBoundary(sigmaVarZ, "trueTriaxial", sigmaVarX, sigmaVarY);
    else if (trueTriaxialType == 2) {
      REAL sigmaX = 0.0, sigmaY = 0.0, sigmaZ = 0.0;
      if (changeDirc == 0) {
        sigmaX = sigmaVar;
        sigmaY = sigmaInit[1];
        sigmaZ = sigmaInit[2];
      } else if (changeDirc == 1) {
        sigmaX = sigmaInit[0];
        sigmaY = sigmaVar;
        sigmaZ = sigmaInit[2];
      } else if (changeDirc == 2) {
        sigmaX = sigmaInit[0];
        sigmaY = sigmaInit[1];
        sigmaZ = sigmaVar;
      }
      dem->updateBoundary(sigmaZ, "trueTriaxial", sigmaX, sigmaY);
    }

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
      dem->printContact(combine(".", "trueTriaxial_contact_", iterSnap, 3));
      ++iterSnap;
    }

    dem->releaseReceivedParticles(); // late release because
                                     // dem->printContact refers to
                                     // received particles
    time1 = MPI_Wtime();
    dem->migrateParticles();
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

    if (trueTriaxialType == 1) {
      if (dem->areBoundaryTractionsEquilibrated(sigmaVarZ, "trueTriaxial", sigmaVarX,
                                     sigmaVarY)) {
        if (dem->getMPIRank() == 0)
          dem->appendToProgressOutputFile(balancedInf, distX, distY, distZ);
        sigmaVarZ += sigmaIncZ;
        sigmaVarX += sigmaIncX;
        sigmaVarY += sigmaIncY;
      }
      if (dem->areBoundaryTractionsEquilibrated(sigmaEndZ, "trueTriaxial", sigmaEndX,
                                     sigmaEndY)) {
        if (dem->getMPIRank() == 0) {
          dem->updateFileNames(iterSnap, ".end");
          dem->writeParticlesToFile(iterSnap);
          dem->printBoundaryContacts();
          dem->printBoundary();
          dem->appendToProgressOutputFile(balancedInf, distX, distY, distZ);
        }
        break;
      }
    } else if (trueTriaxialType == 2) {
      REAL sigmaX = 0.0, sigmaY = 0.0, sigmaZ = 0.0;
      if (changeDirc == 0) {
        sigmaX = sigmaVar;
        sigmaY = sigmaInit[1];
        sigmaZ = sigmaInit[2];
      } else if (changeDirc == 1) {
        sigmaX = sigmaInit[0];
        sigmaY = sigmaVar;
        sigmaZ = sigmaInit[2];
      } else if (changeDirc == 2) {
        sigmaX = sigmaInit[0];
        sigmaY = sigmaInit[1];
        sigmaZ = sigmaVar;
      }
      if (dem->areBoundaryTractionsEquilibrated(sigmaZ, "trueTriaxial", sigmaX, sigmaY)) {
        if (dem->getMPIRank() == 0)
          dem->appendToProgressOutputFile(balancedInf, distX, distY, distZ);
        sigmaVar += sigmaInc;
      }

      if (changeDirc == 0) {
        sigmaX = sigmaEnd;
        sigmaY = sigmaInit[1];
        sigmaZ = sigmaInit[2];
      } else if (changeDirc == 1) {
        sigmaX = sigmaInit[0];
        sigmaY = sigmaEnd;
        sigmaZ = sigmaInit[2];
      } else if (changeDirc == 2) {
        sigmaX = sigmaInit[0];
        sigmaY = sigmaInit[1];
        sigmaZ = sigmaEnd;
      }
      if (dem->areBoundaryTractionsEquilibrated(sigmaZ, "trueTriaxial", sigmaX, sigmaY)) {
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
    dem->updateFileNames(iterSnap, ".end");
    dem->writeParticlesToFile(iterSnap);
    dem->printBoundaryContacts();
    dem->printBoundary();
    dem->appendToProgressOutputFile(progressInf, distX, distY, distZ);
  }

  if (dem->getMPIRank() == 0) {
    dem->closeProgressOutputFile(progressInf);
    dem->closeProgressOutputFile(balancedInf);
  }
}
