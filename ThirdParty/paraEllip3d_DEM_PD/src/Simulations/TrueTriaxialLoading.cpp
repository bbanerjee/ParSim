#include <Simulations/TrueTriaxialLoading.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;
void
TrueTriaxialLoading::execute(Assembly* assembly)
{
  std::ofstream progressInf;
  std::ofstream balancedInf;

  auto trueTriaxialType = util::getParam<std::size_t>("trueTriaxialType");
  if (assembly->getMPIRank() == 0) {
    assembly->readBoundary(
      Parameter::get().datafile["boundaryFile"]);
    assembly->readParticles(
      Parameter::get().datafile["particleFile"]);
    assembly->openCompressProg(progressInf, "trueTriaxial_progress");
    assembly->openCompressProg(balancedInf, "trueTriaxial_balanced");
  }
  assembly->scatterParticle();

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
  if (assembly->getMPIRank() == 0) {

    // Create the output writer in the master process
    // <outputFolder> truetriaxial.pe3d </outputFolder>
    auto folderName =  dem::Parameter::get().datafile["outputFolder"];
    outputFolder = util::createOutputFolder(folderName);
    std::cout << "Output folder = " << outputFolder << "\n";
    assembly->createOutputWriter(outputFolder, iterSnap-1);

    assembly->plotBoundary();
    assembly->plotGrid();
    assembly->plotParticle();
    assembly->printBdryContact();
    assembly->printBoundary();
    assembly->getStartDimension(distX, distY, distZ);
  }
  if (assembly->getMPIRank() == 0)
    debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT"
             << std::setw(OWID) << "migraT" << std::setw(OWID) << "totalT"
             << std::setw(OWID) << "overhead%" << std::endl;

  // Broadcast the output folder to all processes
  broadcast(assembly->getMPIWorld(), outputFolder, 0);

  while (iteration <= endStep) {
    commuT = migraT = gatherT = totalT = 0;
    time0 = MPI_Wtime();
    assembly->commuParticle();
    time2 = MPI_Wtime();
    commuT = time2 - time0;

    assembly->calcTimeStep(); // use values from last step, must call before
                              // findConact
    assembly->findContact();
    if (assembly->isBdryProcess())
      assembly->findBdryContact();

    assembly->clearContactForce();
    assembly->internalForce();
    if (assembly->isBdryProcess())
      assembly->boundaryForce();

    assembly->updateParticle();
    assembly->gatherBdryContact(); // must call before updateBoundary

    if (trueTriaxialType == 1)
      assembly->updateBoundary(sigmaVarZ, "trueTriaxial", sigmaVarX, sigmaVarY);
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
      assembly->updateBoundary(sigmaZ, "trueTriaxial", sigmaX, sigmaY);
    }

    assembly->updateGrid();

    if (iteration % (netStep / netSnap) == 0) {
      time1 = MPI_Wtime();
      assembly->gatherParticle();
      assembly->gatherEnergy();
      time2 = MPI_Wtime();
      gatherT = time2 - time1;

      if (assembly->getMPIRank() == 0) {
        assembly->updateFileNames(iterSnap);
        assembly->plotBoundary();
        assembly->plotGrid();
        assembly->plotParticle();
        assembly->printBdryContact();
        assembly->printBoundary();
        assembly->printCompressProg(progressInf, distX, distY, distZ);
      }
      assembly->printContact(combine(".", "trueTriaxial_contact_", iterSnap, 3));
      ++iterSnap;
    }

    assembly->releaseRecvParticle(); // late release because
                                     // assembly->printContact refers to
                                     // received particles
    time1 = MPI_Wtime();
    assembly->migrateParticle();
    time2 = MPI_Wtime();
    migraT = time2 - time1;
    totalT = time2 - time0;
    if (assembly->getMPIRank() == 0 &&
        (iteration + 1) % (netStep / netSnap) ==
          0) // ignore gather and assembly->print time at this step
      debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT
               << std::setw(OWID) << migraT << std::setw(OWID) << totalT
               << std::setw(OWID) << (commuT + migraT) / totalT * 100
               << std::endl;

    if (trueTriaxialType == 1) {
      if (assembly->tractionErrorTol(sigmaVarZ, "trueTriaxial", sigmaVarX,
                                     sigmaVarY)) {
        if (assembly->getMPIRank() == 0)
          assembly->printCompressProg(balancedInf, distX, distY, distZ);
        sigmaVarZ += sigmaIncZ;
        sigmaVarX += sigmaIncX;
        sigmaVarY += sigmaIncY;
      }
      if (assembly->tractionErrorTol(sigmaEndZ, "trueTriaxial", sigmaEndX,
                                     sigmaEndY)) {
        if (assembly->getMPIRank() == 0) {
          assembly->updateFileNames(iterSnap, ".end");
          assembly->plotParticle();
          assembly->printBdryContact();
          assembly->printBoundary();
          assembly->printCompressProg(balancedInf, distX, distY, distZ);
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
      if (assembly->tractionErrorTol(sigmaZ, "trueTriaxial", sigmaX, sigmaY)) {
        if (assembly->getMPIRank() == 0)
          assembly->printCompressProg(balancedInf, distX, distY, distZ);
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
      if (assembly->tractionErrorTol(sigmaZ, "trueTriaxial", sigmaX, sigmaY)) {
        if (assembly->getMPIRank() == 0) {
          assembly->updateFileNames(iterSnap, ".end");
          assembly->plotParticle();
          assembly->printBdryContact();
          assembly->printBoundary();
          assembly->printCompressProg(balancedInf, distX, distY, distZ);
        }
        break;
      }
    }

    ++iteration;
  }

  if (assembly->getMPIRank() == 0) {
    assembly->updateFileNames(iterSnap, ".end");
    assembly->plotParticle();
    assembly->printBdryContact();
    assembly->printBoundary();
    assembly->printCompressProg(progressInf, distX, distY, distZ);
  }

  if (assembly->getMPIRank() == 0) {
    assembly->closeProg(progressInf);
    assembly->closeProg(balancedInf);
  }
}
