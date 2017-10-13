#include <Simulations/DEM_CFD/CoupledFluidFlow.h>
#include <Core/Util/Utility.h>

using namespace dem;
void
CoupledFluidFlow::execute(DiscreteElements* dem)
{
  Box allContainer = dem->getAllContainer();
  Gradation gradation = dem->getGradation();
  Fluid fluid = dem->getFluid();
  DEMParticlePArray particleVec = dem->getDEMParticleVec();

  std::ofstream progressInf;
  std::ofstream particleInf;

  if (dem->getMPIRank() == 0) {
    dem->readBoundary(
      InputParameter::get().datafile["boundaryFilename"]);
    dem->readParticles(
      InputParameter::get().datafile["particleFilename"]);
    dem->openProgressOutputFile(progressInf, "couple_progress");
    dem->openParticleProg(particleInf, "particle_progress");
    /*1*/ fluid.initParameter(allContainer, gradation);
    /*1*/ fluid.initialize();
  }
  dem->scatterParticles();

  auto startStep = util::getParam<std::size_t>("startStep");
  auto endStep = util::getParam<std::size_t>("endStep");
  auto startSnap = util::getParam<std::size_t>("startSnap");
  auto endSnap = util::getParam<std::size_t>("endSnap");
  std::size_t netStep = endStep - startStep + 1;
  std::size_t netSnap = endSnap - startSnap + 1;
  timeStep = util::getParam<REAL>("timeStep");

  iteration = startStep;
  std::size_t iterSnap = startSnap;
  REAL timeCount = 0;
  timeAccrued = util::getParam<REAL>("timeAccrued");
  REAL timeTotal = timeAccrued + timeStep * netStep;

  std::string outputFolder(".");
  if (dem->getMPIRank() == 0) {

    // Create the output writer in the master process
    // <outputFolder> couple.pe3d </outputFolder>
    auto folderName =  dem::InputParameter::get().datafile["outputFolder"];
    outputFolder = util::createOutputFolder(folderName);
    //std::cout << "Output folder = " << outputFolder << "\n";
    dem->createOutputWriter(outputFolder, iterSnap-1);
    dem->writeBoundaryToFile();
    dem->writePatchGridToFile();
    dem->writeParticlesToFile(iterSnap);
    dem->printBoundaryContacts();
    /*3*/ fluid.plot(util::combine(".", "couple_fluidplot_", iterSnap - 1, 3) + ".dat");
  }
  /*
  if (mpiRank == 0)
    debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT" <<
  std::setw(OWID) << "migraT"
         << std::setw(OWID) << "totalT" << std::setw(OWID) << "overhead%" <<
  std::endl;
  */

  // Broadcast the output folder to all processes
  broadcast(dem->getMPIWorld(), outputFolder, 0);

  while (timeAccrued < timeTotal) {
    // REAL time0 = MPI_Wtime();
    dem->commuParticle();
    // REAL time2 = MPI_Wtime();
    // REAL commuT = time2 - time0;

    dem->calcTimeStep(); // use values from last step, must call before
                              // findConact
    dem->findContact();
    if (dem->isBoundaryProcess())
      dem->findBoundaryContacts();

    dem->clearContactForce();

    /*4*/ fluid.getParticleInfo(particleVec); // not allDEMParticleVec
    /*5*/ fluid.runOneStep();
    /*6*/ fluid.calcParticleForce(particleVec,
                                  particleInf); // not allDEMParticleVec
    /*7*/ fluid.penalize();

    dem->internalForce();
    if (dem->isBoundaryProcess())
      dem->boundaryForce();

    dem->updateParticles();
    dem->updatePatchBoxMaxZ();

    timeCount += timeStep;
    timeAccrued += timeStep;
    if (timeCount >= timeTotal / netSnap) {
      // REAL time1 = MPI_Wtime();
      dem->gatherParticles();
      dem->gatherBoundaryContacts();
      dem->gatherEnergy();
      // time2 = MPI_Wtime();
      // REAL gatherT = time2 - time1;

      if (dem->getMPIRank() == 0) {
        dem->updateFileNames(iterSnap);
        dem->writeBoundaryToFile();
        dem->writePatchGridToFile();
        dem->writeParticlesToFile(iterSnap);
        dem->printBoundaryContacts();
        dem->appendToProgressOutputFile(progressInf);
        /*8*/ fluid.plot(util::combine(".", "couple_fluidplot_", iterSnap, 3) + ".dat");
      }
      dem->printContact(util::combine(".", "couple_contact_", iterSnap, 3));

      timeCount = 0;
      ++iterSnap;
    }

    dem
      ->releaseReceivedParticles(); // late release because printContact refers to
                               // received particles
    // REAL time1 = MPI_Wtime();
    dem->migrateParticles();
    // time2 = MPI_Wtime();
    // REAL migraT = time2 - time1;
    // REAL totalT = time2 - time0;
    /*
    if (mpiRank == 0 && (iteration+1 ) % (netStep / netSnap) == 0) // ignore
  gather and print time at this step
  debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT <<
  std::setw(OWID) << migraT
       << std::setw(OWID) << totalT << std::setw(OWID) << (commuT +
  migraT)/totalT*100 << std::endl;
    */

    ++iteration;
  }

  if (dem->getMPIRank() == 0) {
    dem->closeProgressOutputFile(progressInf);
    dem->closeProgressOutputFile(particleInf);
  }
}
