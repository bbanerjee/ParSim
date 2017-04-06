#include <Simulations/CoupledFluidFlow.h>
#include <Core/Util/Utility.h>

using namespace dem;
void
CoupledFluidFlow::execute(Assembly* assembly)
{
  Box allContainer = assembly->getAllContainer();
  Gradation gradation = assembly->getGradation();
  Fluid fluid = assembly->getFluid();
  ParticlePArray particleVec = assembly->getParticleVec();

  std::ofstream progressInf;
  std::ofstream particleInf;

  if (assembly->getMPIRank() == 0) {
    assembly->readBoundary(
      Parameter::get().datafile["boundaryFile"]);
    assembly->readParticles(
      Parameter::get().datafile["particleFile"]);
    assembly->openDepositProg(progressInf, "couple_progress");
    assembly->openParticleProg(particleInf, "particle_progress");
    /*1*/ fluid.initParameter(allContainer, gradation);
    /*1*/ fluid.initialize();
  }
  assembly->scatterParticle();

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
  if (assembly->getMPIRank() == 0) {
    assembly->plotBoundary(util::combine("couple_bdryplot_", iterSnap - 1, 3) + ".dat");
    assembly->plotGrid(util::combine("couple_gridplot_", iterSnap - 1, 3) + ".dat");
    assembly->printParticle(util::combine("couple_particle_", iterSnap - 1, 3));
    assembly->printBdryContact(util::combine("couple_bdrycntc_", iterSnap - 1, 3));
    /*3*/ fluid.plot(util::combine("couple_fluidplot_", iterSnap - 1, 3) + ".dat");
  }
  /*
  if (mpiRank == 0)
    debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT" <<
  std::setw(OWID) << "migraT"
         << std::setw(OWID) << "totalT" << std::setw(OWID) << "overhead%" <<
  std::endl;
  */
  while (timeAccrued < timeTotal) {
    // REAL time0 = MPI_Wtime();
    assembly->commuParticle();
    // REAL time2 = MPI_Wtime();
    // REAL commuT = time2 - time0;

    assembly->calcTimeStep(); // use values from last step, must call before
                              // findConact
    assembly->findContact();
    if (assembly->isBdryProcess())
      assembly->findBdryContact();

    assembly->clearContactForce();

    /*4*/ fluid.getParticleInfo(particleVec); // not allParticleVec
    /*5*/ fluid.runOneStep();
    /*6*/ fluid.calcParticleForce(particleVec,
                                  particleInf); // not allParticleVec
    /*7*/ fluid.penalize();

    assembly->internalForce();
    if (assembly->isBdryProcess())
      assembly->boundaryForce();

    assembly->updateParticle();
    assembly->updateGridMaxZ();

    timeCount += timeStep;
    timeAccrued += timeStep;
    if (timeCount >= timeTotal / netSnap) {
      // REAL time1 = MPI_Wtime();
      assembly->gatherParticle();
      assembly->gatherBdryContact();
      assembly->gatherEnergy();
      // time2 = MPI_Wtime();
      // REAL gatherT = time2 - time1;

      if (assembly->getMPIRank() == 0) {
        assembly->plotBoundary(util::combine( "couple_bdryplot_", iterSnap, 3) + ".dat");
        assembly->plotGrid(util::combine( "couple_gridplot_", iterSnap, 3) + ".dat");
        assembly->printParticle(util::combine( "couple_particle_", iterSnap, 3));
        assembly->printBdryContact(util::combine( "couple_bdrycntc_", iterSnap, 3));
        assembly->printDepositProg(progressInf);
        /*8*/ fluid.plot(util::combine( "couple_fluidplot_", iterSnap, 3) + ".dat");
      }
      assembly->printContact(util::combine( "couple_contact_", iterSnap, 3));

      timeCount = 0;
      ++iterSnap;
    }

    assembly
      ->releaseRecvParticle(); // late release because printContact refers to
                               // received particles
    // REAL time1 = MPI_Wtime();
    assembly->migrateParticle();
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

  if (assembly->getMPIRank() == 0) {
    assembly->closeProg(progressInf);
    assembly->closeProg(particleInf);
  }
}
