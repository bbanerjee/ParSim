#include <Simulations/OedometerLoading.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;

void
OedometerLoading::execute(Assembly* assembly)
{
  std::ofstream progressInf;
  std::ofstream balancedInf;

  auto odometerType = util::getParam<std::size_t>("odometerType");
  if (assembly->getMPIRank() == 0) {
    assembly->readBoundary(
      Parameter::get().datafile["boundaryFile"]);
    assembly->readParticles(
      Parameter::get().datafile["particleFile"]);
    assembly->openCompressProg(progressInf, "odometer_progress");
    assembly->openCompressProg(balancedInf, "odometer_balanced");
  }
  assembly->scatterParticle();

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
  std::vector<REAL>& sigmaPath = Parameter::get().sigmaPath;
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
  if (assembly->getMPIRank() == 0) {
    assembly->plotBoundary(combine("odometer_bdryplot_", iterSnap - 1, 3) + ".dat");
    assembly->plotGrid(combine("odometer_gridplot_", iterSnap - 1, 3) + ".dat");
    assembly->printParticle(combine("odometer_particle_", iterSnap - 1, 3));
    assembly->printBdryContact(combine("odometer_bdrycntc_", iterSnap - 1, 3));
    assembly->printBoundary(combine("odometer_boundary_", iterSnap - 1, 3));
    assembly->getStartDimension(distX, distY, distZ);
  }
  if (assembly->getMPIRank() == 0)
    debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT"
             << std::setw(OWID) << "migraT" << std::setw(OWID) << "totalT"
             << std::setw(OWID) << "overhead%" << std::endl;
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
    assembly->updateBoundary(sigmaVar, "odometer");
    assembly->updateGrid();

    if (iteration % (netStep / netSnap) == 0) {
      time1 = MPI_Wtime();
      assembly->gatherParticle();
      assembly->gatherEnergy();
      time2 = MPI_Wtime();
      gatherT = time2 - time1;

      if (assembly->getMPIRank() == 0) {
        assembly->plotBoundary(combine( "odometer_bdryplot_", iterSnap, 3) + ".dat");
        assembly->plotGrid(combine( "odometer_gridplot_", iterSnap, 3) + ".dat");
        assembly->printParticle(combine( "odometer_particle_", iterSnap, 3));
        assembly->printBdryContact(combine( "odometer_bdrycntc_", iterSnap, 3));
        assembly->printBoundary(combine( "odometer_boundary_", iterSnap, 3));
        assembly->printCompressProg(progressInf, distX, distY, distZ);
      }
      assembly->printContact(combine( "odometer_contact_", iterSnap, 3));
      ++iterSnap;
    }

    assembly
      ->releaseRecvParticle(); // late release because printContact refers to
                               // received particles
    time1 = MPI_Wtime();
    assembly->migrateParticle();
    time2 = MPI_Wtime();
    migraT = time2 - time1;
    totalT = time2 - time0;
    if (assembly->getMPIRank() == 0 &&
        (iteration + 1) % (netStep / netSnap) ==
          0) // ignore gather and print time at this step
      debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT
               << std::setw(OWID) << migraT << std::setw(OWID) << totalT
               << std::setw(OWID) << (commuT + migraT) / totalT * 100
               << std::endl;

    if (odometerType == 1) {
      if (assembly->tractionErrorTol(sigmaVar, "odometer")) {
        if (assembly->getMPIRank() == 0)
          assembly->printCompressProg(balancedInf, distX, distY, distZ);
        sigmaVar += sigmaInc;
      }
      if (assembly->tractionErrorTol(sigmaEnd, "odometer")) {
        if (assembly->getMPIRank() == 0) {
          assembly->printParticle("odometer_particle_end");
          assembly->printBdryContact("odometer_bdrycntc_end");
          assembly->printBoundary("odometer_boundary_end");
          assembly->printCompressProg(balancedInf, distX, distY, distZ);
        }
        break;
      }
    } else if (odometerType == 2) {
      if (assembly->tractionErrorTol(sigmaVar, "odometer")) {
        if (assembly->getMPIRank() == 0)
          assembly->printCompressProg(balancedInf, distX, distY, distZ);
        sigmaVar += sigmaInc;
        if (sigmaVar == sigmaPath[sigma_i + 1]) {
          sigmaVar = sigmaPath[++sigma_i];
          sigmaInc = (sigmaPath[sigma_i + 1] - sigmaPath[sigma_i]) / sigmaDiv;
        }
      }
      if (assembly->tractionErrorTol(sigmaEnd, "odometer")) {
        if (assembly->getMPIRank() == 0) {
          assembly->printParticle("odometer_particle_end");
          assembly->printBdryContact("odometer_bdrycntc_end");
          assembly->printBoundary("odometer_boundary_end");
          assembly->printCompressProg(balancedInf, distX, distY, distZ);
        }
        break;
      }
    }

    ++iteration;
  }

  if (assembly->getMPIRank() == 0) {
    assembly->closeProg(progressInf);
    assembly->closeProg(balancedInf);
  }
}
