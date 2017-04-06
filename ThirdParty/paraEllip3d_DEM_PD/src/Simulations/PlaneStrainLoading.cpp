#include <Simulations/PlaneStrainLoading.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;
void
PlaneStrainLoading::execute(Assembly* assembly)
{
  std::ofstream progressInf;

  if (assembly->getMPIRank() == 0) {
    assembly->readBoundary(
      dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
    assembly->readParticles(
      dem::Parameter::getSingleton().datafile["particleFile"].c_str());
    assembly->openCompressProg(progressInf, "plnstrn_progress");
  }
  assembly->scatterParticle();

  std::size_t startStep = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["startStep"]);
  std::size_t endStep = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["endStep"]);
  std::size_t startSnap = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["startSnap"]);
  std::size_t endSnap = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["endSnap"]);
  std::size_t netStep = endStep - startStep + 1;
  std::size_t netSnap = endSnap - startSnap + 1;
  REAL sigmaConf = dem::Parameter::getSingleton().parameter["sigmaConf"];
  timeStep = dem::Parameter::getSingleton().parameter["timeStep"];

  REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
  iteration = startStep;
  std::size_t iterSnap = startSnap;
  REAL distX, distY, distZ;
  if (assembly->getMPIRank() == 0) {
    assembly->plotBoundary(combine("plnstrn_bdryplot_", iterSnap - 1, 3) + ".dat");
    assembly->plotGrid(combine("plnstrn_gridplot_", iterSnap - 1, 3) + ".dat");
    assembly->printParticle(combine("plnstrn_particle_", iterSnap - 1, 3));
    assembly->printBdryContact(combine("plnstrn_bdrycntc_", iterSnap - 1, 3));
    assembly->printBoundary(combine("plnstrn_boundary_", iterSnap - 1, 3));
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

    // displacement control relies on constant time step, so do not call
    // calcTimeStep().
    // calcTimeStep(); // use values from last step, must call before findConact
    assembly->findContact();
    if (assembly->isBdryProcess())
      assembly->findBdryContact();

    assembly->clearContactForce();
    assembly->internalForce();
    if (assembly->isBdryProcess())
      assembly->boundaryForce();

    assembly->updateParticle();
    assembly->gatherBdryContact(); // must call before updateBoundary
    assembly->updateBoundary(sigmaConf, "plnstrn");
    assembly->updateGrid();

    if (iteration % (netStep / netSnap) == 0) {
      time1 = MPI_Wtime();
      assembly->gatherParticle();
      assembly->gatherEnergy();
      time2 = MPI_Wtime();
      gatherT = time2 - time1;

      if (assembly->getMPIRank() == 0) {
        assembly->plotBoundary(combine( "plnstrn_bdryplot_", iterSnap, 3) + ".dat");
        assembly->plotGrid(combine( "plnstrn_gridplot_", iterSnap, 3) + ".dat");
        assembly->printParticle(combine( "plnstrn_particle_", iterSnap, 3));
        assembly->printBdryContact(combine( "plnstrn_bdrycntc_", iterSnap, 3));
        assembly->printBoundary(combine( "plnstrn_boundary_", iterSnap, 3));
        // assembly->printCompressProg(progressInf, distX, distY, distZ); //
        // redundant
      }
      assembly->printContact(combine( "plnstrn_contact_", iterSnap, 3));
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

    if (assembly->getMPIRank() == 0 && iteration % 10 == 0)
      assembly->printCompressProg(progressInf, distX, distY, distZ);

    // no break condition, just through top/bottom displacement control
    ++iteration;
  }

  if (assembly->getMPIRank() == 0) {
    assembly->printParticle("plnstrn_particle_end");
    assembly->printBdryContact("plnstrn_bdrycntc_end");
    assembly->printBoundary("plnstrn_boundary_end");
    assembly->printCompressProg(progressInf, distX, distY, distZ);
  }

  if (assembly->getMPIRank() == 0)
    assembly->closeProg(progressInf);
}
