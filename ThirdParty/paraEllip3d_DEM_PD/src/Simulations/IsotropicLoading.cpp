#include <Simulations/IsotropicLoading.h>

using namespace dem;
void
IsotropicLoading::execute(Assembly* assembly)
{
  std::ofstream progressInf;
  std::ofstream balancedInf;

  std::size_t isotropicType = static_cast<std::size_t>(
    dem::Parameter::getSingleton().parameter["isotropicType"]);
  if (assembly->getMPIRank() == 0) {
    assembly->readBoundary(
      dem::Parameter::getSingleton().datafile["boundaryFile"].c_str());
    assembly->readParticle(
      dem::Parameter::getSingleton().datafile["particleFile"].c_str());
    assembly->openCompressProg(progressInf, "isotropic_progress");
    assembly->openCompressProg(balancedInf, "isotropic_balanced");
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
  timeStep = dem::Parameter::getSingleton().parameter["timeStep"];

  REAL sigmaEnd, sigmaInc, sigmaVar;
  std::size_t sigmaDiv;

  sigmaEnd = dem::Parameter::getSingleton().parameter["sigmaEnd"];
  sigmaDiv = dem::Parameter::getSingleton().parameter["sigmaDiv"];
  std::vector<REAL>& sigmaPath = dem::Parameter::getSingleton().sigmaPath;
  std::size_t sigma_i = 0;

  if (isotropicType == 1)
    sigmaVar = sigmaEnd;
  else if (isotropicType == 2) {
    REAL sigmaStart = dem::Parameter::getSingleton().parameter["sigmaStart"];
    sigmaInc = (sigmaEnd - sigmaStart) / sigmaDiv;
    sigmaVar = sigmaStart;
  } else if (isotropicType == 3) {
    sigmaVar = sigmaPath[sigma_i];
    sigmaInc = (sigmaPath[sigma_i + 1] - sigmaPath[sigma_i]) / sigmaDiv;
    sigmaEnd = sigmaPath[sigmaPath.size() - 1];
  }

  REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
  iteration = startStep;
  std::size_t iterSnap = startSnap;
  char cstr0[50];
  REAL distX, distY, distZ;
  if (assembly->getMPIRank() == 0) {
     assembly->plotBoundary(strcat(
      Assembly::combineString(cstr0, "isotropic_bdryplot_", iterSnap - 1, 3), ".dat"));
    assembly->plotGrid(strcat(
      Assembly::combineString(cstr0, "isotropic_gridplot_", iterSnap - 1, 3), ".dat"));
    assembly->printParticle(Assembly::combineString(cstr0, "isotropic_particle_", iterSnap - 1, 3));
    assembly->printBdryContact(
      Assembly::combineString(cstr0, "isotropic_bdrycntc_", iterSnap - 1, 3));
    assembly->printBoundary(Assembly::combineString(cstr0, "isotropic_boundary_", iterSnap - 1, 3));
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

    assembly->calcTimeStep(); // use values from last step, must call before findConact
    assembly->findContact();
    if (assembly->isBdryProcess())
      assembly->findBdryContact();

    assembly->clearContactForce();
    assembly->internalForce();
    if (assembly->isBdryProcess())
      assembly->boundaryForce();

    assembly->updateParticle();
    assembly->gatherBdryContact(); // must call before updateBoundary
    assembly->updateBoundary(sigmaVar, "isotropic");
    assembly->updateGrid();

    if (iteration % (netStep / netSnap) == 0) {
      time1 = MPI_Wtime();
      assembly->gatherParticle();
      assembly->gatherEnergy();
      time2 = MPI_Wtime();
      gatherT = time2 - time1;

      char cstr[50];
      if (assembly->getMPIRank() == 0) {
        assembly->plotBoundary(strcat(
          Assembly::combineString(cstr, "isotropic_bdryplot_", iterSnap, 3), ".dat"));
        assembly->plotGrid(strcat(Assembly::combineString(cstr, "isotropic_gridplot_", iterSnap, 3),
                        ".dat"));
        assembly->printParticle(Assembly::combineString(cstr, "isotropic_particle_", iterSnap, 3));
        assembly->printBdryContact(
          Assembly::combineString(cstr, "isotropic_bdrycntc_", iterSnap, 3));
        assembly->printBoundary(Assembly::combineString(cstr, "isotropic_boundary_", iterSnap, 3));
        assembly->printCompressProg(progressInf, distX, distY, distZ);
      }
      assembly->printContact(Assembly::combineString(cstr, "isotropic_contact_", iterSnap, 3));
      ++iterSnap;
    }

    assembly->releaseRecvParticle(); // late release because printContact refers to
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

    if (isotropicType == 1) {
      if (assembly->tractionErrorTol(sigmaVar, "isotropic")) {
        if (assembly->getMPIRank() == 0) {
          assembly->printParticle("isotropic_particle_end");
          assembly->printBdryContact("isotropic_bdrycntc_end");
          assembly->printBoundary("isotropic_boundary_end");
          assembly->printCompressProg(balancedInf, distX, distY, distZ);
        }
        break;
      }
    } else if (isotropicType == 2) {
      if (assembly->tractionErrorTol(sigmaVar, "isotropic")) {
        if (assembly->getMPIRank() == 0)
          assembly->printCompressProg(balancedInf, distX, distY, distZ);
        sigmaVar += sigmaInc;
      }
      if (assembly->tractionErrorTol(sigmaEnd, "isotropic")) {
        if (assembly->getMPIRank() == 0) {
          assembly->printParticle("isotropic_particle_end");
          assembly->printBdryContact("isotropic_bdrycntc_end");
          assembly->printBoundary("isotropic_boundary_end");
          assembly->printCompressProg(balancedInf, distX, distY, distZ);
        }
        break;
      }
    }

    if (isotropicType == 3) {
      if (assembly->tractionErrorTol(sigmaVar, "isotropic")) {
        if (assembly->getMPIRank() == 0)
          assembly->printCompressProg(balancedInf, distX, distY, distZ);
        sigmaVar += sigmaInc;
        if (sigmaVar == sigmaPath[sigma_i + 1]) {
          sigmaVar = sigmaPath[++sigma_i];
          sigmaInc = (sigmaPath[sigma_i + 1] - sigmaPath[sigma_i]) / sigmaDiv;
        }
      }
      if (assembly->tractionErrorTol(sigmaEnd, "isotropic")) {
        if (assembly->getMPIRank() == 0) {
          assembly->printParticle("isotropic_particle_end");
          assembly->printBdryContact("isotropic_bdrycntc_end");
          assembly->printBoundary("isotropic_boundary_end");
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
