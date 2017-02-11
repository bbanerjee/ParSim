#include <Commands/PeridynamicsPullOutCommand.h>

using namespace dem;
void
PeridynamicsPullOutCommand::execute(Assembly* assembly)
{
  std::ofstream progressInf;
  std::ofstream periProgInf;
  std::ofstream periProgInfHalf;

  REAL distX, distY, distZ;
  if (assembly->getMPIRank() == 0) {
    distX = dem::Parameter::getSingleton().parameter["Xmax"] -
            dem::Parameter::getSingleton().parameter["Xmin"];
    distY = dem::Parameter::getSingleton().parameter["Ymax"] -
            dem::Parameter::getSingleton().parameter["Ymin"];
    distZ = dem::Parameter::getSingleton().parameter["Zmax"] -
            dem::Parameter::getSingleton().parameter["Zmin"];
    REAL x1 = dem::Parameter::getSingleton().parameter["Xmin"];
    REAL x2 = dem::Parameter::getSingleton().parameter["Xmax"];
    REAL y1 = dem::Parameter::getSingleton().parameter["Ymin"];
    REAL y2 = dem::Parameter::getSingleton().parameter["Ymax"];
    REAL z1 = dem::Parameter::getSingleton().parameter["Zmin"];
    REAL z2 = dem::Parameter::getSingleton().parameter["Zmax"];
    assembly->setContainer(Box(x1, y1, z1, x2, y2, z2));
    assembly->setGrid(Box(x1, y1, z1, x2, y2, z2)); // compute grid assumed to be
                                          // the same as container,
                                          // change in scatterParticle()
                                          // if necessary.
    assembly->readParticle(
      dem::Parameter::getSingleton().datafile["particleFile"].c_str());
    assembly->readPeriDynamicsData(
      dem::Parameter::getSingleton().datafile["periFile"].c_str());
    assembly->openCompressProg(progressInf, "rigidInc_progress");
    assembly->openPeriProgress(periProgInf, "rigidInc_peri.dat");
    assembly->openPeriProgress(periProgInfHalf, "rigidInc_peri_half.dat");
    assembly->calcParticleVolume(); // volume and horizon size are related to the mesh,
                          // July 14, 2014
    assembly->calcHorizonSize(); // so these two should be calculated before peri-points
                       // that are within sand particles are deleted
    assembly->removeInsidePeriParticles(); // delete those peri-points that are inside
                                 // sand particles
  }
  assembly->scatterDEMPeriParticle();
  assembly->constructPeriMatrix(); // construct the Matrix members in periParticleVec
  assembly->constructNeighbor();

  // if(assembly->getMPIRank()==0){
  // assembly->printPeriDomain("peridomain_rank_0.dat");
  // assembly->printRecvPeriDomain("peridomain_recv_rank_0.dat");    // bondVec is not
  // empty
  //}
  // if(assembly->getMPIRank()==1){
  // assembly->printPeriDomain("peridomain_rank_1.dat");
  // assembly->printRecvPeriDomain("peridomain_recv_rank_1.dat");    // bondVec is empty
  //}

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

  // REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
  iteration = startStep;
  std::size_t iterSnap = startSnap;
  char cstr0[50];
  if (assembly->getMPIRank() == 0) {
    assembly->plotGrid(strcat(Assembly::combineString(cstr0, "rigidInc_gridplot_", iterSnap - 1, 3),
                    ".dat"));
    assembly->printParticle(Assembly::combineString(cstr0, "rigidInc_particle_", iterSnap - 1, 3));
    assembly->printPeriProgress(periProgInf, 0);
    assembly->printPeriProgressHalf(periProgInfHalf, 0);
  }
  //    if (assembly->getMPIRank() == 0)
  //      debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT"
  //      << std::setw(OWID) << "migraT"
  //           << std::setw(OWID) << "totalT" << std::setw(OWID) << "overhead%"
  //           << std::endl;

  assembly->commuParticle();           // the commuPeriParticle() have to be called after
                             // commuParticle(), since
  assembly->commuPeriParticle();       // rankX1... are calculated in commuParticle()
  assembly->constructRecvPeriMatrix(); // construct the Matrix members in
                             // recvPeriParticle

  assembly->findRecvPeriBonds(); // find the peri-bonds between periParticleVec and
                       // recvPeriParticle
  assembly->calcParticleKinv();
  for (auto & pt : assembly->getPeriParticleVec()) {
    pt->initial();
  }
  for (auto & pt : assembly->getRecvPeriParticleVec()) {
    pt->initial();
  }
  assembly->calcParticleStress();
  assembly->calcParticleAcceleration();
  assembly->releaseRecvParticle();
  //    releaseRecvPeriParticle();
  //    releasePeriBondVec();    // free the bondVec in periParticleVec and
  //    bondVec.clear() before constructNeighbor() again
  assembly->findBoundaryPeriParticles();
  assembly->clearPeriDEMBonds();
  assembly->findPeriDEMBonds();
  while (iteration <= endStep) {
    //      findBoundaryPeriParticles();    // boundary peri-points are
    //      determined by their initial positions, however they are still needed
    //      to be determined
    //      findFixedPeriParticles();        // every step since some
    //      peri-points can pass from one cpu to another
    // commuT = migraT = gatherT = totalT = 0; time0 = MPI_Wtime();
    //      commuParticle(); time2 = MPI_Wtime(); commuT = time2 - time0;

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
    assembly->runFirstHalfStep();

    assembly->applyPeriBoundaryCondition();
    //      commuPeriParticle();
    //      constructRecvPeriMatrix();
    //      constructNeighbor();
    //      checkBondParticleAlive();
    assembly->findRecvPeriBonds();

    assembly->calcParticleStress();
    assembly->calcParticleAcceleration(); // peri-point acceleration is set zero at the
                                // beginning

    assembly->findPeriDEMBonds(); // based on current states of DEM particles and
                        // peri-points
    assembly->applyCoupledForces();
    assembly->runSecondHalfStep();

    assembly->updateParticle();
    assembly->gatherBdryContact(); // must call before updateBoundary
                         //      updateBoundary(sigmaConf, "triaxial");
                         //      updateGrid();
                         //      updatePeriGrid();

    if (iteration % (netStep / netSnap) == 0) {
      // time1 = MPI_Wtime();
      assembly->gatherParticle();
      assembly->gatherPeriParticle();
      assembly->gatherEnergy();
      // time2 = MPI_Wtime(); gatherT = time2 - time1;

      //    checkBondParticleAlive();
      //    findPeriDEMBonds();    // update peri-dem bonds

      char cstr[50];
      if (assembly->getMPIRank() == 0) {
        assembly->plotBoundary(strcat(
          Assembly::combineString(cstr, "rigidInc_bdryplot_", iterSnap, 3), ".dat"));
        assembly->plotGrid(strcat(Assembly::combineString(cstr, "rigidInc_gridplot_", iterSnap, 3),
                        ".dat"));
        assembly->printParticle(Assembly::combineString(cstr, "rigidInc_particle_", iterSnap, 3));
        assembly->printBdryContact(
          Assembly::combineString(cstr, "rigidInc_bdrycntc_", iterSnap, 3));
        assembly->printBoundary(Assembly::combineString(cstr, "rigidInc_boundary_", iterSnap, 3));
        // printCompressProg(progressInf, distX, distY, distZ); // redundant
        assembly->printPeriProgress(periProgInf, iterSnap);
        assembly->printPeriProgressHalf(periProgInfHalf, iterSnap);
      }
      assembly->printContact(Assembly::combineString(cstr, "rigidInc_contact_", iterSnap, 3));
      ++iterSnap;
    }
    //      releaseRecvParticle(); // late release because printContact refers
    //      to received particles
    //      releaseRecvPeriParticle();    // release recvPeriParticleVec and
    //      also remove peri-bonds between recvPeriParticle
    //      releasePeriBondVec();
    // time1 = MPI_Wtime();
    //      migrateParticle();
    //      migratePeriParticle();
    // time2 = MPI_Wtime(); migraT = time2 - time1; totalT = time2 - time0;
    //      if (assembly->getMPIRank() == 0 && (iteration+1 ) % (netStep / netSnap) == 0) //
    //      ignore gather and print time at this step
    //    debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT
    //    << std::setw(OWID) << migraT
    //         << std::setw(OWID) << totalT << std::setw(OWID) << (commuT +
    //         migraT)/totalT*100 << std::endl;

    if (assembly->getMPIRank() == 0 && iteration % 10 == 0)
      assembly->printCompressProg(progressInf, distX, distY, distZ);

    // no break condition, just through top/bottom displacement control
    ++iteration;
  }

  if (assembly->getMPIRank() == 0) {
    assembly->printParticle("rigidInc_particle_end");
    assembly->printBdryContact("rigidInc_bdrycntc_end");
    assembly->printBoundary("rigidInc_boundary_end");
    assembly->printCompressProg(progressInf, distX, distY, distZ);
    //      assembly->printPeriProgress(periProgInf, iterSnap);
  }
  if (assembly->getMPIRank() == 0) {
    assembly->closeProg(progressInf);
    assembly->closeProg(periProgInf);
    assembly->closeProg(periProgInfHalf);
  }
}
