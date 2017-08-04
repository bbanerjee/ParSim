#include <Simulations/DEM_PD/PeridynamicsPullOut.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;

using DiscreteElements = dem::DiscreteElements;
using Peridynamics = pd::Peridynamics;

void
PeridynamicsPullOut::execute(DiscreteElements* dem, Peridynamics* pd)
{
  std::ofstream progressInf;
  std::ofstream periProgInf;
  std::ofstream periProgInfHalf;

  REAL distX, distY, distZ;
  if (dem->getMPIRank() == 0) {
    distX = util::getParam<REAL>("Xmax") -
            util::getParam<REAL>("Xmin");
    distY = util::getParam<REAL>("Ymax") -
            util::getParam<REAL>("Ymin");
    distZ = util::getParam<REAL>("Zmax") -
            util::getParam<REAL>("Zmin");
    REAL x1 = util::getParam<REAL>("Xmin");
    REAL x2 = util::getParam<REAL>("Xmax");
    REAL y1 = util::getParam<REAL>("Ymin");
    REAL y2 = util::getParam<REAL>("Ymax");
    REAL z1 = util::getParam<REAL>("Zmin");
    REAL z2 = util::getParam<REAL>("Zmax");
    dem->setContainer(Box(x1, y1, z1, x2, y2, z2));
    // compute grid assumed to be the same as container,
    // change in scatterParticle()  if necessary.
    dem->setGrid(Box(x1, y1, z1, x2, y2, z2)); 
    pd->setGrid(Box(x1, y1, z1, x2, y2, z2)); 

    dem->readParticles(InputParameter::get().datafile["particleFile"]);
    pd->readPeriDynamicsData(InputParameter::get().datafile["periFile"]);

    dem->openCompressProg(progressInf, "rigidInc_progress");
    pd->openPeriProgress(periProgInf, "rigidInc_peri.dat");
    pd->openPeriProgress(periProgInfHalf, "rigidInc_peri_half.dat");

    // volume and horizon size are related to the mesh, July 14, 2014
    pd->calcParticleVolume(); 
    // so these two should be calculated before  peri-points
    // that are within sand particles are deleted
    pd->calcHorizonSize();    
    // delete those peri-points that are inside sand particles
    pd->removeInsidePeriParticles(dem->getParticleVec()); 
  }

  dem->scatterParticle();
  pd->scatterPeriParticle(dem->getAllContainer());

  // construct the Matrix members in periParticleVec
  pd->constructPeriMatrix(); 
  pd->constructNeighbor();

  // if(dem->getMPIRank()==0){
  // pd->printPeriDomain("peridomain_rank_0.dat");
  // pd->printRecvPeriDomain("peridomain_recv_rank_0.dat");    // bondVec
  // is not
  // empty
  //}
  // if(dem->getMPIRank()==1){
  // pd->printPeriDomain("peridomain_rank_1.dat");
  // pd->printRecvPeriDomain("peridomain_recv_rank_1.dat");    // bondVec
  // is empty
  //}

  auto startStep = util::getParam<std::size_t>("startStep");
  auto endStep = util::getParam<std::size_t>("endStep");
  auto startSnap = util::getParam<std::size_t>("startSnap");
  auto endSnap = util::getParam<std::size_t>("endSnap");
  std::size_t netStep = endStep - startStep + 1;
  std::size_t netSnap = endSnap - startSnap + 1;
  timeStep = util::getParam<REAL>("timeStep");

  // REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
  iteration = startStep;
  std::size_t iterSnap = startSnap;
  std::string outputFolder(".");
  if (dem->getMPIRank() == 0) {

    // Create the output writer in the master process
    // <outputFolder> rigidInc.pe3d </outputFolder>
    auto folderName =  dem::InputParameter::get().datafile["outputFolder"];
    outputFolder = util::createOutputFolder(folderName);
    //std::cout << "Output folder = " << outputFolder << "\n";
    dem->createOutputWriter(outputFolder, iterSnap-1);

    dem->writeGridToFile();
    dem->writeParticlesToFile(iterSnap);
    pd->printPeriProgress(periProgInf, 0);
    pd->printPeriProgressHalf(periProgInfHalf, 0);
  }
  //    if (dem->getMPIRank() == 0)
  //      debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT"
  //      << std::setw(OWID) << "migraT"
  //           << std::setw(OWID) << "totalT" << std::setw(OWID) << "overhead%"
  //           << std::endl;

  // Broadcast the output folder to all processes
  broadcast(dem->getMPIWorld(), outputFolder, 0);

  // the commuPeriParticle() have to be called after  commuParticle()
  dem->commuParticle(); 
  pd->commuPeriParticle(iteration, dem->getGradation().getPtclMaxRadius());

  // construct the Matrix members in recvPeriParticle
  pd->constructRecvPeriMatrix(); 

  // find the peri-bonds between periParticleVec and recvPeriParticle
  pd->findRecvPeriBonds(); 
  pd->calcParticleKinv();
  for (auto& pt : pd->getPeriParticleVec()) {
    pt->initial();
  }
  for (auto& pt : pd->getRecvPeriParticleVec()) {
    pt->initial();
  }
  pd->calcParticleStress();
  pd->calcParticleAcceleration();
  dem->releaseRecvParticle();
  //pd->releaseRecvPeriParticle();
  //    releasePeriBondVec();    // free the bondVec in periParticleVec and
  //    bondVec.clear() before constructNeighbor() again
  pd->findBoundaryPeriParticles();
  pd->clearPeriDEMBonds();
  pd->findPeriDEMBonds(dem->getMergedParticleVec());
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
    dem->findContact();
    if (dem->isBdryProcess())
      dem->findBdryContact();

    dem->clearContactForce();
    dem->internalForce();
    if (dem->isBdryProcess())
      dem->boundaryForce();
    pd->runFirstHalfStep();

    pd->applyPeriBoundaryCondition();
    //      commuPeriParticle();
    //      constructRecvPeriMatrix();
    //      constructNeighbor();
    //      checkBondParticleAlive();
    pd->findRecvPeriBonds();

    pd->calcParticleStress();

    // peri-point acceleration is set zero at the beginning
    pd->calcParticleAcceleration(); 

    // based on current states of DEM particles and peri-points
    pd->findPeriDEMBonds(dem->getMergedParticleVec()); 
    pd->applyCoupledForces();
    pd->runSecondHalfStep();

    dem->updateParticle();
    dem->gatherBdryContact(); // must call before updateBoundary
    //      updateBoundary(sigmaConf, "triaxial");
    //      updateGrid();
    //      updatePeriGrid();

    if (iteration % (netStep / netSnap) == 0) {
      // time1 = MPI_Wtime();
      dem->gatherParticle();
      pd->gatherPeriParticle();
      dem->gatherEnergy();
      // time2 = MPI_Wtime(); gatherT = time2 - time1;

      //    checkBondParticleAlive();
      //    findPeriDEMBonds(dem->getMergedParticleVec());    // update peri-dem bonds

      if (dem->getMPIRank() == 0) {
        dem->updateFileNames(iterSnap);
        dem->writeBoundaryToFile();
        dem->writeGridToFile();
        dem->writeParticlesToFile(iterSnap);
        dem->printBdryContact();
        dem->printBoundary();
        // printCompressProg(progressInf, distX, distY, distZ); // redundant
        pd->printPeriProgress(periProgInf, iterSnap);
        pd->printPeriProgressHalf(periProgInfHalf, iterSnap);
      }
      dem->printContact(combine(".", "rigidInc_contact_", iterSnap, 3));
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
    //      if (dem->getMPIRank() == 0 && (iteration+1 ) % (netStep /
    //      netSnap) == 0) //
    //      ignore gather and print time at this step
    //    debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT
    //    << std::setw(OWID) << migraT
    //         << std::setw(OWID) << totalT << std::setw(OWID) << (commuT +
    //         migraT)/totalT*100 << std::endl;

    if (dem->getMPIRank() == 0 && iteration % 10 == 0)
      dem->printCompressProg(progressInf, distX, distY, distZ);

    // no break condition, just through top/bottom displacement control
    ++iteration;
  }

  if (dem->getMPIRank() == 0) {
    dem->updateFileNames(iterSnap, ".end");
    dem->writeParticlesToFile(iterSnap);
    dem->printBdryContact();
    dem->printBoundary();
    dem->printCompressProg(progressInf, distX, distY, distZ);
    //      dem->printPeriProgress(periProgInf, iterSnap);
  }
  if (dem->getMPIRank() == 0) {
    dem->closeProg(progressInf);
    dem->closeProg(periProgInf);
    dem->closeProg(periProgInfHalf);
  }
}
