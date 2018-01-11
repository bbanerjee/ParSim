#include <Simulations/DEM_PD/PeridynamicsRigidInclusion.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;

using DiscreteElements = dem::DiscreteElements;
using Peridynamics = pd::Peridynamics;

void
PeridynamicsRigidInclusion::execute(DiscreteElements* dem, Peridynamics* pd)
{
  std::ofstream progressInf;
  std::ofstream periProgInf;

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
    dem->setSpatialDomain(Box(x1, y1, z1, x2, y2, z2));

    // compute patchGrid assumed to be
    // the same as domain,
    // change in scatterParticles()
    // if necessary.
    dem->setPatchBox(Box(x1, y1, z1, x2, y2, z2)); 
    pd->setPatchBox(Box(x1, y1, z1, x2, y2, z2)); 

    dem->readParticles(util::getFilename("particleFilename"));
    pd->readPeriDynamicsData(InputParameter::get().datafile["periFilename"]);

    dem->openProgressOutputFile(progressInf, "rigidInc_progress");
    pd->openPeriProgress(periProgInf, "rigidInc_peri.dat");

    // volume and horizon size are related to
    // the mesh,
    // July 14, 2014
    //proc0cout << "**NOTICE** Computing particle volumes\n";
    pd->calcParticleVolume(); 

    // so these two should be dem->calculated
    // before peri-points
    // that are within sand particles are deleted
    //proc0cout << "**NOTICE** Computing particle horizon sizes\n";
    pd->calcHorizonSize(); 

    if (InputParameter::get().param["removePeriParticles"] == 1) {
      // delete those peri-points that are inside
      // sand particles
      //proc0cout << "**NOTICE** Removing inside peri particles\n";
      pd->removeOverlappingParticles(dem->getModifiableAllParticleVec(), true);
    } else {
      // delete those dem particles that contain peri particles
      //proc0cout << "**NOTICE** Removing inside dem particles\n";
      pd->removeOverlappingParticles(dem->getModifiableAllParticleVec(), false);
    }
  }

  //proc0cout << "**NOTICE** Scattering dem particles\n";
  dem->scatterParticles();

  //proc0cout << "**NOTICE** Scattering pd particles\n";
  pd->scatterPeriParticle(dem->getSpatialDomain());

  // dem->construct the Matrix members in
  // periParticleVec
  //proc0cout << "**NOTICE** Constucting peri matrix\n";
  pd->constructPeriMatrix(); 

  //proc0cout << "**NOTICE** Constucting peri neighbor\n";
  pd->constructNeighbor();

  // if(dem->getMPIRank()==0){
  // dem->printPeriDomain("peridomain_rank_0.dat");
  // dem->printRecvPeriDomain("peridomain_recv_rank_0.dat");    // bondVec
  // is not
  // empty
  //}
  // if(dem->getMPIRank()==1){
  // dem->printPeriDomain("peridomain_rank_1.dat");
  // dem->printRecvPeriDomain("peridomain_recv_rank_1.dat");    // bondVec
  // is empty
  //}

  auto startStep = util::getParam<std::size_t>("startStep");
  auto endStep = util::getParam<std::size_t>("endStep");
  auto startSnap = util::getParam<std::size_t>("startSnap");
  auto endSnap = util::getParam<std::size_t>("endSnap");
  std::size_t netStep = endStep - startStep + 1;
  std::size_t netSnap = endSnap - startSnap + 1;
  auto timeStep = util::getParam<REAL>("timeStep");
  auto currentTime = startStep * timeStep;

  // REAL time0, time1, time2, migraT, gatherT, totalT;
  auto iteration = startStep;
  std::size_t iterSnap = startSnap;
  std::string outputFolder(".");
  if (dem->getMPIRank() == 0) {

    // Create the output writer in the master process
    // <outputFolder> couple.pe3d </outputFolder>
    auto folderName =  dem::InputParameter::get().datafile["outputFolder"];
    outputFolder = util::createOutputFolder(folderName);
    //std::cout << "Output folder = " << outputFolder << "\n";
    dem->createOutputWriter(outputFolder, iterSnap-1);
    pd->createOutputWriter(outputFolder, iterSnap-1);

    dem->writePatchGridToFile(currentTime);
    dem->writeParticlesToFile(iterSnap, currentTime);
    pd->writeParticlesToFile(iterSnap, currentTime);
    pd->printPeriProgress(periProgInf, 0);
  }
  //    if (dem->getMPIRank() == 0)
  //      debugInf << std::dem->setw(OWID) << "iter" <<
  //      std::dem->setw(OWID) << "dem->commuT"
  //      << std::dem->setw(OWID) << "migraT"
  //           << std::dem->setw(OWID) << "totalT" <<
  //           std::dem->setw(OWID) << "overhead%"
  //           << std::endl;

  // Broadcast the output folder to all processes
  broadcast(dem->getMPIWorld(), outputFolder, 0);
  dem->printContact();

  //proc0cout << "**NOTICE** Commun dem particle\n";
  dem->communicateGhostParticles();     

  //proc0cout << "**NOTICE** Commun pd particle\n";
  pd->commuPeriParticle(iteration, dem->getGradation().getPtclMaxRadius()); 

  //proc0cout << "**NOTICE** Construct receive pd matrix\n";
  pd->constructRecvPeriMatrix(); 

  // if(dem->getMPIRank()==0){
  // dem->printPeriDomain("peridomain_rank_0.dat");
  // dem->printRecvPeriDomain("peridomain_recv_rank_0.dat");    // bondVec
  // is not
  // empty
  //}
  // if(dem->getMPIRank()==1){
  // dem->printPeriDomain("peridomain_rank_1.dat");
  // dem->printRecvPeriDomain("peridomain_recv_rank_1.dat");    // bondVec
  // is empty
  //}

  //pd->findRecvPeriBonds(); 

  //proc0cout << "**NOTICE** Calc pd Kinv\n";
  pd->calcParticleKinv();

  //proc0cout << "**NOTICE** Initial PD\n";
  for (auto& pt : pd->getPeriParticleVec()) {
    pt->initial();
  }
  for (auto& pt : pd->getRecvPeriParticleVec()) {
    pt->initial();
  }

  //proc0cout << "**NOTICE** PD stress\n";
  pd->calcParticleStress();

  //proc0cout << "**NOTICE** PD acc\n";
  pd->calcParticleAcceleration();

  //proc0cout << "**NOTICE** DEM release\n";
  dem->releaseReceivedParticles();

  //proc0cout << "**NOTICE** PD release\n";
  pd->releaseRecvPeriParticle();

  //proc0cout << "**NOTICE** PD release bonds\n";
  pd->releasePeriBondVec(); // free the bondVec in periParticleVec and
  // bondVec.clear() before dem->constructNeighbor() again

  /*
        std::ofstream ofs_top("topBoundaryInnerPeriPoints.dat");
        if(!ofs_top) {
              //std::cout << "stream error!" << std::endl; exit(-1);
        }

      ofs_top.dem->setf(std::ios::scientific, std::ios::floatfield);
      ofs_top.precision(10);
      for(PeriParticlePArray::const_iterator pt = topBoundaryInnerVec.begin();
  pt!= topBoundaryInnerVec.end(); pt++) {

          ofs_top << std::dem->setw(20) << (*pt)->getInitPosition().x()
  +
  (*pt)->getDisplacement().x()
          << std::dem->setw(20) << (*pt)->getInitPosition().y() +
  (*pt)->getDisplacement().y()
          << std::dem->setw(20) << (*pt)->getInitPosition().z() +
  (*pt)->getDisplacement().z()
          << std::endl;
          ofs_top.flush();
      }
      ofs_top.close();

        std::ofstream ofs_bot("topBoundaryEdgePeriPoints.dat");
        if(!ofs_bot) {
              //std::cout << "stream error!" << std::endl; exit(-1);
        }

      ofs_bot.dem->setf(std::ios::scientific, std::ios::floatfield);
      ofs_bot.precision(10);
      for(PeriParticlePArray::const_iterator pt = topBoundaryEdgeVec.begin();
  pt!= topBoundaryEdgeVec.end(); pt++) {

          ofs_bot << std::dem->setw(20) << (*pt)->getInitPosition().x()
  +
  (*pt)->getDisplacement().x()
          << std::dem->setw(20) << (*pt)->getInitPosition().y() +
  (*pt)->getDisplacement().y()
          << std::dem->setw(20) << (*pt)->getInitPosition().z() +
  (*pt)->getDisplacement().z()
          << std::endl;
          ofs_bot.flush();
      }
      ofs_bot.close();


       std::ofstream ofs_fix("topBoundaryCornerPeriPoints.dat");
        if(!ofs_fix) {
              //std::cout << "stream error!" << std::endl; exit(-1);
        }

      ofs_fix.dem->setf(std::ios::scientific, std::ios::floatfield);
      ofs_fix.precision(10);
      for(PeriParticlePArray::const_iterator pt = topBoundaryCornerVec.begin();
  pt!= topBoundaryCornerVec.end(); pt++) {

          ofs_fix << std::dem->setw(20) << (*pt)->getInitPosition().x()
  +
  (*pt)->getDisplacement().x()
          << std::dem->setw(20) << (*pt)->getInitPosition().y() +
  (*pt)->getDisplacement().y()
          << std::dem->setw(20) << (*pt)->getInitPosition().z() +
  (*pt)->getDisplacement().z()
          << std::endl;
          ofs_fix.flush();
      }
      ofs_fix.close();


       std::ofstream ofs2("bottomBoundaryPoints.dat");
        if(!ofs2) {
              //std::cout << "stream error!" << std::endl; exit(-1);
        }

      ofs2.dem->setf(std::ios::scientific, std::ios::floatfield);
      ofs2.precision(10);
      for(PeriParticlePArray::const_iterator pt = bottomBoundaryVec.begin();
  pt!= bottomBoundaryVec.end(); pt++) {

          ofs2<< std::dem->setw(20) << (*pt)->getInitPosition().x() +
  (*pt)->getDisplacement().x()
          << std::dem->setw(20) << (*pt)->getInitPosition().y() +
  (*pt)->getDisplacement().y()
          << std::dem->setw(20) << (*pt)->getInitPosition().z() +
  (*pt)->getDisplacement().z()
          << std::endl;
          ofs2.flush();
      }
      ofs2.close();


       std::ofstream ofs3("leftxBoundaryPoints.dat");
        if(!ofs3) {
              //std::cout << "stream error!" << std::endl; exit(-1);
        }

      ofs3.dem->setf(std::ios::scientific, std::ios::floatfield);
      ofs3.precision(10);
      for(PeriParticlePArray::const_iterator pt = leftBoundaryVec.begin(); pt!=
  leftBoundaryVec.end(); pt++) {

          ofs3<< std::dem->setw(20) << (*pt)->getInitPosition().x() +
  (*pt)->getDisplacement().x()
          << std::dem->setw(20) << (*pt)->getInitPosition().y() +
  (*pt)->getDisplacement().y()
          << std::dem->setw(20) << (*pt)->getInitPosition().z() +
  (*pt)->getDisplacement().z()
          << std::endl;
          ofs3.flush();
      }
      ofs3.close();


       std::ofstream ofs4("frontyBoundaryPoints.dat");
        if(!ofs4) {
              //std::cout << "stream error!" << std::endl; exit(-1);
        }

      ofs4.dem->setf(std::ios::scientific, std::ios::floatfield);
      ofs4.precision(10);
      for(PeriParticlePArray::const_iterator pt = frontBoundaryVec.begin(); pt!=
  frontBoundaryVec.end(); pt++) {

          ofs4<< std::dem->setw(20) << (*pt)->getInitPosition().x() +
  (*pt)->getDisplacement().x()
          << std::dem->setw(20) << (*pt)->getInitPosition().y() +
  (*pt)->getDisplacement().y()
          << std::dem->setw(20) << (*pt)->getInitPosition().z() +
  (*pt)->getDisplacement().z()
          << std::endl;
          ofs4.flush();
      }
      ofs4.close();

  if(dem->getMPIRank()==0){
       std::ofstream ofs5("fixedPeriPoints_rank_0.dat");
        if(!ofs5) {
              //std::cout << "stream error!" << std::endl; exit(-1);
        }

      ofs5.dem->setf(std::ios::scientific, std::ios::floatfield);
      ofs5.precision(10);
      for(PeriParticlePArray::const_iterator pt = fixedPeriParticleVec.begin();
  pt!= fixedPeriParticleVec.end(); pt++) {

          ofs5<< std::dem->setw(20) << (*pt)->getInitPosition().x() +
  (*pt)->getDisplacement().x()
          << std::dem->setw(20) << (*pt)->getInitPosition().y() +
  (*pt)->getDisplacement().y()
          << std::dem->setw(20) << (*pt)->getInitPosition().z() +
  (*pt)->getDisplacement().z()
          << std::endl;
          ofs5.flush();
      }
      ofs5.close();
  }
  if(dem->getMPIRank()==1){
       std::ofstream ofs5("fixedPeriPoints_rank_1.dat");
        if(!ofs5) {
              //std::cout << "stream error!" << std::endl; exit(-1);
        }

      ofs5.dem->setf(std::ios::scientific, std::ios::floatfield);
      ofs5.precision(10);
      for(PeriParticlePArray::const_iterator pt = fixedPeriParticleVec.begin();
  pt!= fixedPeriParticleVec.end(); pt++) {

          ofs5<< std::dem->setw(20) << (*pt)->getInitPosition().x() +
  (*pt)->getDisplacement().x()
          << std::dem->setw(20) << (*pt)->getInitPosition().y() +
  (*pt)->getDisplacement().y()
          << std::dem->setw(20) << (*pt)->getInitPosition().z() +
  (*pt)->getDisplacement().z()
          << std::endl;
          ofs5.flush();
      }
      ofs5.close();
  }

        std::ofstream ofs6("periParticleVec.dat");
        if(!ofs6) {
              //std::cout << "stream error!" << std::endl; exit(-1);
        }

      ofs6.dem->setf(std::ios::scientific, std::ios::floatfield);
      ofs6.precision(10);
      for(PeriParticlePArray::const_iterator pt = periParticleVec.begin(); pt!=
  periParticleVec.end(); pt++) {

          ofs6 << std::dem->setw(20) << (*pt)->getInitPosition().x() +
  (*pt)->getDisplacement().x()
          << std::dem->setw(20) << (*pt)->getInitPosition().y() +
  (*pt)->getDisplacement().y()
          << std::dem->setw(20) << (*pt)->getInitPosition().z() +
  (*pt)->getDisplacement().z()
          << std::endl;
          ofs6.flush();
      }
      ofs6.close();
  */
  while (iteration <= endStep) {

    proc0cout << "**NOTICE** iteration = " << iteration << "\n";

    // boundary peri-points are determined by
    // their initial positions, however they are
    // still needed to be determined
    //proc0cout << "**NOTICE** Find boundary PD particle\n";
    pd->findBoundaryPeriParticles(); 

    // every step since some peri-points can pass
    // from one cpu to another
    //proc0cout << "**NOTICE** Find fixed PD particle\n";
    pd->findFixedPeriParticles(); 
    // dem->commuT = migraT = gatherT = totalT = 0; time0 = MPI_Wtime();

    //proc0cout << "**NOTICE** DEM commu\n";
    dem->communicateGhostParticles();
    // time2 = MPI_Wtime(); dem->commuT = time2 - time0;

    // displacement control relies on constant time step, so do not call
    // dem->calcTimeStep(timeStep).
    // dem->calcTimeStep(timeStep); // use values from last step, must call before
    // dem->findConact
    //proc0cout << "**NOTICE** DEM findContact\n";
    dem->findContact(iteration);
    if (dem->isBoundaryProcess())
      dem->findBoundaryContacts(iteration);

    dem->initializeForces();

    //proc0cout << "**NOTICE** DEM internal force\n";
    dem->internalForce(timeStep, iteration);
    if (dem->isBoundaryProcess())
      dem->boundaryForce(timeStep, iteration);

    //proc0cout << "**NOTICE** PD first half step\n";
    pd->runFirstHalfStep();

    pd->applyPeriBoundaryCondition();

    //proc0cout << "**NOTICE** PD Commun\n";
    pd->commuPeriParticle(iteration, dem->getGradation().getPtclMaxRadius()); 
    pd->constructRecvPeriMatrix();
    pd->constructNeighbor();
    //pd->findRecvPeriBonds();

    //proc0cout << "**NOTICE** PD stress\n";
    pd->calcParticleStress();
    pd->calcParticleAcceleration(); 

    pd->applyTractionBoundary(iteration - startStep);

    //proc0cout << "**NOTICE** PD second half step\n";
    pd->runSecondHalfStep();

    dem->updateParticles(timeStep, iteration);

    //proc0cout << "**NOTICE** DEM gather boundary contact\n";
    dem->gatherBoundaryContacts(); // must call before updateBoundary
    //      updateBoundary(sigmaConf, "triaxial");
    //      updatePatchBox();
    pd->updatePeriPatchGrid(pd->getPeriParticleVec());

    if (iteration % (netStep / netSnap) == 0) {
      // time1 = MPI_Wtime();
      dem->gatherParticles();
      pd->gatherPeriParticle();
      dem->gatherEnergy();
      // time2 = MPI_Wtime(); gatherT = time2 - time1;

      if (dem->getMPIRank() == 0) {
        dem->updateFileNames(iterSnap);
        pd->updateFileNames(iterSnap);
        dem->writeBoundaryToFile(currentTime);
        dem->writePatchGridToFile(currentTime);
        dem->writeParticlesToFile(iterSnap, currentTime);
        pd->writeParticlesToFile(iterSnap, currentTime);
        dem->printBoundaryContacts();
        dem->printBoundary();
        // dem->appendToProgressOutputFile(progressInf, iteration, timeStep, distX, distY, distZ); //
        // redundant
        pd->printPeriProgress(periProgInf, iterSnap);
      }
      dem->printContact();
      ++iterSnap;
    }
    dem->releaseReceivedParticles();     // late dem->release because
                                         // dem->printContact refers to
                                         // received particles
    pd->releaseRecvPeriParticle(); // dem->release
                                         // recvPeriParticleVec and also
                                         // dem->remove
                                         // peri-bonds between recvPeriParticle
    pd->releasePeriBondVec();
    // time1 = MPI_Wtime();

    //proc0cout << "**NOTICE** DEM migrate\n";
    dem->migrateParticles(iteration);

    //proc0cout << "**NOTICE** PD migrate\n";
    pd->migratePeriParticle(iteration);
    // time2 = MPI_Wtime(); migraT = time2 - time1; totalT = time2 - time0;
    //      if (dem->getMPIRank() == 0 && (iteration+1 ) % (netStep /
    //      netSnap) == 0) //
    //      ignore gather and dem->print time at this step
    //    debugInf << std::dem->setw(OWID) << iteration <<
    //    std::dem->setw(OWID) << dem->commuT
    //    << std::dem->setw(OWID) << migraT
    //         << std::dem->setw(OWID) << totalT <<
    //         std::dem->setw(OWID) << (dem->commuT +
    //         migraT)/totalT*100 << std::endl;

    if (dem->getMPIRank() == 0 && iteration % 10 == 0)
      dem->appendToProgressOutputFile(progressInf, iteration, timeStep, distX, distY, distZ);

    // no break condition, just through top/bottom displacement control
    ++iteration;
    currentTime += timeStep;
  }

  if (dem->getMPIRank() == 0) {
    dem->updateFileNames(iterSnap, ".end");
    dem->writeParticlesToFile(iterSnap, currentTime);
    pd->writeParticlesToFile(iterSnap, currentTime);
    dem->printBoundaryContacts();
    dem->printBoundary();
    dem->appendToProgressOutputFile(progressInf, iteration, timeStep, distX, distY, distZ);
    //      dem->printPeriProgress(periProgInf, iterSnap);
  }
  if (dem->getMPIRank() == 0) {
    dem->closeProgressOutputFile(progressInf);
    dem->closeProgressOutputFile(periProgInf);
  }
}
