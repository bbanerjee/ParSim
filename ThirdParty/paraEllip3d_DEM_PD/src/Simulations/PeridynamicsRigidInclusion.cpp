#include <Simulations/PeridynamicsRigidInclusion.h>
#include <Core/Util/Utility.h>

using namespace dem;
using util::combine;
void
PeridynamicsRigidInclusion::execute(Assembly* assembly)
{
  std::ofstream progressInf;
  std::ofstream periProgInf;

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
    assembly->setGrid(
      Box(x1, y1, z1, x2, y2, z2)); // compute grid assumed to be
                                    // the same as container,
                                    // change in scatterParticle()
                                    // if necessary.
    assembly->readParticles(
      dem::Parameter::getSingleton().datafile["particleFile"].c_str());
    assembly->readPeriDynamicsData(
      dem::Parameter::getSingleton().datafile["periFile"].c_str());
    assembly->openCompressProg(progressInf, "rigidInc_progress");
    assembly->openPeriProgress(periProgInf, "rigidInc_peri.dat");
    assembly->calcParticleVolume(); // volume and horizon size are related to
                                    // the mesh,
                                    // July 14, 2014
    assembly->calcHorizonSize(); // so these two should be assembly->calculated
                                 // before peri-points
                                 // that are within sand particles are deleted
    assembly
      ->removeInsidePeriParticles(); // delete those peri-points that are inside
                                     // sand particles
  }
  assembly->scatterDEMPeriParticle();
  assembly->constructPeriMatrix(); // assembly->construct the Matrix members in
                                   // periParticleVec
  assembly->constructNeighbor();

  // if(assembly->getMPIRank()==0){
  // assembly->printPeriDomain("peridomain_rank_0.dat");
  // assembly->printRecvPeriDomain("peridomain_recv_rank_0.dat");    // bondVec
  // is not
  // empty
  //}
  // if(assembly->getMPIRank()==1){
  // assembly->printPeriDomain("peridomain_rank_1.dat");
  // assembly->printRecvPeriDomain("peridomain_recv_rank_1.dat");    // bondVec
  // is empty
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

  // REAL time0, time1, time2, migraT, gatherT, totalT;
  iteration = startStep;
  std::size_t iterSnap = startSnap;
  if (assembly->getMPIRank() == 0) {
    assembly->plotGrid(combine("rigidInc_gridplot_", iterSnap - 1, 3) + ".dat");
    assembly->printParticle(combine("rigidInc_particle_", iterSnap - 1, 3));
    assembly->printPeriProgress(periProgInf, 0);
  }
  //    if (assembly->getMPIRank() == 0)
  //      debugInf << std::assembly->setw(OWID) << "iter" <<
  //      std::assembly->setw(OWID) << "assembly->commuT"
  //      << std::assembly->setw(OWID) << "migraT"
  //           << std::assembly->setw(OWID) << "totalT" <<
  //           std::assembly->setw(OWID) << "overhead%"
  //           << std::endl;

  assembly->commuParticle();     // the assembly->commuPeriParticle() have to be
                                 // called after
                                 // assembly->commuParticle(), since
  assembly->commuPeriParticle(); // rankX1... are assembly->calculated in
                                 // assembly->commuParticle()
  assembly
    ->constructRecvPeriMatrix(); // assembly->construct the Matrix members in
                                 // recvPeriParticle

  // if(assembly->getMPIRank()==0){
  // assembly->printPeriDomain("peridomain_rank_0.dat");
  // assembly->printRecvPeriDomain("peridomain_recv_rank_0.dat");    // bondVec
  // is not
  // empty
  //}
  // if(assembly->getMPIRank()==1){
  // assembly->printPeriDomain("peridomain_rank_1.dat");
  // assembly->printRecvPeriDomain("peridomain_recv_rank_1.dat");    // bondVec
  // is empty
  //}

  assembly->findRecvPeriBonds(); // assembly->find the peri-bonds between
                                 // periParticleVec and
                                 // recvPeriParticle
  assembly->calcParticleKinv();
  for (auto& pt : assembly->getPeriParticleVec()) {
    pt->initial();
  }
  for (auto& pt : assembly->getRecvPeriParticleVec()) {
    pt->initial();
  }
  assembly->calcParticleStress();
  assembly->calcParticleAcceleration();
  assembly->releaseRecvParticle();
  assembly->releaseRecvPeriParticle();
  assembly->releasePeriBondVec(); // free the bondVec in periParticleVec and
  // bondVec.clear() before assembly->constructNeighbor() again

  /*
        std::ofstream ofs_top("topBoundaryInnerPeriPoints.dat");
        if(!ofs_top) {
              std::cout << "stream error!" << std::endl; exit(-1);
        }

      ofs_top.assembly->setf(std::ios::scientific, std::ios::floatfield);
      ofs_top.precision(10);
      for(PeriParticlePArray::const_iterator pt = topBoundaryInnerVec.begin();
  pt!= topBoundaryInnerVec.end(); pt++) {

          ofs_top << std::assembly->setw(20) << (*pt)->getInitPosition().x()
  +
  (*pt)->getDisplacement().x()
          << std::assembly->setw(20) << (*pt)->getInitPosition().y() +
  (*pt)->getDisplacement().y()
          << std::assembly->setw(20) << (*pt)->getInitPosition().z() +
  (*pt)->getDisplacement().z()
          << std::endl;
          ofs_top.flush();
      }
      ofs_top.close();

        std::ofstream ofs_bot("topBoundaryEdgePeriPoints.dat");
        if(!ofs_bot) {
              std::cout << "stream error!" << std::endl; exit(-1);
        }

      ofs_bot.assembly->setf(std::ios::scientific, std::ios::floatfield);
      ofs_bot.precision(10);
      for(PeriParticlePArray::const_iterator pt = topBoundaryEdgeVec.begin();
  pt!= topBoundaryEdgeVec.end(); pt++) {

          ofs_bot << std::assembly->setw(20) << (*pt)->getInitPosition().x()
  +
  (*pt)->getDisplacement().x()
          << std::assembly->setw(20) << (*pt)->getInitPosition().y() +
  (*pt)->getDisplacement().y()
          << std::assembly->setw(20) << (*pt)->getInitPosition().z() +
  (*pt)->getDisplacement().z()
          << std::endl;
          ofs_bot.flush();
      }
      ofs_bot.close();


       std::ofstream ofs_fix("topBoundaryCornerPeriPoints.dat");
        if(!ofs_fix) {
              std::cout << "stream error!" << std::endl; exit(-1);
        }

      ofs_fix.assembly->setf(std::ios::scientific, std::ios::floatfield);
      ofs_fix.precision(10);
      for(PeriParticlePArray::const_iterator pt = topBoundaryCornerVec.begin();
  pt!= topBoundaryCornerVec.end(); pt++) {

          ofs_fix << std::assembly->setw(20) << (*pt)->getInitPosition().x()
  +
  (*pt)->getDisplacement().x()
          << std::assembly->setw(20) << (*pt)->getInitPosition().y() +
  (*pt)->getDisplacement().y()
          << std::assembly->setw(20) << (*pt)->getInitPosition().z() +
  (*pt)->getDisplacement().z()
          << std::endl;
          ofs_fix.flush();
      }
      ofs_fix.close();


       std::ofstream ofs2("bottomBoundaryPoints.dat");
        if(!ofs2) {
              std::cout << "stream error!" << std::endl; exit(-1);
        }

      ofs2.assembly->setf(std::ios::scientific, std::ios::floatfield);
      ofs2.precision(10);
      for(PeriParticlePArray::const_iterator pt = bottomBoundaryVec.begin();
  pt!= bottomBoundaryVec.end(); pt++) {

          ofs2<< std::assembly->setw(20) << (*pt)->getInitPosition().x() +
  (*pt)->getDisplacement().x()
          << std::assembly->setw(20) << (*pt)->getInitPosition().y() +
  (*pt)->getDisplacement().y()
          << std::assembly->setw(20) << (*pt)->getInitPosition().z() +
  (*pt)->getDisplacement().z()
          << std::endl;
          ofs2.flush();
      }
      ofs2.close();


       std::ofstream ofs3("leftxBoundaryPoints.dat");
        if(!ofs3) {
              std::cout << "stream error!" << std::endl; exit(-1);
        }

      ofs3.assembly->setf(std::ios::scientific, std::ios::floatfield);
      ofs3.precision(10);
      for(PeriParticlePArray::const_iterator pt = leftBoundaryVec.begin(); pt!=
  leftBoundaryVec.end(); pt++) {

          ofs3<< std::assembly->setw(20) << (*pt)->getInitPosition().x() +
  (*pt)->getDisplacement().x()
          << std::assembly->setw(20) << (*pt)->getInitPosition().y() +
  (*pt)->getDisplacement().y()
          << std::assembly->setw(20) << (*pt)->getInitPosition().z() +
  (*pt)->getDisplacement().z()
          << std::endl;
          ofs3.flush();
      }
      ofs3.close();


       std::ofstream ofs4("frontyBoundaryPoints.dat");
        if(!ofs4) {
              std::cout << "stream error!" << std::endl; exit(-1);
        }

      ofs4.assembly->setf(std::ios::scientific, std::ios::floatfield);
      ofs4.precision(10);
      for(PeriParticlePArray::const_iterator pt = frontBoundaryVec.begin(); pt!=
  frontBoundaryVec.end(); pt++) {

          ofs4<< std::assembly->setw(20) << (*pt)->getInitPosition().x() +
  (*pt)->getDisplacement().x()
          << std::assembly->setw(20) << (*pt)->getInitPosition().y() +
  (*pt)->getDisplacement().y()
          << std::assembly->setw(20) << (*pt)->getInitPosition().z() +
  (*pt)->getDisplacement().z()
          << std::endl;
          ofs4.flush();
      }
      ofs4.close();

  if(assembly->getMPIRank()==0){
       std::ofstream ofs5("fixedPeriPoints_rank_0.dat");
        if(!ofs5) {
              std::cout << "stream error!" << std::endl; exit(-1);
        }

      ofs5.assembly->setf(std::ios::scientific, std::ios::floatfield);
      ofs5.precision(10);
      for(PeriParticlePArray::const_iterator pt = fixedPeriParticleVec.begin();
  pt!= fixedPeriParticleVec.end(); pt++) {

          ofs5<< std::assembly->setw(20) << (*pt)->getInitPosition().x() +
  (*pt)->getDisplacement().x()
          << std::assembly->setw(20) << (*pt)->getInitPosition().y() +
  (*pt)->getDisplacement().y()
          << std::assembly->setw(20) << (*pt)->getInitPosition().z() +
  (*pt)->getDisplacement().z()
          << std::endl;
          ofs5.flush();
      }
      ofs5.close();
  }
  if(assembly->getMPIRank()==1){
       std::ofstream ofs5("fixedPeriPoints_rank_1.dat");
        if(!ofs5) {
              std::cout << "stream error!" << std::endl; exit(-1);
        }

      ofs5.assembly->setf(std::ios::scientific, std::ios::floatfield);
      ofs5.precision(10);
      for(PeriParticlePArray::const_iterator pt = fixedPeriParticleVec.begin();
  pt!= fixedPeriParticleVec.end(); pt++) {

          ofs5<< std::assembly->setw(20) << (*pt)->getInitPosition().x() +
  (*pt)->getDisplacement().x()
          << std::assembly->setw(20) << (*pt)->getInitPosition().y() +
  (*pt)->getDisplacement().y()
          << std::assembly->setw(20) << (*pt)->getInitPosition().z() +
  (*pt)->getDisplacement().z()
          << std::endl;
          ofs5.flush();
      }
      ofs5.close();
  }

        std::ofstream ofs6("periParticleVec.dat");
        if(!ofs6) {
              std::cout << "stream error!" << std::endl; exit(-1);
        }

      ofs6.assembly->setf(std::ios::scientific, std::ios::floatfield);
      ofs6.precision(10);
      for(PeriParticlePArray::const_iterator pt = periParticleVec.begin(); pt!=
  periParticleVec.end(); pt++) {

          ofs6 << std::assembly->setw(20) << (*pt)->getInitPosition().x() +
  (*pt)->getDisplacement().x()
          << std::assembly->setw(20) << (*pt)->getInitPosition().y() +
  (*pt)->getDisplacement().y()
          << std::assembly->setw(20) << (*pt)->getInitPosition().z() +
  (*pt)->getDisplacement().z()
          << std::endl;
          ofs6.flush();
      }
      ofs6.close();
  */
  while (iteration <= endStep) {
    assembly
      ->findBoundaryPeriParticles(); // boundary peri-points are determined by
    // their initial positions, however they are
    // still needed to be determined
    assembly
      ->findFixedPeriParticles(); // every step since some peri-points can pass
                                  // from one cpu to another
    // assembly->commuT = migraT = gatherT = totalT = 0; time0 = MPI_Wtime();
    assembly->commuParticle();
    // time2 = MPI_Wtime(); assembly->commuT = time2 - time0;

    // displacement control relies on constant time step, so do not call
    // assembly->calcTimeStep().
    // assembly->calcTimeStep(); // use values from last step, must call before
    // assembly->findConact
    assembly->findContact();
    if (assembly->isBdryProcess())
      assembly->findBdryContact();

    assembly->clearContactForce();
    assembly->internalForce();
    if (assembly->isBdryProcess())
      assembly->boundaryForce();
    assembly->runFirstHalfStep();

    assembly->applyPeriBoundaryCondition();
    assembly->commuPeriParticle();
    assembly->constructRecvPeriMatrix();
    assembly->constructNeighbor();
    assembly->findRecvPeriBonds();

    assembly->calcParticleStress();
    assembly->calcParticleAcceleration(); // peri-point acceleration is
                                          // assembly->set zero at the
                                          // beginning

    assembly->applyTractionBoundary(iteration - startStep);
    assembly->runSecondHalfStep();

    assembly->updateParticle();
    assembly->gatherBdryContact(); // must call before updateBoundary
    //      updateBoundary(sigmaConf, "triaxial");
    //      updateGrid();
    assembly->updatePeriGrid();

    if (iteration % (netStep / netSnap) == 0) {
      // time1 = MPI_Wtime();
      assembly->gatherParticle();
      assembly->gatherPeriParticle();
      assembly->gatherEnergy();
      // time2 = MPI_Wtime(); gatherT = time2 - time1;

      if (assembly->getMPIRank() == 0) {
        assembly->plotBoundary(combine( "rigidInc_bdryplot_", iterSnap, 3) + ".dat");
        assembly->plotGrid(combine( "rigidInc_gridplot_", iterSnap, 3) + ".dat");
        assembly->printParticle(combine( "rigidInc_particle_", iterSnap, 3));
        assembly->printBdryContact(combine( "rigidInc_bdrycntc_", iterSnap, 3));
        assembly->printBoundary(combine( "rigidInc_boundary_", iterSnap, 3));
        // assembly->printCompressProg(progressInf, distX, distY, distZ); //
        // redundant
        assembly->printPeriProgress(periProgInf, iterSnap);
      }
      assembly->printContact(combine( "rigidInc_contact_", iterSnap, 3));
      ++iterSnap;
    }
    assembly->releaseRecvParticle();     // late assembly->release because
                                         // assembly->printContact refers to
                                         // received particles
    assembly->releaseRecvPeriParticle(); // assembly->release
                                         // recvPeriParticleVec and also
                                         // assembly->remove
                                         // peri-bonds between recvPeriParticle
    assembly->releasePeriBondVec();
    // time1 = MPI_Wtime();
    assembly->migrateParticle();
    assembly->migratePeriParticle();
    // time2 = MPI_Wtime(); migraT = time2 - time1; totalT = time2 - time0;
    //      if (assembly->getMPIRank() == 0 && (iteration+1 ) % (netStep /
    //      netSnap) == 0) //
    //      ignore gather and assembly->print time at this step
    //    debugInf << std::assembly->setw(OWID) << iteration <<
    //    std::assembly->setw(OWID) << assembly->commuT
    //    << std::assembly->setw(OWID) << migraT
    //         << std::assembly->setw(OWID) << totalT <<
    //         std::assembly->setw(OWID) << (assembly->commuT +
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
  }
}
