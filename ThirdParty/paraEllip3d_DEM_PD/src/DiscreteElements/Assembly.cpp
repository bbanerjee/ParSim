//
//                                    x1(1)
//                          ---------------------
//                         /                    /|
//                        /                    / |
//                       /                    /  |
//                      /        z2(6)       /   |
//                     /                    /    |
//                    /                    /     |                    z
//                   /                    /      |                    |
//                  |---------------------       |                    |
//            y1(3) |                    | y2(4) |                    |____ y
//                  |                    |       /                   /
//                  |                    |      /                   /
//                  |        x2(2)       |     /                   x
//                  |                    |    /
//                  |                    |   /
//                  |                    |  /
//                  |                    | /
//                  |                    |/
//                  ----------------------
//                             z1(5)
//
//
//    It is preferable to use the description of surface x1, x2, y1, y2, z1, z2,
//    where x1 < x2, y1 < y2, z1 < z2. We also use surface 1, 2, 3, 4, 5, 6
//    accordingly.
//

#include <DiscreteElements/Assembly.h>
#include <DiscreteElements/Patch.h>

#include <Boundary/BoundaryReader.h>
#include <Boundary/CylinderBoundary.h>
#include <Boundary/PlaneBoundary.h>
#include <Core/Const/const.h>
#include <Core/Util/Utility.h>
#include <Core/Math/IntVec.h>
#include <InputOutput/ParticleFileReader.h>
#include <InputOutput/PeriParticleFileReader.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <string>
#include <utility>
#include <array>

//#define BINNING
//#define TIME_PROFILE

// static time_t timeStamp; // for file timestamping
// static struct timeval time_w1, time_w2; // for wall-clock time record
// static struct timeval time_p1, time_p2; // for internal wall-clock time
// profiling, can be used on any piece of code
// static struct timeval time_r1, time_r2; // for internal wall-clock time
// profiling for contact resolution only (excluding space search)

namespace dem {

using pd::PeriParticleP;
using pd::PeriParticlePArray;
using pd::PeriBondP;
using pd::PeriBondPArray;
using pd::PeriDEMBond;
using pd::PeriDEMBondP;
using pd::PeriParticleFileReader;

using util::timediff;
using util::timediffmsec;
using util::timediffsec;
using util::combine;
using Timer = std::chrono::steady_clock;
using Seconds = std::chrono::seconds;

void
Assembly::deposit(const std::string& boundaryFile,
                  const std::string& particleFile)
{
  // The output folder (default name is .)
  std::string outputFolder(".");

  // Read the input data
  if (mpiRank == 0) {
    readBoundary(boundaryFile);
    readParticles(particleFile);
    openDepositProg(progressInf, "deposit_progress");
  }

  //proc0cout << "**NOTICE** Before scatterparticle\n";
  scatterParticle(); // scatter particles only once; also updates grid for the
                     // first time

  auto startStep = util::getParam<std::size_t>("startStep");
  auto endStep = util::getParam<std::size_t>("endStep");
  auto startSnap = util::getParam<std::size_t>("startSnap");
  auto endSnap = util::getParam<std::size_t>("endSnap");

  auto netStep = endStep - startStep + 1;
  auto netSnap = endSnap - startSnap + 1;

  iteration = startStep;
  auto iterSnap = startSnap;
  if (mpiRank == 0) {

    // Create the output writer in the master process
    auto folderName =  dem::Parameter::get().datafile["outputFolder"];
    outputFolder = util::createOutputFolder(folderName);
    //std::cout << "Output folder = " << outputFolder << "\n";
    createOutputWriter(outputFolder, iterSnap-1);

    plotBoundary();
    plotGrid();
    plotParticle();
    printBdryContact();
    debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT"
             << std::setw(OWID) << "migraT" << std::setw(OWID) << "compuT"
             << std::setw(OWID) << "totalT" << std::setw(OWID) << "overhead%"
             << std::endl;
  }

  // Broadcast the output folder to all processes
  broadcast(boostWorld, outputFolder, 0);
  std::cerr << "Proc = " << mpiRank << " outputFolder = " << outputFolder << "\n";

  REAL time0, time1, time2, commuT, migraT, gatherT, totalT;
  REAL timeCount = 0;
  timeStep = util::getParam<REAL>("timeStep");
  timeAccrued = util::getParam<REAL>("timeAccrued");
  REAL timeTotal = timeAccrued + timeStep * netStep;

  while (timeAccrued < timeTotal) {

    bool toCheckTime = (iteration + 1) % (netStep / netSnap) == 0;

    commuT = migraT = gatherT = totalT = 0;
    time0 = MPI_Wtime();

    proc0cout << "**NOTICE** Time = " << timeAccrued << " iteration = " << iteration << "\n";
    commuParticle(iteration);

    if (toCheckTime)
      time2 = MPI_Wtime();
    commuT = time2 - time0;

    //proc0cout << "**NOTICE** Before calcTimeStep\n";
    calcTimeStep(); // use values from last step, must call before
                    // findContact (which clears data)

    //proc0cout << "**NOTICE** Before findContact\n";

    findContact();
    if (isBdryProcess())
      findBdryContact();

    clearContactForce();

    //proc0cout << "**NOTICE** Before internalForce\n";
    internalForce();

    if (isBdryProcess()) {
      //proc0cout << "**NOTICE** Before updateParticle\n";
      boundaryForce();
    }

    //proc0cout << "**NOTICE** Before updateParticle\n";
    updateParticle();
    updateGrid(); // universal; updateGridMaxZ() for deposition only

    /**/ timeCount += timeStep;
    /**/ timeAccrued += timeStep;
    /**/ if (timeCount >= timeTotal / netSnap) {
      // if (iteration % (netStep / netSnap) == 0) {
      if (toCheckTime)
        time1 = MPI_Wtime();

      gatherParticle();

      gatherBdryContact();
      gatherEnergy();
      if (toCheckTime)
        time2 = MPI_Wtime();
      gatherT = time2 - time1;

      if (mpiRank == 0) {
        updateFileNames(iterSnap);
        plotBoundary();
        plotGrid();
        plotParticle();
        printBdryContact();
        printDepositProg(progressInf);
      }
      printContact(combine(outputFolder, "contact_", iterSnap, 5));

      /**/ timeCount = 0;
      ++iterSnap;
    }

    //proc0cout << "**NOTICE** Before releaseRecvParticle\n";
    releaseRecvParticle(); // late release because printContact refers to
                           // received particles
    if (toCheckTime)
      time1 = MPI_Wtime();
    //proc0cout << "**NOTICE** Before migrateParticle\n";
    migrateParticle();
    if (toCheckTime)
      time2 = MPI_Wtime();
    migraT = time2 - time1;
    totalT = time2 - time0;
    if (mpiRank == 0 &&
        toCheckTime) // ignore gather and print time at this step
      debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT
               << std::setw(OWID) << migraT << std::setw(OWID)
               << totalT - commuT - migraT << std::setw(OWID) << totalT
               << std::setw(OWID) << (commuT + migraT) / totalT * 100
               << std::endl;
    ++iteration;
    //proc0cout << "**NOTICE** End of iteration: " << iteration << "\n";
  }

  if (mpiRank == 0)
    closeProg(progressInf);
}

void
Assembly::readBoundary(const std::string& fileName)
{
  // Check whether file is XML or JSON
  std::ifstream ifs(fileName);
  char firstChar;
  ifs >> firstChar;
  ifs.close();
  if (firstChar == '<') { // XML
    // Read the file
    //std::cout << "Using the XML reader\n";
    BoundaryReader reader;
    reader.readXML(fileName, allContainer, grid, boundaryVec);
  } else if (firstChar == '{') { // JSON
    //std::cout << "Using the JSON reader\n";
    BoundaryReader reader;
    reader.readJSON(fileName, allContainer, grid, boundaryVec);
  } else {
    //std::cout << "Using the default text reader\n";
    BoundaryReader reader;
    reader.read(fileName, allContainer, grid, boundaryVec);
  }
}

void
Assembly::readParticles(const std::string& particleFile)
{
  REAL young = util::getParam<REAL>("young");
  REAL poisson = util::getParam<REAL>("poisson");
  bool doInitialize = (util::getParam<int>("toInitParticle") == 1);

  ParticleFileReader reader;
  reader.read(particleFile, young, poisson, doInitialize, allParticleVec,
              gradation);
}

void
Assembly::scatterParticle()
{
  // partition particles and send to each process
  if (mpiRank == 0) { // process 0
    setGrid(Box(grid.getMinCorner().x(), grid.getMinCorner().y(),
                grid.getMinCorner().z(), grid.getMaxCorner().x(),
                grid.getMaxCorner().y(),
                getPtclMaxZ(allParticleVec) + gradation.getPtclMaxRadius()));

    Vec v1 = grid.getMinCorner();
    Vec v2 = grid.getMaxCorner();
    Vec vspan = (v2 - v1) / d_mpiProcs;

    //std::cout << "v1 = " << v1 << "\n";
    //std::cout << "v2 = " << v2 << "\n";
    //std::cout << "d_mpiProcs = " << d_mpiProcs << "\n";
    //std::cout << "vspan = " << vspan << "\n";

    auto reqs = new boost::mpi::request[mpiSize - 1];
    ParticlePArray tmpParticleVec;
    for (int iRank = mpiSize - 1; iRank >= 0; --iRank) {
      tmpParticleVec.clear(); // do not release memory!
      int ndim = 3;
      //int coords[3];
      IntVec coords;
      MPI_Cart_coords(cartComm, iRank, ndim, coords.data());
      //std::cout << "iRank = " << iRank 
      //          << " coords = " << coords << "\n";

      Vec lower = v1 + vspan*coords;
      Vec upper = lower + vspan;
      Box container(lower, upper);
      //std::cout << "lower = " << lower 
      //          << " upper = " << upper << "\n";

      findParticleInBox(container, allParticleVec, tmpParticleVec);

      if (iRank != 0)
        reqs[iRank - 1] = boostWorld.isend(iRank, mpiTag,
                                           tmpParticleVec); // non-blocking send
      if (iRank == 0) {
        particleVec.resize(tmpParticleVec.size());
        for (auto i = 0u; i < particleVec.size(); ++i)
          particleVec[i] = std::make_shared<Particle>(
            *tmpParticleVec[i]); // default synthesized copy constructor
      } // now particleVec do not share memeory with allParticleVec
    }
    boost::mpi::wait_all(reqs, reqs + mpiSize - 1); // for non-blocking send
    delete[] reqs;

  } else { // other processes except 0
    boostWorld.recv(0, mpiTag, particleVec);
  }

  // content of allParticleVec may need to be printed, so do not clear it.
  // if (mpiRank == 0) releaseGatheredParticle();

  // broadcast necessary info
  broadcast(boostWorld, gradation, 0);
  broadcast(boostWorld, boundaryVec, 0);
  broadcast(boostWorld, allContainer, 0);
  broadcast(boostWorld, grid, 0);

  // Create patch for the current process
  REAL ghostWidth = gradation.getPtclMaxRadius() * 2;
  createPatch(0, ghostWidth);
}

void 
Assembly::createPatch(int iteration, const REAL& ghostWidth) 
{
  // determine container of each process
  Vec v1 = grid.getMinCorner();
  Vec v2 = grid.getMaxCorner();
  Vec vspan = (v2 - v1) / d_mpiProcs;
  Vec lower = v1 + vspan * d_mpiCoords;
  Vec upper = lower + vspan;
  d_patchP = std::make_unique<Patch>(cartComm, mpiRank, d_mpiCoords,
                                      lower, upper, ghostWidth, EPS);
  //std::ostringstream out;
  //out << "mpiRank = " << mpiRank 
  //    << " iter = " << iteration << " d_mpiCoords = " << d_mpiCoords 
  //    << " lower = " << lower 
  //    << " upper = " << upper << "\n";
  //std::cout << out.str();
}

void 
Assembly::updatePatch(int iteration, const REAL& ghostWidth)
{
  // determine container of each process
  Vec v1 = grid.getMinCorner();
  Vec v2 = grid.getMaxCorner();
  Vec vspan = (v2 - v1) / d_mpiProcs;
  Vec lower = v1 + vspan * d_mpiCoords;
  Vec upper = lower + vspan;
  d_patchP->update(iteration, lower, upper, ghostWidth);

  //std::ostringstream out;
  //out << "mpiRank = " << mpiRank 
  //    << " iter = " << iteration << " d_mpiCoords = " << d_mpiCoords 
  //    << " lower = " << lower 
  //    << " upper = " << upper << "\n";
  //std::cout << out.str();
}

void
Assembly::commuParticle()
{
  commuParticle(-1);
}

void
Assembly::commuParticle(const int& iteration)
{
  // determine container of each process
  REAL ghostWidth = gradation.getPtclMaxRadius() * 2;
  updatePatch(iteration, ghostWidth);

  // duplicate pointers, pointing to the same memory
  mergeParticleVec.clear();
  mergeParticleVec = particleVec; 

  /*
  std::ostringstream out;
  if (mpiRank == 0) {
    out << "Rank: " << mpiRank << ": in: " << particleVec.size()
        << " merge: " << mergeParticleVec.size();
  }
  */

  d_patchP->sendRecvGhostXMinus(boostWorld, iteration, mergeParticleVec);
  d_patchP->sendRecvGhostXPlus(boostWorld, iteration, mergeParticleVec);
  d_patchP->waitToFinishX(iteration);
  d_patchP->combineReceivedParticlesX(iteration, mergeParticleVec);

  /*
  if (mpiRank == 0) {
    out << " -> " << mergeParticleVec.size();
  }
  */

  d_patchP->sendRecvGhostYMinus(boostWorld, iteration, mergeParticleVec);
  d_patchP->sendRecvGhostYPlus(boostWorld, iteration, mergeParticleVec);
  d_patchP->waitToFinishY(iteration);
  d_patchP->combineReceivedParticlesY(iteration, mergeParticleVec);

  /*
  if (mpiRank == 0) {
    out << " -> " << mergeParticleVec.size();
  }
  */

  d_patchP->sendRecvGhostZMinus(boostWorld, iteration, mergeParticleVec);
  d_patchP->sendRecvGhostZPlus(boostWorld, iteration, mergeParticleVec);
  d_patchP->waitToFinishZ(iteration);
  d_patchP->combineReceivedParticlesZ(iteration, mergeParticleVec);

  /*
  if (mpiRank == 0) {
    out << " -> " << mergeParticleVec.size();
  }
  */

  d_patchP->removeDuplicates(mergeParticleVec);

  /*
  if (mpiRank == 0) {
    out << " -> " << mergeParticleVec.size() << "\n";

    //std::ostringstream out;
    //out << "Rank: " << mpiRank << ": in: " << particleVec.size()
    //    << " recv: " << recvParticleVec.size();
    //out <<  ": out: " << particleVec.size()
    //    << " merge: " << mergeParticleVec.size() << "\n";
    std::cout << out.str();
  }
  */
}

/*
void
Assembly::commuParticle(const int& iteration)
{
*/
/*
  std::ostringstream out;

  // determine container of each process
  Vec v1 = grid.getMinCorner();
  Vec v2 = grid.getMaxCorner();
  Vec vspan = (v2 - v1) / d_mpiProcs;
  Vec lower = v1 + vspan * d_mpiCoords;
  Vec upper = lower + vspan;
  container = Box(lower, upper);
  out << "mpiRank = " << mpiRank 
      << " iter = " << iteration << " d_mpiCoords = " << d_mpiCoords 
      << " lower = " << lower 
      << " upper = " << upper << "\n";

  // find neighboring blocks
  rankX1 = -1;
  rankX2 = -1;
  rankY1 = -1;
  rankY2 = -1;
  rankZ1 = -1;
  rankZ2 = -1;
  rankX1Y1 = -1;
  rankX1Y2 = -1;
  rankX1Z1 = -1;
  rankX1Z2 = -1;
  rankX2Y1 = -1;
  rankX2Y2 = -1;
  rankX2Z1 = -1;
  rankX2Z2 = -1;
  rankY1Z1 = -1;
  rankY1Z2 = -1;
  rankY2Z1 = -1;
  rankY2Z2 = -1;
  rankX1Y1Z1 = -1;
  rankX1Y1Z2 = -1;
  rankX1Y2Z1 = -1;
  rankX1Y2Z2 = -1;
  rankX2Y1Z1 = -1;
  rankX2Y1Z2 = -1;
  rankX2Y2Z1 = -1;
  rankX2Y2Z2 = -1;
  // x1: -x direction
  IntVec neighborCoords = d_mpiCoords;
  --neighborCoords.x();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankX1);
  out << "mpiRank = " << mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankX1 = " << rankX1 << "\n";
  // x2: +x direction
  neighborCoords = d_mpiCoords;
  ++neighborCoords.x();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankX2);
  out << "mpiRank = " << mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankX2 = " << rankX2 << "\n";
  // y1: -y direction
  neighborCoords = d_mpiCoords;
  --neighborCoords.y();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankY1);
  out << "mpiRank = " << mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankY1 = " << rankY1 << "\n";
  // y2: +y direction
  neighborCoords = d_mpiCoords;
  ++neighborCoords.y();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankY2);
  out << "mpiRank = " << mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankY2 = " << rankY2 << "\n";
  // z1: -z direction
  neighborCoords = d_mpiCoords;
  --neighborCoords.z();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankZ1);
  out << "mpiRank = " << mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankZ1 = " << rankZ1 << "\n";
  // z2: +z direction
  neighborCoords = d_mpiCoords;
  ++neighborCoords.z();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankZ2);
  out << "mpiRank = " << mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankZ2 = " << rankZ2 << "\n";
  // x1y1
  neighborCoords = d_mpiCoords;
  --neighborCoords.x();
  --neighborCoords.y();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankX1Y1);
  out << "mpiRank = " << mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankX1Y1 = " << rankX1Y1 << "\n";
  // x1y2
  neighborCoords = d_mpiCoords;
  --neighborCoords.x();
  ++neighborCoords.y();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankX1Y2);
  out << "mpiRank = " << mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankX1Y2 = " << rankX1Y2 << "\n";
  // x1z1
  neighborCoords = d_mpiCoords;
  --neighborCoords.x();
  --neighborCoords.z();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankX1Z1);
  out << "mpiRank = " << mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankX1Z1 = " << rankX1Z1 << "\n";
  // x1z2
  neighborCoords = d_mpiCoords;
  --neighborCoords.x();
  ++neighborCoords.z();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankX1Z2);
  out << "mpiRank = " << mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankX1Z2 = " << rankX1Z2 << "\n";
  // x2y1
  neighborCoords = d_mpiCoords;
  ++neighborCoords.x();
  --neighborCoords.y();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankX2Y1);
  // x2y2
  neighborCoords = d_mpiCoords;
  ++neighborCoords.x();
  ++neighborCoords.y();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankX2Y2);
  // x2z1
  neighborCoords = d_mpiCoords;
  ++neighborCoords.x();
  --neighborCoords.z();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankX2Z1);
  // x2z2
  neighborCoords = d_mpiCoords;
  ++neighborCoords.x();
  ++neighborCoords.z();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankX2Z2);
  // y1z1
  neighborCoords = d_mpiCoords;
  --neighborCoords.y();
  --neighborCoords.z();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankY1Z1);
  // y1z2
  neighborCoords = d_mpiCoords;
  --neighborCoords.y();
  ++neighborCoords.z();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankY1Z2);
  // y2z1
  neighborCoords = d_mpiCoords;
  ++neighborCoords.y();
  --neighborCoords.z();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankY2Z1);
  // y2z2
  neighborCoords = d_mpiCoords;
  ++neighborCoords.y();
  ++neighborCoords.z();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankY2Z2);
  // x1y1z1
  neighborCoords = d_mpiCoords;
  --neighborCoords.x();
  --neighborCoords.y();
  --neighborCoords.z();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankX1Y1Z1);
  // x1y1z2
  neighborCoords = d_mpiCoords;
  --neighborCoords.x();
  --neighborCoords.y();
  ++neighborCoords.z();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankX1Y1Z2);
  // x1y2z1
  neighborCoords = d_mpiCoords;
  --neighborCoords.x();
  ++neighborCoords.y();
  --neighborCoords.z();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankX1Y2Z1);
  // x1y2z2
  neighborCoords = d_mpiCoords;
  --neighborCoords.x();
  ++neighborCoords.y();
  ++neighborCoords.z();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankX1Y2Z2);
  // x2y1z1
  neighborCoords = d_mpiCoords;
  ++neighborCoords.x();
  --neighborCoords.y();
  --neighborCoords.z();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankX2Y1Z1);
  // x2y1z2
  neighborCoords = d_mpiCoords;
  ++neighborCoords.x();
  --neighborCoords.y();
  ++neighborCoords.z();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankX2Y1Z2);
  // x2y2z1
  neighborCoords = d_mpiCoords;
  ++neighborCoords.x();
  ++neighborCoords.y();
  --neighborCoords.z();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankX2Y2Z1);
  // x2y2z2
  neighborCoords = d_mpiCoords;
  ++neighborCoords.x();
  ++neighborCoords.y();
  ++neighborCoords.z();
  MPI_Cart_rank(cartComm, neighborCoords.data(), &rankX2Y2Z2);

  // if found, communicate with neighboring blocks
  ParticlePArray particleX1, particleX2;
  ParticlePArray particleY1, particleY2;
  ParticlePArray particleZ1, particleZ2;
  ParticlePArray particleX1Y1, particleX1Y2, particleX1Z1, particleX1Z2;
  ParticlePArray particleX2Y1, particleX2Y2, particleX2Z1, particleX2Z2;
  ParticlePArray particleY1Z1, particleY1Z2, particleY2Z1, particleY2Z2;
  ParticlePArray particleX1Y1Z1, particleX1Y1Z2, particleX1Y2Z1, particleX1Y2Z2;
  ParticlePArray particleX2Y1Z1, particleX2Y1Z2, particleX2Y2Z1, particleX2Y2Z2;
  boost::mpi::request reqX1[2], reqX2[2];
  boost::mpi::request reqY1[2], reqY2[2];
  boost::mpi::request reqZ1[2], reqZ2[2];
  boost::mpi::request reqX1Y1[2], reqX1Y2[2], reqX1Z1[2], reqX1Z2[2];
  boost::mpi::request reqX2Y1[2], reqX2Y2[2], reqX2Z1[2], reqX2Z2[2];
  boost::mpi::request reqY1Z1[2], reqY1Z2[2], reqY2Z1[2], reqY2Z2[2];
  boost::mpi::request reqX1Y1Z1[2], reqX1Y1Z2[2], reqX1Y2Z1[2], reqX1Y2Z2[2];
  boost::mpi::request reqX2Y1Z1[2], reqX2Y1Z2[2], reqX2Y2Z1[2], reqX2Y2Z2[2];
  v1 = container.getMinCorner(); // redefine v1, v2 in terms of process
  v2 = container.getMaxCorner();
  */
  /*
  debugInf << "rank=" << mpiRank
     << ' ' << v1.x() << ' ' << v1.y() << ' ' << v1.z()
     << ' ' << v2.x() << ' ' << v2.y() << ' ' << v2.z()
     << std::endl;
  */
  /*
  REAL cellSize = gradation.getPtclMaxRadius() * 2;
  // 6 surfaces
  if (rankX1 >= 0) { // surface x1
    Box containerX1(v1.x(), v1.y(), v1.z(), v1.x() + cellSize, v2.y(), v2.z());
    findParticleInBox(containerX1, particleVec, particleX1);
    reqX1[0] = boostWorld.isend(rankX1, mpiTag, particleX1);
    reqX1[1] = boostWorld.irecv(rankX1, mpiTag, rParticleX1);
    out << "rankX1: mpiRank = " << mpiRank 
        << " iter = " << iteration << "\n"
        << " box = " << containerX1
        << " num particles sent =  " << particleX1.size()
        << " num particles recv =  " << rParticleX1.size() << "\n";
  }
  if (rankX2 >= 0) { // surface x2
    Box containerX2(v2.x() - cellSize, v1.y(), v1.z(), v2.x(), v2.y(), v2.z());
    findParticleInBox(containerX2, particleVec, particleX2);
    reqX2[0] = boostWorld.isend(rankX2, mpiTag, particleX2);
    reqX2[1] = boostWorld.irecv(rankX2, mpiTag, rParticleX2);
    out << "rankX2: mpiRank = " << mpiRank 
        << " iter = " << iteration << "\n"
        << " box = " << containerX2
        << " num particles sent =  " << particleX2.size()
        << " num particles recv =  " << rParticleX2.size() << "\n";
  }
  if (rankY1 >= 0) { // surface y1
    Box containerY1(v1.x(), v1.y(), v1.z(), v2.x(), v1.y() + cellSize, v2.z());
    findParticleInBox(containerY1, particleVec, particleY1);
    reqY1[0] = boostWorld.isend(rankY1, mpiTag, particleY1);
    reqY1[1] = boostWorld.irecv(rankY1, mpiTag, rParticleY1);
  }
  if (rankY2 >= 0) { // surface y2
    Box containerY2(v1.x(), v2.y() - cellSize, v1.z(), v2.x(), v2.y(), v2.z());
    findParticleInBox(containerY2, particleVec, particleY2);
    reqY2[0] = boostWorld.isend(rankY2, mpiTag, particleY2);
    reqY2[1] = boostWorld.irecv(rankY2, mpiTag, rParticleY2);
  }
  if (rankZ1 >= 0) { // surface z1
    Box containerZ1(v1.x(), v1.y(), v1.z(), v2.x(), v2.y(), v1.z() + cellSize);
    findParticleInBox(containerZ1, particleVec, particleZ1);
    reqZ1[0] = boostWorld.isend(rankZ1, mpiTag, particleZ1);
    reqZ1[1] = boostWorld.irecv(rankZ1, mpiTag, rParticleZ1);
  }
  if (rankZ2 >= 0) { // surface z2
    Box containerZ2(v1.x(), v1.y(), v2.z() - cellSize, v2.x(), v2.y(), v2.z());
    findParticleInBox(containerZ2, particleVec, particleZ2);
    reqZ2[0] = boostWorld.isend(rankZ2, mpiTag, particleZ2);
    reqZ2[1] = boostWorld.irecv(rankZ2, mpiTag, rParticleZ2);
  }
  // 12 edges
  if (rankX1Y1 >= 0) { // edge x1y1
    Box containerX1Y1(v1.x(), v1.y(), v1.z(), v1.x() + cellSize,
                      v1.y() + cellSize, v2.z());
    findParticleInBox(containerX1Y1, particleVec, particleX1Y1);
    reqX1Y1[0] = boostWorld.isend(rankX1Y1, mpiTag, particleX1Y1);
    reqX1Y1[1] = boostWorld.irecv(rankX1Y1, mpiTag, rParticleX1Y1);
  }
  if (rankX1Y2 >= 0) { // edge x1y2
    Box containerX1Y2(v1.x(), v2.y() - cellSize, v1.z(), v1.x() + cellSize,
                      v2.y(), v2.z());
    findParticleInBox(containerX1Y2, particleVec, particleX1Y2);
    reqX1Y2[0] = boostWorld.isend(rankX1Y2, mpiTag, particleX1Y2);
    reqX1Y2[1] = boostWorld.irecv(rankX1Y2, mpiTag, rParticleX1Y2);
  }
  if (rankX1Z1 >= 0) { // edge x1z1
    Box containerX1Z1(v1.x(), v1.y(), v1.z(), v1.x() + cellSize, v2.y(),
                      v1.z() + cellSize);
    findParticleInBox(containerX1Z1, particleVec, particleX1Z1);
    reqX1Z1[0] = boostWorld.isend(rankX1Z1, mpiTag, particleX1Z1);
    reqX1Z1[1] = boostWorld.irecv(rankX1Z1, mpiTag, rParticleX1Z1);
  }
  if (rankX1Z2 >= 0) { // edge x1z2
    Box containerX1Z2(v1.x(), v1.y(), v2.z() - cellSize, v1.x() + cellSize,
                      v2.y(), v2.z());
    findParticleInBox(containerX1Z2, particleVec, particleX1Z2);
    reqX1Z2[0] = boostWorld.isend(rankX1Z2, mpiTag, particleX1Z2);
    reqX1Z2[1] = boostWorld.irecv(rankX1Z2, mpiTag, rParticleX1Z2);
  }
  if (rankX2Y1 >= 0) { // edge x2y1
    Box containerX2Y1(v2.x() - cellSize, v1.y(), v1.z(), v2.x(),
                      v1.y() + cellSize, v2.z());
    findParticleInBox(containerX2Y1, particleVec, particleX2Y1);
    reqX2Y1[0] = boostWorld.isend(rankX2Y1, mpiTag, particleX2Y1);
    reqX2Y1[1] = boostWorld.irecv(rankX2Y1, mpiTag, rParticleX2Y1);
  }
  if (rankX2Y2 >= 0) { // edge x2y2
    Box containerX2Y2(v2.x() - cellSize, v2.y() - cellSize, v1.z(), v2.x(),
                      v2.y(), v2.z());
    findParticleInBox(containerX2Y2, particleVec, particleX2Y2);
    reqX2Y2[0] = boostWorld.isend(rankX2Y2, mpiTag, particleX2Y2);
    reqX2Y2[1] = boostWorld.irecv(rankX2Y2, mpiTag, rParticleX2Y2);
  }
  if (rankX2Z1 >= 0) { // edge x2z1
    Box containerX2Z1(v2.x() - cellSize, v1.y(), v1.z(), v2.x(), v2.y(),
                      v1.z() + cellSize);
    findParticleInBox(containerX2Z1, particleVec, particleX2Z1);
    reqX2Z1[0] = boostWorld.isend(rankX2Z1, mpiTag, particleX2Z1);
    reqX2Z1[1] = boostWorld.irecv(rankX2Z1, mpiTag, rParticleX2Z1);
  }
  if (rankX2Z2 >= 0) { // edge x2z2
    Box containerX2Z2(v2.x() - cellSize, v1.y(), v2.z() - cellSize, v2.x(),
                      v2.y(), v2.z());
    findParticleInBox(containerX2Z2, particleVec, particleX2Z2);
    reqX2Z2[0] = boostWorld.isend(rankX2Z2, mpiTag, particleX2Z2);
    reqX2Z2[1] = boostWorld.irecv(rankX2Z2, mpiTag, rParticleX2Z2);
  }
  if (rankY1Z1 >= 0) { // edge y1z1
    Box containerY1Z1(v1.x(), v1.y(), v1.z(), v2.x(), v1.y() + cellSize,
                      v1.z() + cellSize);
    findParticleInBox(containerY1Z1, particleVec, particleY1Z1);
    reqY1Z1[0] = boostWorld.isend(rankY1Z1, mpiTag, particleY1Z1);
    reqY1Z1[1] = boostWorld.irecv(rankY1Z1, mpiTag, rParticleY1Z1);
  }
  if (rankY1Z2 >= 0) { // edge y1z2
    Box containerY1Z2(v1.x(), v1.y(), v2.z() - cellSize, v2.x(),
                      v1.y() + cellSize, v2.z());
    findParticleInBox(containerY1Z2, particleVec, particleY1Z2);
    reqY1Z2[0] = boostWorld.isend(rankY1Z2, mpiTag, particleY1Z2);
    reqY1Z2[1] = boostWorld.irecv(rankY1Z2, mpiTag, rParticleY1Z2);
  }
  if (rankY2Z1 >= 0) { // edge y2z1
    Box containerY2Z1(v1.x(), v2.y() - cellSize, v1.z(), v2.x(), v2.y(),
                      v1.z() + cellSize);
    findParticleInBox(containerY2Z1, particleVec, particleY2Z1);
    reqY2Z1[0] = boostWorld.isend(rankY2Z1, mpiTag, particleY2Z1);
    reqY2Z1[1] = boostWorld.irecv(rankY2Z1, mpiTag, rParticleY2Z1);
  }
  if (rankY2Z2 >= 0) { // edge y2z2
    Box containerY2Z2(v1.x(), v2.y() - cellSize, v2.z() - cellSize, v2.x(),
                      v2.y(), v2.z());
    findParticleInBox(containerY2Z2, particleVec, particleY2Z2);
    reqY2Z2[0] = boostWorld.isend(rankY2Z2, mpiTag, particleY2Z2);
    reqY2Z2[1] = boostWorld.irecv(rankY2Z2, mpiTag, rParticleY2Z2);
  }
  // 8 vertices
  if (rankX1Y1Z1 >= 0) { // edge x1y1z1
    Box containerX1Y1Z1(v1.x(), v1.y(), v1.z(), v1.x() + cellSize,
                        v1.y() + cellSize, v1.z() + cellSize);
    findParticleInBox(containerX1Y1Z1, particleVec, particleX1Y1Z1);
    reqX1Y1Z1[0] = boostWorld.isend(rankX1Y1Z1, mpiTag, particleX1Y1Z1);
    reqX1Y1Z1[1] = boostWorld.irecv(rankX1Y1Z1, mpiTag, rParticleX1Y1Z1);
  }
  if (rankX1Y1Z2 >= 0) { // edge x1y1z2
    Box containerX1Y1Z2(v1.x(), v1.y(), v2.z() - cellSize, v1.x() + cellSize,
                        v1.y() + cellSize, v2.z());
    findParticleInBox(containerX1Y1Z2, particleVec, particleX1Y1Z2);
    reqX1Y1Z2[0] = boostWorld.isend(rankX1Y1Z2, mpiTag, particleX1Y1Z2);
    reqX1Y1Z2[1] = boostWorld.irecv(rankX1Y1Z2, mpiTag, rParticleX1Y1Z2);
  }
  if (rankX1Y2Z1 >= 0) { // edge x1y2z1
    Box containerX1Y2Z1(v1.x(), v2.y() - cellSize, v1.z(), v1.x() + cellSize,
                        v2.y(), v1.z() + cellSize);
    findParticleInBox(containerX1Y2Z1, particleVec, particleX1Y2Z1);
    reqX1Y2Z1[0] = boostWorld.isend(rankX1Y2Z1, mpiTag, particleX1Y2Z1);
    reqX1Y2Z1[1] = boostWorld.irecv(rankX1Y2Z1, mpiTag, rParticleX1Y2Z1);
  }
  if (rankX1Y2Z2 >= 0) { // edge x1y2z2
    Box containerX1Y2Z2(v1.x(), v2.y() - cellSize, v2.z() - cellSize,
                        v1.x() + cellSize, v2.y() + cellSize, v2.z());
    findParticleInBox(containerX1Y2Z2, particleVec, particleX1Y2Z2);
    reqX1Y2Z2[0] = boostWorld.isend(rankX1Y2Z2, mpiTag, particleX1Y2Z2);
    reqX1Y2Z2[1] = boostWorld.irecv(rankX1Y2Z2, mpiTag, rParticleX1Y2Z2);
  }
  if (rankX2Y1Z1 >= 0) { // edge x2y1z1
    Box containerX2Y1Z1(v2.x() - cellSize, v1.y(), v1.z(), v2.x(),
                        v1.y() + cellSize, v1.z() + cellSize);
    findParticleInBox(containerX2Y1Z1, particleVec, particleX2Y1Z1);
    reqX2Y1Z1[0] = boostWorld.isend(rankX2Y1Z1, mpiTag, particleX2Y1Z1);
    reqX2Y1Z1[1] = boostWorld.irecv(rankX2Y1Z1, mpiTag, rParticleX2Y1Z1);
  }
  if (rankX2Y1Z2 >= 0) { // edge x2y1z2
    Box containerX2Y1Z2(v2.x() - cellSize, v1.y(), v2.z() - cellSize, v2.x(),
                        v1.y() + cellSize, v2.z());
    findParticleInBox(containerX2Y1Z2, particleVec, particleX2Y1Z2);
    reqX2Y1Z2[0] = boostWorld.isend(rankX2Y1Z2, mpiTag, particleX2Y1Z2);
    reqX2Y1Z2[1] = boostWorld.irecv(rankX2Y1Z2, mpiTag, rParticleX2Y1Z2);
  }
  if (rankX2Y2Z1 >= 0) { // edge x2y2z1
    Box containerX2Y2Z1(v2.x() - cellSize, v2.y() - cellSize, v1.z(), v2.x(),
                        v2.y(), v1.z() + cellSize);
    findParticleInBox(containerX2Y2Z1, particleVec, particleX2Y2Z1);
    reqX2Y2Z1[0] = boostWorld.isend(rankX2Y2Z1, mpiTag, particleX2Y2Z1);
    reqX2Y2Z1[1] = boostWorld.irecv(rankX2Y2Z1, mpiTag, rParticleX2Y2Z1);
  }
  if (rankX2Y2Z2 >= 0) { // edge x2y2z2
    Box containerX2Y2Z2(v2.x() - cellSize, v2.y() - cellSize, v2.z() - cellSize,
                        v2.x(), v2.y(), v2.z());
    findParticleInBox(containerX2Y2Z2, particleVec, particleX2Y2Z2);
    reqX2Y2Z2[0] = boostWorld.isend(rankX2Y2Z2, mpiTag, particleX2Y2Z2);
    reqX2Y2Z2[1] = boostWorld.irecv(rankX2Y2Z2, mpiTag, rParticleX2Y2Z2);
  }


  // 6 surfaces
  if (rankX1 >= 0) {
    boost::mpi::wait_all(reqX1, reqX1 + 2);
    out << "rankX1 after wait: mpiRank = " << mpiRank 
        << " iter = " << iteration << "\n"
        << " num particles sent =  " << particleX1.size()
        << " num particles recv =  " << rParticleX1.size() << "\n";
  }
  if (rankX2 >= 0) {
    boost::mpi::wait_all(reqX2, reqX2 + 2);
    out << "rankX2 after wait: mpiRank = " << mpiRank 
        << " iter = " << iteration << "\n"
        << " num particles sent =  " << particleX2.size()
        << " num particles recv =  " << rParticleX2.size() << "\n";
  }
  if (rankY1 >= 0)
    boost::mpi::wait_all(reqY1, reqY1 + 2);
  if (rankY2 >= 0)
    boost::mpi::wait_all(reqY2, reqY2 + 2);
  if (rankZ1 >= 0)
    boost::mpi::wait_all(reqZ1, reqZ1 + 2);
  if (rankZ2 >= 0)
    boost::mpi::wait_all(reqZ2, reqZ2 + 2);
  // 12 edges
  if (rankX1Y1 >= 0)
    boost::mpi::wait_all(reqX1Y1, reqX1Y1 + 2);
  if (rankX1Y2 >= 0)
    boost::mpi::wait_all(reqX1Y2, reqX1Y2 + 2);
  if (rankX1Z1 >= 0)
    boost::mpi::wait_all(reqX1Z1, reqX1Z1 + 2);
  if (rankX1Z2 >= 0)
    boost::mpi::wait_all(reqX1Z2, reqX1Z2 + 2);
  if (rankX2Y1 >= 0)
    boost::mpi::wait_all(reqX2Y1, reqX2Y1 + 2);
  if (rankX2Y2 >= 0)
    boost::mpi::wait_all(reqX2Y2, reqX2Y2 + 2);
  if (rankX2Z1 >= 0)
    boost::mpi::wait_all(reqX2Z1, reqX2Z1 + 2);
  if (rankX2Z2 >= 0)
    boost::mpi::wait_all(reqX2Z2, reqX2Z2 + 2);
  if (rankY1Z1 >= 0)
    boost::mpi::wait_all(reqY1Z1, reqY1Z1 + 2);
  if (rankY1Z2 >= 0)
    boost::mpi::wait_all(reqY1Z2, reqY1Z2 + 2);
  if (rankY2Z1 >= 0)
    boost::mpi::wait_all(reqY2Z1, reqY2Z1 + 2);
  if (rankY2Z2 >= 0)
    boost::mpi::wait_all(reqY2Z2, reqY2Z2 + 2);
  // 8 vertices
  if (rankX1Y1Z1 >= 0)
    boost::mpi::wait_all(reqX1Y1Z1, reqX1Y1Z1 + 2);
  if (rankX1Y1Z2 >= 0)
    boost::mpi::wait_all(reqX1Y1Z2, reqX1Y1Z2 + 2);
  if (rankX1Y2Z1 >= 0)
    boost::mpi::wait_all(reqX1Y2Z1, reqX1Y2Z1 + 2);
  if (rankX1Y2Z2 >= 0)
    boost::mpi::wait_all(reqX1Y2Z2, reqX1Y2Z2 + 2);
  if (rankX2Y1Z1 >= 0)
    boost::mpi::wait_all(reqX2Y1Z1, reqX2Y1Z1 + 2);
  if (rankX2Y1Z2 >= 0)
    boost::mpi::wait_all(reqX2Y1Z2, reqX2Y1Z2 + 2);
  if (rankX2Y2Z1 >= 0)
    boost::mpi::wait_all(reqX2Y2Z1, reqX2Y2Z1 + 2);
  if (rankX2Y2Z2 >= 0)
    boost::mpi::wait_all(reqX2Y2Z2, reqX2Y2Z2 + 2);

  std::cout << out.str();

  // merge: particles inside container (at front) + particles from neighoring
  // blocks (at end)
  recvParticleVec.clear();
  // 6 surfaces
  if (rankX1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX1.begin(),
                           rParticleX1.end());
  if (rankX2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX2.begin(),
                           rParticleX2.end());
  if (rankY1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleY1.begin(),
                           rParticleY1.end());
  if (rankY2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleY2.begin(),
                           rParticleY2.end());
  if (rankZ1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleZ1.begin(),
                           rParticleZ1.end());
  if (rankZ2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleZ2.begin(),
                           rParticleZ2.end());
  // 12 edges
  if (rankX1Y1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1.begin(),
                           rParticleX1Y1.end());
  if (rankX1Y2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2.begin(),
                           rParticleX1Y2.end());
  if (rankX1Z1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX1Z1.begin(),
                           rParticleX1Z1.end());
  if (rankX1Z2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX1Z2.begin(),
                           rParticleX1Z2.end());
  if (rankX2Y1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1.begin(),
                           rParticleX2Y1.end());
  if (rankX2Y2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2.begin(),
                           rParticleX2Y2.end());
  if (rankX2Z1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX2Z1.begin(),
                           rParticleX2Z1.end());
  if (rankX2Z2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX2Z2.begin(),
                           rParticleX2Z2.end());
  if (rankY1Z1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleY1Z1.begin(),
                           rParticleY1Z1.end());
  if (rankY1Z2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleY1Z2.begin(),
                           rParticleY1Z2.end());
  if (rankY2Z1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleY2Z1.begin(),
                           rParticleY2Z1.end());
  if (rankY2Z2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleY2Z2.begin(),
                           rParticleY2Z2.end());
  // 8 vertices
  if (rankX1Y1Z1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1Z1.begin(),
                           rParticleX1Y1Z1.end());
  if (rankX1Y1Z2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1Z2.begin(),
                           rParticleX1Y1Z2.end());
  if (rankX1Y2Z1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2Z1.begin(),
                           rParticleX1Y2Z1.end());
  if (rankX1Y2Z2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2Z2.begin(),
                           rParticleX1Y2Z2.end());
  if (rankX2Y1Z1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1Z1.begin(),
                           rParticleX2Y1Z1.end());
  if (rankX2Y1Z2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1Z2.begin(),
                           rParticleX2Y1Z2.end());
  if (rankX2Y2Z1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2Z1.begin(),
                           rParticleX2Y2Z1.end());
  if (rankX2Y2Z2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2Z2.begin(),
                           rParticleX2Y2Z2.end());

  mergeParticleVec.clear();
  mergeParticleVec =
    particleVec; // duplicate pointers, pointing to the same memory
  mergeParticleVec.insert(mergeParticleVec.end(), recvParticleVec.begin(),
                          recvParticleVec.end());
  */
  /*
    ParticlePArray testParticleVec;
    testParticleVec.insert(testParticleVec.end(), rParticleX1.begin(),
    rParticleX1.end());
    testParticleVec.insert(testParticleVec.end(), rParticleX2.begin(),
    rParticleX2.end());
    testParticleVec.insert(testParticleVec.end(), rParticleY1.begin(),
    rParticleY1.end());
    testParticleVec.insert(testParticleVec.end(), rParticleY2.begin(),
    rParticleY2.end());
    testParticleVec.insert(testParticleVec.end(), rParticleZ1.begin(),
    rParticleZ1.end());
    testParticleVec.insert(testParticleVec.end(), rParticleZ2.begin(),
    rParticleZ2.end());
    debugInf << "iter=" << std::setw(4) << iteration << " rank=" << std::setw(4)
    << mpiRank
    << " ptclNum=" << std::setw(4) << particleVec.size()
    << " surface="
    << std::setw(4) << particleX1.size()  << std::setw(4) << particleX2.size()
    << std::setw(4) << particleY1.size()  << std::setw(4) << particleY2.size()
    << std::setw(4) << particleZ1.size()  << std::setw(4) << particleZ2.size()
    << " recv="
    << std::setw(4) << rParticleX1.size() << std::setw(4) << rParticleX2.size()
    << std::setw(4) << rParticleY1.size() << std::setw(4) << rParticleY2.size()
    << std::setw(4) << rParticleZ1.size() << std::setw(4) << rParticleZ2.size()
    << " rNum="
    << std::setw(4) << recvParticleVec.size() << ": ";

    for (ParticlePArray::const_iterator it = testParticleVec.begin(); it !=
    testParticleVec.end();++it)
    debugInf << (*it)->getId() << ' ';
    debugInf << std::endl;
    testParticleVec.clear();
  */
/*
}
*/

void
Assembly::calcTimeStep()
{
  calcVibraTimeStep();
  calcImpactTimeStep();
  calcContactNum();

  REAL CFL = 0.5;
  std::valarray<REAL> dt(3);
  dt[0] = util::getParam<REAL>("timeStep");
  dt[1] = CFL * vibraTimeStep;
  dt[2] = CFL * impactTimeStep;

  timeStep = dt.min();
  //if (mpiRank == 0) {
    //std::cout << "Timestep = " << timeStep
    //          << " Vibration timestep = " << vibraTimeStep
    //          << " Impact timestep = " << impactTimeStep << "\n";
  //}
}

void
Assembly::calcVibraTimeStep()
{
  REAL pTimeStep = 1 / EPS;
  if (contactVec.size() == 0) {
    pTimeStep = 1 / EPS;
  } else {
    auto it = contactVec.cbegin();
    pTimeStep = it->getVibraTimeStep();
    for (++it; it != contactVec.cend(); ++it) {
      REAL val = it->getVibraTimeStep();
      pTimeStep = val < pTimeStep ? val : pTimeStep;
    }
  }

  MPI_Allreduce(&pTimeStep, &vibraTimeStep, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);
}

void
Assembly::calcImpactTimeStep()
{
  REAL pTimeStep = 1 / EPS;
  if (contactVec.size() == 0)
    pTimeStep = 1 / EPS;
  else {
    auto it = contactVec.cbegin();
    pTimeStep = it->getImpactTimeStep();
    for (++it; it != contactVec.cend(); ++it) {
      REAL val = it->getImpactTimeStep();
      pTimeStep = val < pTimeStep ? val : pTimeStep;
    }
  }

  MPI_Allreduce(&pTimeStep, &impactTimeStep, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);
}

void
Assembly::calcContactNum()
{
  std::size_t pContactNum = contactVec.size();
  MPI_Reduce(&pContactNum, &allContactNum, 1, MPI_INT, MPI_SUM, 0, mpiWorld);
}

void
Assembly::findContact()
{ // various implementations
  int ompThreads = util::getParam<int>("ompThreads");

  if (ompThreads == 1) { // non-openmp single-thread version, time complexity
                         // bigO(n x n), n is the number of particles.
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    //std::cout << "\t FindContact: MPI Cart rank = " << mpiRank
    //          << " World rank = " << world_rank << "\n";
    findContactSingleThread();
  } else if (ompThreads > 1) { // openmp implementation: various loop scheduling
                               // - (static), (static,1), (dynamic), (dynamic,1)
    findContactMultiThread(ompThreads);

  } // end of openmp implementation
}

void
Assembly::findContactSingleThread()
{

  contactVec.clear();

#ifdef TIME_PROFILE
  Timer::time_point startOuter, endOuter;
  Timer::time_point startInner, endInner, durationInner = 0;
  startOuter = Timer::now();
#endif

  auto num1 = particleVec.size();      // particles inside container
  auto num2 = mergeParticleVec.size(); // particles inside container
                                       // (at front) + particles from
                                       // neighboring blocks (at end)
  // NOT (num1 - 1), in parallel situation where one particle
  // could contact received particles!
  for (auto i = 0; i < num1; ++i) {

    auto particle = particleVec[i];
    Vec u = particle->currentPos();
    auto particleType = particle->getType();
    auto particleRad = particle->getA();

    for (auto j = i + 1; j < num2; ++j) {

      auto mergeParticle = mergeParticleVec[j];
      Vec v = mergeParticle->currentPos();
      auto mergeParticleType = mergeParticle->getType();
      auto mergeParticleRad = mergeParticle->getA();

      if ((particle->getId() == 2 && mergeParticle->getId() == 94) ||
          (particle->getId() == 94 && mergeParticle->getId() == 2)) {
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        //std::cout << " MPI_Rank = " << world_rank 
        //          << " : particle = (" << particle->getId()
        //          << ", " << particleType
        //          << "), mergeParticle = (" << mergeParticle->getId()
        //          << ", " << mergeParticleType << ")\n";
      } 

      if ((vfabs(v - u) < particleRad + mergeParticleRad) &&
          // not both are fixed particles
          (particleType != 1 || mergeParticleType != 1) &&
          // not both are free boundary particles
          (particleType != 5 || mergeParticleType != 5) &&
          // not both are ghost particles
          (particleType != 10 || mergeParticleType != 10)) {

        Contact tmpContact(particle.get(), mergeParticle.get());

#ifdef TIME_PROFILE
        startInner = Timer::now();
#endif
        if (tmpContact.isOverlapped()) {
          contactVec.push_back(tmpContact); // containers use value
                                            // semantics, so a "copy" is
                                            // pushed back.
        }
#ifdef TIME_PROFILE
        endInner = Timer::now();
        durationInner += (endInner - startInner);
#endif
      }
    }
  }
#ifdef TIME_PROFILE
  endOuter = Timer::now();
  debugInf << std::setw(OWID) << "findContact=" << std::setw(OWID)
           << std::chrono::duration_cast<Seconds>(endOuter - startOuter).count()
           << std::setw(OWID) << "isOverlapped=" << std::setw(OWID)
           << std::chrono::duration_cast<Seconds>(durationInner).count();
#endif
}

void
Assembly::findContactMultiThread(int ompThreads)
{

  contactVec.clear();

#ifdef TIME_PROFILE
  Timer::time_point startOuter, endOuter;
  startOuter = Timer::now();
#endif

  std::size_t i, j, particleType, mergeParticleType;
  Vec u, v;
  auto num1 = particleVec.size();      // particles inside container
  auto num2 = mergeParticleVec.size(); // particles inside container
                                       // (at front) + particles from
                                       // neighboring blocks (at end)
#pragma omp parallel for num_threads(ompThreads) private(                      \
  i, j, u, v, particleType, mergeParticleType)                                 \
    shared(num1, num2) schedule(dynamic)
  for (i = 0; i < num1; ++i) {
    u = particleVec[i]->currentPos();
    particleType = particleVec[i]->getType();
    for (j = i + 1; j < num2; ++j) {
      Vec v = mergeParticleVec[j]->currentPos();
      mergeParticleType = mergeParticleVec[j]->getType();
      if ((vfabs(v - u) <
           particleVec[i]->getA() + mergeParticleVec[j]->getA()) &&
          // not both are fixed particles
          (particleType != 1 || mergeParticleType != 1) &&
          // not both are free boundary particles
          (particleType != 5 || mergeParticleType != 5) &&
          // not both are ghost particles
          (particleType != 10 || mergeParticleType != 10)) {

        Contact tmpContact(particleVec[i].get(), mergeParticleVec[j].get());

        if (tmpContact.isOverlapped()) {
#pragma omp critical
          contactVec.push_back(tmpContact); // containers use value
                                            // semantics, so a "copy" is
                                            // pushed back.
        }
      }
    }
  }
#ifdef TIME_PROFILE
  endOuter = Timer::now();
  debugInf
    << std::setw(OWID) << "findContact=" << std::setw(OWID)
    << std::chrono::duration_cast<Seconds>(endOuter - startOuter).count();
#endif
}

void
Assembly::findBdryContact()
{
  for (auto& boundary : boundaryVec) {
    boundary->findBdryContact(particleVec);
  }
}

void
Assembly::clearContactForce()
{
  for (auto& particle : particleVec) {
    particle->clearContactForce();
  }
}

void
Assembly::internalForce()
{
  /*
    std::ostringstream msg;
    msg << "iteration = " << iteration << " proc = " << mpiRank;
    msg << " : in internalForce : " << std::endl;
    cerr << msg.str();
  */

  REAL pAvg[3], sumAvg[3];
  for (auto i = 0u; i < 3; ++i) {
    pAvg[i] = 0;
    sumAvg[i] = 0;
  }

  // checkin previous tangential force and displacment
  for (auto& contact : contactVec) {
    contact.checkinPrevTgt(contactTgtVec);
  }

#ifdef TIME_PROFILE
  auto start = Timer::now();
#endif

  // contactTgtVec must be cleared before filling in new values.
  contactTgtVec.clear();

  for (auto& contact : contactVec) {
    // cannot be parallelized as it may change a
    // particle's force simultaneously.
    contact.contactForce();

    // checkout current tangential force and displacment
    contact.checkoutTgt(contactTgtVec);

    pAvg[0] += contact.getNormalForce();
    pAvg[1] += contact.getTgtForce();
    pAvg[2] += contact.getPenetration();
  }

  if (contactVec.size() != 0) {
    for (double& pVal : pAvg) {
      pVal /= contactVec.size();
    }
  }

#ifdef TIME_PROFILE
  auto end = Timer::now();
  debugInf << std::setw(OWID) << "internalForce=" << std::setw(OWID)
           << std::chrono::duration_cast<Seconds>(end - start).count()
           << std::endl;
#endif

  MPI_Reduce(pAvg, sumAvg, 3, MPI_DOUBLE, MPI_SUM, 0, mpiWorld);
  avgNormal = sumAvg[0] / mpiSize;
  avgShear = sumAvg[1] / mpiSize;
  avgPenetr = sumAvg[2] / mpiSize;
}

void
Assembly::boundaryForce()
{
  for (auto& boundary : boundaryVec) {
    boundary->boundaryForce(boundaryTgtMap);
  }
}

void
Assembly::updateParticle()
{
  for (auto& particle : particleVec)
    particle->update();
}

void
Assembly::updateGrid()
{
  updateGridMinX();
  updateGridMaxX();
  updateGridMinY();
  updateGridMaxY();
  updateGridMinZ();
  updateGridMaxZ();
}

void
Assembly::updateGridMinX()
{
  REAL pMinX = getPtclMinX(particleVec);
  REAL minX = 0;
  MPI_Allreduce(&pMinX, &minX, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);

  setGrid(Box(minX - gradation.getPtclMaxRadius(), grid.getMinCorner().y(),
              grid.getMinCorner().z(), grid.getMaxCorner().x(),
              grid.getMaxCorner().y(), grid.getMaxCorner().z()));
}

void
Assembly::updateGridMaxX()
{
  REAL pMaxX = getPtclMaxX(particleVec);
  REAL maxX = 0;
  MPI_Allreduce(&pMaxX, &maxX, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);

  setGrid(Box(grid.getMinCorner().x(), grid.getMinCorner().y(),
              grid.getMinCorner().z(), maxX + gradation.getPtclMaxRadius(),
              grid.getMaxCorner().y(), grid.getMaxCorner().z()));
}

void
Assembly::updateGridMinY()
{
  REAL pMinY = getPtclMinY(particleVec);
  REAL minY = 0;
  MPI_Allreduce(&pMinY, &minY, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);

  setGrid(Box(grid.getMinCorner().x(), minY - gradation.getPtclMaxRadius(),
              grid.getMinCorner().z(), grid.getMaxCorner().x(),
              grid.getMaxCorner().y(), grid.getMaxCorner().z()));
}

void
Assembly::updateGridMaxY()
{
  REAL pMaxY = getPtclMaxY(particleVec);
  REAL maxY = 0;
  MPI_Allreduce(&pMaxY, &maxY, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);

  setGrid(Box(grid.getMinCorner().x(), grid.getMinCorner().y(),
              grid.getMinCorner().z(), grid.getMaxCorner().x(),
              maxY + gradation.getPtclMaxRadius(), grid.getMaxCorner().z()));
}

void
Assembly::updateGridMinZ()
{
  REAL pMinZ = getPtclMinZ(particleVec);
  REAL minZ = 0;
  MPI_Allreduce(&pMinZ, &minZ, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);

  setGrid(Box(grid.getMinCorner().x(), grid.getMinCorner().y(),
              minZ - gradation.getPtclMaxRadius(), grid.getMaxCorner().x(),
              grid.getMaxCorner().y(), grid.getMaxCorner().z()));
}

void
Assembly::updateGridMaxZ()
{
  // update compute grids adaptively due to particle motion
  REAL pMaxZ = getPtclMaxZ(particleVec);
  REAL maxZ = 0;
  MPI_Allreduce(&pMaxZ, &maxZ, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);

  // no need to broadcast grid as it is updated in each process
  setGrid(Box(grid.getMinCorner().x(), grid.getMinCorner().y(),
              grid.getMinCorner().z(), grid.getMaxCorner().x(),
              grid.getMaxCorner().y(), maxZ + gradation.getPtclMaxRadius()));
}

void
Assembly::findParticleInBox(const Box& container,
                            const ParticlePArray& allParticles,
                            ParticlePArray& insideParticles)
{
  insideParticles.clear();
  for (const auto& particle : allParticles) {
    // it is critical to use EPS
    if (container.inside(particle->currentPos(), EPS)) {
      insideParticles.push_back(particle);
    }
  }
}

/*
void
Assembly::findParticleInBox(const Box& container,
                            const ParticlePArray& allParticles,
                            ParticlePArray& insideParticles)
{
  Vec  v1 = container.getMinCorner();
  Vec  v2 = container.getMaxCorner();
  REAL x1 = v1.x();
  REAL y1 = v1.y();
  REAL z1 = v1.z();
  REAL x2 = v2.x();
  REAL y2 = v2.y();
  REAL z2 = v2.z();
  for (const auto& particle : allParticles) {
    Vec center = particle->currentPos();
    // it is critical to use EPS
    if (center.x() - x1 >= -EPS && center.x() - x2 < -EPS &&
        center.y() - y1 >= -EPS && center.y() - y2 < -EPS &&
        center.z() - z1 >= -EPS && center.z() - z2 < -EPS) { 
      insideParticles.push_back(particle);
    }
  }
}
*/

void
Assembly::plotBoundary() const
{
  d_writer->writeDomain(&allContainer);
}

void
Assembly::printBoundary() const
{
  std::ofstream ofs(d_writer->getBoundaryFileName());
  if (!ofs) {
    debugInf << "stream error: printBoundary" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);

  Vec v1 = allContainer.getMinCorner();
  Vec v2 = allContainer.getMaxCorner();
  REAL x1 = v1.x();
  REAL y1 = v1.y();
  REAL z1 = v1.z();
  REAL x2 = v2.x();
  REAL y2 = v2.y();
  REAL z2 = v2.z();

  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z1
      << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z2
      << std::endl
      << std::endl
      << std::setw(OWID) << boundaryVec.size() << std::endl;

  for (const auto& it : boundaryVec) {
    it->print(ofs);
  }

  ofs.close();
}

void
Assembly::printBdryContact() const
{
  std::ofstream ofs(d_writer->getBdryContactFileName());
  if (!ofs) {
    debugInf << "stream error: printBdryContact" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  for (const auto& it : mergeBoundaryVec) {
    it->printContactInfo(ofs);
  }

  ofs.close();
}

void
Assembly::plotGrid() const
{
  d_writer->setMPIComm(cartComm);
  d_writer->setMPIProc(d_mpiProcs.x(), d_mpiProcs.y(), d_mpiProcs.z());
  d_writer->writeGrid(&grid);
}

void
Assembly::plotParticle() const
{
  d_writer->writeParticles(&allParticleVec);
  d_writer->writeSieves(&gradation);
}

void
Assembly::plotParticle(ParticlePArray& particles) const
{
  d_writer->writeParticles(&particles);
}

void
Assembly::printParticle(const std::string& fileName) const
{
  OutputTecplot writer(".", 0);
  writer.setParticleFileName(fileName);
  writer.writeParticles(&allParticleVec);
}

void
Assembly::printParticle(const std::string& fileName, ParticlePArray& particles) const
{
  OutputTecplot writer(".", 0);
  writer.setParticleFileName(fileName);
  writer.writeParticles(&particles);
}

void
Assembly::openDepositProg(std::ofstream& ofs, const std::string& str)
{
  ofs.open(str);
  if (!ofs) {
    debugInf << "stream error: openDepositProg" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  ofs << std::setw(OWID) << "iteration" << std::setw(OWID) << "normal_x1"
      << std::setw(OWID) << "normal_x2" << std::setw(OWID) << "normal_y1"
      << std::setw(OWID) << "normal_y2" << std::setw(OWID) << "normal_z1"
      << std::setw(OWID) << "normal_z2"

      << std::setw(OWID) << "contact_x1" << std::setw(OWID) << "contact_x2"
      << std::setw(OWID) << "contact_y1" << std::setw(OWID) << "contact_y2"
      << std::setw(OWID) << "contact_z1" << std::setw(OWID) << "contact_z2"
      << std::setw(OWID) << "contact_inside"

      << std::setw(OWID) << "penetr_x1" << std::setw(OWID) << "penetr_x2"
      << std::setw(OWID) << "penetr_y1" << std::setw(OWID) << "penetr_y2"
      << std::setw(OWID) << "penetr_z1" << std::setw(OWID) << "penetr_z2"

      << std::setw(OWID) << "avgNormal" << std::setw(OWID) << "avgShear"
      << std::setw(OWID) << "avgPenetr"

      << std::setw(OWID) << "transEnergy" << std::setw(OWID) << "rotatEnergy"
      << std::setw(OWID) << "kinetEnergy" << std::setw(OWID) << "graviEnergy"
      << std::setw(OWID) << "mechaEnergy"

      << std::setw(OWID) << "vibra_est_dt" << std::setw(OWID) << "impact_est_dt"
      << std::setw(OWID) << "actual_dt" << std::setw(OWID) << "accruedTime"

      << std::endl;
}

void
Assembly::printDepositProg(std::ofstream& ofs)
{
  REAL var[6];

  // normalForce
  for (double& i : var)
    i = 0;
  for (BoundaryPArray::const_iterator it = mergeBoundaryVec.begin();
       it != mergeBoundaryVec.end(); ++it) {
    std::size_t id = (*it)->getId();
    Vec normal = (*it)->getNormalForce();
    //std::cout << "normal force = " << std::setprecision(16) << normal << "\n";

    switch (id) {
      case 1:
        var[0] = fabs(normal.x());
        break;
      case 2:
        var[1] = normal.x();
        break;
      case 3:
        var[2] = fabs(normal.y());
        break;
      case 4:
        var[3] = normal.y();
        break;
      case 5:
        var[4] = fabs(normal.z());
        break;
      case 6:
        var[5] = normal.z();
        break;
    }
  }
  ofs << std::setw(OWID) << iteration;
  for (double i : var)
    ofs << std::setw(OWID) << i;

  // contactNum
  for (double& i : var)
    i = 0;
  for (BoundaryPArray::const_iterator it = mergeBoundaryVec.begin();
       it != mergeBoundaryVec.end(); ++it) {
    std::size_t id = (*it)->getId();
    var[id - 1] = (*it)->getContactNum();
  }
  for (double i : var)
    ofs << std::setw(OWID) << static_cast<std::size_t>(i);
  ofs << std::setw(OWID) << allContactNum;

  // avgPenetr
  for (double& i : var)
    i = 0;
  for (BoundaryPArray::const_iterator it = mergeBoundaryVec.begin();
       it != mergeBoundaryVec.end(); ++it) {
    std::size_t id = (*it)->getId();
    var[id - 1] = (*it)->getAvgPenetr();
  }
  for (double i : var)
    ofs << std::setw(OWID) << i;

  // average data
  ofs << std::setw(OWID) << avgNormal << std::setw(OWID) << avgShear
      << std::setw(OWID) << avgPenetr;

  // energy
  ofs << std::setw(OWID) << transEnergy << std::setw(OWID) << rotatEnergy
      << std::setw(OWID) << kinetEnergy << std::setw(OWID) << graviEnergy
      << std::setw(OWID) << mechaEnergy;

  // time
  ofs << std::setw(OWID) << vibraTimeStep << std::setw(OWID) << impactTimeStep
      << std::setw(OWID) << timeStep << std::setw(OWID) << timeAccrued;

  ofs << std::endl;
}

bool
Assembly::tractionErrorTol(REAL sigma, std::string type, REAL sigmaX,
                           REAL sigmaY)
{
  // sigma implies sigmaZ
  REAL tol = util::getParam<REAL>("tractionErrorTol");

  std::map<std::string, REAL> normalForce;
  REAL x1, x2, y1, y2, z1, z2;
  // do not use mergeBoundaryVec because each process calls this function.
  for (const auto& it : boundaryVec) {
    std::size_t id = it->getId();
    Vec normal = it->getNormalForce();
    Vec point = it->getPoint();
    switch (id) {
      case 1:
        normalForce["x1"] = fabs(normal.x());
        x1 = point.x();
        break;
      case 2:
        normalForce["x2"] = normal.x();
        x2 = point.x();
        break;
      case 3:
        normalForce["y1"] = fabs(normal.y());
        y1 = point.y();
        break;
      case 4:
        normalForce["y2"] = normal.y();
        y2 = point.y();
        break;
      case 5:
        normalForce["z1"] = fabs(normal.z());
        z1 = point.z();
        break;
      case 6:
        normalForce["z2"] = normal.z();
        z2 = point.z();
        break;
    }
  }
  REAL areaX = (y2 - y1) * (z2 - z1);
  REAL areaY = (z2 - z1) * (x2 - x1);
  REAL areaZ = (x2 - x1) * (y2 - y1);

  if (type.compare("isotropic") == 0)
    return (fabs(normalForce["x1"] / areaX - sigma) / sigma <= tol &&
            fabs(normalForce["x2"] / areaX - sigma) / sigma <= tol &&
            fabs(normalForce["y1"] / areaY - sigma) / sigma <= tol &&
            fabs(normalForce["y2"] / areaY - sigma) / sigma <= tol &&
            fabs(normalForce["z1"] / areaZ - sigma) / sigma <= tol &&
            fabs(normalForce["z2"] / areaZ - sigma) / sigma <= tol);

  else if (type.compare("odometer") == 0)
    return (fabs(normalForce["z1"] / areaZ - sigma) / sigma <= tol &&
            fabs(normalForce["z2"] / areaZ - sigma) / sigma <= tol);

  else if (type.compare("triaxial") == 0)
    return true; // always near equilibrium

  else if (type.compare("trueTriaxial") == 0)
    return (fabs(normalForce["x1"] / areaX - sigmaX) / sigmaX <= tol &&
            fabs(normalForce["x2"] / areaX - sigmaX) / sigmaX <= tol &&
            fabs(normalForce["y1"] / areaY - sigmaY) / sigmaY <= tol &&
            fabs(normalForce["y2"] / areaY - sigmaY) / sigmaY <= tol &&
            fabs(normalForce["z1"] / areaZ - sigma) / sigma <= tol &&
            fabs(normalForce["z2"] / areaZ - sigma) / sigma <= tol);

  return false;
}

// particleLayers:
// 0 - one free particle
// 1 - a horizontal layer of free particles
// 2 - multiple layers of free particles
void
Assembly::generateParticle(std::size_t particleLayers,
                           const std::string& genParticle)
{
  REAL young = util::getParam<REAL>("young");
  REAL poisson = util::getParam<REAL>("poisson");

  REAL x, y, z;
  std::size_t particleNum = 0;
  REAL diameter = gradation.getPtclMaxRadius() * 2.0;

  REAL offset = 0;
  REAL edge = diameter;
  if (gradation.getSize().size() == 1 && gradation.getPtclRatioBA() == 1.0 &&
      gradation.getPtclRatioCA() == 1.0) {
    edge = diameter * 2.0;
    offset = diameter * 0.25;
  }

  REAL x1 = allContainer.getMinCorner().x() + edge;
  REAL y1 = allContainer.getMinCorner().y() + edge;
  REAL z1 = allContainer.getMinCorner().z() + diameter;
  REAL x2 = allContainer.getMaxCorner().x() - edge;
  REAL y2 = allContainer.getMaxCorner().y() - edge;
  // REAL z2 = allContainer.getMaxCorner().z() - diameter;
  REAL z2 = util::getParam<REAL>("floatMaxZ") - diameter;
  REAL x0 = allContainer.getCenter().x();
  REAL y0 = allContainer.getCenter().y();
  REAL z0 = allContainer.getCenter().z();

  if (particleLayers == 0) { // just one free particle
    ParticleP newptcl = std::make_shared<Particle>(
      particleNum + 1, 0, Vec(x0, y0, z0), gradation, young, poisson);
    allParticleVec.push_back(newptcl);
    particleNum++;
  } else if (particleLayers == 1) { // a horizontal layer of free particles
    for (x = x1; x - x2 < EPS; x += diameter)
      for (y = y1; y - y2 < EPS; y += diameter) {
        ParticleP newptcl = std::make_shared<Particle>(
          particleNum + 1, 0, Vec(x, y, z0), gradation, young, poisson);
        allParticleVec.push_back(newptcl);
        particleNum++;
      }
  } else if (particleLayers == 2) { // multiple layers of free particles
    for (z = z1; z - z2 < EPS; z += diameter) {
      for (x = x1 + offset; x - x2 < EPS; x += diameter)
        for (y = y1 + offset; y - y2 < EPS; y += diameter) {
          ParticleP newptcl = std::make_shared<Particle>(
            particleNum + 1, 0, Vec(x, y, z), gradation, young, poisson);
          allParticleVec.push_back(newptcl);
          particleNum++;
        }
      offset *= -1;
    }
  }

  printParticle(genParticle);
}

void
Assembly::trim(bool toRebuild, const std::string& inputParticle,
               const std::string& trmParticle)
{
  if (toRebuild)
    readParticles(inputParticle);
  trimHistoryNum = allParticleVec.size();

  Vec v1 = allContainer.getMinCorner();
  Vec v2 = allContainer.getMaxCorner();
  REAL x1 = v1.x();
  REAL y1 = v1.y();
  REAL z1 = v1.z();
  REAL x2 = v2.x();
  REAL y2 = v2.y();
  REAL z2 = v2.z();
  REAL maxR = gradation.getPtclMaxRadius();

  // BB: Feb 2, 2017:
  // Not an efficient operation
  // Better approach may be to use a list if random access of vector
  // members is not needed
  allParticleVec.erase(
    std::remove_if(allParticleVec.begin(), allParticleVec.end(),
                   [&x1, &y1, &z1, &x2, &y2, &z2, &maxR](ParticleP particle) {
                     Vec center = particle->currentPos();
                     if (center.x() < x1 || center.x() > x2 ||
                         center.y() < y1 || center.y() > y2 ||
                         center.z() < z1 || center.z() + maxR > z2) {
                       return true;
                     }
                     return false;
                   }),
    allParticleVec.end());

  /*
  ParticlePArray::iterator itr;
  Vec center;
  for (auto itr = allParticleVec.begin(); itr != allParticleVec.end(); ) {
    center=(*itr)->currentPos();
    if(center.x() < x1 || center.x() > x2 ||
   center.y() < y1 || center.y() > y2 ||
   center.z() < z1 || center.z() + maxR > z2)
  {
    delete (*itr); // release memory
    itr = allParticleVec.erase(itr);
  }
    else
  ++itr;
  }
  */

  printParticle(trmParticle);
}

void
Assembly::findPeriParticleInBox(const Box& container,
                                const PeriParticlePArray& inputParticle,
                                PeriParticlePArray& foundParticle)
{
  Vec v1 = container.getMinCorner();
  Vec v2 = container.getMaxCorner();
  REAL x1 = v1.x();
  REAL y1 = v1.y();
  REAL z1 = v1.z();
  REAL x2 = v2.x();
  REAL y2 = v2.y();
  REAL z2 = v2.z();
  for (const auto& pt : inputParticle) {
    Vec center = pt->currentPosition();
    // it is critical to use EPS
    if (center.x() - x1 >= -EPS && center.x() - x2 < -EPS &&
        center.y() - y1 >= -EPS && center.y() - y2 < -EPS &&
        center.z() - z1 >= -EPS && center.z() - z2 < -EPS)
      foundParticle.push_back(pt);
  }
}

void
Assembly::removeParticleOutBox()
{
  Vec v1 = container.getMinCorner();
  Vec v2 = container.getMaxCorner();
  REAL x1 = v1.x();
  REAL y1 = v1.y();
  REAL z1 = v1.z();
  REAL x2 = v2.x();
  REAL y2 = v2.y();
  REAL z2 = v2.z();

  // BB: Feb 2, 2017:
  // Not an efficient operation
  // Better approach may be to use a list if random access of vector
  // members is not needed
  //std::cout << "MPI Cart rank = " << mpiRank << "\n";
  REAL epsilon = EPS;
  particleVec.erase(
    std::remove_if(
      particleVec.begin(), particleVec.end(),
      [&x1, &y1, &z1, &x2, &y2, &z2, &epsilon](ParticleP particle) {
        Vec center = particle->currentPos();
        if (!(center.x() - x1 >= -epsilon && center.x() - x2 < -epsilon &&
              center.y() - y1 >= -epsilon && center.y() - y2 < -epsilon &&
              center.z() - z1 >= -epsilon && center.z() - z2 < -epsilon)) {
          /*
          if (particle->getId() == 2 || particle->getId() == 94) {
            //std::cout << "**WARNING** Removing particle " << particle->getId()
                      << " from container with \n "
                      << " x : " << center.x() << " not in [" << x1 << ","  << x2 << "]\n"
                      << " y : " << center.y() << " not in [" << y1 << ","  << y2 << "]\n"
                      << " z : " << center.z() << " not in [" << z1 << ","  << z2 << "]\n";
          }
          */
          return true;
        }
        return false;
      }),
    particleVec.end());

  /*
  ParticlePArray::iterator itr;
  Vec center;
  //std::size_t flag = 0;

  for (itr = particleVec.begin(); itr != particleVec.end(); ) {
    center=(*itr)->currentPos();
    // it is critical to use EPS
    if ( !(center.x() - x1 >= -EPS && center.x() - x2 < -EPS &&
       center.y() - y1 >= -EPS && center.y() - y2 < -EPS &&
       center.z() - z1 >= -EPS && center.z() - z2 < -EPS) )
  {
    // debugInf << "iter=" << std::setw(8) << iteration << " rank=" <<
  std::setw(2) << mpiRank
    // << " removed=" << std::setw(3) << (*itr)->getId();
    // flag = 1;
    delete (*itr); // release memory
    itr = particleVec.erase(itr);
  }
    else
  ++itr;
  }
  */

  /*
    if (flag == 1) {
    debugInf << " now " << particleVec.size() << ": ";
    for (ParticlePArray::const_iterator it = particleVec.begin(); it !=
    particleVec.end(); ++it)
    debugInf << std::setw(3) << (*it)->getId();
    debugInf << std::endl;
    }
  */
}

void
Assembly::removePeriParticleOutBox()
{
  Vec v1 = container.getMinCorner();
  Vec v2 = container.getMaxCorner();
  REAL x1 = v1.x();
  REAL y1 = v1.y();
  REAL z1 = v1.z();
  REAL x2 = v2.x();
  REAL y2 = v2.y();
  REAL z2 = v2.z();

  // BB: Feb 3, 2017:
  // Not an efficient operation
  // Better approach may be to use a list if random access of vector
  // members is not needed
  REAL epsilon = EPS;
  periParticleVec.erase(
    std::remove_if(
      periParticleVec.begin(), periParticleVec.end(),
      [&x1, &y1, &z1, &x2, &y2, &z2, &epsilon](PeriParticleP particle) {
        Vec center = particle->currentPosition();
        if (!(center.x() - x1 >= -epsilon && center.x() - x2 < -epsilon &&
              center.y() - y1 >= -epsilon && center.y() - y2 < -epsilon &&
              center.z() - z1 >= -epsilon && center.z() - z2 < -epsilon)) {
          return true;
        }
        return false;
      }),
    periParticleVec.end());

  /*
  PeriParticlePArray::iterator itr;
  Vec center;
  //std::size_t flag = 0;

  for (itr = periParticleVec.begin(); itr != periParticleVec.end(); ) {
    center=(*itr)->currentPosition();
    // it is critical to use EPS
    if ( !(center.x() - x1 >= -EPS && center.x() - x2 < -EPS &&
       center.y() - y1 >= -EPS && center.y() - y2 < -EPS &&
       center.z() - z1 >= -EPS && center.z() - z2 < -EPS) )
  {
    // debugInf << "iter=" << std::setw(8) << iteration << " rank=" <<
  std::setw(2) << mpiRank
    // << " removed=" << std::setw(3) << (*itr)->getId();
    // flag = 1;
    delete (*itr); // release memory
    itr = periParticleVec.erase(itr);
  }
    else
  ++itr;
  }
  */
  /*
    if (flag == 1) {
    debugInf << " now " << particleVec.size() << ": ";
    for (ParticlePArray::const_iterator it = particleVec.begin(); it !=
    particleVec.end(); ++it)
    debugInf << std::setw(3) << (*it)->getId();
    debugInf << std::endl;
    }
  */
}

REAL
Assembly::getPtclMaxX(const ParticlePArray& inputParticle) const
{
  if (inputParticle.size() == 0)
    return -1 / EPS;

  auto it = inputParticle.cbegin();
  REAL x0 = (*it)->currentPos().x();
  for (; it != inputParticle.cend(); ++it) {
    if ((*it)->currentPos().x() > x0)
      x0 = (*it)->currentPos().x();
  }
  return x0;
}

REAL
Assembly::getPtclMinX(const ParticlePArray& inputParticle) const
{
  if (inputParticle.size() == 0)
    return 1 / EPS;

  auto it = inputParticle.cbegin();
  REAL x0 = (*it)->currentPos().x();
  for (; it != inputParticle.cend(); ++it) {
    if ((*it)->currentPos().x() < x0)
      x0 = (*it)->currentPos().x();
  }
  return x0;
}

REAL
Assembly::getPtclMaxY(const ParticlePArray& inputParticle) const
{
  if (inputParticle.size() == 0)
    return -1 / EPS;

  auto it = inputParticle.cbegin();
  REAL y0 = (*it)->currentPos().y();
  for (; it != inputParticle.cend(); ++it) {
    if ((*it)->currentPos().y() > y0)
      y0 = (*it)->currentPos().y();
  }
  return y0;
}

REAL
Assembly::getPtclMinY(const ParticlePArray& inputParticle) const
{
  if (inputParticle.size() == 0)
    return 1 / EPS;

  auto it = inputParticle.cbegin();
  REAL y0 = (*it)->currentPos().y();
  for (; it != inputParticle.cend(); ++it) {
    if ((*it)->currentPos().y() < y0)
      y0 = (*it)->currentPos().y();
  }
  return y0;
}

REAL
Assembly::getPtclMaxZ(const ParticlePArray& inputParticle) const
{
  if (inputParticle.size() == 0)
    return -1 / EPS;

  auto it = inputParticle.cbegin();
  REAL z0 = (*it)->currentPos().z();
  for (; it != inputParticle.cend(); ++it) {
    if ((*it)->currentPos().z() > z0)
      z0 = (*it)->currentPos().z();
  }
  return z0;
}

REAL
Assembly::getPtclMinZ(const ParticlePArray& inputParticle) const
{
  if (inputParticle.size() == 0)
    return 1 / EPS;

  auto it = inputParticle.cbegin();
  REAL z0 = (*it)->currentPos().z();
  for (; it != inputParticle.cend(); ++it) {
    if ((*it)->currentPos().z() < z0)
      z0 = (*it)->currentPos().z();
  }
  return z0;
}

void
Assembly::setCommunicator(boost::mpi::communicator& comm)
{
  boostWorld = comm;
  mpiWorld = MPI_Comm(comm);
  int mpiProcX = util::getParam<int>("mpiProcX");
  int mpiProcY = util::getParam<int>("mpiProcY");
  int mpiProcZ = util::getParam<int>("mpiProcZ");
  d_mpiProcs = {{ mpiProcX, mpiProcY, mpiProcZ }};

  // create Cartesian virtual topology (unavailable in boost.mpi)
  int ndim = 3;
  int periods[3] = { 0, 0, 0 };
  int reorder = 0; // mpiRank not reordered
  MPI_Cart_create(mpiWorld, ndim, d_mpiProcs.data(), periods, reorder, &cartComm);
  MPI_Comm_rank(cartComm, &mpiRank);
  MPI_Comm_size(cartComm, &mpiSize);
  MPI_Cart_coords(cartComm, mpiRank, ndim, d_mpiCoords.data());
  mpiTag = 0;
  assert(mpiRank == boostWorld.rank());
  // debugInf << mpiRank << " " << d_mpiCoords[0] << " " << d_mpiCoords[1] << " " <<
  // d_mpiCoords[2] << std::endl;

  for (int iRank = 0; iRank < mpiSize; ++iRank) {
    int ndim = 3;
    int coords[3];
    MPI_Cart_coords(cartComm, iRank, ndim, coords);
    if (coords[0] == 0 || coords[0] == mpiProcX - 1 || coords[1] == 0 ||
        coords[1] == mpiProcY - 1 || coords[2] == 0 ||
        coords[2] == mpiProcZ - 1)
      bdryProcess.push_back(iRank);
  }
}

// this is to scatter the dem and sph particle
// two point: (1) Partition the peri-points before any bonds constructed, all
// peri-points will be treated as free peri-points.
//            This is also what the coupling model shows, i.e. even the
//            peri-points that are bonded to DEM particles still have properties
//            of free peri-points
//            (2) Before partition there are no any bonds in PeriParticle, after
//            partition, constructNeighbor() is called, then these peri-bonds,
//            boundary-bonds,
//            and peri-DEM-bonds will be generated
void
Assembly::scatterDEMPeriParticle()
{
  // partition particles and send to each process
  if (mpiRank == 0) { // process 0
    setGrid(Box(allContainer.getMinCorner().x() - point_interval * 0.2,
                allContainer.getMinCorner().y() - point_interval * 0.2,
                allContainer.getMinCorner().z() - point_interval * 0.2,
                allContainer.getMaxCorner().x() + point_interval * 0.2,
                allContainer.getMaxCorner().y() + point_interval * 0.2,
                allContainer.getMaxCorner().z() + point_interval * 0.2));

    Vec v1 = grid.getMinCorner();
    Vec v2 = grid.getMaxCorner();
    Vec vspan = v2 - v1;

    auto reqs = new boost::mpi::request[mpiSize - 1];
    ParticlePArray tmpParticleVec;
    for (int iRank = mpiSize - 1; iRank >= 0; --iRank) {
      tmpParticleVec.clear(); // do not release memory!
      int ndim = 3;
      int coords[3];
      MPI_Cart_coords(cartComm, iRank, ndim, coords);
      Box container(v1.x() + vspan.x() / d_mpiProcs.x() * coords[0],
                    v1.y() + vspan.y() / d_mpiProcs.y() * coords[1],
                    v1.z() + vspan.z() / d_mpiProcs.z() * coords[2],
                    v1.x() + vspan.x() / d_mpiProcs.x() * (coords[0] + 1),
                    v1.y() + vspan.y() / d_mpiProcs.y() * (coords[1] + 1),
                    v1.z() + vspan.z() / d_mpiProcs.z() * (coords[2] + 1));
      findParticleInBox(container, allParticleVec, tmpParticleVec);
      if (iRank != 0)
        reqs[iRank - 1] = boostWorld.isend(iRank, mpiTag,
                                           tmpParticleVec); // non-blocking send
      // before send, the SPHParticle.demParticle == NULL, since NULL is
      // assigned when SPHParticle is created
      if (iRank == 0) {
        particleVec.resize(tmpParticleVec.size());
        for (auto i = 0u; i < particleVec.size(); ++i)
          particleVec[i] = std::make_shared<Particle>(
            *tmpParticleVec[i]); // default synthesized copy constructor
      } // now particleVec do not share memeory with allParticleVec
    }
    boost::mpi::wait_all(reqs, reqs + mpiSize - 1); // for non-blocking send
    delete[] reqs;

  } else { // other processes except 0
    boostWorld.recv(0, mpiTag, particleVec);
  }

  // content of allParticleVec may need to be printed, so do not clear it.
  // if (mpiRank == 0) releaseGatheredParticle();

  ///////////////////////////////////////////////////////////////////////////////
  // partition peri-points and send to each process
  if (mpiRank == 0) { // process 0
    setGrid(Box(allContainer.getMinCorner().x() - point_interval * 0.2,
                allContainer.getMinCorner().y() - point_interval * 0.2,
                allContainer.getMinCorner().z() - point_interval * 0.2,
                allContainer.getMaxCorner().x() + point_interval * 0.2,
                allContainer.getMaxCorner().y() + point_interval * 0.2,
                allContainer.getMaxCorner().z() + point_interval * 0.2));

    Vec v1 = grid.getMinCorner();
    Vec v2 = grid.getMaxCorner();
    Vec vspan = v2 - v1;
    auto reqs = new boost::mpi::request[mpiSize - 1];
    PeriParticlePArray tmpPeriParticleVec;
    for (int iRank = mpiSize - 1; iRank >= 0; --iRank) {
      tmpPeriParticleVec.clear(); // do not release memory!
      int ndim = 3;
      int coords[3];
      MPI_Cart_coords(cartComm, iRank, ndim, coords);
      Box container(v1.x() + vspan.x() / d_mpiProcs.x() * coords[0],
                    v1.y() + vspan.y() / d_mpiProcs.y() * coords[1],
                    v1.z() + vspan.z() / d_mpiProcs.z() * coords[2],
                    v1.x() + vspan.x() / d_mpiProcs.x() * (coords[0] + 1),
                    v1.y() + vspan.y() / d_mpiProcs.y() * (coords[1] + 1),
                    v1.z() + vspan.z() / d_mpiProcs.z() * (coords[2] + 1));
      findPeriParticleInBox(container, allPeriParticleVec, tmpPeriParticleVec);
      if (iRank != 0)
        reqs[iRank - 1] = boostWorld.isend(
          iRank, mpiTag, tmpPeriParticleVec); // non-blocking send
      if (iRank == 0) {
        periParticleVec.resize(tmpPeriParticleVec.size());
        for (auto i = 0u; i < periParticleVec.size(); ++i)
          periParticleVec[i] = std::make_shared<pd::PeriParticle>(
            *tmpPeriParticleVec[i]); // default synthesized copy constructor
      } // now particleVec do not share memeory with allParticleVec
    }
    boost::mpi::wait_all(reqs, reqs + mpiSize - 1); // for non-blocking send
    delete[] reqs;

  } else { // other processes except 0
    boostWorld.recv(0, mpiTag, periParticleVec);
  }

  // broadcast necessary info
  broadcast(boostWorld, gradation, 0);
  broadcast(boostWorld, boundaryVec, 0);
  broadcast(boostWorld, allContainer, 0);
  broadcast(boostWorld, grid, 0);
  broadcast(boostWorld, point_interval, 0);
  broadcast(boostWorld, maxHorizonSize, 0);
} // scatterDEMPeriParticle

bool
Assembly::isBdryProcess()
{
  return (d_mpiCoords.x() == 0 || d_mpiCoords.x() == d_mpiProcs.x() - 1 ||
          d_mpiCoords.y() == 0 || d_mpiCoords.y() == d_mpiProcs.y() - 1 ||
          d_mpiCoords.z() == 0 || d_mpiCoords.z() == d_mpiProcs.z() - 1);
}

// before the transfer to peri-points, their bondVec should be cleared
void
Assembly::commuPeriParticle()
{
  // determine container of each process
  Vec v1 = grid.getMinCorner();
  Vec v2 = grid.getMaxCorner();
  Vec vspan = v2 - v1;
  container = Box(v1.x() + vspan.x() / d_mpiProcs.x() * d_mpiCoords.x(),
                  v1.y() + vspan.y() / d_mpiProcs.y() * d_mpiCoords.y(),
                  v1.z() + vspan.z() / d_mpiProcs.z() * d_mpiCoords.z(),
                  v1.x() + vspan.x() / d_mpiProcs.x() * (d_mpiCoords.x() + 1),
                  v1.y() + vspan.y() / d_mpiProcs.y() * (d_mpiCoords.y() + 1),
                  v1.z() + vspan.z() / d_mpiProcs.z() * (d_mpiCoords.z() + 1));

  // if found, communicate with neighboring blocks
  PeriParticlePArray periParticleX1, periParticleX2;
  PeriParticlePArray periParticleY1, periParticleY2;
  PeriParticlePArray periParticleZ1, periParticleZ2;
  PeriParticlePArray periParticleX1Y1, periParticleX1Y2, periParticleX1Z1,
    periParticleX1Z2;
  PeriParticlePArray periParticleX2Y1, periParticleX2Y2, periParticleX2Z1,
    periParticleX2Z2;
  PeriParticlePArray periParticleY1Z1, periParticleY1Z2, periParticleY2Z1,
    periParticleY2Z2;
  PeriParticlePArray periParticleX1Y1Z1, periParticleX1Y1Z2, periParticleX1Y2Z1,
    periParticleX1Y2Z2;
  PeriParticlePArray periParticleX2Y1Z1, periParticleX2Y1Z2, periParticleX2Y2Z1,
    periParticleX2Y2Z2;
  boost::mpi::request reqX1[2], reqX2[2];
  boost::mpi::request reqY1[2], reqY2[2];
  boost::mpi::request reqZ1[2], reqZ2[2];
  boost::mpi::request reqX1Y1[2], reqX1Y2[2], reqX1Z1[2], reqX1Z2[2];
  boost::mpi::request reqX2Y1[2], reqX2Y2[2], reqX2Z1[2], reqX2Z2[2];
  boost::mpi::request reqY1Z1[2], reqY1Z2[2], reqY2Z1[2], reqY2Z2[2];
  boost::mpi::request reqX1Y1Z1[2], reqX1Y1Z2[2], reqX1Y2Z1[2], reqX1Y2Z2[2];
  boost::mpi::request reqX2Y1Z1[2], reqX2Y1Z2[2], reqX2Y2Z1[2], reqX2Y2Z2[2];
  v1 = container.getMinCorner(); // redefine v1, v2 in terms of process
  v2 = container.getMaxCorner();
  // debugInf << "rank=" << mpiRank << ' ' << v1.x() << ' ' << v1.y() << '
  // ' << v1.z() << ' '  << v2.x() << ' ' << v2.y() << ' ' << v2.z()
  // << std::endl;
  REAL cellSize = std::max(4 * maxHorizonSize,
                           gradation.getPtclMaxRadius() +
                             3 * point_interval); // constructNeighbor() is
                                                  // based on 2*horizonSize,
  // refer to PeriParitcle.calcAcceleration(), in this function the
  // deformationGradient and sigma of the other peri-points in the bond are
  // needed,
  // this means that the deformationGradient, sigma and Kinv of the other
  // peri-point should also be calculated even if it is in recvParticleVec, thus
  // we need to transfer 2*cellSize peri-points, the peri-points in the outer
  // cell are used to calculate the deformationGradient, sigma and Kinv
  // of the peri-points in inner cell
  // 6 surfaces
  if (rankX1 >= 0) { // surface x1
    Box containerX1(v1.x(), v1.y(), v1.z(), v1.x() + cellSize, v2.y(), v2.z());
    findPeriParticleInBox(containerX1, periParticleVec, periParticleX1);
    reqX1[0] = boostWorld.isend(rankX1, mpiTag, periParticleX1);
    reqX1[1] = boostWorld.irecv(rankX1, mpiTag, rperiParticleX1);
  }
  if (rankX2 >= 0) { // surface x2
    Box containerX2(v2.x() - cellSize, v1.y(), v1.z(), v2.x(), v2.y(), v2.z());
    findPeriParticleInBox(containerX2, periParticleVec, periParticleX2);
    reqX2[0] = boostWorld.isend(rankX2, mpiTag, periParticleX2);
    reqX2[1] = boostWorld.irecv(rankX2, mpiTag, rperiParticleX2);
  }
  if (rankY1 >= 0) { // surface y1
    Box containerY1(v1.x(), v1.y(), v1.z(), v2.x(), v1.y() + cellSize, v2.z());
    findPeriParticleInBox(containerY1, periParticleVec, periParticleY1);
    reqY1[0] = boostWorld.isend(rankY1, mpiTag, periParticleY1);
    reqY1[1] = boostWorld.irecv(rankY1, mpiTag, rperiParticleY1);
  }
  if (rankY2 >= 0) { // surface y2
    Box containerY2(v1.x(), v2.y() - cellSize, v1.z(), v2.x(), v2.y(), v2.z());
    findPeriParticleInBox(containerY2, periParticleVec, periParticleY2);
    reqY2[0] = boostWorld.isend(rankY2, mpiTag, periParticleY2);
    reqY2[1] = boostWorld.irecv(rankY2, mpiTag, rperiParticleY2);
  }
  if (rankZ1 >= 0) { // surface z1
    Box containerZ1(v1.x(), v1.y(), v1.z(), v2.x(), v2.y(), v1.z() + cellSize);
    findPeriParticleInBox(containerZ1, periParticleVec, periParticleZ1);
    reqZ1[0] = boostWorld.isend(rankZ1, mpiTag, periParticleZ1);
    reqZ1[1] = boostWorld.irecv(rankZ1, mpiTag, rperiParticleZ1);
  }
  if (rankZ2 >= 0) { // surface z2
    Box containerZ2(v1.x(), v1.y(), v2.z() - cellSize, v2.x(), v2.y(), v2.z());
    findPeriParticleInBox(containerZ2, periParticleVec, periParticleZ2);
    reqZ2[0] = boostWorld.isend(rankZ2, mpiTag, periParticleZ2);
    reqZ2[1] = boostWorld.irecv(rankZ2, mpiTag, rperiParticleZ2);
  }
  // 12 edges
  if (rankX1Y1 >= 0) { // edge x1y1
    Box containerX1Y1(v1.x(), v1.y(), v1.z(), v1.x() + cellSize,
                      v1.y() + cellSize, v2.z());
    findPeriParticleInBox(containerX1Y1, periParticleVec, periParticleX1Y1);
    reqX1Y1[0] = boostWorld.isend(rankX1Y1, mpiTag, periParticleX1Y1);
    reqX1Y1[1] = boostWorld.irecv(rankX1Y1, mpiTag, rperiParticleX1Y1);
  }
  if (rankX1Y2 >= 0) { // edge x1y2
    Box containerX1Y2(v1.x(), v2.y() - cellSize, v1.z(), v1.x() + cellSize,
                      v2.y(), v2.z());
    findPeriParticleInBox(containerX1Y2, periParticleVec, periParticleX1Y2);
    reqX1Y2[0] = boostWorld.isend(rankX1Y2, mpiTag, periParticleX1Y2);
    reqX1Y2[1] = boostWorld.irecv(rankX1Y2, mpiTag, rperiParticleX1Y2);
  }
  if (rankX1Z1 >= 0) { // edge x1z1
    Box containerX1Z1(v1.x(), v1.y(), v1.z(), v1.x() + cellSize, v2.y(),
                      v1.z() + cellSize);
    findPeriParticleInBox(containerX1Z1, periParticleVec, periParticleX1Z1);
    reqX1Z1[0] = boostWorld.isend(rankX1Z1, mpiTag, periParticleX1Z1);
    reqX1Z1[1] = boostWorld.irecv(rankX1Z1, mpiTag, rperiParticleX1Z1);
  }
  if (rankX1Z2 >= 0) { // edge x1z2
    Box containerX1Z2(v1.x(), v1.y(), v2.z() - cellSize, v1.x() + cellSize,
                      v2.y(), v2.z());
    findPeriParticleInBox(containerX1Z2, periParticleVec, periParticleX1Z2);
    reqX1Z2[0] = boostWorld.isend(rankX1Z2, mpiTag, periParticleX1Z2);
    reqX1Z2[1] = boostWorld.irecv(rankX1Z2, mpiTag, rperiParticleX1Z2);
  }
  if (rankX2Y1 >= 0) { // edge x2y1
    Box containerX2Y1(v2.x() - cellSize, v1.y(), v1.z(), v2.x(),
                      v1.y() + cellSize, v2.z());
    findPeriParticleInBox(containerX2Y1, periParticleVec, periParticleX2Y1);
    reqX2Y1[0] = boostWorld.isend(rankX2Y1, mpiTag, periParticleX2Y1);
    reqX2Y1[1] = boostWorld.irecv(rankX2Y1, mpiTag, rperiParticleX2Y1);
  }
  if (rankX2Y2 >= 0) { // edge x2y2
    Box containerX2Y2(v2.x() - cellSize, v2.y() - cellSize, v1.z(), v2.x(),
                      v2.y(), v2.z());
    findPeriParticleInBox(containerX2Y2, periParticleVec, periParticleX2Y2);
    reqX2Y2[0] = boostWorld.isend(rankX2Y2, mpiTag, periParticleX2Y2);
    reqX2Y2[1] = boostWorld.irecv(rankX2Y2, mpiTag, rperiParticleX2Y2);
  }
  if (rankX2Z1 >= 0) { // edge x2z1
    Box containerX2Z1(v2.x() - cellSize, v1.y(), v1.z(), v2.x(), v2.y(),
                      v1.z() + cellSize);
    findPeriParticleInBox(containerX2Z1, periParticleVec, periParticleX2Z1);
    reqX2Z1[0] = boostWorld.isend(rankX2Z1, mpiTag, periParticleX2Z1);
    reqX2Z1[1] = boostWorld.irecv(rankX2Z1, mpiTag, rperiParticleX2Z1);
  }
  if (rankX2Z2 >= 0) { // edge x2z2
    Box containerX2Z2(v2.x() - cellSize, v1.y(), v2.z() - cellSize, v2.x(),
                      v2.y(), v2.z());
    findPeriParticleInBox(containerX2Z2, periParticleVec, periParticleX2Z2);
    reqX2Z2[0] = boostWorld.isend(rankX2Z2, mpiTag, periParticleX2Z2);
    reqX2Z2[1] = boostWorld.irecv(rankX2Z2, mpiTag, rperiParticleX2Z2);
  }
  if (rankY1Z1 >= 0) { // edge y1z1
    Box containerY1Z1(v1.x(), v1.y(), v1.z(), v2.x(), v1.y() + cellSize,
                      v1.z() + cellSize);
    findPeriParticleInBox(containerY1Z1, periParticleVec, periParticleY1Z1);
    reqY1Z1[0] = boostWorld.isend(rankY1Z1, mpiTag, periParticleY1Z1);
    reqY1Z1[1] = boostWorld.irecv(rankY1Z1, mpiTag, rperiParticleY1Z1);
  }
  if (rankY1Z2 >= 0) { // edge y1z2
    Box containerY1Z2(v1.x(), v1.y(), v2.z() - cellSize, v2.x(),
                      v1.y() + cellSize, v2.z());
    findPeriParticleInBox(containerY1Z2, periParticleVec, periParticleY1Z2);
    reqY1Z2[0] = boostWorld.isend(rankY1Z2, mpiTag, periParticleY1Z2);
    reqY1Z2[1] = boostWorld.irecv(rankY1Z2, mpiTag, rperiParticleY1Z2);
  }
  if (rankY2Z1 >= 0) { // edge y2z1
    Box containerY2Z1(v1.x(), v2.y() - cellSize, v1.z(), v2.x(), v2.y(),
                      v1.z() + cellSize);
    findPeriParticleInBox(containerY2Z1, periParticleVec, periParticleY2Z1);
    reqY2Z1[0] = boostWorld.isend(rankY2Z1, mpiTag, periParticleY2Z1);
    reqY2Z1[1] = boostWorld.irecv(rankY2Z1, mpiTag, rperiParticleY2Z1);
  }
  if (rankY2Z2 >= 0) { // edge y2z2
    Box containerY2Z2(v1.x(), v2.y() - cellSize, v2.z() - cellSize, v2.x(),
                      v2.y(), v2.z());
    findPeriParticleInBox(containerY2Z2, periParticleVec, periParticleY2Z2);
    reqY2Z2[0] = boostWorld.isend(rankY2Z2, mpiTag, periParticleY2Z2);
    reqY2Z2[1] = boostWorld.irecv(rankY2Z2, mpiTag, rperiParticleY2Z2);
  }
  // 8 vertices
  if (rankX1Y1Z1 >= 0) { // edge x1y1z1
    Box containerX1Y1Z1(v1.x(), v1.y(), v1.z(), v1.x() + cellSize,
                        v1.y() + cellSize, v1.z() + cellSize);
    findPeriParticleInBox(containerX1Y1Z1, periParticleVec, periParticleX1Y1Z1);
    reqX1Y1Z1[0] = boostWorld.isend(rankX1Y1Z1, mpiTag, periParticleX1Y1Z1);
    reqX1Y1Z1[1] = boostWorld.irecv(rankX1Y1Z1, mpiTag, rperiParticleX1Y1Z1);
  }
  if (rankX1Y1Z2 >= 0) { // edge x1y1z2
    Box containerX1Y1Z2(v1.x(), v1.y(), v2.z() - cellSize, v1.x() + cellSize,
                        v1.y() + cellSize, v2.z());
    findPeriParticleInBox(containerX1Y1Z2, periParticleVec, periParticleX1Y1Z2);
    reqX1Y1Z2[0] = boostWorld.isend(rankX1Y1Z2, mpiTag, periParticleX1Y1Z2);
    reqX1Y1Z2[1] = boostWorld.irecv(rankX1Y1Z2, mpiTag, rperiParticleX1Y1Z2);
  }
  if (rankX1Y2Z1 >= 0) { // edge x1y2z1
    Box containerX1Y2Z1(v1.x(), v2.y() - cellSize, v1.z(), v1.x() + cellSize,
                        v2.y(), v1.z() + cellSize);
    findPeriParticleInBox(containerX1Y2Z1, periParticleVec, periParticleX1Y2Z1);
    reqX1Y2Z1[0] = boostWorld.isend(rankX1Y2Z1, mpiTag, periParticleX1Y2Z1);
    reqX1Y2Z1[1] = boostWorld.irecv(rankX1Y2Z1, mpiTag, rperiParticleX1Y2Z1);
  }
  if (rankX1Y2Z2 >= 0) { // edge x1y2z2
    Box containerX1Y2Z2(v1.x(), v2.y() - cellSize, v2.z() - cellSize,
                        v1.x() + cellSize, v2.y() + cellSize, v2.z());
    findPeriParticleInBox(containerX1Y2Z2, periParticleVec, periParticleX1Y2Z2);
    reqX1Y2Z2[0] = boostWorld.isend(rankX1Y2Z2, mpiTag, periParticleX1Y2Z2);
    reqX1Y2Z2[1] = boostWorld.irecv(rankX1Y2Z2, mpiTag, rperiParticleX1Y2Z2);
  }
  if (rankX2Y1Z1 >= 0) { // edge x2y1z1
    Box containerX2Y1Z1(v2.x() - cellSize, v1.y(), v1.z(), v2.x(),
                        v1.y() + cellSize, v1.z() + cellSize);
    findPeriParticleInBox(containerX2Y1Z1, periParticleVec, periParticleX2Y1Z1);
    reqX2Y1Z1[0] = boostWorld.isend(rankX2Y1Z1, mpiTag, periParticleX2Y1Z1);
    reqX2Y1Z1[1] = boostWorld.irecv(rankX2Y1Z1, mpiTag, rperiParticleX2Y1Z1);
  }
  if (rankX2Y1Z2 >= 0) { // edge x2y1z2
    Box containerX2Y1Z2(v2.x() - cellSize, v1.y(), v2.z() - cellSize, v2.x(),
                        v1.y() + cellSize, v2.z());
    findPeriParticleInBox(containerX2Y1Z2, periParticleVec, periParticleX2Y1Z2);
    reqX2Y1Z2[0] = boostWorld.isend(rankX2Y1Z2, mpiTag, periParticleX2Y1Z2);
    reqX2Y1Z2[1] = boostWorld.irecv(rankX2Y1Z2, mpiTag, rperiParticleX2Y1Z2);
  }
  if (rankX2Y2Z1 >= 0) { // edge x2y2z1
    Box containerX2Y2Z1(v2.x() - cellSize, v2.y() - cellSize, v1.z(), v2.x(),
                        v2.y(), v1.z() + cellSize);
    findPeriParticleInBox(containerX2Y2Z1, periParticleVec, periParticleX2Y2Z1);
    reqX2Y2Z1[0] = boostWorld.isend(rankX2Y2Z1, mpiTag, periParticleX2Y2Z1);
    reqX2Y2Z1[1] = boostWorld.irecv(rankX2Y2Z1, mpiTag, rperiParticleX2Y2Z1);
  }
  if (rankX2Y2Z2 >= 0) { // edge x2y2z2
    Box containerX2Y2Z2(v2.x() - cellSize, v2.y() - cellSize, v2.z() - cellSize,
                        v2.x(), v2.y(), v2.z());
    findPeriParticleInBox(containerX2Y2Z2, periParticleVec, periParticleX2Y2Z2);
    reqX2Y2Z2[0] = boostWorld.isend(rankX2Y2Z2, mpiTag, periParticleX2Y2Z2);
    reqX2Y2Z2[1] = boostWorld.irecv(rankX2Y2Z2, mpiTag, rperiParticleX2Y2Z2);
  }

  // 6 surfaces
  if (rankX1 >= 0)
    boost::mpi::wait_all(reqX1, reqX1 + 2);
  if (rankX2 >= 0)
    boost::mpi::wait_all(reqX2, reqX2 + 2);
  if (rankY1 >= 0)
    boost::mpi::wait_all(reqY1, reqY1 + 2);
  if (rankY2 >= 0)
    boost::mpi::wait_all(reqY2, reqY2 + 2);
  if (rankZ1 >= 0)
    boost::mpi::wait_all(reqZ1, reqZ1 + 2);
  if (rankZ2 >= 0)
    boost::mpi::wait_all(reqZ2, reqZ2 + 2);
  // 12 edges
  if (rankX1Y1 >= 0)
    boost::mpi::wait_all(reqX1Y1, reqX1Y1 + 2);
  if (rankX1Y2 >= 0)
    boost::mpi::wait_all(reqX1Y2, reqX1Y2 + 2);
  if (rankX1Z1 >= 0)
    boost::mpi::wait_all(reqX1Z1, reqX1Z1 + 2);
  if (rankX1Z2 >= 0)
    boost::mpi::wait_all(reqX1Z2, reqX1Z2 + 2);
  if (rankX2Y1 >= 0)
    boost::mpi::wait_all(reqX2Y1, reqX2Y1 + 2);
  if (rankX2Y2 >= 0)
    boost::mpi::wait_all(reqX2Y2, reqX2Y2 + 2);
  if (rankX2Z1 >= 0)
    boost::mpi::wait_all(reqX2Z1, reqX2Z1 + 2);
  if (rankX2Z2 >= 0)
    boost::mpi::wait_all(reqX2Z2, reqX2Z2 + 2);
  if (rankY1Z1 >= 0)
    boost::mpi::wait_all(reqY1Z1, reqY1Z1 + 2);
  if (rankY1Z2 >= 0)
    boost::mpi::wait_all(reqY1Z2, reqY1Z2 + 2);
  if (rankY2Z1 >= 0)
    boost::mpi::wait_all(reqY2Z1, reqY2Z1 + 2);
  if (rankY2Z2 >= 0)
    boost::mpi::wait_all(reqY2Z2, reqY2Z2 + 2);
  // 8 vertices
  if (rankX1Y1Z1 >= 0)
    boost::mpi::wait_all(reqX1Y1Z1, reqX1Y1Z1 + 2);
  if (rankX1Y1Z2 >= 0)
    boost::mpi::wait_all(reqX1Y1Z2, reqX1Y1Z2 + 2);
  if (rankX1Y2Z1 >= 0)
    boost::mpi::wait_all(reqX1Y2Z1, reqX1Y2Z1 + 2);
  if (rankX1Y2Z2 >= 0)
    boost::mpi::wait_all(reqX1Y2Z2, reqX1Y2Z2 + 2);
  if (rankX2Y1Z1 >= 0)
    boost::mpi::wait_all(reqX2Y1Z1, reqX2Y1Z1 + 2);
  if (rankX2Y1Z2 >= 0)
    boost::mpi::wait_all(reqX2Y1Z2, reqX2Y1Z2 + 2);
  if (rankX2Y2Z1 >= 0)
    boost::mpi::wait_all(reqX2Y2Z1, reqX2Y2Z1 + 2);
  if (rankX2Y2Z2 >= 0)
    boost::mpi::wait_all(reqX2Y2Z2, reqX2Y2Z2 + 2);

  // merge: periParticles inside container (at front) + periParticles from
  // neighoring blocks (at end)
  recvPeriParticleVec.clear();
  // 6 surfaces
  if (rankX1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX1.begin(), rperiParticleX1.end());
  if (rankX2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX2.begin(), rperiParticleX2.end());
  if (rankY1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleY1.begin(), rperiParticleY1.end());
  if (rankY2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleY2.begin(), rperiParticleY2.end());
  if (rankZ1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleZ1.begin(), rperiParticleZ1.end());
  if (rankZ2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleZ2.begin(), rperiParticleZ2.end());
  // 12 edges
  if (rankX1Y1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX1Y1.begin(),
                               rperiParticleX1Y1.end());
  if (rankX1Y2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX1Y2.begin(),
                               rperiParticleX1Y2.end());
  if (rankX1Z1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX1Z1.begin(),
                               rperiParticleX1Z1.end());
  if (rankX1Z2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX1Z2.begin(),
                               rperiParticleX1Z2.end());
  if (rankX2Y1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX2Y1.begin(),
                               rperiParticleX2Y1.end());
  if (rankX2Y2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX2Y2.begin(),
                               rperiParticleX2Y2.end());
  if (rankX2Z1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX2Z1.begin(),
                               rperiParticleX2Z1.end());
  if (rankX2Z2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX2Z2.begin(),
                               rperiParticleX2Z2.end());
  if (rankY1Z1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleY1Z1.begin(),
                               rperiParticleY1Z1.end());
  if (rankY1Z2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleY1Z2.begin(),
                               rperiParticleY1Z2.end());
  if (rankY2Z1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleY2Z1.begin(),
                               rperiParticleY2Z1.end());
  if (rankY2Z2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleY2Z2.begin(),
                               rperiParticleY2Z2.end());
  // 8 vertices
  if (rankX1Y1Z1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX1Y1Z1.begin(),
                               rperiParticleX1Y1Z1.end());
  if (rankX1Y1Z2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX1Y1Z2.begin(),
                               rperiParticleX1Y1Z2.end());
  if (rankX1Y2Z1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX1Y2Z1.begin(),
                               rperiParticleX1Y2Z1.end());
  if (rankX1Y2Z2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX1Y2Z2.begin(),
                               rperiParticleX1Y2Z2.end());
  if (rankX2Y1Z1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX2Y1Z1.begin(),
                               rperiParticleX2Y1Z1.end());
  if (rankX2Y1Z2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX2Y1Z2.begin(),
                               rperiParticleX2Y1Z2.end());
  if (rankX2Y2Z1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX2Y2Z1.begin(),
                               rperiParticleX2Y2Z1.end());
  if (rankX2Y2Z2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX2Y2Z2.begin(),
                               rperiParticleX2Y2Z2.end());

  mergePeriParticleVec.clear();
  mergePeriParticleVec =
    periParticleVec; // duplicate pointers, pointing to the same memory
  mergePeriParticleVec.insert(mergePeriParticleVec.end(),
                              recvPeriParticleVec.begin(),
                              recvPeriParticleVec.end());

  /*
    ParticlePArray testParticleVec;
    testParticleVec.insert(testParticleVec.end(), rParticleX1.begin(),
    rParticleX1.end());
    testParticleVec.insert(testParticleVec.end(), rParticleX2.begin(),
    rParticleX2.end());
    testParticleVec.insert(testParticleVec.end(), rParticleY1.begin(),
    rParticleY1.end());
    testParticleVec.insert(testParticleVec.end(), rParticleY2.begin(),
    rParticleY2.end());
    testParticleVec.insert(testParticleVec.end(), rParticleZ1.begin(),
    rParticleZ1.end());
    testParticleVec.insert(testParticleVec.end(), rParticleZ2.begin(),
    rParticleZ2.end());
    debugInf << "iter=" << std::setw(4) << iteration << " rank=" << std::setw(4)
    << mpiRank
    << " ptclNum=" << std::setw(4) << particleVec.size()
    << " surface="
    << std::setw(4) << particleX1.size()  << std::setw(4) << particleX2.size()
    << std::setw(4) << particleY1.size()  << std::setw(4) << particleY2.size()
    << std::setw(4) << particleZ1.size()  << std::setw(4) << particleZ2.size()
    << " recv="
    << std::setw(4) << rParticleX1.size() << std::setw(4) << rParticleX2.size()
    << std::setw(4) << rParticleY1.size() << std::setw(4) << rParticleY2.size()
    << std::setw(4) << rParticleZ1.size() << std::setw(4) << rParticleZ2.size()
    << " rNum="
    << std::setw(4) << recvParticleVec.size() << ": ";

    for (ParticlePArray::const_iterator it = testParticleVec.begin(); it !=
    testParticleVec.end();++it)
    debugInf << (*it)->getId() << ' ';
    debugInf << std::endl;
    testParticleVec.clear();
  */
}

void
Assembly::releaseRecvParticle()
{
  // release memory of received particles
  /*
  for (ParticlePArray::iterator it = recvParticleVec.begin(); it !=
  recvParticleVec.end(); ++it)
    delete (*it);
  */
  recvParticleVec.clear();
  // 6 surfaces
  rParticleX1.clear();
  rParticleX2.clear();
  rParticleY1.clear();
  rParticleY2.clear();
  rParticleZ1.clear();
  rParticleZ2.clear();
  // 12 edges
  rParticleX1Y1.clear();
  rParticleX1Y2.clear();
  rParticleX1Z1.clear();
  rParticleX1Z2.clear();
  rParticleX2Y1.clear();
  rParticleX2Y2.clear();
  rParticleX2Z1.clear();
  rParticleX2Z2.clear();
  rParticleY1Z1.clear();
  rParticleY1Z2.clear();
  rParticleY2Z1.clear();
  rParticleY2Z2.clear();
  // 8 vertices
  rParticleX1Y1Z1.clear();
  rParticleX1Y1Z2.clear();
  rParticleX1Y2Z1.clear();
  rParticleX1Y2Z2.clear();
  rParticleX2Y1Z1.clear();
  rParticleX2Y1Z2.clear();
  rParticleX2Y2Z1.clear();
  rParticleX2Y2Z2.clear();
}

void
Assembly::releaseRecvPeriParticle()
{
  //    // delete those periBonds between recvPeriParticle
  //    for (PeriParticlePArray::iterator it = periParticleVec.begin(); it !=
  //    periParticleVec.end(); ++it){
  //    (*it)->eraseRecvPeriBonds();
  //    }
  // for(auto pt=recvPeriBondVec.begin(); pt!=recvPeriBondVec.end(); pt++){
  //   delete (*pt);
  //}
  recvPeriBondVec.clear();
  PeriBondPArray().swap(recvPeriBondVec); // actual memory release

  // release memory of received particles
  // for (auto it = recvPeriParticleVec.begin(); it !=
  // recvPeriParticleVec.end(); ++it){
  //  delete (*it);
  //}
  recvPeriParticleVec.clear();
  PeriParticlePArray().swap(recvPeriParticleVec); // actual memory release
  // 6 surfaces
  rperiParticleX1.clear();
  rperiParticleX2.clear();
  rperiParticleY1.clear();
  rperiParticleY2.clear();
  rperiParticleZ1.clear();
  rperiParticleZ2.clear();
  // 12 edges
  rperiParticleX1Y1.clear();
  rperiParticleX1Y2.clear();
  rperiParticleX1Z1.clear();
  rperiParticleX1Z2.clear();
  rperiParticleX2Y1.clear();
  rperiParticleX2Y2.clear();
  rperiParticleX2Z1.clear();
  rperiParticleX2Z2.clear();
  rperiParticleY1Z1.clear();
  rperiParticleY1Z2.clear();
  rperiParticleY2Z1.clear();
  rperiParticleY2Z2.clear();
  // 8 vertices
  rperiParticleX1Y1Z1.clear();
  rperiParticleX1Y1Z2.clear();
  rperiParticleX1Y2Z1.clear();
  rperiParticleX1Y2Z2.clear();
  rperiParticleX2Y1Z1.clear();
  rperiParticleX2Y1Z2.clear();
  rperiParticleX2Y2Z1.clear();
  rperiParticleX2Y2Z2.clear();
}

void
Assembly::releasePeriBondVec()
{
  // for(auto pt=periBondVec.begin(); pt!=periBondVec.end(); pt++){
  //  delete (*pt);
  //}
  periBondVec.clear();
  PeriBondPArray().swap(periBondVec); // actual memory release
}

void
Assembly::updatePeriGrid()
{
  PeriParticlePArray::const_iterator pit = periParticleVec.begin();
  REAL pMaxX = (*pit)->currentPosition().x();
  REAL pMinX = (*pit)->currentPosition().x();
  REAL pMaxY = (*pit)->currentPosition().y();
  REAL pMinY = (*pit)->currentPosition().y();
  REAL pMaxZ = (*pit)->currentPosition().z();
  REAL pMinZ = (*pit)->currentPosition().z();

  REAL tmpx, tmpy, tmpz;
  dem::Vec tmp_xyz;
  for (pit = periParticleVec.begin(); pit != periParticleVec.end(); pit++) {
    tmp_xyz = (*pit)->currentPosition();
    tmpx = tmp_xyz.x();
    tmpy = tmp_xyz.y();
    tmpz = tmp_xyz.z();
    if (pMaxX < tmpx)
      pMaxX = tmpx;
    if (pMinX > tmpx)
      pMinX = tmpx;
    if (pMaxY < tmpy)
      pMaxY = tmpy;
    if (pMinY > tmpy)
      pMinY = tmpy;
    if (pMaxZ < tmpz)
      pMaxZ = tmpz;
    if (pMinZ > tmpz)
      pMinZ = tmpz;
  }
  REAL maxX = 0;
  REAL maxY = 0;
  REAL maxZ = 0;
  REAL minX = 0;
  REAL minY = 0;
  REAL minZ = 0;
  MPI_Allreduce(&pMaxX, &maxX, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);
  MPI_Allreduce(&pMinX, &minX, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);
  MPI_Allreduce(&pMaxY, &maxY, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);
  MPI_Allreduce(&pMinY, &minY, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);
  MPI_Allreduce(&pMaxZ, &maxZ, 1, MPI_DOUBLE, MPI_MAX, mpiWorld);
  MPI_Allreduce(&pMinZ, &minZ, 1, MPI_DOUBLE, MPI_MIN, mpiWorld);

  // no need to broadcast grid as it is updated in each process
  setGrid(Box(minX - point_interval * 0.2, minY - point_interval * 0.2,
              minZ - point_interval * 0.2, maxX + point_interval * 0.2,
              maxY + point_interval * 0.2, maxZ + point_interval * 0.2));
}

void
Assembly::migrateParticle()
{
  //std::ostringstream out;
  //out << "Migrate: Rank: " << mpiRank << ": in: " << particleVec.size();

  Vec vspan = grid.getMaxCorner() - grid.getMinCorner();
  Vec width = vspan / d_mpiProcs;

  sentParticleVec.clear();
  recvParticleVec.clear();
  d_patchP->sendRecvMigrateXMinus(boostWorld, iteration, width, particleVec);
  d_patchP->sendRecvMigrateXPlus(boostWorld, iteration, width, particleVec);
  d_patchP->waitToFinishX(iteration);
  d_patchP->combineSentParticlesX(iteration, sentParticleVec);
  d_patchP->combineReceivedParticlesX(iteration, recvParticleVec);
  d_patchP->deleteSentParticles(iteration, sentParticleVec, particleVec);
  d_patchP->addReceivedParticles(iteration, recvParticleVec, particleVec);
  //out << " sentX : " << sentParticleVec.size()
  //    << " recvX : " << recvParticleVec.size();

  sentParticleVec.clear();
  recvParticleVec.clear();
  d_patchP->sendRecvMigrateYMinus(boostWorld, iteration, width, particleVec);
  d_patchP->sendRecvMigrateYPlus(boostWorld, iteration, width, particleVec);
  d_patchP->waitToFinishY(iteration);
  d_patchP->combineSentParticlesY(iteration, sentParticleVec);
  d_patchP->combineReceivedParticlesY(iteration, recvParticleVec);
  d_patchP->deleteSentParticles(iteration, sentParticleVec, particleVec);
  d_patchP->addReceivedParticles(iteration, recvParticleVec, particleVec);
  //out << " sentY : " << sentParticleVec.size()
  //    << " recvY : " << recvParticleVec.size();

  sentParticleVec.clear();
  recvParticleVec.clear();
  d_patchP->sendRecvMigrateZMinus(boostWorld, iteration, width, particleVec);
  d_patchP->sendRecvMigrateZPlus(boostWorld, iteration, width, particleVec);
  d_patchP->waitToFinishZ(iteration);
  d_patchP->combineSentParticlesZ(iteration, sentParticleVec);
  d_patchP->combineReceivedParticlesZ(iteration, recvParticleVec);
  d_patchP->deleteSentParticles(iteration, sentParticleVec, particleVec);
  d_patchP->addReceivedParticles(iteration, recvParticleVec, particleVec);
  //out << " sentZ : " << sentParticleVec.size()
  //    << " recvZ : " << recvParticleVec.size();

  // delete outgoing particles
  d_patchP->removeParticlesOutsidePatch(particleVec);
  //out << " outside: " << particleVec.size();

  //d_patchP->removeDuplicates(particleVec);
  //out <<  ": dup out: " << particleVec.size() << "\n";
  //std::cout << out.str();

}

/*
void
Assembly::migrateParticle()
{
*/
/*
  Vec vspan = grid.getMaxCorner() - grid.getMinCorner();
  Vec seg = vspan / d_mpiProcs;
  REAL segX = seg.x();
  REAL segY = seg.y();
  REAL segZ = seg.z();
  Vec v1 = container.getMinCorner(); // v1, v2 in terms of process
  Vec v2 = container.getMaxCorner();

  // if a neighbor exists, transfer particles crossing the boundary in between.
  ParticlePArray particleX1, particleX2;
  ParticlePArray particleY1, particleY2;
  ParticlePArray particleZ1, particleZ2;
  ParticlePArray particleX1Y1, particleX1Y2, particleX1Z1, particleX1Z2;
  ParticlePArray particleX2Y1, particleX2Y2, particleX2Z1, particleX2Z2;
  ParticlePArray particleY1Z1, particleY1Z2, particleY2Z1, particleY2Z2;
  ParticlePArray particleX1Y1Z1, particleX1Y1Z2, particleX1Y2Z1, particleX1Y2Z2;
  ParticlePArray particleX2Y1Z1, particleX2Y1Z2, particleX2Y2Z1, particleX2Y2Z2;
  boost::mpi::request reqX1[2], reqX2[2];
  boost::mpi::request reqY1[2], reqY2[2];
  boost::mpi::request reqZ1[2], reqZ2[2];
  boost::mpi::request reqX1Y1[2], reqX1Y2[2], reqX1Z1[2], reqX1Z2[2];
  boost::mpi::request reqX2Y1[2], reqX2Y2[2], reqX2Z1[2], reqX2Z2[2];
  boost::mpi::request reqY1Z1[2], reqY1Z2[2], reqY2Z1[2], reqY2Z2[2];
  boost::mpi::request reqX1Y1Z1[2], reqX1Y1Z2[2], reqX1Y2Z1[2], reqX1Y2Z2[2];
  boost::mpi::request reqX2Y1Z1[2], reqX2Y1Z2[2], reqX2Y2Z1[2], reqX2Y2Z2[2];

  // 6 surfaces
  if (rankX1 >= 0) { // surface x1
    Box containerX1(v1.x() - segX, v1.y(), v1.z(), v1.x(), v2.y(), v2.z());
    findParticleInBox(containerX1, particleVec, particleX1);
    reqX1[0] = boostWorld.isend(rankX1, mpiTag, particleX1);
    reqX1[1] = boostWorld.irecv(rankX1, mpiTag, rParticleX1);
  }
  if (rankX2 >= 0) { // surface x2
    Box containerX2(v2.x(), v1.y(), v1.z(), v2.x() + segX, v2.y(), v2.z());
    findParticleInBox(containerX2, particleVec, particleX2);
    reqX2[0] = boostWorld.isend(rankX2, mpiTag, particleX2);
    reqX2[1] = boostWorld.irecv(rankX2, mpiTag, rParticleX2);
  }
  if (rankY1 >= 0) { // surface y1
    Box containerY1(v1.x(), v1.y() - segY, v1.z(), v2.x(), v1.y(), v2.z());
    findParticleInBox(containerY1, particleVec, particleY1);
    reqY1[0] = boostWorld.isend(rankY1, mpiTag, particleY1);
    reqY1[1] = boostWorld.irecv(rankY1, mpiTag, rParticleY1);
  }
  if (rankY2 >= 0) { // surface y2
    Box containerY2(v1.x(), v2.y(), v1.z(), v2.x(), v2.y() + segY, v2.z());
    findParticleInBox(containerY2, particleVec, particleY2);
    reqY2[0] = boostWorld.isend(rankY2, mpiTag, particleY2);
    reqY2[1] = boostWorld.irecv(rankY2, mpiTag, rParticleY2);
  }
  if (rankZ1 >= 0) { // surface z1
    Box containerZ1(v1.x(), v1.y(), v1.z() - segZ, v2.x(), v2.y(), v1.z());
    findParticleInBox(containerZ1, particleVec, particleZ1);
    reqZ1[0] = boostWorld.isend(rankZ1, mpiTag, particleZ1);
    reqZ1[1] = boostWorld.irecv(rankZ1, mpiTag, rParticleZ1);
  }
  if (rankZ2 >= 0) { // surface z2
    Box containerZ2(v1.x(), v1.y(), v2.z(), v2.x(), v2.y(), v2.z() + segZ);
    findParticleInBox(containerZ2, particleVec, particleZ2);
    reqZ2[0] = boostWorld.isend(rankZ2, mpiTag, particleZ2);
    reqZ2[1] = boostWorld.irecv(rankZ2, mpiTag, rParticleZ2);
  }
  // 12 edges
  if (rankX1Y1 >= 0) { // edge x1y1
    Box containerX1Y1(v1.x() - segX, v1.y() - segY, v1.z(), v1.x(), v1.y(),
                      v2.z());
    findParticleInBox(containerX1Y1, particleVec, particleX1Y1);
    reqX1Y1[0] = boostWorld.isend(rankX1Y1, mpiTag, particleX1Y1);
    reqX1Y1[1] = boostWorld.irecv(rankX1Y1, mpiTag, rParticleX1Y1);
  }
  if (rankX1Y2 >= 0) { // edge x1y2
    Box containerX1Y2(v1.x() - segX, v2.y(), v1.z(), v1.x(), v2.y() + segY,
                      v2.z());
    findParticleInBox(containerX1Y2, particleVec, particleX1Y2);
    reqX1Y2[0] = boostWorld.isend(rankX1Y2, mpiTag, particleX1Y2);
    reqX1Y2[1] = boostWorld.irecv(rankX1Y2, mpiTag, rParticleX1Y2);
  }
  if (rankX1Z1 >= 0) { // edge x1z1
    Box containerX1Z1(v1.x() - segX, v1.y(), v1.z() - segZ, v1.x(), v2.y(),
                      v1.z());
    findParticleInBox(containerX1Z1, particleVec, particleX1Z1);
    reqX1Z1[0] = boostWorld.isend(rankX1Z1, mpiTag, particleX1Z1);
    reqX1Z1[1] = boostWorld.irecv(rankX1Z1, mpiTag, rParticleX1Z1);
  }
  if (rankX1Z2 >= 0) { // edge x1z2
    Box containerX1Z2(v1.x() - segX, v1.y(), v2.z(), v1.x(), v2.y(),
                      v2.z() + segZ);
    findParticleInBox(containerX1Z2, particleVec, particleX1Z2);
    reqX1Z2[0] = boostWorld.isend(rankX1Z2, mpiTag, particleX1Z2);
    reqX1Z2[1] = boostWorld.irecv(rankX1Z2, mpiTag, rParticleX1Z2);
  }
  if (rankX2Y1 >= 0) { // edge x2y1
    Box containerX2Y1(v2.x(), v1.y() - segY, v1.z(), v2.x() + segX, v1.y(),
                      v2.z());
    findParticleInBox(containerX2Y1, particleVec, particleX2Y1);
    reqX2Y1[0] = boostWorld.isend(rankX2Y1, mpiTag, particleX2Y1);
    reqX2Y1[1] = boostWorld.irecv(rankX2Y1, mpiTag, rParticleX2Y1);
  }
  if (rankX2Y2 >= 0) { // edge x2y2
    Box containerX2Y2(v2.x(), v2.y(), v1.z(), v2.x() + segX, v2.y() + segY,
                      v2.z());
    findParticleInBox(containerX2Y2, particleVec, particleX2Y2);
    reqX2Y2[0] = boostWorld.isend(rankX2Y2, mpiTag, particleX2Y2);
    reqX2Y2[1] = boostWorld.irecv(rankX2Y2, mpiTag, rParticleX2Y2);
  }
  if (rankX2Z1 >= 0) { // edge x2z1
    Box containerX2Z1(v2.x(), v1.y(), v1.z() - segZ, v2.x() + segX, v2.y(),
                      v1.z());
    findParticleInBox(containerX2Z1, particleVec, particleX2Z1);
    reqX2Z1[0] = boostWorld.isend(rankX2Z1, mpiTag, particleX2Z1);
    reqX2Z1[1] = boostWorld.irecv(rankX2Z1, mpiTag, rParticleX2Z1);
  }
  if (rankX2Z2 >= 0) { // edge x2z2
    Box containerX2Z2(v2.x(), v1.y(), v2.z(), v2.x() + segX, v2.y(),
                      v2.z() + segZ);
    findParticleInBox(containerX2Z2, particleVec, particleX2Z2);
    reqX2Z2[0] = boostWorld.isend(rankX2Z2, mpiTag, particleX2Z2);
    reqX2Z2[1] = boostWorld.irecv(rankX2Z2, mpiTag, rParticleX2Z2);
  }
  if (rankY1Z1 >= 0) { // edge y1z1
    Box containerY1Z1(v1.x(), v1.y() - segY, v1.z() - segZ, v2.x(), v1.y(),
                      v1.z());
    findParticleInBox(containerY1Z1, particleVec, particleY1Z1);
    reqY1Z1[0] = boostWorld.isend(rankY1Z1, mpiTag, particleY1Z1);
    reqY1Z1[1] = boostWorld.irecv(rankY1Z1, mpiTag, rParticleY1Z1);
  }
  if (rankY1Z2 >= 0) { // edge y1z2
    Box containerY1Z2(v1.x(), v1.y() - segY, v2.z(), v2.x(), v1.y(),
                      v2.z() + segZ);
    findParticleInBox(containerY1Z2, particleVec, particleY1Z2);
    reqY1Z2[0] = boostWorld.isend(rankY1Z2, mpiTag, particleY1Z2);
    reqY1Z2[1] = boostWorld.irecv(rankY1Z2, mpiTag, rParticleY1Z2);
  }
  if (rankY2Z1 >= 0) { // edge y2z1
    Box containerY2Z1(v1.x(), v2.y(), v1.z() - segZ, v2.x(), v2.y() + segY,
                      v1.z());
    findParticleInBox(containerY2Z1, particleVec, particleY2Z1);
    reqY2Z1[0] = boostWorld.isend(rankY2Z1, mpiTag, particleY2Z1);
    reqY2Z1[1] = boostWorld.irecv(rankY2Z1, mpiTag, rParticleY2Z1);
  }
  if (rankY2Z2 >= 0) { // edge y2z2
    Box containerY2Z2(v1.x(), v2.y(), v2.z(), v2.x(), v2.y() + segY,
                      v2.z() + segZ);
    findParticleInBox(containerY2Z2, particleVec, particleY2Z2);
    reqY2Z2[0] = boostWorld.isend(rankY2Z2, mpiTag, particleY2Z2);
    reqY2Z2[1] = boostWorld.irecv(rankY2Z2, mpiTag, rParticleY2Z2);
  }
  // 8 vertices
  if (rankX1Y1Z1 >= 0) { // edge x1y1z1
    Box containerX1Y1Z1(v1.x() - segX, v1.y() - segY, v1.z() - segZ, v1.x(),
                        v1.y(), v1.z());
    findParticleInBox(containerX1Y1Z1, particleVec, particleX1Y1Z1);
    reqX1Y1Z1[0] = boostWorld.isend(rankX1Y1Z1, mpiTag, particleX1Y1Z1);
    reqX1Y1Z1[1] = boostWorld.irecv(rankX1Y1Z1, mpiTag, rParticleX1Y1Z1);
  }
  if (rankX1Y1Z2 >= 0) { // edge x1y1z2
    Box containerX1Y1Z2(v1.x() - segX, v1.y() - segY, v2.z(), v1.x(), v1.y(),
                        v2.z() + segZ);
    findParticleInBox(containerX1Y1Z2, particleVec, particleX1Y1Z2);
    reqX1Y1Z2[0] = boostWorld.isend(rankX1Y1Z2, mpiTag, particleX1Y1Z2);
    reqX1Y1Z2[1] = boostWorld.irecv(rankX1Y1Z2, mpiTag, rParticleX1Y1Z2);
  }
  if (rankX1Y2Z1 >= 0) { // edge x1y2z1
    Box containerX1Y2Z1(v1.x() - segX, v2.y(), v1.z() - segZ, v1.x(),
                        v2.y() + segY, v1.z());
    findParticleInBox(containerX1Y2Z1, particleVec, particleX1Y2Z1);
    reqX1Y2Z1[0] = boostWorld.isend(rankX1Y2Z1, mpiTag, particleX1Y2Z1);
    reqX1Y2Z1[1] = boostWorld.irecv(rankX1Y2Z1, mpiTag, rParticleX1Y2Z1);
  }
  if (rankX1Y2Z2 >= 0) { // edge x1y2z2
    Box containerX1Y2Z2(v1.x() - segX, v2.y(), v2.z(), v1.x(), v2.y() + segY,
                        v2.z() + segZ);
    findParticleInBox(containerX1Y2Z2, particleVec, particleX1Y2Z2);
    reqX1Y2Z2[0] = boostWorld.isend(rankX1Y2Z2, mpiTag, particleX1Y2Z2);
    reqX1Y2Z2[1] = boostWorld.irecv(rankX1Y2Z2, mpiTag, rParticleX1Y2Z2);
  }
  if (rankX2Y1Z1 >= 0) { // edge x2y1z1
    Box containerX2Y1Z1(v2.x(), v1.y() - segY, v1.z() - segZ, v2.x() + segX,
                        v1.y(), v1.z());
    findParticleInBox(containerX2Y1Z1, particleVec, particleX2Y1Z1);
    reqX2Y1Z1[0] = boostWorld.isend(rankX2Y1Z1, mpiTag, particleX2Y1Z1);
    reqX2Y1Z1[1] = boostWorld.irecv(rankX2Y1Z1, mpiTag, rParticleX2Y1Z1);
  }
  if (rankX2Y1Z2 >= 0) { // edge x2y1z2
    Box containerX2Y1Z2(v2.x(), v1.y() - segY, v2.z(), v2.x() + segX, v1.y(),
                        v2.z() + segZ);
    findParticleInBox(containerX2Y1Z2, particleVec, particleX2Y1Z2);
    reqX2Y1Z2[0] = boostWorld.isend(rankX2Y1Z2, mpiTag, particleX2Y1Z2);
    reqX2Y1Z2[1] = boostWorld.irecv(rankX2Y1Z2, mpiTag, rParticleX2Y1Z2);
  }
  if (rankX2Y2Z1 >= 0) { // edge x2y2z1
    Box containerX2Y2Z1(v2.x(), v2.y(), v1.z() - segZ, v2.x() + segX,
                        v2.y() + segY, v1.z());
    findParticleInBox(containerX2Y2Z1, particleVec, particleX2Y2Z1);
    reqX2Y2Z1[0] = boostWorld.isend(rankX2Y2Z1, mpiTag, particleX2Y2Z1);
    reqX2Y2Z1[1] = boostWorld.irecv(rankX2Y2Z1, mpiTag, rParticleX2Y2Z1);
  }
  if (rankX2Y2Z2 >= 0) { // edge x2y2z2
    Box containerX2Y2Z2(v2.x(), v2.y(), v2.z(), v2.x() + segX, v2.y() + segY,
                        v2.z() + segZ);
    findParticleInBox(containerX2Y2Z2, particleVec, particleX2Y2Z2);
    reqX2Y2Z2[0] = boostWorld.isend(rankX2Y2Z2, mpiTag, particleX2Y2Z2);
    reqX2Y2Z2[1] = boostWorld.irecv(rankX2Y2Z2, mpiTag, rParticleX2Y2Z2);
  }
  // 6 surfaces
  if (rankX1 >= 0)
    boost::mpi::wait_all(reqX1, reqX1 + 2);
  if (rankX2 >= 0)
    boost::mpi::wait_all(reqX2, reqX2 + 2);
  if (rankY1 >= 0)
    boost::mpi::wait_all(reqY1, reqY1 + 2);
  if (rankY2 >= 0)
    boost::mpi::wait_all(reqY2, reqY2 + 2);
  if (rankZ1 >= 0)
    boost::mpi::wait_all(reqZ1, reqZ1 + 2);
  if (rankZ2 >= 0)
    boost::mpi::wait_all(reqZ2, reqZ2 + 2);
  // 12 edges
  if (rankX1Y1 >= 0)
    boost::mpi::wait_all(reqX1Y1, reqX1Y1 + 2);
  if (rankX1Y2 >= 0)
    boost::mpi::wait_all(reqX1Y2, reqX1Y2 + 2);
  if (rankX1Z1 >= 0)
    boost::mpi::wait_all(reqX1Z1, reqX1Z1 + 2);
  if (rankX1Z2 >= 0)
    boost::mpi::wait_all(reqX1Z2, reqX1Z2 + 2);
  if (rankX2Y1 >= 0)
    boost::mpi::wait_all(reqX2Y1, reqX2Y1 + 2);
  if (rankX2Y2 >= 0)
    boost::mpi::wait_all(reqX2Y2, reqX2Y2 + 2);
  if (rankX2Z1 >= 0)
    boost::mpi::wait_all(reqX2Z1, reqX2Z1 + 2);
  if (rankX2Z2 >= 0)
    boost::mpi::wait_all(reqX2Z2, reqX2Z2 + 2);
  if (rankY1Z1 >= 0)
    boost::mpi::wait_all(reqY1Z1, reqY1Z1 + 2);
  if (rankY1Z2 >= 0)
    boost::mpi::wait_all(reqY1Z2, reqY1Z2 + 2);
  if (rankY2Z1 >= 0)
    boost::mpi::wait_all(reqY2Z1, reqY2Z1 + 2);
  if (rankY2Z2 >= 0)
    boost::mpi::wait_all(reqY2Z2, reqY2Z2 + 2);
  // 8 vertices
  if (rankX1Y1Z1 >= 0)
    boost::mpi::wait_all(reqX1Y1Z1, reqX1Y1Z1 + 2);
  if (rankX1Y1Z2 >= 0)
    boost::mpi::wait_all(reqX1Y1Z2, reqX1Y1Z2 + 2);
  if (rankX1Y2Z1 >= 0)
    boost::mpi::wait_all(reqX1Y2Z1, reqX1Y2Z1 + 2);
  if (rankX1Y2Z2 >= 0)
    boost::mpi::wait_all(reqX1Y2Z2, reqX1Y2Z2 + 2);
  if (rankX2Y1Z1 >= 0)
    boost::mpi::wait_all(reqX2Y1Z1, reqX2Y1Z1 + 2);
  if (rankX2Y1Z2 >= 0)
    boost::mpi::wait_all(reqX2Y1Z2, reqX2Y1Z2 + 2);
  if (rankX2Y2Z1 >= 0)
    boost::mpi::wait_all(reqX2Y2Z1, reqX2Y2Z1 + 2);
  if (rankX2Y2Z2 >= 0)
    boost::mpi::wait_all(reqX2Y2Z2, reqX2Y2Z2 + 2);

  // delete outgoing particles
  removeParticleOutBox();

  // add incoming particles
  recvParticleVec.clear(); // new use of recvParticleVec
  // 6 surfaces
  if (rankX1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX1.begin(),
                           rParticleX1.end());
  if (rankX2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX2.begin(),
                           rParticleX2.end());
  if (rankY1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleY1.begin(),
                           rParticleY1.end());
  if (rankY2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleY2.begin(),
                           rParticleY2.end());
  if (rankZ1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleZ1.begin(),
                           rParticleZ1.end());
  if (rankZ2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleZ2.begin(),
                           rParticleZ2.end());
  // 12 edges
  if (rankX1Y1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1.begin(),
                           rParticleX1Y1.end());
  if (rankX1Y2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2.begin(),
                           rParticleX1Y2.end());
  if (rankX1Z1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX1Z1.begin(),
                           rParticleX1Z1.end());
  if (rankX1Z2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX1Z2.begin(),
                           rParticleX1Z2.end());
  if (rankX2Y1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1.begin(),
                           rParticleX2Y1.end());
  if (rankX2Y2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2.begin(),
                           rParticleX2Y2.end());
  if (rankX2Z1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX2Z1.begin(),
                           rParticleX2Z1.end());
  if (rankX2Z2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX2Z2.begin(),
                           rParticleX2Z2.end());
  if (rankY1Z1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleY1Z1.begin(),
                           rParticleY1Z1.end());
  if (rankY1Z2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleY1Z2.begin(),
                           rParticleY1Z2.end());
  if (rankY2Z1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleY2Z1.begin(),
                           rParticleY2Z1.end());
  if (rankY2Z2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleY2Z2.begin(),
                           rParticleY2Z2.end());
  // 8 vertices
  if (rankX1Y1Z1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1Z1.begin(),
                           rParticleX1Y1Z1.end());
  if (rankX1Y1Z2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y1Z2.begin(),
                           rParticleX1Y1Z2.end());
  if (rankX1Y2Z1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2Z1.begin(),
                           rParticleX1Y2Z1.end());
  if (rankX1Y2Z2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX1Y2Z2.begin(),
                           rParticleX1Y2Z2.end());
  if (rankX2Y1Z1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1Z1.begin(),
                           rParticleX2Y1Z1.end());
  if (rankX2Y1Z2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y1Z2.begin(),
                           rParticleX2Y1Z2.end());
  if (rankX2Y2Z1 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2Z1.begin(),
                           rParticleX2Y2Z1.end());
  if (rankX2Y2Z2 >= 0)
    recvParticleVec.insert(recvParticleVec.end(), rParticleX2Y2Z2.begin(),
                           rParticleX2Y2Z2.end());

  particleVec.insert(particleVec.end(), recvParticleVec.begin(),
                     recvParticleVec.end());

  */
  /*
    if (recvParticleVec.size() > 0) {
    debugInf << "iter=" << std::setw(8) << iteration << " rank=" << std::setw(2)
    << mpiRank
    << "   added=";
    for (ParticlePArray::const_iterator it = recvParticleVec.begin(); it !=
    recvParticleVec.end(); ++it)
    debugInf << std::setw(3) << (*it)->getId();
    debugInf << " now " << particleVec.size() << ": ";
    for (ParticlePArray::const_iterator it = particleVec.begin(); it !=
    particleVec.end(); ++it)
    debugInf << std::setw(3) << (*it)->getId();
    debugInf << std::endl;
    }
  */

  /*
  // do not release memory of received particles because they are part of and
  // managed by particleVec
  // 6 surfaces
  rParticleX1.clear();
  rParticleX2.clear();
  rParticleY1.clear();
  rParticleY2.clear();
  rParticleZ1.clear();
  rParticleZ2.clear();
  // 12 edges
  rParticleX1Y1.clear();
  rParticleX1Y2.clear();
  rParticleX1Z1.clear();
  rParticleX1Z2.clear();
  rParticleX2Y1.clear();
  rParticleX2Y2.clear();
  rParticleX2Z1.clear();
  rParticleX2Z2.clear();
  rParticleY1Z1.clear();
  rParticleY1Z2.clear();
  rParticleY2Z1.clear();
  rParticleY2Z2.clear();
  // 8 vertices
  rParticleX1Y1Z1.clear();
  rParticleX1Y1Z2.clear();
  rParticleX1Y2Z1.clear();
  rParticleX1Y2Z2.clear();
  rParticleX2Y1Z1.clear();
  rParticleX2Y1Z2.clear();
  rParticleX2Y2Z1.clear();
  rParticleX2Y2Z2.clear();

  recvParticleVec.clear();
}
*/

void
Assembly::migratePeriParticle()
{
  Vec vspan = grid.getMaxCorner() - grid.getMinCorner();
  REAL segX = vspan.x() / d_mpiProcs.x();
  REAL segY = vspan.y() / d_mpiProcs.y();
  REAL segZ = vspan.z() / d_mpiProcs.z();
  Vec v1 = container.getMinCorner(); // v1, v2 in terms of process
  Vec v2 = container.getMaxCorner();

  // if a neighbor exists, transfer particles crossing the boundary in between.
  PeriParticlePArray periParticleX1, periParticleX2;
  PeriParticlePArray periParticleY1, periParticleY2;
  PeriParticlePArray periParticleZ1, periParticleZ2;
  PeriParticlePArray periParticleX1Y1, periParticleX1Y2, periParticleX1Z1,
    periParticleX1Z2;
  PeriParticlePArray periParticleX2Y1, periParticleX2Y2, periParticleX2Z1,
    periParticleX2Z2;
  PeriParticlePArray periParticleY1Z1, periParticleY1Z2, periParticleY2Z1,
    periParticleY2Z2;
  PeriParticlePArray periParticleX1Y1Z1, periParticleX1Y1Z2, periParticleX1Y2Z1,
    periParticleX1Y2Z2;
  PeriParticlePArray periParticleX2Y1Z1, periParticleX2Y1Z2, periParticleX2Y2Z1,
    periParticleX2Y2Z2;
  boost::mpi::request reqX1[2], reqX2[2];
  boost::mpi::request reqY1[2], reqY2[2];
  boost::mpi::request reqZ1[2], reqZ2[2];
  boost::mpi::request reqX1Y1[2], reqX1Y2[2], reqX1Z1[2], reqX1Z2[2];
  boost::mpi::request reqX2Y1[2], reqX2Y2[2], reqX2Z1[2], reqX2Z2[2];
  boost::mpi::request reqY1Z1[2], reqY1Z2[2], reqY2Z1[2], reqY2Z2[2];
  boost::mpi::request reqX1Y1Z1[2], reqX1Y1Z2[2], reqX1Y2Z1[2], reqX1Y2Z2[2];
  boost::mpi::request reqX2Y1Z1[2], reqX2Y1Z2[2], reqX2Y2Z1[2], reqX2Y2Z2[2];

  // 6 surfaces
  if (rankX1 >= 0) { // surface x1
    Box containerX1(v1.x() - segX, v1.y(), v1.z(), v1.x(), v2.y(), v2.z());
    findPeriParticleInBox(containerX1, periParticleVec, periParticleX1);
    reqX1[0] = boostWorld.isend(rankX1, mpiTag, periParticleX1);
    reqX1[1] = boostWorld.irecv(rankX1, mpiTag, rperiParticleX1);
  }
  if (rankX2 >= 0) { // surface x2
    Box containerX2(v2.x(), v1.y(), v1.z(), v2.x() + segX, v2.y(), v2.z());
    findPeriParticleInBox(containerX2, periParticleVec, periParticleX2);
    reqX2[0] = boostWorld.isend(rankX2, mpiTag, periParticleX2);
    reqX2[1] = boostWorld.irecv(rankX2, mpiTag, rperiParticleX2);
  }
  if (rankY1 >= 0) { // surface y1
    Box containerY1(v1.x(), v1.y() - segY, v1.z(), v2.x(), v1.y(), v2.z());
    findPeriParticleInBox(containerY1, periParticleVec, periParticleY1);
    reqY1[0] = boostWorld.isend(rankY1, mpiTag, periParticleY1);
    reqY1[1] = boostWorld.irecv(rankY1, mpiTag, rperiParticleY1);
  }
  if (rankY2 >= 0) { // surface y2
    Box containerY2(v1.x(), v2.y(), v1.z(), v2.x(), v2.y() + segY, v2.z());
    findPeriParticleInBox(containerY2, periParticleVec, periParticleY2);
    reqY2[0] = boostWorld.isend(rankY2, mpiTag, periParticleY2);
    reqY2[1] = boostWorld.irecv(rankY2, mpiTag, rperiParticleY2);
  }
  if (rankZ1 >= 0) { // surface z1
    Box containerZ1(v1.x(), v1.y(), v1.z() - segZ, v2.x(), v2.y(), v1.z());
    findPeriParticleInBox(containerZ1, periParticleVec, periParticleZ1);
    reqZ1[0] = boostWorld.isend(rankZ1, mpiTag, periParticleZ1);
    reqZ1[1] = boostWorld.irecv(rankZ1, mpiTag, rperiParticleZ1);
  }
  if (rankZ2 >= 0) { // surface z2
    Box containerZ2(v1.x(), v1.y(), v2.z(), v2.x(), v2.y(), v2.z() + segZ);
    findPeriParticleInBox(containerZ2, periParticleVec, periParticleZ2);
    reqZ2[0] = boostWorld.isend(rankZ2, mpiTag, periParticleZ2);
    reqZ2[1] = boostWorld.irecv(rankZ2, mpiTag, rperiParticleZ2);
  }
  // 12 edges
  if (rankX1Y1 >= 0) { // edge x1y1
    Box containerX1Y1(v1.x() - segX, v1.y() - segY, v1.z(), v1.x(), v1.y(),
                      v2.z());
    findPeriParticleInBox(containerX1Y1, periParticleVec, periParticleX1Y1);
    reqX1Y1[0] = boostWorld.isend(rankX1Y1, mpiTag, periParticleX1Y1);
    reqX1Y1[1] = boostWorld.irecv(rankX1Y1, mpiTag, rperiParticleX1Y1);
  }
  if (rankX1Y2 >= 0) { // edge x1y2
    Box containerX1Y2(v1.x() - segX, v2.y(), v1.z(), v1.x(), v2.y() + segY,
                      v2.z());
    findPeriParticleInBox(containerX1Y2, periParticleVec, periParticleX1Y2);
    reqX1Y2[0] = boostWorld.isend(rankX1Y2, mpiTag, periParticleX1Y2);
    reqX1Y2[1] = boostWorld.irecv(rankX1Y2, mpiTag, rperiParticleX1Y2);
  }
  if (rankX1Z1 >= 0) { // edge x1z1
    Box containerX1Z1(v1.x() - segX, v1.y(), v1.z() - segZ, v1.x(), v2.y(),
                      v1.z());
    findPeriParticleInBox(containerX1Z1, periParticleVec, periParticleX1Z1);
    reqX1Z1[0] = boostWorld.isend(rankX1Z1, mpiTag, periParticleX1Z1);
    reqX1Z1[1] = boostWorld.irecv(rankX1Z1, mpiTag, rperiParticleX1Z1);
  }
  if (rankX1Z2 >= 0) { // edge x1z2
    Box containerX1Z2(v1.x() - segX, v1.y(), v2.z(), v1.x(), v2.y(),
                      v2.z() + segZ);
    findPeriParticleInBox(containerX1Z2, periParticleVec, periParticleX1Z2);
    reqX1Z2[0] = boostWorld.isend(rankX1Z2, mpiTag, periParticleX1Z2);
    reqX1Z2[1] = boostWorld.irecv(rankX1Z2, mpiTag, rperiParticleX1Z2);
  }
  if (rankX2Y1 >= 0) { // edge x2y1
    Box containerX2Y1(v2.x(), v1.y() - segY, v1.z(), v2.x() + segX, v1.y(),
                      v2.z());
    findPeriParticleInBox(containerX2Y1, periParticleVec, periParticleX2Y1);
    reqX2Y1[0] = boostWorld.isend(rankX2Y1, mpiTag, periParticleX2Y1);
    reqX2Y1[1] = boostWorld.irecv(rankX2Y1, mpiTag, rperiParticleX2Y1);
  }
  if (rankX2Y2 >= 0) { // edge x2y2
    Box containerX2Y2(v2.x(), v2.y(), v1.z(), v2.x() + segX, v2.y() + segY,
                      v2.z());
    findPeriParticleInBox(containerX2Y2, periParticleVec, periParticleX2Y2);
    reqX2Y2[0] = boostWorld.isend(rankX2Y2, mpiTag, periParticleX2Y2);
    reqX2Y2[1] = boostWorld.irecv(rankX2Y2, mpiTag, rperiParticleX2Y2);
  }
  if (rankX2Z1 >= 0) { // edge x2z1
    Box containerX2Z1(v2.x(), v1.y(), v1.z() - segZ, v2.x() + segX, v2.y(),
                      v1.z());
    findPeriParticleInBox(containerX2Z1, periParticleVec, periParticleX2Z1);
    reqX2Z1[0] = boostWorld.isend(rankX2Z1, mpiTag, periParticleX2Z1);
    reqX2Z1[1] = boostWorld.irecv(rankX2Z1, mpiTag, rperiParticleX2Z1);
  }
  if (rankX2Z2 >= 0) { // edge x2z2
    Box containerX2Z2(v2.x(), v1.y(), v2.z(), v2.x() + segX, v2.y(),
                      v2.z() + segZ);
    findPeriParticleInBox(containerX2Z2, periParticleVec, periParticleX2Z2);
    reqX2Z2[0] = boostWorld.isend(rankX2Z2, mpiTag, periParticleX2Z2);
    reqX2Z2[1] = boostWorld.irecv(rankX2Z2, mpiTag, rperiParticleX2Z2);
  }
  if (rankY1Z1 >= 0) { // edge y1z1
    Box containerY1Z1(v1.x(), v1.y() - segY, v1.z() - segZ, v2.x(), v1.y(),
                      v1.z());
    findPeriParticleInBox(containerY1Z1, periParticleVec, periParticleY1Z1);
    reqY1Z1[0] = boostWorld.isend(rankY1Z1, mpiTag, periParticleY1Z1);
    reqY1Z1[1] = boostWorld.irecv(rankY1Z1, mpiTag, rperiParticleY1Z1);
  }
  if (rankY1Z2 >= 0) { // edge y1z2
    Box containerY1Z2(v1.x(), v1.y() - segY, v2.z(), v2.x(), v1.y(),
                      v2.z() + segZ);
    findPeriParticleInBox(containerY1Z2, periParticleVec, periParticleY1Z2);
    reqY1Z2[0] = boostWorld.isend(rankY1Z2, mpiTag, periParticleY1Z2);
    reqY1Z2[1] = boostWorld.irecv(rankY1Z2, mpiTag, rperiParticleY1Z2);
  }
  if (rankY2Z1 >= 0) { // edge y2z1
    Box containerY2Z1(v1.x(), v2.y(), v1.z() - segZ, v2.x(), v2.y() + segY,
                      v1.z());
    findPeriParticleInBox(containerY2Z1, periParticleVec, periParticleY2Z1);
    reqY2Z1[0] = boostWorld.isend(rankY2Z1, mpiTag, periParticleY2Z1);
    reqY2Z1[1] = boostWorld.irecv(rankY2Z1, mpiTag, rperiParticleY2Z1);
  }
  if (rankY2Z2 >= 0) { // edge y2z2
    Box containerY2Z2(v1.x(), v2.y(), v2.z(), v2.x(), v2.y() + segY,
                      v2.z() + segZ);
    findPeriParticleInBox(containerY2Z2, periParticleVec, periParticleY2Z2);
    reqY2Z2[0] = boostWorld.isend(rankY2Z2, mpiTag, periParticleY2Z2);
    reqY2Z2[1] = boostWorld.irecv(rankY2Z2, mpiTag, rperiParticleY2Z2);
  }
  // 8 vertices
  if (rankX1Y1Z1 >= 0) { // edge x1y1z1
    Box containerX1Y1Z1(v1.x() - segX, v1.y() - segY, v1.z() - segZ, v1.x(),
                        v1.y(), v1.z());
    findPeriParticleInBox(containerX1Y1Z1, periParticleVec, periParticleX1Y1Z1);
    reqX1Y1Z1[0] = boostWorld.isend(rankX1Y1Z1, mpiTag, periParticleX1Y1Z1);
    reqX1Y1Z1[1] = boostWorld.irecv(rankX1Y1Z1, mpiTag, rperiParticleX1Y1Z1);
  }
  if (rankX1Y1Z2 >= 0) { // edge x1y1z2
    Box containerX1Y1Z2(v1.x() - segX, v1.y() - segY, v2.z(), v1.x(), v1.y(),
                        v2.z() + segZ);
    findPeriParticleInBox(containerX1Y1Z2, periParticleVec, periParticleX1Y1Z2);
    reqX1Y1Z2[0] = boostWorld.isend(rankX1Y1Z2, mpiTag, periParticleX1Y1Z2);
    reqX1Y1Z2[1] = boostWorld.irecv(rankX1Y1Z2, mpiTag, rperiParticleX1Y1Z2);
  }
  if (rankX1Y2Z1 >= 0) { // edge x1y2z1
    Box containerX1Y2Z1(v1.x() - segX, v2.y(), v1.z() - segZ, v1.x(),
                        v2.y() + segY, v1.z());
    findPeriParticleInBox(containerX1Y2Z1, periParticleVec, periParticleX1Y2Z1);
    reqX1Y2Z1[0] = boostWorld.isend(rankX1Y2Z1, mpiTag, periParticleX1Y2Z1);
    reqX1Y2Z1[1] = boostWorld.irecv(rankX1Y2Z1, mpiTag, rperiParticleX1Y2Z1);
  }
  if (rankX1Y2Z2 >= 0) { // edge x1y2z2
    Box containerX1Y2Z2(v1.x() - segX, v2.y(), v2.z(), v1.x(), v2.y() + segY,
                        v2.z() + segZ);
    findPeriParticleInBox(containerX1Y2Z2, periParticleVec, periParticleX1Y2Z2);
    reqX1Y2Z2[0] = boostWorld.isend(rankX1Y2Z2, mpiTag, periParticleX1Y2Z2);
    reqX1Y2Z2[1] = boostWorld.irecv(rankX1Y2Z2, mpiTag, rperiParticleX1Y2Z2);
  }
  if (rankX2Y1Z1 >= 0) { // edge x2y1z1
    Box containerX2Y1Z1(v2.x(), v1.y() - segY, v1.z() - segZ, v2.x() + segX,
                        v1.y(), v1.z());
    findPeriParticleInBox(containerX2Y1Z1, periParticleVec, periParticleX2Y1Z1);
    reqX2Y1Z1[0] = boostWorld.isend(rankX2Y1Z1, mpiTag, periParticleX2Y1Z1);
    reqX2Y1Z1[1] = boostWorld.irecv(rankX2Y1Z1, mpiTag, rperiParticleX2Y1Z1);
  }
  if (rankX2Y1Z2 >= 0) { // edge x2y1z2
    Box containerX2Y1Z2(v2.x(), v1.y() - segY, v2.z(), v2.x() + segX, v1.y(),
                        v2.z() + segZ);
    findPeriParticleInBox(containerX2Y1Z2, periParticleVec, periParticleX2Y1Z2);
    reqX2Y1Z2[0] = boostWorld.isend(rankX2Y1Z2, mpiTag, periParticleX2Y1Z2);
    reqX2Y1Z2[1] = boostWorld.irecv(rankX2Y1Z2, mpiTag, rperiParticleX2Y1Z2);
  }
  if (rankX2Y2Z1 >= 0) { // edge x2y2z1
    Box containerX2Y2Z1(v2.x(), v2.y(), v1.z() - segZ, v2.x() + segX,
                        v2.y() + segY, v1.z());
    findPeriParticleInBox(containerX2Y2Z1, periParticleVec, periParticleX2Y2Z1);
    reqX2Y2Z1[0] = boostWorld.isend(rankX2Y2Z1, mpiTag, periParticleX2Y2Z1);
    reqX2Y2Z1[1] = boostWorld.irecv(rankX2Y2Z1, mpiTag, rperiParticleX2Y2Z1);
  }
  if (rankX2Y2Z2 >= 0) { // edge x2y2z2
    Box containerX2Y2Z2(v2.x(), v2.y(), v2.z(), v2.x() + segX, v2.y() + segY,
                        v2.z() + segZ);
    findPeriParticleInBox(containerX2Y2Z2, periParticleVec, periParticleX2Y2Z2);
    reqX2Y2Z2[0] = boostWorld.isend(rankX2Y2Z2, mpiTag, periParticleX2Y2Z2);
    reqX2Y2Z2[1] = boostWorld.irecv(rankX2Y2Z2, mpiTag, rperiParticleX2Y2Z2);
  }
  // 6 surfaces
  if (rankX1 >= 0)
    boost::mpi::wait_all(reqX1, reqX1 + 2);
  if (rankX2 >= 0)
    boost::mpi::wait_all(reqX2, reqX2 + 2);
  if (rankY1 >= 0)
    boost::mpi::wait_all(reqY1, reqY1 + 2);
  if (rankY2 >= 0)
    boost::mpi::wait_all(reqY2, reqY2 + 2);
  if (rankZ1 >= 0)
    boost::mpi::wait_all(reqZ1, reqZ1 + 2);
  if (rankZ2 >= 0)
    boost::mpi::wait_all(reqZ2, reqZ2 + 2);
  // 12 edges
  if (rankX1Y1 >= 0)
    boost::mpi::wait_all(reqX1Y1, reqX1Y1 + 2);
  if (rankX1Y2 >= 0)
    boost::mpi::wait_all(reqX1Y2, reqX1Y2 + 2);
  if (rankX1Z1 >= 0)
    boost::mpi::wait_all(reqX1Z1, reqX1Z1 + 2);
  if (rankX1Z2 >= 0)
    boost::mpi::wait_all(reqX1Z2, reqX1Z2 + 2);
  if (rankX2Y1 >= 0)
    boost::mpi::wait_all(reqX2Y1, reqX2Y1 + 2);
  if (rankX2Y2 >= 0)
    boost::mpi::wait_all(reqX2Y2, reqX2Y2 + 2);
  if (rankX2Z1 >= 0)
    boost::mpi::wait_all(reqX2Z1, reqX2Z1 + 2);
  if (rankX2Z2 >= 0)
    boost::mpi::wait_all(reqX2Z2, reqX2Z2 + 2);
  if (rankY1Z1 >= 0)
    boost::mpi::wait_all(reqY1Z1, reqY1Z1 + 2);
  if (rankY1Z2 >= 0)
    boost::mpi::wait_all(reqY1Z2, reqY1Z2 + 2);
  if (rankY2Z1 >= 0)
    boost::mpi::wait_all(reqY2Z1, reqY2Z1 + 2);
  if (rankY2Z2 >= 0)
    boost::mpi::wait_all(reqY2Z2, reqY2Z2 + 2);
  // 8 vertices
  if (rankX1Y1Z1 >= 0)
    boost::mpi::wait_all(reqX1Y1Z1, reqX1Y1Z1 + 2);
  if (rankX1Y1Z2 >= 0)
    boost::mpi::wait_all(reqX1Y1Z2, reqX1Y1Z2 + 2);
  if (rankX1Y2Z1 >= 0)
    boost::mpi::wait_all(reqX1Y2Z1, reqX1Y2Z1 + 2);
  if (rankX1Y2Z2 >= 0)
    boost::mpi::wait_all(reqX1Y2Z2, reqX1Y2Z2 + 2);
  if (rankX2Y1Z1 >= 0)
    boost::mpi::wait_all(reqX2Y1Z1, reqX2Y1Z1 + 2);
  if (rankX2Y1Z2 >= 0)
    boost::mpi::wait_all(reqX2Y1Z2, reqX2Y1Z2 + 2);
  if (rankX2Y2Z1 >= 0)
    boost::mpi::wait_all(reqX2Y2Z1, reqX2Y2Z1 + 2);
  if (rankX2Y2Z2 >= 0)
    boost::mpi::wait_all(reqX2Y2Z2, reqX2Y2Z2 + 2);

  // delete outgoing particles
  removePeriParticleOutBox();

  // add incoming particles
  recvPeriParticleVec.clear();
  // 6 surfaces
  if (rankX1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX1.begin(), rperiParticleX1.end());
  if (rankX2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX2.begin(), rperiParticleX2.end());
  if (rankY1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleY1.begin(), rperiParticleY1.end());
  if (rankY2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleY2.begin(), rperiParticleY2.end());
  if (rankZ1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleZ1.begin(), rperiParticleZ1.end());
  if (rankZ2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleZ2.begin(), rperiParticleZ2.end());
  // 12 edges
  if (rankX1Y1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX1Y1.begin(),
                               rperiParticleX1Y1.end());
  if (rankX1Y2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX1Y2.begin(),
                               rperiParticleX1Y2.end());
  if (rankX1Z1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX1Z1.begin(),
                               rperiParticleX1Z1.end());
  if (rankX1Z2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX1Z2.begin(),
                               rperiParticleX1Z2.end());
  if (rankX2Y1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX2Y1.begin(),
                               rperiParticleX2Y1.end());
  if (rankX2Y2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX2Y2.begin(),
                               rperiParticleX2Y2.end());
  if (rankX2Z1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX2Z1.begin(),
                               rperiParticleX2Z1.end());
  if (rankX2Z2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX2Z2.begin(),
                               rperiParticleX2Z2.end());
  if (rankY1Z1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleY1Z1.begin(),
                               rperiParticleY1Z1.end());
  if (rankY1Z2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleY1Z2.begin(),
                               rperiParticleY1Z2.end());
  if (rankY2Z1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleY2Z1.begin(),
                               rperiParticleY2Z1.end());
  if (rankY2Z2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleY2Z2.begin(),
                               rperiParticleY2Z2.end());
  // 8 vertices
  if (rankX1Y1Z1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX1Y1Z1.begin(),
                               rperiParticleX1Y1Z1.end());
  if (rankX1Y1Z2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX1Y1Z2.begin(),
                               rperiParticleX1Y1Z2.end());
  if (rankX1Y2Z1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX1Y2Z1.begin(),
                               rperiParticleX1Y2Z1.end());
  if (rankX1Y2Z2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX1Y2Z2.begin(),
                               rperiParticleX1Y2Z2.end());
  if (rankX2Y1Z1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX2Y1Z1.begin(),
                               rperiParticleX2Y1Z1.end());
  if (rankX2Y1Z2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX2Y1Z2.begin(),
                               rperiParticleX2Y1Z2.end());
  if (rankX2Y2Z1 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX2Y2Z1.begin(),
                               rperiParticleX2Y2Z1.end());
  if (rankX2Y2Z2 >= 0)
    recvPeriParticleVec.insert(recvPeriParticleVec.end(),
                               rperiParticleX2Y2Z2.begin(),
                               rperiParticleX2Y2Z2.end());

  for (auto& pit : recvPeriParticleVec)
    pit->constructMatrixMember();

  periParticleVec.insert(periParticleVec.end(), recvPeriParticleVec.begin(),
                         recvPeriParticleVec.end());

  /*
    if (recvParticleVec.size() > 0) {
    debugInf << "iter=" << std::setw(8) << iteration << " rank=" << std::setw(2)
    << mpiRank
    << "   added=";
    for (ParticlePArray::const_iterator it = recvParticleVec.begin(); it !=
    recvParticleVec.end(); ++it)
    debugInf << std::setw(3) << (*it)->getId();
    debugInf << " now " << particleVec.size() << ": ";
    for (ParticlePArray::const_iterator it = particleVec.begin(); it !=
    particleVec.end(); ++it)
    debugInf << std::setw(3) << (*it)->getId();
    debugInf << std::endl;
    }
  */

  // do not release memory of received particles because they are part of and
  // managed by particleVec
  // 6 surfaces
  rParticleX1.clear();
  rParticleX2.clear();
  rParticleY1.clear();
  rParticleY2.clear();
  rParticleZ1.clear();
  rParticleZ2.clear();
  // 12 edges
  rParticleX1Y1.clear();
  rParticleX1Y2.clear();
  rParticleX1Z1.clear();
  rParticleX1Z2.clear();
  rParticleX2Y1.clear();
  rParticleX2Y2.clear();
  rParticleX2Z1.clear();
  rParticleX2Z2.clear();
  rParticleY1Z1.clear();
  rParticleY1Z2.clear();
  rParticleY2Z1.clear();
  rParticleY2Z2.clear();
  // 8 vertices
  rParticleX1Y1Z1.clear();
  rParticleX1Y1Z2.clear();
  rParticleX1Y2Z1.clear();
  rParticleX1Y2Z2.clear();
  rParticleX2Y1Z1.clear();
  rParticleX2Y1Z2.clear();
  rParticleX2Y2Z1.clear();
  rParticleX2Y2Z2.clear();

  recvPeriParticleVec.clear();
}

void
Assembly::gatherParticle()
{
  // update allParticleVec: process 0 collects all updated particles from each
  // other process
  if (mpiRank != 0) { // each process except 0
    boostWorld.send(0, mpiTag, particleVec);
  } else { // process 0
    // allParticleVec is cleared before filling with new data
    releaseGatheredParticle();

    // duplicate particleVec so that it is not destroyed by allParticleVec in
    // next iteration,
    // otherwise it causes memory error.
    ParticlePArray dupParticleVec(particleVec.size());
    for (std::size_t i = 0; i < dupParticleVec.size(); ++i)
      dupParticleVec[i] = std::make_shared<Particle>(*particleVec[i]);

    // fill allParticleVec with dupParticleVec and received particles
    allParticleVec.insert(allParticleVec.end(), dupParticleVec.begin(),
                          dupParticleVec.end());

    ParticlePArray tmpParticleVec;
    long gatherRam = 0;
    for (int iRank = 1; iRank < mpiSize; ++iRank) {
      tmpParticleVec.clear(); // do not destroy particles!
      boostWorld.recv(iRank, mpiTag, tmpParticleVec);
      allParticleVec.insert(allParticleVec.end(), tmpParticleVec.begin(),
                            tmpParticleVec.end());
      gatherRam += tmpParticleVec.size();
    }
    // debugInf << "gather: particleNum = " << gatherRam <<  " particleRam = "
    // << gatherRam * sizeof(Particle) << std::endl;
  }
}

void
Assembly::releaseGatheredParticle()
{
  // clear allParticleVec, avoid long time memory footprint.
  /*
  for (ParticlePArray::iterator it = allParticleVec.begin(); it !=
  allParticleVec.end(); ++it)
    delete (*it);
  */
  allParticleVec.clear();
  ParticlePArray().swap(allParticleVec); // actual memory release
}

void
Assembly::gatherPeriParticle()
{
  // update allPeriParticleVec: process 0 collects all updated particles from
  // each other process
  for (auto& it : periParticleVec) {
    it->assignSigma(); // store Matrix sigma to each value
  }
  if (mpiRank != 0) { // each process except 0
    boostWorld.send(0, mpiTag, periParticleVec);
  } else { // process 0
    // allPeriParticleVec is cleared before filling with new data
    releaseGatheredPeriParticle();

    // duplicate PeriParticleVec so that it is not destroyed by
    // allPeriParticleVec in next iteration,
    // otherwise it causes memory error.
    PeriParticlePArray dupPeriParticleVec(periParticleVec.size());
    for (std::size_t i = 0; i < dupPeriParticleVec.size(); ++i) {
      dupPeriParticleVec[i] =
        std::make_shared<pd::PeriParticle>(*periParticleVec[i]);
      dupPeriParticleVec[i]->releaseBondVec();
    }

    // fill allParticleVec with dupParticleVec and received particles
    allPeriParticleVec.insert(allPeriParticleVec.end(),
                              dupPeriParticleVec.begin(),
                              dupPeriParticleVec.end());
    PeriParticlePArray tmpPeriParticleVec;
    long gatherRam = 0;
    for (int iRank = 1; iRank < mpiSize; ++iRank) {
      tmpPeriParticleVec.clear(); // do not destroy particles!
      boostWorld.recv(iRank, mpiTag, tmpPeriParticleVec);
      allPeriParticleVec.insert(allPeriParticleVec.end(),
                                tmpPeriParticleVec.begin(),
                                tmpPeriParticleVec.end());
      gatherRam += tmpPeriParticleVec.size();
    }
    // debugInf << "gather: particleNum = " << gatherRam <<  " particleRam = "
    // << gatherRam * sizeof(Particle) << std::endl;
  }
}

void
Assembly::releaseGatheredPeriParticle()
{
  // clear allPeriParticleVec, avoid long time memory footprint.
  // for (auto it = allPeriParticleVec.begin(); it != allPeriParticleVec.end();
  // ++it)
  //  delete (*it);
  allPeriParticleVec.clear();
  PeriParticlePArray().swap(allPeriParticleVec); // actual memory release
}

void
Assembly::gatherBdryContact()
{
  if (isBdryProcess()) {
    if (mpiRank != 0)
      boostWorld.send(0, mpiTag, boundaryVec);
  }

  if (mpiRank == 0) {
    mergeBoundaryVec.clear();
    BoundaryPArray().swap(mergeBoundaryVec); // actual memory release
    mergeBoundaryVec = boundaryVec;

    BoundaryPArray tmpBoundaryVec;
    for (unsigned long bdryProces : bdryProcess) {
      if (bdryProces != 0) {    // not root process
        tmpBoundaryVec.clear(); // do not destroy particles!
        boostWorld.recv(bdryProces, mpiTag, tmpBoundaryVec);
        // merge tmpBoundaryVec into mergeBoundaryVec
        assert(tmpBoundaryVec.size() == mergeBoundaryVec.size());
        for (std::size_t jt = 0; jt < tmpBoundaryVec.size(); ++jt)
          mergeBoundaryVec[jt]->getContactInfo().insert(
            mergeBoundaryVec[jt]->getContactInfo().end(),
            tmpBoundaryVec[jt]->getContactInfo().begin(),
            tmpBoundaryVec[jt]->getContactInfo().end());
      }
    }

    // must update after collecting all boundary contact info
    for (auto& it : mergeBoundaryVec)
      it->updateStatForce();
  }
}

void
Assembly::gatherEnergy()
{
  calcTransEnergy();
  calcRotatEnergy();
  calcKinetEnergy();
  calcGraviEnergy(allContainer.getMinCorner().z());
  calcMechaEnergy();
}

void
Assembly::closeProg(std::ofstream& ofs)
{
  ofs.close();
}

void
Assembly::getStartDimension(REAL& distX, REAL& distY, REAL& distZ)
{
  REAL x1, x2, y1, y2, z1, z2;
  // use boundaryVec
  for (BoundaryPArray::const_iterator it = boundaryVec.begin();
       it != boundaryVec.end(); ++it) {
    switch ((*it)->getId()) {
      case 1:
        x1 = (*it)->getPoint().x();
        break;
      case 2:
        x2 = (*it)->getPoint().x();
        break;
      case 3:
        y1 = (*it)->getPoint().y();
        break;
      case 4:
        y2 = (*it)->getPoint().y();
        break;
      case 5:
        z1 = (*it)->getPoint().z();
        break;
      case 6:
        z2 = (*it)->getPoint().z();
        break;
    }
  }
  distX = x2 - x1;
  distY = y2 - y1;
  distZ = z2 - z1;
}

void
Assembly::openCompressProg(std::ofstream& ofs, const std::string& str)
{
  ofs.open(str);
  if (!ofs) {
    debugInf << "stream error: openCompressProg" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  ofs << std::setw(OWID) << "iteration" << std::setw(OWID) << "traction_x1"
      << std::setw(OWID) << "traction_x2" << std::setw(OWID) << "traction_y1"
      << std::setw(OWID) << "traction_y2" << std::setw(OWID) << "traction_z1"
      << std::setw(OWID) << "traction_z2" << std::setw(OWID) << "mean_stress"

      << std::setw(OWID) << "bulk_volume" << std::setw(OWID) << "density"
      << std::setw(OWID) << "epsilon_x" << std::setw(OWID) << "epsilon_y"
      << std::setw(OWID) << "epsilon_z" << std::setw(OWID) << "epsilon_v"
      << std::setw(OWID) << "void_ratio" << std::setw(OWID) << "porosity"

      << std::setw(OWID) << "velocity_x1" << std::setw(OWID) << "velocity_x2"
      << std::setw(OWID) << "velocity_y1" << std::setw(OWID) << "velocity_y2"
      << std::setw(OWID) << "velocity_z1" << std::setw(OWID) << "velocity_z2"

      << std::setw(OWID) << "contact_x1" << std::setw(OWID) << "contact_x2"
      << std::setw(OWID) << "contact_y1" << std::setw(OWID) << "contact_y2"
      << std::setw(OWID) << "contact_z1" << std::setw(OWID) << "contact_z2"
      << std::setw(OWID) << "contact_inside"

      << std::setw(OWID) << "penetr_x1" << std::setw(OWID) << "penetr_x2"
      << std::setw(OWID) << "penetr_y1" << std::setw(OWID) << "penetr_y2"
      << std::setw(OWID) << "penetr_z1" << std::setw(OWID) << "penetr_z2"

      << std::setw(OWID) << "avgNormal" << std::setw(OWID) << "avgShear"
      << std::setw(OWID) << "avgPenetr"

      << std::setw(OWID) << "transEnergy" << std::setw(OWID) << "rotatEnergy"
      << std::setw(OWID) << "kinetEnergy" << std::setw(OWID) << "graviEnergy"
      << std::setw(OWID) << "mechaEnergy"

      << std::setw(OWID) << "vibra_est_dt" << std::setw(OWID) << "impact_est_dt"
      << std::setw(OWID) << "actual_dt"

      << std::endl;
}

void
Assembly::printCompressProg(std::ofstream& ofs, REAL distX, REAL distY,
                            REAL distZ)
{
  REAL x1, x2, y1, y2, z1, z2;
  for (BoundaryPArray::const_iterator it = mergeBoundaryVec.begin();
       it != mergeBoundaryVec.end(); ++it) {
    switch ((*it)->getId()) {
      case 1:
        x1 = (*it)->getPoint().x();
        break;
      case 2:
        x2 = (*it)->getPoint().x();
        break;
      case 3:
        y1 = (*it)->getPoint().y();
        break;
      case 4:
        y2 = (*it)->getPoint().y();
        break;
      case 5:
        z1 = (*it)->getPoint().z();
        break;
      case 6:
        z2 = (*it)->getPoint().z();
        break;
    }
  }
  REAL areaX = (y2 - y1) * (z2 - z1);
  REAL areaY = (z2 - z1) * (x2 - x1);
  REAL areaZ = (x2 - x1) * (y2 - y1);
  REAL bulkVolume = (x2 - x1) * (y2 - y1) * (z2 - z1);
  REAL voidRatio = bulkVolume / getParticleVolume() - 1;

  REAL var[6], vel[6];
  // normalForce
  for (std::size_t i = 0; i < 6; ++i) {
    var[i] = 0;
    vel[i] = 0;
  }
  for (BoundaryPArray::const_iterator it = mergeBoundaryVec.begin();
       it != mergeBoundaryVec.end(); ++it) {
    std::size_t id = (*it)->getId();
    Vec normal = (*it)->getNormalForce();
    Vec veloc = (*it)->getVeloc();
    switch (id) {
      case 1:
        var[0] = fabs(normal.x()) / areaX;
        vel[0] = veloc.x();
        break;
      case 2:
        var[1] = normal.x() / areaX;
        vel[1] = veloc.x();
        break;
      case 3:
        var[2] = fabs(normal.y()) / areaY;
        vel[2] = veloc.y();
        break;
      case 4:
        var[3] = normal.y() / areaY;
        vel[3] = veloc.y();
        break;
      case 5:
        var[4] = fabs(normal.z()) / areaZ;
        vel[4] = veloc.z();
        break;
      case 6:
        var[5] = normal.z() / areaZ;
        vel[5] = veloc.z();
        break;
    }
  }
  ofs << std::setw(OWID) << iteration;
  REAL avg = 0;
  for (double i : var) {
    ofs << std::setw(OWID) << i;
    avg += i;
  }
  ofs << std::setw(OWID) << avg / 6;

  // volume
  ofs << std::setw(OWID) << bulkVolume << std::setw(OWID)
      << getMass() / bulkVolume << std::setw(OWID) << 1 - (x2 - x1) / distX
      << std::setw(OWID) << 1 - (y2 - y1) / distY << std::setw(OWID)
      << 1 - (z2 - z1) / distZ << std::setw(OWID)
      << 3 - (x2 - x1) / distX - (y2 - y1) / distY - (z2 - z1) / distZ
      << std::setw(OWID) << voidRatio << std::setw(OWID)
      << voidRatio / (1 + voidRatio);

  // velocity
  for (double i : vel)
    ofs << std::setw(OWID) << i;

  // contactNum
  for (double& i : var)
    i = 0;
  for (BoundaryPArray::const_iterator it = mergeBoundaryVec.begin();
       it != mergeBoundaryVec.end(); ++it) {
    std::size_t id = (*it)->getId();
    var[id - 1] = (*it)->getContactNum();
  }
  for (double i : var)
    ofs << std::setw(OWID) << static_cast<std::size_t>(i);
  ofs << std::setw(OWID) << allContactNum;

  // avgPenetr
  for (double& i : var)
    i = 0;
  for (BoundaryPArray::const_iterator it = mergeBoundaryVec.begin();
       it != mergeBoundaryVec.end(); ++it) {
    std::size_t id = (*it)->getId();
    var[id - 1] = (*it)->getAvgPenetr();
  }
  for (double i : var)
    ofs << std::setw(OWID) << i;

  // average data
  ofs << std::setw(OWID) << avgNormal << std::setw(OWID) << avgShear
      << std::setw(OWID) << avgPenetr;

  // energy
  ofs << std::setw(OWID) << transEnergy << std::setw(OWID) << rotatEnergy
      << std::setw(OWID) << kinetEnergy << std::setw(OWID) << graviEnergy
      << std::setw(OWID) << mechaEnergy;

  // time
  ofs << std::setw(OWID) << vibraTimeStep << std::setw(OWID) << impactTimeStep
      << std::setw(OWID) << timeStep;

  ofs << std::endl;
}

void
Assembly::openPeriProgress(std::ofstream& ofs, const std::string& str)
{
  ofs.open(str);
  if (!ofs) {
    debugInf << "stream error: openPeriProgress" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);
  ofs << "Title = \"Particle Information\"" << std::endl;
  ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\" \"Vy\" "
         "\"Vz\" \"KE\" \"P\" \"Mises\" \"Volume\" \"horizonSize\""
      << std::endl;
}

void
Assembly::printPeriProgress(std::ofstream& ofs, const int iframe) const
{
  ofs << "ZONE T =\" " << iframe << "-th Load Step\" " << std::endl;
  // Output the coordinates and the array information
  REAL pressure, vonMisesStress;
  REAL sigma11, sigma12; // sigma13;
  REAL /*sigma21,*/ sigma22, sigma23;
  REAL sigma31, /*sigma32,*/ sigma33;
  for (const auto& pt : allPeriParticleVec) {
    //        if( (*pt)->getInitPosition().x()>
    //        0.5*(util::getParam<REAL>("Xmin")+util::getParam<REAL>("Xmax"))
    //        )
    //        continue;
    sigma11 = pt->getSigma11();
    sigma12 = pt->getSigma12();
    // sigma13=(*pt)->getSigma13();
    // sigma21=(*pt)->getSigma21();
    sigma22 = pt->getSigma22();
    sigma23 = pt->getSigma23();
    sigma31 = pt->getSigma31();
    // sigma32=(*pt)->getSigma32();
    sigma33 = pt->getSigma33();
    pressure = sigma11 + sigma22 + sigma33;
    vonMisesStress =
      sqrt(0.5 * ((sigma11 - sigma22) * (sigma11 - sigma22) +
                  (sigma22 - sigma33) * (sigma22 - sigma33) +
                  (sigma11 - sigma33) * (sigma11 - sigma33)) +
           3 * (sigma12 * sigma12 + sigma23 * sigma23 + sigma31 * sigma31));
    ofs << std::setw(20)
        << pt->getInitPosition().x() + pt->getDisplacement().x()
        << std::setw(20)
        << pt->getInitPosition().y() + pt->getDisplacement().y()
        << std::setw(20)
        << pt->getInitPosition().z() + pt->getDisplacement().z()
        << std::setw(20) << pt->getDisplacement().x() << std::setw(20)
        << pt->getDisplacement().y() << std::setw(20)
        << pt->getDisplacement().z() << std::setw(20) << pt->getVelocity().x()
        << std::setw(20) << pt->getVelocity().y() << std::setw(20)
        << pt->getVelocity().z() << std::setw(20) << vfabs(pt->getVelocity())
        << std::setw(20) << pressure << std::setw(20) << vonMisesStress
        << std::setw(20) << pt->getParticleVolume() << std::setw(20)
        << pt->getHorizonSize() << std::endl;
    ofs.flush();
  }
}

void
Assembly::printPeriProgressHalf(std::ofstream& ofs, const int iframe) const
{
  ofs << "ZONE T =\" " << iframe << "-th Load Step\" " << std::endl;
  // Output the coordinates and the array information
  REAL pressure, vonMisesStress;
  REAL sigma11, sigma12; // sigma13;
  REAL /*sigma21,*/ sigma22, sigma23;
  REAL sigma31, /*sigma32,*/ sigma33;
  for (const auto& pt : allPeriParticleVec) {
    if (pt->getInitPosition().x() >
        0.5 * (util::getParam<REAL>("Xmin") + util::getParam<REAL>("Xmax")))
      continue;
    sigma11 = pt->getSigma11();
    sigma12 = pt->getSigma12();
    // sigma13=(*pt)->getSigma13();
    // sigma21=(*pt)->getSigma21();
    sigma22 = pt->getSigma22();
    sigma23 = pt->getSigma23();
    sigma31 = pt->getSigma31();
    // sigma32=(*pt)->getSigma32();
    sigma33 = pt->getSigma33();
    pressure = sigma11 + sigma22 + sigma33;
    vonMisesStress =
      sqrt(0.5 * ((sigma11 - sigma22) * (sigma11 - sigma22) +
                  (sigma22 - sigma33) * (sigma22 - sigma33) +
                  (sigma11 - sigma33) * (sigma11 - sigma33)) +
           3 * (sigma12 * sigma12 + sigma23 * sigma23 + sigma31 * sigma31));
    ofs << std::setw(20)
        << pt->getInitPosition().x() + pt->getDisplacement().x()
        << std::setw(20)
        << pt->getInitPosition().y() + pt->getDisplacement().y()
        << std::setw(20)
        << pt->getInitPosition().z() + pt->getDisplacement().z()
        << std::setw(20) << pt->getDisplacement().x() << std::setw(20)
        << pt->getDisplacement().y() << std::setw(20)
        << pt->getDisplacement().z() << std::setw(20) << pt->getVelocity().x()
        << std::setw(20) << pt->getVelocity().y() << std::setw(20)
        << pt->getVelocity().z() << std::setw(20) << vfabs(pt->getVelocity())
        << std::setw(20) << pressure << std::setw(20) << vonMisesStress
        << std::setw(20) << pt->getParticleVolume() << std::setw(20)
        << pt->getHorizonSize() << std::endl;
    ofs.flush();
  }
}

void
Assembly::openParticleProg(std::ofstream& ofs, const std::string& str)
{
  ofs.open(str);
  if (!ofs) {
    debugInf << "stream error: openParticleProg" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  ofs << std::setw(OWID) << "iteration" << std::setw(OWID) << "accruedTime"
      << std::setw(OWID) << "penalFx" << std::setw(OWID) << "penalFy"
      << std::setw(OWID) << "penalFz" << std::setw(OWID) << "pressureFx"
      << std::setw(OWID) << "pressureFy" << std::setw(OWID) << "pressureFz"
      << std::setw(OWID) << "penalMx" << std::setw(OWID) << "penalMy"
      << std::setw(OWID) << "penalMz" << std::setw(OWID) << "pressureMx"
      << std::setw(OWID) << "pressureMy" << std::setw(OWID) << "pressureMz"
      << std::setw(OWID) << "accelX" << std::setw(OWID) << "accelY"
      << std::setw(OWID) << "accelZ" << std::setw(OWID) << "velocX"
      << std::setw(OWID) << "velocY" << std::setw(OWID) << "velocZ"
      << std::endl;
}

void
Assembly::updateBoundary(REAL sigma, std::string type, REAL sigmaX, REAL sigmaY)
{
  if (mpiRank == 0) {
    REAL x1, x2, y1, y2, z1, z2;
    for (BoundaryPArray::const_iterator it = mergeBoundaryVec.begin();
         it != mergeBoundaryVec.end(); ++it) {
      switch ((*it)->getId()) {
        case 1:
          x1 = (*it)->getPoint().x();
          break;
        case 2:
          x2 = (*it)->getPoint().x();
          break;
        case 3:
          y1 = (*it)->getPoint().y();
          break;
        case 4:
          y2 = (*it)->getPoint().y();
          break;
        case 5:
          z1 = (*it)->getPoint().z();
          break;
        case 6:
          z2 = (*it)->getPoint().z();
          break;
      }
    }
    REAL areaX = (y2 - y1) * (z2 - z1);
    REAL areaY = (z2 - z1) * (x2 - x1);
    REAL areaZ = (x2 - x1) * (y2 - y1);

    if (type.compare("isotropic") == 0) {
      for (auto& it : mergeBoundaryVec)
        it->updateIsotropic(sigma, areaX, areaY, areaZ);
    } else if (type.compare("odometer") == 0) {
      for (auto& it : mergeBoundaryVec)
        it->updateOdometer(sigma, areaX, areaY, areaZ);
    } else if (type.compare("triaxial") == 0) {
      for (auto& it : mergeBoundaryVec)
        it->updateTriaxial(sigma, areaX, areaY, areaZ);
    } else if (type.compare("plnstrn") == 0) {
      for (auto& it : mergeBoundaryVec)
        it->updatePlaneStrain(sigma, areaX, areaY, areaZ);
    } else if (type.compare("trueTriaxial") == 0) {
      for (auto& it : mergeBoundaryVec)
        it->updateTrueTriaxial(sigma, areaX, areaY, areaZ, sigmaX, sigmaY);
    }

    // update boundaryVec from mergeBoundaryVec and remove contactInfo to reduce
    // MPI transmission
    boundaryVec = mergeBoundaryVec;
    for (auto& it : boundaryVec)
      it->clearContactInfo();

    // update allContainer
    for (BoundaryPArray::const_iterator it = boundaryVec.begin();
         it != boundaryVec.end(); ++it) {
      switch ((*it)->getId()) {
        case 1:
          x1 = (*it)->getPoint().x();
          break;
        case 2:
          x2 = (*it)->getPoint().x();
          break;
        case 3:
          y1 = (*it)->getPoint().y();
          break;
        case 4:
          y2 = (*it)->getPoint().y();
          break;
        case 5:
          z1 = (*it)->getPoint().z();
          break;
        case 6:
          z2 = (*it)->getPoint().z();
          break;
      }
    }
    setContainer(Box(x1, y1, z1, x2, y2, z2));
  }

  broadcast(boostWorld, boundaryVec, 0);
}

void
Assembly::printContact(const std::string& str) const
{
  // There are two implementions of printContact
  // implementation 1: parallel IO, each process prints to a data file using a
  // shared pointer.
  //                   and use post-processing tool to remove redundant info.
  MPI_Status status;
  MPI_File contactFile;
  MPI_File_open(mpiWorld, str.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &contactFile);
  if (boostWorld.rank() == 0 && !contactFile) {
    std::cerr << "stream error: printContact" << std::endl;
    exit(-1);
  }

  std::stringstream inf;
  inf.setf(std::ios::scientific, std::ios::floatfield);

  for (const auto& it : contactVec)
    inf << std::setw(OWID) << it.getP1()->getId() << std::setw(OWID)
        << it.getP2()->getId() << std::setw(OWID) << it.getPoint1().x()
        << std::setw(OWID) << it.getPoint1().y() << std::setw(OWID)
        << it.getPoint1().z() << std::setw(OWID) << it.getPoint2().x()
        << std::setw(OWID) << it.getPoint2().y() << std::setw(OWID)
        << it.getPoint2().z() << std::setw(OWID) << it.getRadius1()
        << std::setw(OWID) << it.getRadius2() << std::setw(OWID)
        << it.getPenetration() << std::setw(OWID) << it.getTgtDisp()
        << std::setw(OWID) << it.getContactRadius() << std::setw(OWID)
        << it.getR0() << std::setw(OWID) << it.getE0() << std::setw(OWID)
        << it.getNormalForce() << std::setw(OWID) << it.getTgtForce()
        << std::setw(OWID) << (it.getPoint1().x() + it.getPoint2().x()) / 2
        << std::setw(OWID) << (it.getPoint1().y() + it.getPoint2().y()) / 2
        << std::setw(OWID) << (it.getPoint1().z() + it.getPoint2().z()) / 2
        << std::setw(OWID) << it.normalForceVec().x() << std::setw(OWID)
        << it.normalForceVec().y() << std::setw(OWID) << it.normalForceVec().z()
        << std::setw(OWID) << it.tgtForceVec().x() << std::setw(OWID)
        << it.tgtForceVec().y() << std::setw(OWID) << it.tgtForceVec().z()
        << std::setw(OWID) << it.getVibraTimeStep() << std::setw(OWID)
        << it.getImpactTimeStep() << std::endl;

  int length = (OWID * 28 + 1) * contactVec.size();
  // write a file at a location specified by a shared file pointer (blocking,
  // collective)
  // note MPI_File_write_shared is non-collective
  MPI_File_write_ordered(contactFile, const_cast<char*>(inf.str().c_str()),
                         length, MPI_CHAR, &status);
  MPI_File_close(&contactFile);

  // implementation 2: each process prints to an individual file.
  //                   use post-processing tool to merge files and remove
  //                   redundance.
  /*
  char csuf[10];
  combine(csuf, ".p", mpiRank, 5);
  strcat(str, csuf);

  std::ofstream ofs(str);
  if(!ofs) { debugInf << "stream error: printContact" << std::endl; exit(-1); }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  ofs << std::setw(OWID) << contactVec.size() << std::endl;
  ofs << std::setw(OWID) << "ptcl_1"
  << std::setw(OWID) << "ptcl_2"
  << std::setw(OWID) << "point1_x"
  << std::setw(OWID) << "point1_y"
  << std::setw(OWID) << "point1_z"
  << std::setw(OWID) << "point2_x"
  << std::setw(OWID) << "point2_y"
  << std::setw(OWID) << "point2_z"
  << std::setw(OWID) << "radius_1"
  << std::setw(OWID) << "radius_2"
  << std::setw(OWID) << "penetration"
  << std::setw(OWID) << "tangt_disp"
  << std::setw(OWID) << "contact_radius"
  << std::setw(OWID) << "R0"
  << std::setw(OWID) << "E0"
  << std::setw(OWID) << "normal_force"
  << std::setw(OWID) << "tangt_force"
  << std::setw(OWID) << "contact_x"
  << std::setw(OWID) << "contact_y"
  << std::setw(OWID) << "contact_z"
  << std::setw(OWID) << "normal_x"
  << std::setw(OWID) << "normal_y"
  << std::setw(OWID) << "normal_z"
  << std::setw(OWID) << "tangt_x"
  << std::setw(OWID) << "tangt_y"
  << std::setw(OWID) << "tangt_z"
  << std::setw(OWID) << "vibra_t_step"
  << std::setw(OWID) << "impact_t_step"
  << std::endl;

  ContactArray::const_iterator it;
  for (it = contactVec.begin(); it != contactVec.end(); ++it)
    ofs << std::setw(OWID) << it->getP1()->getId()
    << std::setw(OWID) << it->getP2()->getId()
    << std::setw(OWID) << it->getPoint1().x()
    << std::setw(OWID) << it->getPoint1().y()
    << std::setw(OWID) << it->getPoint1().z()
    << std::setw(OWID) << it->getPoint2().x()
    << std::setw(OWID) << it->getPoint2().y()
    << std::setw(OWID) << it->getPoint2().z()
    << std::setw(OWID) << it->getRadius1()
    << std::setw(OWID) << it->getRadius2()
    << std::setw(OWID) << it->getPenetration()
    << std::setw(OWID) << it->getTgtDisp()
    << std::setw(OWID) << it->getContactRadius()
    << std::setw(OWID) << it->getR0()
    << std::setw(OWID) << it->getE0()
    << std::setw(OWID) << it->getNormalForce()
    << std::setw(OWID) << it->getTgtForce()
    << std::setw(OWID) << ( it->getPoint1().x() + it->getPoint2().x() )/2
    << std::setw(OWID) << ( it->getPoint1().y() + it->getPoint2().y() )/2
    << std::setw(OWID) << ( it->getPoint1().z() + it->getPoint2().z() )/2
    << std::setw(OWID) << it->normalForceVec().x()
    << std::setw(OWID) << it->normalForceVec().y()
    << std::setw(OWID) << it->normalForceVec().z()
    << std::setw(OWID) << it->tgtForceVec().x()
    << std::setw(OWID) << it->tgtForceVec().y()
    << std::setw(OWID) << it->tgtForceVec().z()
    << std::setw(OWID) << it->getVibraTimeStep()
    << std::setw(OWID) << it->getImpactTimeStep()
    << std::endl;
  ofs.close();
  */
}

void
Assembly::calcTransEnergy()
{
  REAL pEngy = 0;
  ParticlePArray::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it) {
    if ((*it)->getType() == 0)
      pEngy += (*it)->getTransEnergy();
  }
  MPI_Reduce(&pEngy, &transEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, mpiWorld);
}

void
Assembly::calcRotatEnergy()
{
  REAL pEngy = 0;
  ParticlePArray::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it) {
    if ((*it)->getType() == 0)
      pEngy += (*it)->getRotatEnergy();
  }
  MPI_Reduce(&pEngy, &rotatEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, mpiWorld);
}

void
Assembly::calcKinetEnergy()
{
  REAL pEngy = 0;
  ParticlePArray::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it) {
    if ((*it)->getType() == 0)
      pEngy += (*it)->getKinetEnergy();
  }
  MPI_Reduce(&pEngy, &kinetEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, mpiWorld);
}

void
Assembly::calcGraviEnergy(REAL ref)
{
  REAL pEngy = 0;
  ParticlePArray::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it) {
    if ((*it)->getType() == 0)
      pEngy += (*it)->getPotenEnergy(ref);
  }
  MPI_Reduce(&pEngy, &graviEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, mpiWorld);
}

void
Assembly::calcMechaEnergy()
{
  mechaEnergy = kinetEnergy + graviEnergy;
}

REAL
Assembly::getMass() const
{
  REAL var = 0;
  for (const auto& it : allParticleVec)
    var += it->getMass();
  return var;
}

REAL
Assembly::getParticleVolume() const
{
  REAL var = 0;
  for (const auto& it : allParticleVec)
    if (it->getType() == 0)
      var += it->getVolume();
  return var;
}

REAL
Assembly::getAvgTransVelocity() const
{
  REAL avgv = 0;
  std::size_t count = 0;
  ParticlePArray::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it)
    if ((*it)->getType() == 0) {
      avgv += vfabs((*it)->currentVel());
      ++count;
    }
  return avgv /= count;
}

REAL
Assembly::getAvgRotatVelocity() const
{
  REAL avgv = 0;
  std::size_t count = 0;
  ParticlePArray::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it)
    if ((*it)->getType() == 0) {
      avgv += vfabs((*it)->currentOmega());
      ++count;
    }
  return avgv /= count;
}

REAL
Assembly::getAvgForce() const
{
  REAL avgv = 0;
  std::size_t count = 0;
  ParticlePArray::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it)
    if ((*it)->getType() == 0) {
      avgv += vfabs((*it)->getForce());
      ++count;
    }
  return avgv / count;
}

REAL
Assembly::getAvgMoment() const
{
  REAL avgv = 0;
  std::size_t count = 0;
  ParticlePArray::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it)
    if ((*it)->getType() == 0) {
      avgv += vfabs((*it)->getMoment());
      ++count;
    }
  return avgv /= count;
}

void
Assembly::buildBoundary(std::size_t boundaryNum,
                        const std::string& boundaryFile)
{
  std::ofstream ofs(boundaryFile);
  if (!ofs) {
    debugInf << "stream error: buildBoundary" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);

  Vec v1 = allContainer.getMinCorner();
  Vec v2 = allContainer.getMaxCorner();
  Vec v0 = allContainer.getCenter();
  REAL x1 = v1.x();
  REAL y1 = v1.y();
  REAL z1 = v1.z();
  REAL x2 = v2.x();
  REAL y2 = v2.y();
  REAL z2 = v2.z();
  REAL x0 = v0.x();
  REAL y0 = v0.y();
  REAL z0 = v0.z();

  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z1
      << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z2
      << std::endl
      << std::endl
      << std::setw(OWID) << boundaryNum << std::endl
      << std::endl;

  if (boundaryNum == 1) { // only a bottom boundary, i.e., boundary 5
    ofs << std::setw(OWID) << 1 << std::setw(OWID) << 0 << std::endl

        << std::setw(OWID) << 5 << std::setw(OWID) << 0 << std::setw(OWID) << 0
        << std::setw(OWID) << -1 << std::setw(OWID) << x0 << std::setw(OWID)
        << y0 << std::setw(OWID) << z1 << std::endl
        << std::endl;

  } else if (boundaryNum == 5) { // no top boundary, i.e., no boundary 6
    // boundary 1
    ofs << std::setw(OWID) << 1 << std::setw(OWID) << 1 << std::endl

        << std::setw(OWID) << 1 << std::setw(OWID) << -1 << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::setw(OWID) << x1 << std::setw(OWID)
        << y0 << std::setw(OWID) << z0 << std::endl

        << std::setw(OWID) << " " << std::setw(OWID) << 0 << std::setw(OWID)
        << 0 << std::setw(OWID) << 1 << std::setw(OWID) << x0 << std::setw(OWID)
        << y0 << std::setw(OWID) << z2 << std::endl
        << std::endl

        // boundary 2
        << std::setw(OWID) << 1 << std::setw(OWID) << 1 << std::endl

        << std::setw(OWID) << 2 << std::setw(OWID) << 1 << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::setw(OWID) << x2 << std::setw(OWID)
        << y0 << std::setw(OWID) << z0 << std::endl

        << std::setw(OWID) << " " << std::setw(OWID) << 0 << std::setw(OWID)
        << 0 << std::setw(OWID) << 1 << std::setw(OWID) << x0 << std::setw(OWID)
        << y0 << std::setw(OWID) << z2 << std::endl
        << std::endl

        // boundary 3
        << std::setw(OWID) << 1 << std::setw(OWID) << 1 << std::endl

        << std::setw(OWID) << 3 << std::setw(OWID) << 0 << std::setw(OWID) << -1
        << std::setw(OWID) << 0 << std::setw(OWID) << x0 << std::setw(OWID)
        << y1 << std::setw(OWID) << z0 << std::endl

        << std::setw(OWID) << " " << std::setw(OWID) << 0 << std::setw(OWID)
        << 0 << std::setw(OWID) << 1 << std::setw(OWID) << x0 << std::setw(OWID)
        << y0 << std::setw(OWID) << z2 << std::endl
        << std::endl

        // boundary 4
        << std::setw(OWID) << 1 << std::setw(OWID) << 1 << std::endl

        << std::setw(OWID) << 4 << std::setw(OWID) << 0 << std::setw(OWID) << 1
        << std::setw(OWID) << 0 << std::setw(OWID) << x0 << std::setw(OWID)
        << y2 << std::setw(OWID) << z0 << std::endl

        << std::setw(OWID) << " " << std::setw(OWID) << 0 << std::setw(OWID)
        << 0 << std::setw(OWID) << 1 << std::setw(OWID) << x0 << std::setw(OWID)
        << y0 << std::setw(OWID) << z2 << std::endl
        << std::endl

        // boundary 5
        << std::setw(OWID) << 1 << std::setw(OWID) << 0 << std::endl

        << std::setw(OWID) << 5 << std::setw(OWID) << 0 << std::setw(OWID) << 0
        << std::setw(OWID) << -1 << std::setw(OWID) << x0 << std::setw(OWID)
        << y0 << std::setw(OWID) << z1 << std::endl
        << std::endl;
  } else if (boundaryNum == 6) { // all 6 boundaries
    // boundary 1
    ofs << std::setw(OWID) << 1 << std::setw(OWID) << 0 << std::endl

        << std::setw(OWID) << 1 << std::setw(OWID) << -1 << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::setw(OWID) << x1 << std::setw(OWID)
        << y0 << std::setw(OWID) << z0 << std::endl
        << std::endl

        // boundary 2
        << std::setw(OWID) << 1 << std::setw(OWID) << 0 << std::endl

        << std::setw(OWID) << 2 << std::setw(OWID) << 1 << std::setw(OWID) << 0
        << std::setw(OWID) << 0 << std::setw(OWID) << x2 << std::setw(OWID)
        << y0 << std::setw(OWID) << z0 << std::endl
        << std::endl

        // boundary 3
        << std::setw(OWID) << 1 << std::setw(OWID) << 0 << std::endl

        << std::setw(OWID) << 3 << std::setw(OWID) << 0 << std::setw(OWID) << -1
        << std::setw(OWID) << 0 << std::setw(OWID) << x0 << std::setw(OWID)
        << y1 << std::setw(OWID) << z0 << std::endl
        << std::endl

        // boundary 4
        << std::setw(OWID) << 1 << std::setw(OWID) << 0 << std::endl

        << std::setw(OWID) << 4 << std::setw(OWID) << 0 << std::setw(OWID) << 1
        << std::setw(OWID) << 0 << std::setw(OWID) << x0 << std::setw(OWID)
        << y2 << std::setw(OWID) << z0 << std::endl
        << std::endl

        // boundary 5
        << std::setw(OWID) << 1 << std::setw(OWID) << 0 << std::endl

        << std::setw(OWID) << 5 << std::setw(OWID) << 0 << std::setw(OWID) << 0
        << std::setw(OWID) << -1 << std::setw(OWID) << x0 << std::setw(OWID)
        << y0 << std::setw(OWID) << z1 << std::endl
        << std::endl

        // boundary 6
        << std::setw(OWID) << 1 << std::setw(OWID) << 0 << std::endl

        << std::setw(OWID) << 6 << std::setw(OWID) << 0 << std::setw(OWID) << 0
        << std::setw(OWID) << 1 << std::setw(OWID) << x0 << std::setw(OWID)
        << y0 << std::setw(OWID) << z2 << std::endl
        << std::endl;
  }

  ofs.close();
}

//////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// pd part
/////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void
Assembly::initialPeriDynamics(const std::string& inputFile)
{
  /* // not used in DEM-PD coupling code, i.e. rigidInclusion()
      //std::cout <<
     "------------------------------------------------------------------------------"
     << std::endl;
      //std::cout << "Problem Initilization " << std::endl;
      //std::cout <<
     "------------------------------------------------------------------------------"
     << std::endl;
      //std::cout << "Read data file ..." << std::endl;
      readPeriDynamicsData(inputFile);
      //std::cout << "Calculate particle volume ..." << std::endl;
      calcParticleVolume();
      // writeMeshCheckVolume("checkv.dat"); exit(1);
      //std::cout << "Calculate horizon size ..." << std::endl;
      calcHorizonSize();
      //std::cout << "Construct neighor list ..." << std::endl;
      constructNeighbor();
      //std::cout << "Calculate Kinv ..." << std::endl;
      calcParticleKinv();
      for(PeriParticlePArray::iterator pt=periParticleVec.begin();
     pt!=periParticleVec.end(); pt++){
          (*pt)->initial();
      }

      // prescrible the essential boundary condition
      //std::cout << "Prescribe the boundary condition ..." << std::endl;
      prescribeEssentialBoundaryCondition(0);

      // calculate the stress at each particle
      //std::cout << "Calculate particle stress ..." << std::endl;
      calcParticleStress();

      // calculate the acceleration for each particle
      //std::cout << "Calculate particle acceleration ..." << std::endl;
      calcParticleAcceleration();

      // traction boundary
      ApplyExternalForce(0);

      // apply initial velocity boundary condition, in this case give the
     cubicTopBoundaryVec particles initial velocities
      //for(std::vector<PeriParticle*>::iterator pt=cubicTopBoundaryVec.begin();
     pt!=cubicTopBoundaryVec.end(); pt++){
      //    (*pt)->setInitVelocity(Vec(0.0, 0.0, 1.0));
      //}
  */
} // end initialPeriDynamics

void
Assembly::prescribeEssentialBoundaryCondition(const int istep)
{
  //    for(PeriParticlePArray::iterator pt=bottomBoundaryVec.begin();
  //    pt!=bottomBoundaryVec.end(); pt++){
  //        (*pt)->prescribeBottomDisplacement(0.0);    // fix z displacement in
  //        the z direction
  //    }

  // REAL dispz;
  // if(istep <= 200) {
  //    dispz = 1.8*0.05*double(istep)/200.;
  //}else {
  //    dispz = 1.8*0.05;
  //}

  // for(PeriParticlePArray::iterator pt=topBoundaryVec.begin();
  // pt!=topBoundaryVec.end(); pt++){
  //    (*pt)->prescribeTopDisplacement(dispz);    // fix z displacement in the
  //    z direction
  //}

} // end prescribeEssentialBoundaryCondition()

void
Assembly::solve(const std::string& outputFile)
{
  /* // not used in the DEM-PD coupling code, i.e. rigidInclusion()
      // open the tecplot file for output
      std::ofstream ofs(outputFile);
      int iframe = 0;
      writeParticleTecplot(ofs,iframe);
      //std::cout <<
  "------------------------------------------------------------------------------"
  << std::endl;
      //std::cout << "Start of the time loop " << std::endl;
      //std::cout <<
  "------------------------------------------------------------------------------"
  << std::endl;
      std::ofstream datafile("uxyz.dat");
      datafile.setf(std::ios::scientific, std::ios::floatfield);
      datafile.precision(10);
      datafile << "VARIABLES = \"Time step\", \"UX\", \"UY\", \"UZ\"" <<
  std::endl;

      int nsteps = 10000;    // is not true
      for(int istep = 1; istep <= nsteps; istep++) {

          runFirstHalfStep();

          prescribeEssentialBoundaryCondition(istep);

          checkBondParticleAlive();

          calcParticleStress();

          calcParticleAcceleration();

          ApplyExternalForce(istep);

          runSecondHalfStep();

  //        if( istep % printInterval == 0) {
  //        //std::cout << "*** current time step is    " << istep << std::endl;
  //        iframe++;
  //        writeParticleTecplot(ofs,iframe);
  //        datafile << istep
  //        << std::setw(20) << periParticleVec[568]->getDisplacement().x()
  //        << std::setw(20) << periParticleVec[568]->getDisplacement().y()
  //        << std::setw(20) << periParticleVec[568]->getDisplacement().z()
  << std::endl;
  //        }
  //        if( istep % 200 == 0) {
  //        writeDisplacementData("ux.dat","uy.dat","uz.dat");
  //        }

      } // time loop
      ofs.close();
      datafile.close();
      //std::cout <<
  "------------------------------------------------------------------------------"
  << std::endl;
      //std::cout << "Simulation Finished !" << std::endl;
      //std::cout <<
  "------------------------------------------------------------------------------"
  << std::endl;
      writeDisplacementData("ux.dat","uy.dat","uz.dat");
  */
} // end solve()

void
Assembly::writeDisplacementData(const std::string& outputFilex,
                                const std::string& outputFiley,
                                const std::string& outputFilez)
{
  /* // not used in the DEM-PD coupling code, i.e. rigidInclusion()
      // displacment along the x axis
      std::ofstream ofs(outputFilex);
      ofs.setf(std::ios::scientific, std::ios::floatfield);
      ofs.precision(10);
      ofs << "VARIABLES = \"X\", \"UX\"" << std::endl;
  //    for(int index = 0; index < 5; index++){
  //        int node = Uxindex[index];
  //        ofs << std::setw(20) <<
  periParticleVec[node]->getInitPosition().x()
  //        << std::setw(20) << periParticleVec[node]->getDisplacement().x()
  << std::endl;
  //    }
      ofs.flush();
      ofs.close();

      //dispalcement along the y axis
      ofs.open(outputFiley);
      ofs.setf(std::ios::scientific, std::ios::floatfield);
      ofs.precision(10);
      ofs << "VARIABLES = \"Y\", \"UY\"" << std::endl;
  //    for(int index = 0; index < 5; index++){
  //        int node = Uyindex[index];
  //        ofs << std::setw(20) <<
  periParticleVec[node]->getInitPosition().y()
  //        << std::setw(20) << periParticleVec[node]->getDisplacement().y()
  << std::endl;
  //    }
      ofs.flush();
      ofs.close();

      //dispalcement along the z axis
      ofs.open(outputFilez);
      ofs.setf(std::ios::scientific, std::ios::floatfield);
      ofs.precision(10);
      ofs << "VARIABLES = \"Z\", \"UZ\"" << std::endl;
  //    for(int index = 0; index < 39; index++){
  //        int node = Uzindex[index];
  //        ofs << std::setw(20) <<
  periParticleVec[node]->getInitPosition().z()
  //        << std::setw(20) << periParticleVec[node]->getDisplacement().z()
  << std::endl;
  //    }
      ofs.flush();
      ofs.close();
  */
} // end writeDisplacementData

void
Assembly::runFirstHalfStep()
{
  int ompThreads = util::getParam<int>("ompThreads");
  int num; // number of peri-points
  int i;
  num = periParticleVec.size();
#pragma omp parallel for num_threads(ompThreads) private(i) shared(num)        \
  schedule(dynamic)
  for (i = 0; i < num; i++) {
    periParticleVec[i]->updateDisplacement();
  }
} // end runFirstHalfStep()

void
Assembly::runSecondHalfStep()
{
  int ompThreads = util::getParam<int>("ompThreads");
  int num; // number of peri-points
  int i;
  num = periParticleVec.size();
#pragma omp parallel for num_threads(ompThreads) private(i) shared(num)        \
  schedule(dynamic)
  for (i = 0; i < num; i++) {
    periParticleVec[i]->updateVelocity();
  }
} // end runSecondHalfStep()

void
Assembly::constructNeighbor()
{ // this function should be called after
  // scattering in each cpu to construct
  // peri-bonds
  // neighbor - searches and constructs the particle neighborlists
  // construct neighborlist for all particles ...
  // compute the weighting function for all particles ...
  if (periParticleVec.empty())
    return;
  for (auto& i_nt : periParticleVec) {
    i_nt->clearPeriBonds(); // bondVec should be empty at this time
  }
  periBondVec.clear();
  std::size_t num = periParticleVec.size();
  int ompThreads = util::getParam<int>("ompThreads");
  std::size_t i_nt, j_nt;
#pragma omp parallel for num_threads(ompThreads) private(i_nt, j_nt) shared(   \
  num) schedule(dynamic)
  for (i_nt = 0; i_nt < num - 1; i_nt++) {
    Vec coord0_i = periParticleVec[i_nt]->getInitPosition();
    REAL horizonSize_i = periParticleVec[i_nt]->getHorizonSize();
    for (j_nt = i_nt + 1; j_nt < num; j_nt++) {
      Vec coord0_j = periParticleVec[j_nt]->getInitPosition();
      REAL tmp_length = vfabs(coord0_i - coord0_j);

      REAL horizonSize_j = periParticleVec[j_nt]->getHorizonSize();
      REAL horizonSize_ij =
        (horizonSize_i + horizonSize_j) *
        0.5; // This will lead to the fact that horizion is not a sphere!!!

      REAL ratio = tmp_length / horizonSize_ij;

      // establish the neighbor list
      if (ratio <= 2.0) {
        // create bond
        PeriBondP bond_pt = std::make_shared<pd::PeriBond>(
          tmp_length, periParticleVec[i_nt], periParticleVec[j_nt]);
#pragma omp critical
        {
          periParticleVec[i_nt]->pushBackBondVec(bond_pt);
          periParticleVec[j_nt]->pushBackBondVec(bond_pt);
          periBondVec.push_back(bond_pt);
        }
        REAL factor =
          3.0 / (2.0 * Pi * horizonSize_ij * horizonSize_ij *
                 horizonSize_ij); // for the factor of 3d window function

        // weighting function (influence function)
        if (ratio < 1.0) {
          bond_pt->setWeight(
            factor * (2.0 / 3.0 - ratio * ratio + 0.5 * ratio * ratio * ratio));
        } else {
          bond_pt->setWeight(factor * (2.0 - ratio) * (2.0 - ratio) *
                             (2.0 - ratio) / 6.0);
        }
      } // if(ratio<2.0)

    } // end j_nt
  }   // end i_nt

} // end constNeighbor()

void
Assembly::findRecvPeriBonds()
{ // this function should be called after
  // commuPeriParticle() in each cpu to
  // construct peri-bonds
  // neighbor - searches and constructs the particle neighborlists
  // construct neighborlist for all particles ...
  // compute the weighting function for all particles ...
  if (recvPeriParticleVec.empty())
    return;
  for (auto i_nt = recvPeriParticleVec.begin();
       i_nt < recvPeriParticleVec.end(); i_nt++) {
    (*i_nt)->clearPeriBonds(); // bondVec should be empty at this time
  }
  recvPeriBondVec.clear();
  for (auto i_nt = recvPeriParticleVec.begin();
       i_nt < recvPeriParticleVec.end(); i_nt++) {
    Vec coord0_i = (*i_nt)->getInitPosition();
    REAL horizonSize_i = (*i_nt)->getHorizonSize();
    for (auto j_nt = periParticleVec.begin(); j_nt < periParticleVec.end();
         j_nt++) {
      Vec coord0_j = (*j_nt)->getInitPosition();
      REAL tmp_length = vfabs(coord0_i - coord0_j);

      REAL horizonSize_j = (*j_nt)->getHorizonSize();
      REAL horizonSize_ij =
        (horizonSize_i + horizonSize_j) *
        0.5; // This will lead to the fact that horizion is not a sphere!!!

      REAL ratio = tmp_length / horizonSize_ij;

      // establish the neighbor list
      if (ratio <= 2.0) {
        // create bond
        PeriBondP bond_pt =
          std::make_shared<pd::PeriBond>(tmp_length, *i_nt, *j_nt);
        (*j_nt)->pushBackBondVec(bond_pt);
        (*i_nt)->pushBackBondVec(bond_pt); // i_nt is in recvPeriParticleVec,
                                           // this is to calculate the
                                           // deformationGradient, sigma and
                                           // Kinv
        // for the peri-points in inner cell of recvPeriParticleVec, refer to
        // commuPeriParticle()
        //            bond_pt->setIsRecv();    // bond_pt->isRecv=true;
        recvPeriBondVec.push_back(bond_pt);

        REAL factor =
          3.0 / (2.0 * Pi * horizonSize_ij * horizonSize_ij *
                 horizonSize_ij); // for the factor of 3d window function

        // weighting function (influence function)
        if (ratio < 1.0) {
          bond_pt->setWeight(
            factor * (2.0 / 3.0 - ratio * ratio + 0.5 * ratio * ratio * ratio));
        } else {
          bond_pt->setWeight(factor * (2.0 - ratio) * (2.0 - ratio) *
                             (2.0 - ratio) / 6.0);
        }
      } // if(ratio<2.0)

    } // end j_nt
  }   // end i_nt

  // since the calculation of PeriParticle.calcAcceleration() needs to know the
  // deformationGradient, sigma, and Kinv of the peri-points in peri-bonds
  // even these peri-points are in recvPeriParticleVec, thus commuPeriParticle()
  // communicate two cellSize peri-points, the deformationGradient, sigma and
  // Kinv
  // of the inner peri-points needs to be calculated exactly, thus the
  // peri-bonds between the recvPeriParticles should also be constructed
  for (auto i_nt = recvPeriParticleVec.begin();
       i_nt < recvPeriParticleVec.end() - 1; i_nt++) {
    Vec coord0_i = (*i_nt)->getInitPosition();
    REAL horizonSize_i = (*i_nt)->getHorizonSize();
    for (auto j_nt = i_nt + 1; j_nt < recvPeriParticleVec.end(); j_nt++) {
      Vec coord0_j = (*j_nt)->getInitPosition();    // need to use "<" since if
                                                    // recvPeriParticleVec is
                                                    // empty, then j_nt
      REAL tmp_length = vfabs(coord0_i - coord0_j); // will exceed the limit of
                                                    // recvPeriParticleVec,
                                                    // then segmentational
                                                    // fault

      REAL horizonSize_j = (*j_nt)->getHorizonSize();
      REAL horizonSize_ij =
        (horizonSize_i + horizonSize_j) *
        0.5; // This will lead to the fact that horizion is not a sphere!!!

      REAL ratio = tmp_length / horizonSize_ij;
      // establish the neighbor list
      if (ratio <= 2.0) {
        // create bond
        PeriBondP bond_pt =
          std::make_shared<pd::PeriBond>(tmp_length, *i_nt, *j_nt);
        (*i_nt)->pushBackBondVec(bond_pt);
        (*j_nt)->pushBackBondVec(bond_pt);
        //            bond_pt->setIsRecv();    // bond_pt->isRecv=true;
        recvPeriBondVec.push_back(bond_pt);

        REAL factor =
          3.0 / (2.0 * Pi * horizonSize_ij * horizonSize_ij *
                 horizonSize_ij); // for the factor of 3d window function

        // weighting function (influence function)
        if (ratio < 1.0) {
          bond_pt->setWeight(
            factor * (2.0 / 3.0 - ratio * ratio + 0.5 * ratio * ratio * ratio));
        } else {
          bond_pt->setWeight(factor * (2.0 - ratio) * (2.0 - ratio) *
                             (2.0 - ratio) / 6.0);
        }
      } // if(ratio<2.0)

    } // end j_nt
  }   // end i_nt

} // end findRecvPeriBonds()

void
Assembly::readPeriDynamicsData(const std::string& inputFile)
{ // should only be called by master cpu
  // readData - reads controlling parameters, particle positions and mesh
  // connectivities
  // @param std::string&  - reference of the input file name

  PeriParticleFileReader reader;
  reader.read(inputFile, allPeriParticleVecInitial, connectivity);
  point_interval =
    vfabs(allPeriParticleVecInitial[connectivity[0][0]]->getInitPosition() -
          allPeriParticleVecInitial[connectivity[0][1]]->getInitPosition());

} // readPeriDynamicsData()

void
Assembly::writeMesh(const std::string& outputFile)
{
  /* // not used in DEM-PD coupling code, i.e. rigidInclustion()
      std::ofstream ofs(outputFile);
      ofs.setf(std::ios::scientific, std::ios::floatfield);
      ofs.precision(10);
      ofs << "Title = \"Mesh Checking\"" << std::endl;
      ofs << "VARIABLES = \"X\", \"Y\",\"Z\"" << std::endl;
      ofs << "ZONE N = " << nPeriParticle << " E = " << nele << ", F = FEPOINT
     ET = BRICK" << std::endl;
      for(int node = 0; node < nPeriParticle; node++){
          ofs << std::setw(20) <<
     periParticleVec[node]->getInitPosition().x()
          << std::setw(20) << periParticleVec[node]->getInitPosition().y()
          << std::setw(20) << periParticleVec[node]->getInitPosition().z() <<
     std::endl;
      }
      for(int iel = 0; iel < nele; iel++){
          for(int node = 0; node < 8; node++){
          ofs << std::setw(10) << connectivity[iel][node];
          }
          ofs << std::endl;
      }
      ofs.close();
  */
} // end writeMesh

void
Assembly::writeMeshCheckVolume(const std::string& outputFile)
{
  /* // not used in DEM-PD coupling code, i.e. rigidInclusion()
      std::ofstream ofs(outputFile);
      ofs.setf(std::ios::scientific, std::ios::floatfield);
      ofs.precision(10);
      ofs << "Title = \"Volume Checking\"" << std::endl;
      ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"V\"" << std::endl;
      ofs << "ZONE N = " << nPeriParticle << " E = " << nele << ", F = FEPOINT
     ET = BRICK" << std::endl;
      // Output the coordinates and the array information
      for(PeriParticlePArray::iterator pt = periParticleVecInitial.begin(); pt!=
     periParticleVecInitial.end(); pt++) {
          ofs << std::setw(20) << (*pt)->getInitPosition().x()
          << std::setw(20) << (*pt)->getInitPosition().y()
          << std::setw(20) << (*pt)->getInitPosition().z()
          << std::setw(20) << (*pt)->getParticleVolume()
          << std::endl;
      }
      for(int iel = 0; iel < nele; iel++){
          for(int node = 0; node < 8; node++){
          ofs << std::setw(10) << connectivity[iel][node];
          }
          ofs << std::endl;
      }
      ofs.close();
  */

} // end writeMeshCheckVolume

void
Assembly::writeParticleTecplot(std::ofstream& ofs, const int iframe) const
{
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(10);
  if (iframe == 0) {
    ofs << "Title = \"Particle Information\"" << std::endl;
    ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\" \"Vy\" "
           "\"Vz\" \"KE\" \"P\" \"Mises\""
        << std::endl;
  }
  ofs << "ZONE T =\" " << iframe << "-th Load Step\" " << std::endl;
  // Output the coordinates and the array information
  REAL pressure, vonMisesStress;
  Matrix sigma;
  for (const auto& pt : allPeriParticleVec) {
    sigma = pt->getSigma();
    pressure = sigma(1, 1) + sigma(2, 2) + sigma(3, 3);
    vonMisesStress =
      sqrt(0.5 * ((sigma(1, 1) - sigma(2, 2)) * (sigma(1, 1) - sigma(2, 2)) +
                  (sigma(2, 2) - sigma(3, 3)) * (sigma(2, 2) - sigma(3, 3)) +
                  (sigma(1, 1) - sigma(3, 3)) * (sigma(1, 1) - sigma(3, 3))) +
           3 * (sigma(1, 2) * sigma(1, 2) + sigma(2, 3) * sigma(2, 3) +
                sigma(3, 1) * sigma(3, 1)));
    ofs << std::setw(20)
        << pt->getInitPosition().x() + pt->getDisplacement().x()
        << std::setw(20)
        << pt->getInitPosition().y() + pt->getDisplacement().y()
        << std::setw(20)
        << pt->getInitPosition().z() + pt->getDisplacement().z()
        << std::setw(20) << pt->getDisplacement().x() << std::setw(20)
        << pt->getDisplacement().y() << std::setw(20)
        << pt->getDisplacement().z() << std::setw(20) << pt->getVelocity().x()
        << std::setw(20) << pt->getVelocity().y() << std::setw(20)
        << pt->getVelocity().z() << std::setw(20) << vfabs(pt->getVelocity())
        << std::setw(20) << pressure << std::setw(20) << vonMisesStress
        << std::endl;
    ofs.flush();
  }

} // end writeparticleTecplot

// this function has some problems, allPeriParticleVecInitial does not exist
// anymore after the first gather
void
Assembly::printPeriDomain(const std::string& str) const
{
  std::ofstream ofs(str);
  if (!ofs) {
    //std::cout << "stream error!" << std::endl;
    exit(-1);
  }

  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(10);
  ofs << "Title = \"Particle Information\"" << std::endl;
  ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\" \"Vy\" "
         "\"Vz\" \"Ax\" \"Ay\" \"Az\" \"KE\" \"P\" \"Mises\" \"Volume\" "
         "\"horizonSize\" \"bondsNum\" \"Kinv_det\" \"deformation\""
      << std::endl;
  //    ofs << "ZONE T = \"periDomain\", DATAPACKING=POINT, NODES=" <<
  //    nPeriParticle << ", ELEMENTS=" << nele << ", ZONETYPE=FEBRICK" <<
  //    std::endl;
  // Output the coordinates and the array information
  REAL pressure, vonMisesStress;
  Matrix sigma;
  Matrix Kinv_tmp;
  Matrix deformationG;
  for (const auto& pt : periParticleVec) {
    sigma = pt->getSigma();
    pressure = sigma(1, 1) + sigma(2, 2) + sigma(3, 3);
    vonMisesStress =
      sqrt(0.5 * ((sigma(1, 1) - sigma(2, 2)) * (sigma(1, 1) - sigma(2, 2)) +
                  (sigma(2, 2) - sigma(3, 3)) * (sigma(2, 2) - sigma(3, 3)) +
                  (sigma(1, 1) - sigma(3, 3)) * (sigma(1, 1) - sigma(3, 3))) +
           3 * (sigma(1, 2) * sigma(1, 2) + sigma(2, 3) * sigma(2, 3) +
                sigma(3, 1) * sigma(3, 1)));
    Kinv_tmp = pt->getParticleKinv();
    deformationG = pt->getDeformationGradient();
    ofs << std::setw(20)
        << pt->getInitPosition().x() + pt->getDisplacement().x()
        << std::setw(20)
        << pt->getInitPosition().y() + pt->getDisplacement().y()
        << std::setw(20)
        << pt->getInitPosition().z() + pt->getDisplacement().z()
        << std::setw(20) << pt->getDisplacement().x() << std::setw(20)
        << pt->getDisplacement().y() << std::setw(20)
        << pt->getDisplacement().z() << std::setw(20) << pt->getVelocity().x()
        << std::setw(20) << pt->getVelocity().y() << std::setw(20)
        << pt->getVelocity().z() << std::setw(20) << pt->getAcceleration().x()
        << std::setw(20) << pt->getAcceleration().y() << std::setw(20)
        << pt->getAcceleration().z() << std::setw(20)
        << vfabs(pt->getVelocity()) << std::setw(20) << pressure
        << std::setw(20) << vonMisesStress << std::setw(20)
        << pt->getParticleVolume() << std::setw(20) << pt->getHorizonSize()
        << std::setw(20) << pt->getBondsNumber() << std::setw(20)
        << dem::det(Kinv_tmp) << std::setw(20) << dem::det(deformationG)
        << std::endl;
    ofs.flush();
  }
  //    for(int iel = 0; iel < nele; iel++){
  //        for(int node = 0; node < 8; node++){
  //        ofs << std::setw(10) << connectivity[iel][node];
  //        }
  //        ofs << std::endl;
  //    }

  ofs.close();

} // end printPeriDomain

// this function has some problems, allPeriParticleVecInitial does not exist
// anymore after the first gather
void
Assembly::printRecvPeriDomain(const std::string& str) const
{
  std::ofstream ofs(str);
  if (!ofs) {
    //std::cout << "stream error!" << std::endl;
    exit(-1);
  }

  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(10);
  ofs << "Title = \"Particle Information\"" << std::endl;
  ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\" \"Vy\" "
         "\"Vz\" \"Ax\" \"Ay\" \"Az\" \"KE\" \"P\" \"Mises\" \"Volume\" "
         "\"horizonSize\" \"bondsNum\" \"Kinv_det\" \"deformation\""
      << std::endl;
  //    ofs << "ZONE T = \"periDomain\", DATAPACKING=POINT, NODES=" <<
  //    nPeriParticle << ", ELEMENTS=" << nele << ", ZONETYPE=FEBRICK" <<
  //    std::endl;
  // Output the coordinates and the array information
  REAL pressure, vonMisesStress;
  Matrix sigma;
  Matrix Kinv_tmp;
  Matrix deformationG;
  for (const auto& pt : recvPeriParticleVec) {
    sigma = pt->getSigma();
    pressure = sigma(1, 1) + sigma(2, 2) + sigma(3, 3);
    vonMisesStress =
      sqrt(0.5 * ((sigma(1, 1) - sigma(2, 2)) * (sigma(1, 1) - sigma(2, 2)) +
                  (sigma(2, 2) - sigma(3, 3)) * (sigma(2, 2) - sigma(3, 3)) +
                  (sigma(1, 1) - sigma(3, 3)) * (sigma(1, 1) - sigma(3, 3))) +
           3 * (sigma(1, 2) * sigma(1, 2) + sigma(2, 3) * sigma(2, 3) +
                sigma(3, 1) * sigma(3, 1)));
    Kinv_tmp = pt->getParticleKinv();
    deformationG = pt->getDeformationGradient();
    ofs << std::setw(20)
        << pt->getInitPosition().x() + pt->getDisplacement().x()
        << std::setw(20)
        << pt->getInitPosition().y() + pt->getDisplacement().y()
        << std::setw(20)
        << pt->getInitPosition().z() + pt->getDisplacement().z()
        << std::setw(20) << pt->getDisplacement().x() << std::setw(20)
        << pt->getDisplacement().y() << std::setw(20)
        << pt->getDisplacement().z() << std::setw(20) << pt->getVelocity().x()
        << std::setw(20) << pt->getVelocity().y() << std::setw(20)
        << pt->getVelocity().z() << std::setw(20) << pt->getAcceleration().x()
        << std::setw(20) << pt->getAcceleration().y() << std::setw(20)
        << pt->getAcceleration().z() << std::setw(20)
        << vfabs(pt->getVelocity()) << std::setw(20) << pressure
        << std::setw(20) << vonMisesStress << std::setw(20)
        << pt->getParticleVolume() << std::setw(20) << pt->getHorizonSize()
        << std::setw(20) << pt->getBondsNumber() << std::setw(20)
        << dem::det(Kinv_tmp) << std::setw(20) << dem::det(deformationG)
        << std::endl;
    ofs.flush();
  }
  //    for(int iel = 0; iel < nele; iel++){
  //        for(int node = 0; node < 8; node++){
  //        ofs << std::setw(10) << connectivity[iel][node];
  //        }
  //        ofs << std::endl;
  //    }

  ofs.close();

} // end printRecvPeriDomain

void
Assembly::printPeriParticle(const std::string& str) const
{
  std::ofstream ofs(str);
  if (!ofs) {
    //std::cout << "stream error!" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  ofs << std::setw(OWID) << "Node ID" << std::setw(OWID) << "U1"
      << std::setw(OWID) << "U2" << std::setw(OWID) << "U3" << std::setw(OWID)
      << "S.Mises" << std::setw(OWID) << "S11" << std::setw(OWID) << "S22"
      << std::setw(OWID) << "S33" << std::setw(OWID) << "S12" << std::setw(OWID)
      << "S13" << std::setw(OWID) << "S23" << std::setw(OWID) << "X"
      << std::setw(OWID) << "Y" << std::setw(OWID) << "Z" << std::endl;

  int id = 0;
  Matrix sigma;
  REAL vonMisesStress;
  for (const auto& pt : allPeriParticleVec) {
    sigma = pt->getSigma();
    id++;
    vonMisesStress =
      sqrt(0.5 * ((sigma(1, 1) - sigma(2, 2)) * (sigma(1, 1) - sigma(2, 2)) +
                  (sigma(2, 2) - sigma(3, 3)) * (sigma(2, 2) - sigma(3, 3)) +
                  (sigma(1, 1) - sigma(3, 3)) * (sigma(1, 1) - sigma(3, 3))) +
           3 * (sigma(1, 2) * sigma(1, 2) + sigma(2, 3) * sigma(2, 3) +
                sigma(3, 1) * sigma(3, 1)));
    ofs << std::setw(OWID) << id << std::setw(OWID) << pt->getDisplacement().x()
        << std::setw(OWID) << pt->getDisplacement().y() << std::setw(OWID)
        << pt->getDisplacement().z() << std::setw(OWID) << vonMisesStress
        << std::setw(OWID) << sigma(1, 1) << std::setw(OWID) << sigma(2, 2)
        << std::setw(OWID) << sigma(3, 3) << std::setw(OWID) << sigma(1, 2)
        << std::setw(OWID) << sigma(1, 3) << std::setw(OWID) << sigma(2, 3)
        << std::setw(OWID) << pt->getInitPosition().x() << std::setw(OWID)
        << pt->getInitPosition().y() << std::setw(OWID)
        << pt->getInitPosition().z() << std::endl;
  }

  ofs.close();
}

void
Assembly::printPeriDomainSphere(const std::string& str) const
{
  /* // not used in DEM-PD coupling code, i.e. rigidInclusion()
        std::ofstream ofs(str);
        if(!ofs) {
              //std::cout << "stream error!" << endl; exit(-1);
        }

      ofs.setf(std::ios::scientific, std::ios::floatfield);
      ofs.precision(10);
      ofs << "Title = \"Particle Information\"" << std::endl;
      ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"Sigma_rr\" \"Sigma_tt\"
     \"Tao_rt\" " << std::endl;
      ofs << "ZONE T = \"periDomain\" "<< std::endl;
      // Output the coordinates and the array information
      REAL pressure, vonMisesStress;
      Matrix sigma;
      // this is not to print all the peri-points, it just print the peri-points
     in the interface between the sphere and the box
      for(PeriParticlePArray::const_iterator pt =
     interfacePeriParticleVec.begin(); pt!= interfacePeriParticleVec.end();
     pt++) {
          sigma = (*pt)->getSigma();
          Vec currPosition = (*pt)->currentPos();
          REAL R = vfabs(currPosition);
          REAL theta = acos(currPosition.z()/R);

          // for phi should notice that x can be zero
          REAL phi;
          if( abs(currPosition.x())<EPS ){
          if( currPosition.y()>0 ){
              phi = Pi*0.5;
          }
          else{
              phi = Pi*1.5;
          }
          }
          else {
              phi = atan(currPosition.y()/currPosition.x());
          }

          Matrix Qtrans(3,3);    // the transfor tensor from cartesian to
     spherical coordinate
                  // refer to
     http://www.brown.edu/Departments/Engineering/Courses/En221/Notes/Polar_Coords/Polar_Coords.htm
          Qtrans(1,1) = sin(theta)*cos(phi); Qtrans(1,2) = sin(theta)*sin(phi);
     Qtrans(1,3) = cos(theta);
          Qtrans(2,1) = cos(theta)*cos(phi); Qtrans(2,2) = cos(theta)*sin(phi);
     Qtrans(2,3) = -sin(theta);
          Qtrans(3,1) = -sin(phi);           Qtrans(3,2) = cos(phi);
     Qtrans(3,3) = 0;

          Matrix sigma_sphere = Qtrans*sigma*trans(Qtrans);

          ofs << std::setw(20) << currPosition.x()
          << std::setw(20) << currPosition.y()
          << std::setw(20) << currPosition.z()
          << std::setw(20) << sigma_sphere(1,1)
          << std::setw(20) << sigma_sphere(2,2)
          << std::setw(20) << sigma_sphere(1,2)
          << std::endl;
          ofs.flush();
      }

      ofs.close();
  */
} // end printPeriDomainSphere

void
Assembly::calcDeformationGradient()
{
  // calcDeformationGradient - calculates the deformation gradient for all
  // peri-particles
}

void
Assembly::calcHorizonSize()
{
  maxHorizonSize = 0;
  for (int iel = 0; iel < nele; iel++) {
    // get initial positions of the particles in this connectivity
    int n1 = connectivity[iel][0]; // index of first node
    int n2 = connectivity[iel][1];
    int n4 = connectivity[iel][3];
    int n5 = connectivity[iel][4];

    Vec coord1 =
      allPeriParticleVecInitial[n1 - 1]
        ->getInitPosition(); // the index of input file is starting from 1
    Vec coord2 = allPeriParticleVecInitial[n2 - 1]->getInitPosition();
    Vec coord4 = allPeriParticleVecInitial[n4 - 1]->getInitPosition();
    Vec coord5 = allPeriParticleVecInitial[n5 - 1]->getInitPosition();

    REAL tmp1, tmp2, tmp3;
    tmp1 = vfabs(coord2 - coord1);
    tmp2 = vfabs(coord4 - coord1);
    tmp3 = vfabs(coord5 - coord1);

    REAL tmpmax; // max number in tmp1, tmp2 and tmp3
    tmpmax = std::max(std::max(tmp1, tmp2), std::max(tmp2, tmp3));
    if (1.5075 * tmpmax > maxHorizonSize)
      maxHorizonSize = 1.5075 * tmpmax;

    for (int node = 0; node < 8; node++) {
      allPeriParticleVecInitial[connectivity[iel][node] - 1]
        ->replaceHorizonSizeIfLarger(1.5075 * tmpmax);
    }
  }
  //std::cout << "maxHorizonSize: " << maxHorizonSize << std::endl;

} // end calcHorizonSize()

void
Assembly::setInitIsv()
{
  REAL isv_tmp;
  if (util::getParam<int>("typeConstitutive") ==
      1) { // 1---implicit, 2---explicit
    isv_tmp = util::getParam<REAL>("Chi");
  } else {
    isv_tmp = util::getParam<REAL>("c");
  }

  for (auto& pt : allPeriParticleVecInitial) {
    pt->setInitIsv(isv_tmp);
  }

} // setInitIsv()

void
Assembly::calcParticleVolume()
{
  int* numofpieces;
  REAL* particleVolume;
  REAL xi, eta, zeta;
  numofpieces = new int[nPeriParticle];
  for (int node = 0; node < nPeriParticle; numofpieces[node] = 0, node++)
    ; // initialize
  int nip = 2;
  Matrix gp_loc3D;
  Matrix gp_weight3D;
  pd::gauss3D(nip, gp_loc3D, gp_weight3D);
  particleVolume = new REAL[nPeriParticle];
  for (int node = 0; node < nPeriParticle; particleVolume[node] = 0.0, node++)
    ; // initialize
  Matrix xl(3, 8);
  Matrix shp;

  for (int iel = 0; iel < nele; iel++) {
    /*        // check if this element contains the peri-point that is in the
       sphere (fixedPeriParticleVec)
            int isContain = 0;    // if contains

            // check the first node
            for(PeriParticlePArray::iterator pt=fixedPeriParticleVec.begin();
       pt!=fixedPeriParticleVec.end(); pt++){
            if(periParticleVecInitial[connectivity[iel][0]-1] == (*pt) ||
       periParticleVecInitial[connectivity[iel][1]-1] == (*pt)
            || periParticleVecInitial[connectivity[iel][2]-1] == (*pt) ||
       periParticleVecInitial[connectivity[iel][3]-1] == (*pt)
            || periParticleVecInitial[connectivity[iel][4]-1] == (*pt) ||
       periParticleVecInitial[connectivity[iel][5]-1] == (*pt)
            || periParticleVecInitial[connectivity[iel][6]-1] == (*pt) ||
       periParticleVecInitial[connectivity[iel][7]-1] == (*pt) ){
                isContain = 1;
                break;
            }
            }
            if(isContain == 1){
            continue;
            }
    */

    for (int node = 0; node < 8; node++) {
      int nposition = connectivity[iel][node] - 1;
      xl(1, node + 1) =
        allPeriParticleVecInitial[nposition]->getInitPosition().x();
      xl(2, node + 1) =
        allPeriParticleVecInitial[nposition]->getInitPosition().y();
      xl(3, node + 1) =
        allPeriParticleVecInitial[nposition]->getInitPosition().z();
      numofpieces[nposition] += 1;
    }
    for (int ik = 0; ik < 8; ik++) {
      xi = gp_loc3D(1, ik + 1);
      eta = gp_loc3D(2, ik + 1);
      zeta = gp_loc3D(3, ik + 1);
      // call fem function to get the volume
      REAL xsj = 0.0;
      pd::shp3d(xi, eta, zeta, xl, shp, xsj);
      for (int node = 0; node < 8; node++) {
        int nposition = connectivity[iel][node] - 1;
        particleVolume[nposition] =
          particleVolume[nposition] + gp_weight3D(1, ik + 1) * xsj / 8.0;
      }
    }
  }

  // Commented to compare result with the fortran code
  for (int node = 0; node < nPeriParticle; node++) {
    particleVolume[node] =
      particleVolume[node] * 8.0 / (REAL(numofpieces[node]));
  }

  // store the particle volume into the object PeriParticle
  for (int node = 0; node < nPeriParticle; node++) {
    allPeriParticleVecInitial[node]->setParticleVolume(particleVolume[node]);
  }

  delete[] numofpieces;
  delete[] particleVolume;
} // end calcParticleVolume

/*    void Assembly::checkBondParticleAlive(){
    // compute the bond length and check for whether a bond or particle is alive
   or not
    for(PeriParticlePArray::iterator pt=periParticleVec.begin();
   pt!=periParticleVec.end(); pt++){
        if( (*pt)->getIsAlive() ){ // particle not alive, then go to next
   particle
        (*pt)->checkParticleAlive(util::getParam<REAL>("bondStretchLimit"));
        }

    } // end particle

    } // end checkBondParticleAlive()
*/

void
Assembly::checkBondParticleAlive()
{
  int ompThreads = util::getParam<int>("ompThreads");
  int num; // number of peri-points
  int pt;
  num = periBondVec.size();
#pragma omp parallel for num_threads(ompThreads) private(pt) shared(num)       \
  schedule(dynamic)
  for (pt = 0; pt < num; pt++) {
    periBondVec[pt]->checkIfAlive();
  }

  // compute the bond length and check for whether a bond or particle is alive
  // or not
  num = periParticleVec.size();
#pragma omp parallel for num_threads(ompThreads) private(pt) shared(num)       \
  schedule(dynamic)
  for (pt = 0; pt < num; pt++) {
    if (periParticleVec[pt]
          ->getIsAlive()) { // particle not alive, then go to next particle
      periParticleVec[pt]->checkParticleAlive();
    }
  } // end particle

} // end checkBondParticleAlive()

void
Assembly::calcParticleKinv()
{ // this should be called after
  // commuDEMPeriParticle(), and find
  // peri-bonds between communicated
  // periParticles
  // Compute the inverse of the shape tensor K
  for (auto& pt : periParticleVec) {
    pt->calcParticleKinv();
  }
  for (auto& pt : recvPeriParticleVec) {
    pt->calcParticleKinv();
  }

} // end calcParticleKinv()

void
Assembly::calcParticleStress()
{
  int ompThreads = util::getParam<int>("ompThreads");
  int num; // number of peri-points
  int i;
  num = periParticleVec.size();
#pragma omp parallel for num_threads(ompThreads) private(i) shared(num)        \
  schedule(dynamic)
  for (i = 0; i < num; i++) {
    periParticleVec[i]->calcParticleStress();
  }
  for (auto& pt : recvPeriParticleVec) {
    pt->calcParticleStress();
  }
} // calcParticleStress()

void
Assembly::calcParticleAcceleration()
{
  int ompThreads = util::getParam<int>("ompThreads");
  int num; // number of peri-points
  int i;
  num = periParticleVec.size();
#pragma omp parallel for num_threads(ompThreads) private(i) shared(num)        \
  schedule(dynamic)
  for (i = 0; i < num; i++) {
    periParticleVec[i]->calcParticleAcceleration();
  }
} // calcParticleAcceleration()

void
Assembly::calcRecvParticleKinv()
{ // this should be called after
  // commuDEMPeriParticle(), and find
  // peri-bonds between communicated
  // periParticles
  // Compute the inverse of the shape tensor K
  for (auto& pt : recvPeriParticleVec) {
    pt->calcParticleKinv();
  }

} // end calcParticleKinv()

void
Assembly::ApplyExternalForce(int istep)
{
  // deal with the external force, applied at the top of the boundary
  REAL factor = 0.0;
  REAL rampStep = util::getParam<REAL>("rampStep");
  if (istep <= rampStep) {
    factor = REAL(istep) / REAL(rampStep);
  } else {
    factor = REAL(1.0);
  }
  //    for(PeriParticlePArray::iterator pt=topBoundaryVec.begin();
  //    pt!=topBoundaryVec.end(); pt++){
  //        (*pt)->setAcceleration(Vec(0.0,0.0,factor*util::getParam<REAL>("bodyDensity")));
  //    }
} // end ApplyExternalForce

// delete those peri-points that are inside sand particles
// actually what are deleted are the peri-points that are inside
// a enlarged sand particle to avoid peri-points go to inside the
// sand after a few steps. Here the enlarged sand particle is a little
// larger than the sand particle, the delta = 0.5*point_interval and
// delta is the difference of principle lengths
void
Assembly::removeInsidePeriParticles()
{
  bool is_inside;
  allPeriParticleVec.clear();
  for (auto& pt : allPeriParticleVecInitial) {
    Vec coord_pt = pt->currentPosition();
    is_inside = false; // if this peri-point is inside sand particles
                       //        if(coord_pt.x()==0 && coord_pt.y()==0 &&
                       //        coord_pt.z()==1.5){
                       //        is_inside = true;
                       //        }

    // remove the inside peri-points that are in the box mesh
    for (auto& dem_pt : allParticleVec) {
      REAL a = dem_pt->getA() + 0.5 * point_interval; // enlarged sand particle
      REAL b = dem_pt->getB() + 0.5 * point_interval;
      REAL c = dem_pt->getC() + 0.5 * point_interval;

      Vec coord_pt_tmp = dem_pt->globalToLocal(
        coord_pt -
        dem_pt->currentPos()); // this is important to get local coordinate

      REAL x_pt = coord_pt_tmp.x();
      REAL y_pt = coord_pt_tmp.y();
      REAL z_pt = coord_pt_tmp.z();

      if ((x_pt / a) * (x_pt / a) + (y_pt / b) * (y_pt / b) +
            (z_pt / c) * (z_pt / c) <
          1) { // this peri-points is inside sand particle
        is_inside = true;
        break; // do not need to check other sand particles
      }
    }

    if (is_inside ==
        false) { // this peri-point is not within any sand particles
      allPeriParticleVec.push_back(pt);
    }
  }

  /*
  // remove the inside peri-points that are in the box mesh
  //    bool is_inside;
  //    periParticleVec.clear();
      for(PeriParticlePArray::iterator pt=periParticleVecInitial.begin();
  pt!=periParticleVecInitial.end(); pt++){
          Vec coord_pt = (*pt)->currentPos();
          is_inside = false;    // if this peri-point is inside sand particles
          for(std::vector<particle*>::iterator dem_pt=ParticleVec.begin();
  dem_pt!=ParticleVec.end(); dem_pt++){
  //        REAL a = (*dem_pt)->getA()+0.5*point_interval;    // enlarged sand
  particle
  //        REAL b = (*dem_pt)->getB()+0.5*point_interval;
  //        REAL c = (*dem_pt)->getC()+0.5*point_interval;

          // since we know the model, the sphere's radius is 0.1m, then we can
  delete points within radius as 0.095m
          REAL a = 0.099;    // enlarged sand particle
          REAL b = 0.099;
          REAL c = 0.099;

          coord_pt = (*dem_pt)->localVec(coord_pt-(*dem_pt)->currentPos());
  // this is important to get local coordinate

          REAL x_pt = coord_pt.x();
          REAL y_pt = coord_pt.y();
          REAL z_pt = coord_pt.z();

          if( (x_pt/a)*(x_pt/a) + (y_pt/b)*(y_pt/b) + (z_pt/c)*(z_pt/c) < 1 ){
  // this peri-points is inside sand particle
              is_inside = true;
              break;    // do not need to check other sand particles
          }
          }


          if(is_inside == false){ // this peri-point is not within any sand
  particles
          periParticleVec.push_back(*pt);
          }
      }
  */

} // end removeInsidePeriParticles()

void
Assembly::clearPeriDEMBonds()
{
  // for (auto bt = periDEMBondVec.begin(); bt != periDEMBondVec.end(); bt++) {
  //  delete (*bt);
  //}
  periDEMBondVec.clear();
  PeriDEMBondPArray().swap(periDEMBondVec); // actual memory release
  plan_gravity = util::getParam<REAL>("periDensity") * point_interval *
                 point_interval * 9.8; // rho*l^2*g, used in PeriDEMBond.cpp
}

void
Assembly::eraseBrokenPeriDEMBonds()
{
  // BB: Feb 3, 2017:
  // Not an efficient operation
  // Better approach may be to use a list if random access of vector
  // members is not needed
  periDEMBondVec.erase(std::remove_if(periDEMBondVec.begin(),
                                      periDEMBondVec.end(),
                                      [](PeriDEMBondP bond) {
                                        if (bond->getIsAlive()) {
                                          return false;
                                        }
                                        return true;
                                      }),
                       periDEMBondVec.end());

  /*
  for (PeriDEMBondPArray::iterator bt = periDEMBondVec.begin();
       bt != periDEMBondVec.end();) {
    if ((*bt)->getIsAlive()) {  // still alive
      bt++;
    } else {
      delete (*bt);
      bt = periDEMBondVec.erase(bt);  // release memory
    }
  }
  */
}

// if a peri-point is first time (i.e. previous step peri-point is not bonded to
// dem particle) entering the influence domain of a DEM particle,
// then construct peri-DEM bond between this peri-point and this DEM particle
void
Assembly::findPeriDEMBonds()
{
  eraseBrokenPeriDEMBonds();
  // construct sand peri-bonds
  // here mergePeriParticleVec is used, since if we have a DEM particle that is
  // crossing the boundary between two cpus, then sand-peri bonds between
  // the communicated peri-points can provide force on the DEM particle from the
  // peri-points in neighboring cpus
  // mergeParticleVec is used, since if a DEM particle is crossing the boundary
  // between two cpus, then these sand-peri bonds between the communicated
  // DEM particles can provide force on peri-points from the DEM particles in
  // neighboring cpus
  REAL delta = point_interval * 3.2; // 3 layers of peri-points

  int ompThreads = util::getParam<int>("ompThreads");
  int num; // number of peri-points
  int peri_pt;
  num = mergePeriParticleVec.size();
#pragma omp parallel for num_threads(ompThreads) private(peri_pt)              \
  shared(num, delta) schedule(dynamic)
  for (peri_pt = 0; peri_pt < num; peri_pt++) {
    Vec xyz_peri = mergePeriParticleVec[peri_pt]->currentPosition();

    for (auto& dem_pt : mergeParticleVec) {
      // check and construct the periDEMBondVec in this particle
      REAL ra = dem_pt->getA();
      REAL rb = dem_pt->getB();
      REAL rc = dem_pt->getC();
      Vec xyz_peri_tmp = dem_pt->globalToLocal(
        xyz_peri - dem_pt->currentPos()); // this is very important, since
                                          // all calculations below for
                                          // ellipsoid
      REAL x_peri =
        xyz_peri_tmp.x(); // are based on the local coordinate of the ellipsoid
      REAL y_peri = xyz_peri_tmp.y();
      REAL z_peri = xyz_peri_tmp.z();
      REAL x_peri_de =
        x_peri /
        (ra + delta); // are based on the local coordinate of the ellipsoid
      REAL y_peri_de = y_peri / (rb + delta);
      REAL z_peri_de = z_peri / (rc + delta);
      if (x_peri_de * x_peri_de + y_peri_de * y_peri_de +
            z_peri_de * z_peri_de <
          1) {
        // check if (*peri_pt) is first time entering the influence domain of
        // (*dem_pt)
        if (mergePeriParticleVec[peri_pt]->isBonded(dem_pt->getId()))
          continue; // already bonded

        // calcualte the local coordinates of the projector on the surface of
        // the ellipsoid particle
        Vec projector_local;

        // get the projector of this peri-point, *peri_pt, on the surface of
        // this ellipsoid, *dem_pt

        // kd should be positive, since (x1,y1,z1) and (x0,y0,z0) are in the
        // same octant
        REAL kd = sqrt(ra * ra * rb * rb * rc * rc /
                       (x_peri * x_peri * rb * rb * rc * rc +
                        y_peri * y_peri * ra * ra * rc * rc +
                        z_peri * z_peri * ra * ra * rb * rb));
        REAL dd =
          sqrt(x_peri * x_peri + y_peri * y_peri + z_peri * z_peri) * (1 - kd);

        REAL n1 = 2 * x_peri / (ra + dd) / (ra + dd);
        REAL n2 = 2 * y_peri / (rb + dd) / (rb + dd);
        REAL n3 = 2 * z_peri / (rc + dd) / (rc + dd);

        REAL aa = n1 * n1 * rb * rb * rc * rc + n2 * n2 * ra * ra * rc * rc +
                  n3 * n3 * ra * ra * rb * rb;
        REAL bb = 2 * (n1 * x_peri * rb * rb * rc * rc +
                       n2 * y_peri * ra * ra * rc * rc +
                       n3 * z_peri * ra * ra * rb * rb);
        REAL cc = rb * rb * rc * rc * x_peri * x_peri +
                  ra * ra * rc * rc * y_peri * y_peri +
                  ra * ra * rb * rb * z_peri * z_peri -
                  ra * ra * rb * rb * rc * rc;

        REAL kn = (-bb + sqrt(bb * bb - 4 * aa * cc)) * 0.5 / aa;

        // //std::cout << "kn: " << kn << std::endl;

        projector_local.setX(x_peri + kn * n1);
        projector_local.setY(y_peri + kn * n2);
        projector_local.setZ(z_peri + kn * n3);

        //            // this is to use the cross point instead of the normal
        //            point as the projector point
        //            projector_local.setX(x_peri*kd);
        //            projector_local.setY(y_peri*kd);
        //            projector_local.setZ(z_peri*kd);

        PeriDEMBondP bond_tmp = std::make_shared<PeriDEMBond>(
          projector_local, dem_pt.get(), mergePeriParticleVec[peri_pt].get());

//            // this is used to test the coupled force model, October 10, 2014
//            // in this test model, the sand-peri-points will move along the
//            dem-particle
//            PeriDEMBondP bond_tmp =
//            std::make_shared<PeriDEMBond>(projector_local,
//               (*dem_pt).get(), (*peri_pt).get());

#pragma omp critical
        periDEMBondVec.push_back(bond_tmp);
      }
    }
  }

} // findPeriDEMBonds()

// construct boundary bonds, sand bonds, July 15, 2014
void
Assembly::constructBoundarySandPeriBonds()
{
  /*
      // delete previous bonds vector
          for(std::vector<dem::PeriBoundaryBond*>::iterator
  pt=topBoundaryBondVec.begin(); pt!=topBoundaryBondVec.end(); pt++){
          delete (*pt);
          }
          for(std::vector<dem::PeriBoundaryBond*>::iterator
  pt=bottomBoundaryBondVec.begin(); pt!=bottomBoundaryBondVec.end(); pt++){
          delete (*pt);
          }
          for(std::vector<dem::PeriBoundaryBond*>::iterator
  pt=leftBoundaryBondVec.begin(); pt!=leftBoundaryBondVec.end(); pt++){
          delete (*pt);
          }
          for(std::vector<dem::PeriBoundaryBond*>::iterator
  pt=rightBoundaryBondVec.begin(); pt!=rightBoundaryBondVec.end(); pt++){
          delete (*pt);
          }
          for(std::vector<dem::PeriBoundaryBond*>::iterator
  pt=frontBoundaryBondVec.begin(); pt!=frontBoundaryBondVec.end(); pt++){
          delete (*pt);
          }
          for(std::vector<dem::PeriBoundaryBond*>::iterator
  pt=backBoundaryBondVec.begin(); pt!=backBoundaryBondVec.end(); pt++){
          delete (*pt);
          }

          topBoundaryBondVec.clear();
          bottomBoundaryBondVec.clear();
          leftBoundaryBondVec.clear();
          rightBoundaryBondVec.clear();
          frontBoundaryBondVec.clear();
          backBoundaryBondVec.clear();

          for(std::vector<dem::PeriDEMBond*>::iterator
  pt=periDEMBondVec.begin(); pt!=periDEMBondVec.end(); pt++){
          delete (*pt);
          }
          periDEMBondVec.clear();


      REAL delta = point_interval*3;    // 3 layers of peri-points

      // boundary coordinates
      REAL x_min = getApt(3).x();
       REAL x_max = getApt(1).x();
      REAL y_min = getApt(4).y();
      REAL y_max = getApt(2).y();
      REAL z_min = getApt(6).z();
      REAL z_max = getApt(5).z();

  // no boundaries in pile penetration
  //    // construct boundary peri-bonds
  //    PeriBoundaryBond* bond_tmp;
  //    for(PeriParticlePArray::iterator pt=periParticleVec.begin();
  pt!=periParticleVec.end(); pt++){ // overlap over peri-points
  //        Vec xyz_pt = (*pt)->currentPos();
  //        REAL x_pt = xyz_pt.x();
  //        REAL y_pt = xyz_pt.y();
  //        REAL z_pt = xyz_pt.z();
  //
  //        if( x_pt-x_min<delta ){    // back boundary points
  //        bond_tmp = new PeriBoundaryBond(Vec(x_min,y_pt,z_pt), (*pt));
  //        backBoundaryBondVec.push_back(bond_tmp);
  //        }
  //        if( x_max-x_pt<delta ){    // front boundary points
  //        bond_tmp = new PeriBoundaryBond(Vec(x_max,y_pt,z_pt), (*pt));
  //        frontBoundaryBondVec.push_back(bond_tmp);
  //        }
  //        if( y_pt-y_min<delta ){    // left boundary points
  //        bond_tmp = new PeriBoundaryBond(Vec(x_pt,y_min,z_pt), (*pt));
  //        leftBoundaryBondVec.push_back(bond_tmp);
  //        }
  //        if( y_max-y_pt<delta ){    // right boundary points
  //        bond_tmp = new PeriBoundaryBond(Vec(x_pt,y_max,z_pt), (*pt));
  //        rightBoundaryBondVec.push_back(bond_tmp);
  //        }
  //        if( z_pt-z_min<delta ){    // bottom boundary points
  //        bond_tmp = new PeriBoundaryBond(Vec(x_pt,y_pt,z_min), (*pt));
  //        bottomBoundaryBondVec.push_back(bond_tmp);
  //        }
  //        if( z_max-z_pt<delta ){    // top boundary points
  //        bond_tmp = new PeriBoundaryBond(Vec(x_pt,y_pt,z_max), (*pt));
  //        topBoundaryBondVec.push_back(bond_tmp);
  //        }
  //
  //    }
  //

      // construct sand peri-bonds
      for(PeriParticlePArray::iterator peri_pt=interfacePeriParticleVec.begin();
  peri_pt!=interfacePeriParticleVec.end(); peri_pt++){
          Vec xyz_peri = (*peri_pt)->currentPos();

          for(ParticlePArray::iterator dem_pt=ParticleVec.begin();
  dem_pt!=ParticleVec.end(); dem_pt++){
          // check and construct the periDEMBondVec in this particle
          REAL ra = (*dem_pt)->getA();
          REAL rb = (*dem_pt)->getB();
          REAL rc = (*dem_pt)->getC();
              xyz_peri = (*dem_pt)->localVec(xyz_peri-(*dem_pt)->currentPos());
  // this is very important, since all calculations below for ellipsoid
              REAL x_peri = xyz_peri.x();                        // are based
  on the local coordinate of the ellipsoid
              REAL y_peri = xyz_peri.y();
              REAL z_peri = xyz_peri.z();
  //        if(
  (x_peri/(ra+delta))*(x_peri/(ra+delta))+(y_peri/(rb+delta))*(y_peri/(rb+delta))+(z_peri/(rc+delta))*(z_peri/(rc+delta))
  < 1 ){
                  // calcualte the local coordinates of the projector on the
  surface of the ellipsoid particle
                  Vec projector_local;

                  // get the projector of this peri-point, *peri_pt, on the
  surface of this ellipsoid, *dem_pt

              // kd should be positive, since (x1,y1,z1) and (x0,y0,z0) are in
  the same octant
              REAL kd =
  sqrt(ra*ra*rb*rb*rc*rc/(x_peri*x_peri*rb*rb*rc*rc+y_peri*y_peri*ra*ra*rc*rc+z_peri*z_peri*ra*ra*rb*rb));
  //            REAL dd =
  sqrt(x_peri*x_peri+y_peri*y_peri+z_peri*z_peri)*(1-kd);
  //
  //            REAL n1 = 2*x_peri/(ra+dd)/(ra+dd);
  //            REAL n2 = 2*y_peri/(rb+dd)/(rb+dd);
  //            REAL n3 = 2*z_peri/(rc+dd)/(rc+dd);
  //
  //            REAL aa = n1*n1*rb*rb*rc*rc+n2*n2*ra*ra*rc*rc+n3*n3*ra*ra*rb*rb;
  //            REAL bb =
  2*(n1*x_peri*rb*rb*rc*rc+n2*y_peri*ra*ra*rc*rc+n3*z_peri*ra*ra*rb*rb);
  //            REAL cc =
  rb*rb*rc*rc*x_peri*x_peri+ra*ra*rc*rc*y_peri*y_peri+ra*ra*rb*rb*z_peri*z_peri-ra*ra*rb*rb*rc*rc;
  //
  //            REAL kn = (-bb+sqrt(bb*bb-4*aa*cc))*0.5/aa;
  //
  ////std::cout << "kn: " << kn << std::endl;
  //
  //
  //            projector_local.setX(x_peri+kn*n1);
  //            projector_local.setY(y_peri+kn*n2);
  //            projector_local.setZ(z_peri+kn*n3);
  //
  //
              // this is to use the cross point instead of the normal point as
  the projector point
              projector_local.setX(x_peri*kd);
              projector_local.setY(y_peri*kd);
              projector_local.setZ(z_peri*kd);

              PeriDEMBond* bond_tmp = new PeriDEMBond(projector_local, *dem_pt,
  *peri_pt);

  //            // this is used to test the coupled force model, October 10,
  2014
  //            // in this test model, the sand-peri-points will move along the
  dem-particle
  //            PeriDEMBond* bond_tmp = new PeriDEMBond(projector_local,
  *dem_pt, *peri_pt);

              periDEMBondVec.push_back(bond_tmp);

  //        }
          }
      }

  */
} // end constructBoundarySandPeriBonds

void
Assembly::applyCoupledForces()
{
  /* no boundaries in the pile penetration
      // calculate the coupled forces between boundaries
      // calculate current length of the bonds
      // boundary coordinates
      REAL x_min = getApt(3).x();
       REAL x_max = getApt(1).x();
      REAL y_min = getApt(4).y();
      REAL y_max = getApt(2).y();
      REAL z_min = getApt(6).z();
      REAL z_max = getApt(5).z();
      for(PeriBoundaryBondPArray::iterator
     bond_pt=bottomBoundaryBondVec.begin();
     bond_pt!=bottomBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(z_min, 6);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=topBoundaryBondVec.begin();
     bond_pt!=topBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(z_max, 5);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=backBoundaryBondVec.begin();
     bond_pt!=backBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(x_min, 3);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=frontBoundaryBondVec.begin();
     bond_pt!=frontBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(x_max, 1);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=leftBoundaryBondVec.begin();
     bond_pt!=leftBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(y_min, 4);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=rightBoundaryBondVec.begin();
     bond_pt!=rightBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(y_max, 2);
          // pass the boundary coordinate and the number of the boundary
      }
  */
  // calculate coupled forces between sand particles and peri-points
  for (auto& bond_pt : periDEMBondVec) {
    bond_pt->applyBondForce();
  }
} // end applyCoupledForces()

// this is used to test the coupling model. In this test model, the
// sand-peri-points will move along the
// dem-particles. And the peri-points have no influence on the dem-particles.
// October 10, 2014
void
Assembly::applyCoupledBoundary()
{
  /* no boundaries in the pile penetration
      // calculate the coupled forces between boundaries
      // calculate current length of the bonds
      // boundary coordinates
      REAL x_min = getApt(3).x();
       REAL x_max = getApt(1).x();
      REAL y_min = getApt(4).y();
      REAL y_max = getApt(2).y();
      REAL z_min = getApt(6).z();
      REAL z_max = getApt(5).z();
      for(PeriBoundaryBondPArray::iterator
     bond_pt=bottomBoundaryBondVec.begin();
     bond_pt!=bottomBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(z_min, 6);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=topBoundaryBondVec.begin();
     bond_pt!=topBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(z_max, 5);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=backBoundaryBondVec.begin();
     bond_pt!=backBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(x_min, 3);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=frontBoundaryBondVec.begin();
     bond_pt!=frontBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(x_max, 1);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=leftBoundaryBondVec.begin();
     bond_pt!=leftBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(y_min, 4);
          // pass the boundary coordinate and the number of the boundary
      }

      for(PeriBoundaryBondPArray::iterator bond_pt=rightBoundaryBondVec.begin();
     bond_pt!=rightBoundaryBondVec.end(); bond_pt++){
          (*bond_pt)->applyBondForce(y_max, 2);
          // pass the boundary coordinate and the number of the boundary
      }
  */

  // calculate coupled forces between sand particles and peri-points
  for (auto& bond_pt : periDEMBondVec) {
    bond_pt->applyBondBoundary();
  }

} // end applyCoupledBoundary()

void
Assembly::constructPeriMatrix()
{ // construct the Matrix members in
  // periParticleVec,
  // since currently the transfer of pointer array in Matrix is not implemented
  // well
  for (auto& pt : periParticleVec) {
    pt->constructMatrixMember(); // construct these Matrix members
  }

} // constructPeriMatrix()

void
Assembly::constructRecvPeriMatrix()
{ // construct the Matrix members in
  // periParticleVec,
  // since currently the transfer of pointer array in Matrix is not implemented
  // well
  for (auto& pt : recvPeriParticleVec) {
    pt->constructMatrixMember(); // construct these Matrix members
  }

} // constructRecvPeriMatrix()

void
Assembly::findBoundaryPeriParticles()
{ // it does not matter since this
  // will be only called once
  bottomBoundaryVec.clear();
  frontBoundaryVec.clear();
  leftBoundaryVec.clear();
  topBoundaryInnerVec.clear();
  topBoundaryEdgeVec.clear();
  topBoundaryCornerVec.clear();

  REAL x1 = util::getParam<REAL>("Xmin");
  REAL x2 = util::getParam<REAL>("Xmax");
  REAL y1 = util::getParam<REAL>("Ymin");
  REAL y2 = util::getParam<REAL>("Ymax");
  REAL z1 = util::getParam<REAL>("Zmin");
  REAL z2 = util::getParam<REAL>("Zmax");
  for (auto& peri_pt : periParticleVec) {
    REAL x_peri = peri_pt->getInitPosition().x();
    REAL y_peri = peri_pt->getInitPosition().y();
    REAL z_peri = peri_pt->getInitPosition().z();
    if (z_peri == z2) { // top boundary
      //        if( (x_peri==x2&&y_peri==y2) || x_peri==x1&&y_peri==y2 ||
      //        x_peri==x2&&y_peri==y1 || x_peri==x1&&y_peri==y1){    // corner
      //        points
      //        topBoundaryCornerVec.push_back(*peri_pt);
      //        }
      //        else if(x_peri==x1 || y_peri==y1 || x_peri==x2 || y_peri==y2){
      //        // edge points
      //        topBoundaryEdgeVec.push_back(*peri_pt);
      //        }
      //        else {    // inner points
      //            topBoundaryInnerVec.push_back(*peri_pt);
      //        }
    } else if (z_peri < z1 + 4.5 * point_interval) { // bottom boundary
      bottomBoundaryVec.push_back(peri_pt);
    }

    if (x_peri < x1 + 4.5 * point_interval ||
        x_peri > x2 - 4.5 * point_interval) {
      leftBoundaryVec.push_back(peri_pt);
    }
    if (y_peri < y1 + 4.5 * point_interval ||
        y_peri > y2 - 4.5 * point_interval) {
      frontBoundaryVec.push_back(peri_pt);
    }
  }

} // findBoundaryPeriParticles()

void
Assembly::findFixedPeriParticles()
{ // it does not matter since this will
  // be only called once
  fixedPeriParticleVec.clear();
  REAL radius = util::getParam<REAL>("fixRadius");
  radius = radius +
           point_interval *
             0.2; // in order to find all points on the spherical surface
  REAL x0 = util::getParam<REAL>("periFixCentroidX");
  REAL y0 = util::getParam<REAL>("periFixCentroidY");
  REAL z0 = util::getParam<REAL>("periFixCentroidZ");
  for (auto& peri_pt : periParticleVec) {
    REAL x_peri = peri_pt->getInitPosition().x();
    REAL y_peri = peri_pt->getInitPosition().y();
    REAL z_peri = peri_pt->getInitPosition().z();
    if ((x_peri - x0) * (x_peri - x0) + (y_peri - y0) * (y_peri - y0) +
          (z_peri - z0) * (z_peri - z0) <
        radius * radius + EPS) {
      fixedPeriParticleVec.push_back(peri_pt);
    }
  }

} // findFixedPeriParticles()

void
Assembly::applyPeriBoundaryCondition()
{
  dem::Vec zero_vec = dem::Vec(0, 0, 0);
  //    for(PeriParticlePArray::iterator peri_pt=fixedPeriParticleVec.begin();
  //    peri_pt!=fixedPeriParticleVec.end(); peri_pt++){
  //        (*peri_pt)->prescribeDisplacement(zero_vec);
  //    }

  int ompThreads = util::getParam<int>("ompThreads");
  int num; // number of peri-points
  int peri_pt;
  num = bottomBoundaryVec.size();
  //#pragma omp parallel for num_threads(ompThreads) private(peri_pt)
  // shared(num) schedule(dynamic)
  //    for(peri_pt=0; peri_pt<num; peri_pt++){
  //        bottomBoundaryVec[peri_pt]->prescribeBottomDisplacement(0);
  //    }
  num = leftBoundaryVec.size();
#pragma omp parallel for num_threads(ompThreads) private(peri_pt) shared(num)  \
  schedule(dynamic)
  for (peri_pt = 0; peri_pt < num; peri_pt++) {
    leftBoundaryVec[peri_pt]->prescribeDisplacement(zero_vec);
  }
  num = frontBoundaryVec.size();
#pragma omp parallel for num_threads(ompThreads) private(peri_pt) shared(num)  \
  schedule(dynamic)
  for (peri_pt = 0; peri_pt < num; peri_pt++) {
    frontBoundaryVec[peri_pt]->prescribeDisplacement(zero_vec);
  }

} // applyPeriBoundaryCondition()

void
Assembly::applyTractionBoundary(int g_iteration)
{
  // set traction boundary
  REAL force = util::getParam<REAL>("periForce");
  REAL rampStep = util::getParam<REAL>("rampStep");
  REAL framp = 0;
  if (g_iteration <= rampStep) {
    framp = g_iteration / rampStep;
  } else {
    framp = 1;
  }
  framp = framp * force;
  dem::Vec tmp_vec = dem::Vec(0, 0, framp);
  for (auto& peri_pt : topBoundaryInnerVec) {
    peri_pt->addAccelerationByForce(tmp_vec);
  }
  for (auto& peri_pt : topBoundaryEdgeVec) {
    peri_pt->addAccelerationByForce(tmp_vec * 0.5);
  }
  for (auto& peri_pt : topBoundaryCornerVec) {
    peri_pt->addAccelerationByForce(tmp_vec * 0.25);
  }
  for (auto& peri_pt : bottomBoundaryVec) {
    peri_pt->addAccelerationByForce(-tmp_vec);
  }

} // applyTractionBoundary()

void
Assembly::rigidInclusion()
{

} // rigidInclusion

void
Assembly::pullOutPeri()
{

} // pullOutPeri

} // namespace dem ends
