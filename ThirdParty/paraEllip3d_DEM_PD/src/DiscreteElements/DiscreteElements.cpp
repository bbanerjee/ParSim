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

#include <DiscreteElements/DiscreteElements.h>
#include <Core/Parallel/Patch.h>

#include <Boundary/BoundaryReader.h>
#include <Boundary/CylinderBoundary.h>
#include <Boundary/PlaneBoundary.h>
#include <Core/Const/const.h>
#include <Core/Util/Utility.h>
#include <Core/Math/IntVec.h>
#include <InputOutput/DEMParticleFileReader.h>
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

// Initialize static MPI comm variables
MPI_Comm DiscreteElements::s_mpiWorld = MPI_Comm(-1);
MPI_Comm DiscreteElements::s_cartComm = MPI_Comm(-1);
int DiscreteElements::s_mpiRank = -1;
int DiscreteElements::s_mpiSize = -1;
IntVec DiscreteElements::s_mpiProcs = IntVec(-1,-1,-1);
IntVec DiscreteElements::s_mpiCoords = IntVec(-1,-1,-1);


void
DiscreteElements::deposit(const std::string& boundaryFile,
                  const std::string& particleFile)
{
  // The output folder (default name is .)
  std::string outputFolder(".");

  // Read the input data
  if (s_mpiRank == 0) {
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
  if (s_mpiRank == 0) {

    // Create the output writer in the master process
    auto folderName =  dem::InputParameter::get().datafile["outputFolder"];
    outputFolder = util::createOutputFolder(folderName);
    //std::cout << "Output folder = " << outputFolder << "\n";
    createOutputWriter(outputFolder, iterSnap-1);

    writeBoundaryToFile();
    writeGridToFile();
    writeParticlesToFile(iterSnap);
    printBdryContact();
    debugInf << std::setw(OWID) << "iter" << std::setw(OWID) << "commuT"
             << std::setw(OWID) << "migraT" << std::setw(OWID) << "compuT"
             << std::setw(OWID) << "totalT" << std::setw(OWID) << "overhead%"
             << std::endl;
  }

  // Broadcast the output folder to all processes
  broadcast(boostWorld, outputFolder, 0);
  std::cerr << "Proc = " << s_mpiRank << " outputFolder = " << outputFolder << "\n";

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

      if (s_mpiRank == 0) {
        updateFileNames(iterSnap);
        writeBoundaryToFile();
        writeGridToFile();
        writeParticlesToFile(iterSnap);
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
    if (s_mpiRank == 0 &&
        toCheckTime) // ignore gather and print time at this step
      debugInf << std::setw(OWID) << iteration << std::setw(OWID) << commuT
               << std::setw(OWID) << migraT << std::setw(OWID)
               << totalT - commuT - migraT << std::setw(OWID) << totalT
               << std::setw(OWID) << (commuT + migraT) / totalT * 100
               << std::endl;
    ++iteration;
    //proc0cout << "**NOTICE** End of iteration: " << iteration << "\n";
  }

  if (s_mpiRank == 0)
    closeProg(progressInf);
}

void
DiscreteElements::readBoundary(const std::string& fileName)
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
DiscreteElements::readParticles(const std::string& particleFile)
{
  REAL young = util::getParam<REAL>("young");
  REAL poisson = util::getParam<REAL>("poisson");
  bool doInitialize = (util::getParam<int>("toInitParticle") == 1);

  DEMParticleFileReader reader;
  reader.read(particleFile, young, poisson, doInitialize, allParticleVec,
              gradation);
}

void
DiscreteElements::scatterParticle()
{
  // partition particles and send to each process
  if (s_mpiRank == 0) { // process 0
    setGrid(Box(grid.getMinCorner().x(), grid.getMinCorner().y(),
                grid.getMinCorner().z(), grid.getMaxCorner().x(),
                grid.getMaxCorner().y(),
                getPtclMaxZ(allParticleVec) + gradation.getPtclMaxRadius()));

    Vec v1 = grid.getMinCorner();
    Vec v2 = grid.getMaxCorner();
    Vec vspan = (v2 - v1) / s_mpiProcs;

    //std::cout << "v1 = " << v1 << "\n";
    //std::cout << "v2 = " << v2 << "\n";
    //std::cout << "s_mpiProcs = " << s_mpiProcs << "\n";
    //std::cout << "vspan = " << vspan << "\n";

    auto reqs = new boost::mpi::request[s_mpiSize - 1];
    DEMParticlePArray tmpParticleVec;
    for (int iRank = s_mpiSize - 1; iRank >= 0; --iRank) {
      tmpParticleVec.clear(); // do not release memory!
      int ndim = 3;
      //int coords[3];
      IntVec coords;
      MPI_Cart_coords(s_cartComm, iRank, ndim, coords.data());
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
          particleVec[i] = std::make_shared<DEMParticle>(
            *tmpParticleVec[i]); // default synthesized copy constructor
      } // now particleVec do not share memeory with allParticleVec
    }
    boost::mpi::wait_all(reqs, reqs + s_mpiSize - 1); // for non-blocking send
    delete[] reqs;

  } else { // other processes except 0
    boostWorld.recv(0, mpiTag, particleVec);
  }

  // content of allParticleVec may need to be printed, so do not clear it.
  // if (s_mpiRank == 0) releaseGatheredParticle();
  //proc0cout << "DEM::scatter:: num particles = " << particleVec.size() << "\n";

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
DiscreteElements::createPatch(int iteration, const REAL& ghostWidth) 
{
  // determine container of each process
  Vec v1 = grid.getMinCorner();
  Vec v2 = grid.getMaxCorner();
  Vec vspan = (v2 - v1) / s_mpiProcs;
  Vec lower = v1 + vspan * s_mpiCoords;
  Vec upper = lower + vspan;
  d_patchP = std::make_unique<Patch<DEMParticlePArray>>(s_cartComm, s_mpiRank, s_mpiCoords,
                                      lower, upper, ghostWidth, EPS);
  //std::ostringstream out;
  //out << "s_mpiRank = " << s_mpiRank 
  //    << " iter = " << iteration << " s_mpiCoords = " << s_mpiCoords 
  //    << " lower = " << lower 
  //    << " upper = " << upper << "\n";
  //std::cout << out.str();
}

void 
DiscreteElements::updatePatch(int iteration, const REAL& ghostWidth)
{
  // determine container of each process
  Vec v1 = grid.getMinCorner();
  Vec v2 = grid.getMaxCorner();
  Vec vspan = (v2 - v1) / s_mpiProcs;
  Vec lower = v1 + vspan * s_mpiCoords;
  Vec upper = lower + vspan;
  d_patchP->update(iteration, lower, upper, ghostWidth);

  //std::ostringstream out;
  //out << "s_mpiRank = " << s_mpiRank 
  //    << " iter = " << iteration << " s_mpiCoords = " << s_mpiCoords 
  //    << " lower = " << lower 
  //    << " upper = " << upper << "\n";
  //std::cout << out.str();
}

void
DiscreteElements::commuParticle()
{
  commuParticle(-1);
}

void
DiscreteElements::commuParticle(const int& iteration)
{
  // determine container of each process
  REAL ghostWidth = gradation.getPtclMaxRadius() * 2;
  updatePatch(iteration, ghostWidth);

  // duplicate pointers, pointing to the same memory
  mergeParticleVec.clear();
  mergeParticleVec = particleVec; 

  /*
  std::ostringstream out;
  if (s_mpiRank == 0) {
    out << "Rank: " << s_mpiRank << ": in: " << particleVec.size()
        << " merge: " << mergeParticleVec.size();
  }
  */

  d_patchP->sendRecvGhostXMinus(boostWorld, iteration, mergeParticleVec);
  d_patchP->sendRecvGhostXPlus(boostWorld, iteration, mergeParticleVec);
  d_patchP->waitToFinishX(iteration);
  d_patchP->combineReceivedParticlesX(iteration, mergeParticleVec);

  /*
  if (s_mpiRank == 0) {
    out << " -> " << mergeParticleVec.size();
  }
  */

  d_patchP->sendRecvGhostYMinus(boostWorld, iteration, mergeParticleVec);
  d_patchP->sendRecvGhostYPlus(boostWorld, iteration, mergeParticleVec);
  d_patchP->waitToFinishY(iteration);
  d_patchP->combineReceivedParticlesY(iteration, mergeParticleVec);

  /*
  if (s_mpiRank == 0) {
    out << " -> " << mergeParticleVec.size();
  }
  */

  d_patchP->sendRecvGhostZMinus(boostWorld, iteration, mergeParticleVec);
  d_patchP->sendRecvGhostZPlus(boostWorld, iteration, mergeParticleVec);
  d_patchP->waitToFinishZ(iteration);
  d_patchP->combineReceivedParticlesZ(iteration, mergeParticleVec);

  /*
  if (s_mpiRank == 0) {
    out << " -> " << mergeParticleVec.size();
  }
  */

  //d_patchP->removeDuplicates(mergeParticleVec);

  /*
  if (s_mpiRank == 0) {
    out << " -> " << mergeParticleVec.size() << "\n";

    //std::ostringstream out;
    //out << "Rank: " << s_mpiRank << ": in: " << particleVec.size()
    //    << " recv: " << recvParticleVec.size();
    //out <<  ": out: " << particleVec.size()
    //    << " merge: " << mergeParticleVec.size() << "\n";
    std::cout << out.str();
  }
  */
}

/*
void
DiscreteElements::commuParticle(const int& iteration)
{
*/
/*
  std::ostringstream out;

  // determine container of each process
  Vec v1 = grid.getMinCorner();
  Vec v2 = grid.getMaxCorner();
  Vec vspan = (v2 - v1) / s_mpiProcs;
  Vec lower = v1 + vspan * s_mpiCoords;
  Vec upper = lower + vspan;
  container = Box(lower, upper);
  out << "s_mpiRank = " << s_mpiRank 
      << " iter = " << iteration << " s_mpiCoords = " << s_mpiCoords 
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
  IntVec neighborCoords = s_mpiCoords;
  --neighborCoords.x();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankX1);
  out << "s_mpiRank = " << s_mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankX1 = " << rankX1 << "\n";
  // x2: +x direction
  neighborCoords = s_mpiCoords;
  ++neighborCoords.x();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankX2);
  out << "s_mpiRank = " << s_mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankX2 = " << rankX2 << "\n";
  // y1: -y direction
  neighborCoords = s_mpiCoords;
  --neighborCoords.y();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankY1);
  out << "s_mpiRank = " << s_mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankY1 = " << rankY1 << "\n";
  // y2: +y direction
  neighborCoords = s_mpiCoords;
  ++neighborCoords.y();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankY2);
  out << "s_mpiRank = " << s_mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankY2 = " << rankY2 << "\n";
  // z1: -z direction
  neighborCoords = s_mpiCoords;
  --neighborCoords.z();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankZ1);
  out << "s_mpiRank = " << s_mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankZ1 = " << rankZ1 << "\n";
  // z2: +z direction
  neighborCoords = s_mpiCoords;
  ++neighborCoords.z();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankZ2);
  out << "s_mpiRank = " << s_mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankZ2 = " << rankZ2 << "\n";
  // x1y1
  neighborCoords = s_mpiCoords;
  --neighborCoords.x();
  --neighborCoords.y();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankX1Y1);
  out << "s_mpiRank = " << s_mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankX1Y1 = " << rankX1Y1 << "\n";
  // x1y2
  neighborCoords = s_mpiCoords;
  --neighborCoords.x();
  ++neighborCoords.y();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankX1Y2);
  out << "s_mpiRank = " << s_mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankX1Y2 = " << rankX1Y2 << "\n";
  // x1z1
  neighborCoords = s_mpiCoords;
  --neighborCoords.x();
  --neighborCoords.z();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankX1Z1);
  out << "s_mpiRank = " << s_mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankX1Z1 = " << rankX1Z1 << "\n";
  // x1z2
  neighborCoords = s_mpiCoords;
  --neighborCoords.x();
  ++neighborCoords.z();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankX1Z2);
  out << "s_mpiRank = " << s_mpiRank 
      << " iter = " << iteration << " neighbor = " << neighborCoords 
      << " rankX1Z2 = " << rankX1Z2 << "\n";
  // x2y1
  neighborCoords = s_mpiCoords;
  ++neighborCoords.x();
  --neighborCoords.y();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankX2Y1);
  // x2y2
  neighborCoords = s_mpiCoords;
  ++neighborCoords.x();
  ++neighborCoords.y();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankX2Y2);
  // x2z1
  neighborCoords = s_mpiCoords;
  ++neighborCoords.x();
  --neighborCoords.z();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankX2Z1);
  // x2z2
  neighborCoords = s_mpiCoords;
  ++neighborCoords.x();
  ++neighborCoords.z();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankX2Z2);
  // y1z1
  neighborCoords = s_mpiCoords;
  --neighborCoords.y();
  --neighborCoords.z();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankY1Z1);
  // y1z2
  neighborCoords = s_mpiCoords;
  --neighborCoords.y();
  ++neighborCoords.z();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankY1Z2);
  // y2z1
  neighborCoords = s_mpiCoords;
  ++neighborCoords.y();
  --neighborCoords.z();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankY2Z1);
  // y2z2
  neighborCoords = s_mpiCoords;
  ++neighborCoords.y();
  ++neighborCoords.z();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankY2Z2);
  // x1y1z1
  neighborCoords = s_mpiCoords;
  --neighborCoords.x();
  --neighborCoords.y();
  --neighborCoords.z();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankX1Y1Z1);
  // x1y1z2
  neighborCoords = s_mpiCoords;
  --neighborCoords.x();
  --neighborCoords.y();
  ++neighborCoords.z();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankX1Y1Z2);
  // x1y2z1
  neighborCoords = s_mpiCoords;
  --neighborCoords.x();
  ++neighborCoords.y();
  --neighborCoords.z();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankX1Y2Z1);
  // x1y2z2
  neighborCoords = s_mpiCoords;
  --neighborCoords.x();
  ++neighborCoords.y();
  ++neighborCoords.z();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankX1Y2Z2);
  // x2y1z1
  neighborCoords = s_mpiCoords;
  ++neighborCoords.x();
  --neighborCoords.y();
  --neighborCoords.z();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankX2Y1Z1);
  // x2y1z2
  neighborCoords = s_mpiCoords;
  ++neighborCoords.x();
  --neighborCoords.y();
  ++neighborCoords.z();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankX2Y1Z2);
  // x2y2z1
  neighborCoords = s_mpiCoords;
  ++neighborCoords.x();
  ++neighborCoords.y();
  --neighborCoords.z();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankX2Y2Z1);
  // x2y2z2
  neighborCoords = s_mpiCoords;
  ++neighborCoords.x();
  ++neighborCoords.y();
  ++neighborCoords.z();
  MPI_Cart_rank(s_cartComm, neighborCoords.data(), &rankX2Y2Z2);

  // if found, communicate with neighboring blocks
  DEMParticlePArray particleX1, particleX2;
  DEMParticlePArray particleY1, particleY2;
  DEMParticlePArray particleZ1, particleZ2;
  DEMParticlePArray particleX1Y1, particleX1Y2, particleX1Z1, particleX1Z2;
  DEMParticlePArray particleX2Y1, particleX2Y2, particleX2Z1, particleX2Z2;
  DEMParticlePArray particleY1Z1, particleY1Z2, particleY2Z1, particleY2Z2;
  DEMParticlePArray particleX1Y1Z1, particleX1Y1Z2, particleX1Y2Z1, particleX1Y2Z2;
  DEMParticlePArray particleX2Y1Z1, particleX2Y1Z2, particleX2Y2Z1, particleX2Y2Z2;
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
  debugInf << "rank=" << s_mpiRank
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
    out << "rankX1: s_mpiRank = " << s_mpiRank 
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
    out << "rankX2: s_mpiRank = " << s_mpiRank 
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
    out << "rankX1 after wait: s_mpiRank = " << s_mpiRank 
        << " iter = " << iteration << "\n"
        << " num particles sent =  " << particleX1.size()
        << " num particles recv =  " << rParticleX1.size() << "\n";
  }
  if (rankX2 >= 0) {
    boost::mpi::wait_all(reqX2, reqX2 + 2);
    out << "rankX2 after wait: s_mpiRank = " << s_mpiRank 
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
    DEMParticlePArray testParticleVec;
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
    << s_mpiRank
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

    for (DEMParticlePArray::const_iterator it = testParticleVec.begin(); it !=
    testParticleVec.end();++it)
    debugInf << (*it)->getId() << ' ';
    debugInf << std::endl;
    testParticleVec.clear();
  */
/*
}
*/

void
DiscreteElements::calcTimeStep()
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
  //if (s_mpiRank == 0) {
    //std::cout << "Timestep = " << timeStep
    //          << " Vibration timestep = " << vibraTimeStep
    //          << " Impact timestep = " << impactTimeStep << "\n";
  //}
}

void
DiscreteElements::calcVibraTimeStep()
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

  MPI_Allreduce(&pTimeStep, &vibraTimeStep, 1, MPI_DOUBLE, MPI_MIN, s_mpiWorld);
}

void
DiscreteElements::calcImpactTimeStep()
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

  MPI_Allreduce(&pTimeStep, &impactTimeStep, 1, MPI_DOUBLE, MPI_MIN, s_mpiWorld);
}

void
DiscreteElements::calcContactNum()
{
  std::size_t pContactNum = contactVec.size();
  MPI_Reduce(&pContactNum, &allContactNum, 1, MPI_INT, MPI_SUM, 0, s_mpiWorld);
}

void
DiscreteElements::findContact()
{ // various implementations
  int ompThreads = util::getParam<int>("ompThreads");

  if (ompThreads == 1) { // non-openmp single-thread version, time complexity
                         // bigO(n x n), n is the number of particles.
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    //std::cout << "\t FindContact: MPI Cart rank = " << s_mpiRank
    //          << " World rank = " << world_rank << "\n";
    findContactSingleThread();
  } else if (ompThreads > 1) { // openmp implementation: various loop scheduling
                               // - (static), (static,1), (dynamic), (dynamic,1)
    findContactMultiThread(ompThreads);

  } // end of openmp implementation
}

void
DiscreteElements::findContactSingleThread()
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
    Vec u = particle->currentPosition();
    auto particleType = particle->getType();
    auto particleRad = particle->getA();

    for (auto j = i + 1; j < num2; ++j) {

      auto mergeParticle = mergeParticleVec[j];
      Vec v = mergeParticle->currentPosition();
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

      if ((vnormL2(v - u) < particleRad + mergeParticleRad) &&
          // not both are fixed particles
          (particleType != 1 || mergeParticleType != 1) &&
          // not both are free boundary particles
          (particleType != 5 || mergeParticleType != 5) &&
          // not both are ghost particles
          (particleType != 10 || mergeParticleType != 10)) {

        DEMContact tmpContact(particle.get(), mergeParticle.get());

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
DiscreteElements::findContactMultiThread(int ompThreads)
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
    u = particleVec[i]->currentPosition();
    particleType = particleVec[i]->getType();
    for (j = i + 1; j < num2; ++j) {
      Vec v = mergeParticleVec[j]->currentPosition();
      mergeParticleType = mergeParticleVec[j]->getType();
      if ((vnormL2(v - u) <
           particleVec[i]->getA() + mergeParticleVec[j]->getA()) &&
          // not both are fixed particles
          (particleType != 1 || mergeParticleType != 1) &&
          // not both are free boundary particles
          (particleType != 5 || mergeParticleType != 5) &&
          // not both are ghost particles
          (particleType != 10 || mergeParticleType != 10)) {

        DEMContact tmpContact(particleVec[i].get(), mergeParticleVec[j].get());

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
DiscreteElements::findBdryContact()
{
  for (auto& boundary : boundaryVec) {
    boundary->findBdryContact(particleVec);
  }
}

void
DiscreteElements::clearContactForce()
{
  for (auto& particle : particleVec) {
    particle->clearContactForce();
  }
}

void
DiscreteElements::internalForce()
{
  /*
    std::ostringstream msg;
    msg << "iteration = " << iteration << " proc = " << s_mpiRank;
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

  MPI_Reduce(pAvg, sumAvg, 3, MPI_DOUBLE, MPI_SUM, 0, s_mpiWorld);
  avgNormal = sumAvg[0] / s_mpiSize;
  avgShear = sumAvg[1] / s_mpiSize;
  avgPenetr = sumAvg[2] / s_mpiSize;
}

void
DiscreteElements::boundaryForce()
{
  for (auto& boundary : boundaryVec) {
    boundary->boundaryForce(boundaryTgtMap);
  }
}

void
DiscreteElements::updateParticle()
{
  //proc0cout << "Num DEM particles = " << particleVec.size() << "\n";
  for (auto& particle : particleVec)
    particle->update();
}

void
DiscreteElements::updateGrid()
{
  updateGridMinX();
  updateGridMaxX();
  updateGridMinY();
  updateGridMaxY();
  updateGridMinZ();
  updateGridMaxZ();
}

void
DiscreteElements::updateGridMinX()
{
  REAL pMinX = getPtclMinX(particleVec);
  REAL minX = 0;
  MPI_Allreduce(&pMinX, &minX, 1, MPI_DOUBLE, MPI_MIN, s_mpiWorld);

  setGrid(Box(minX - gradation.getPtclMaxRadius(), grid.getMinCorner().y(),
              grid.getMinCorner().z(), grid.getMaxCorner().x(),
              grid.getMaxCorner().y(), grid.getMaxCorner().z()));
}

void
DiscreteElements::updateGridMaxX()
{
  REAL pMaxX = getPtclMaxX(particleVec);
  REAL maxX = 0;
  MPI_Allreduce(&pMaxX, &maxX, 1, MPI_DOUBLE, MPI_MAX, s_mpiWorld);

  setGrid(Box(grid.getMinCorner().x(), grid.getMinCorner().y(),
              grid.getMinCorner().z(), maxX + gradation.getPtclMaxRadius(),
              grid.getMaxCorner().y(), grid.getMaxCorner().z()));
}

void
DiscreteElements::updateGridMinY()
{
  REAL pMinY = getPtclMinY(particleVec);
  REAL minY = 0;
  MPI_Allreduce(&pMinY, &minY, 1, MPI_DOUBLE, MPI_MIN, s_mpiWorld);

  setGrid(Box(grid.getMinCorner().x(), minY - gradation.getPtclMaxRadius(),
              grid.getMinCorner().z(), grid.getMaxCorner().x(),
              grid.getMaxCorner().y(), grid.getMaxCorner().z()));
}

void
DiscreteElements::updateGridMaxY()
{
  REAL pMaxY = getPtclMaxY(particleVec);
  REAL maxY = 0;
  MPI_Allreduce(&pMaxY, &maxY, 1, MPI_DOUBLE, MPI_MAX, s_mpiWorld);

  setGrid(Box(grid.getMinCorner().x(), grid.getMinCorner().y(),
              grid.getMinCorner().z(), grid.getMaxCorner().x(),
              maxY + gradation.getPtclMaxRadius(), grid.getMaxCorner().z()));
}

void
DiscreteElements::updateGridMinZ()
{
  REAL pMinZ = getPtclMinZ(particleVec);
  REAL minZ = 0;
  MPI_Allreduce(&pMinZ, &minZ, 1, MPI_DOUBLE, MPI_MIN, s_mpiWorld);

  setGrid(Box(grid.getMinCorner().x(), grid.getMinCorner().y(),
              minZ - gradation.getPtclMaxRadius(), grid.getMaxCorner().x(),
              grid.getMaxCorner().y(), grid.getMaxCorner().z()));
}

void
DiscreteElements::updateGridMaxZ()
{
  // update compute grids adaptively due to particle motion
  REAL pMaxZ = getPtclMaxZ(particleVec);
  REAL maxZ = 0;
  MPI_Allreduce(&pMaxZ, &maxZ, 1, MPI_DOUBLE, MPI_MAX, s_mpiWorld);

  // no need to broadcast grid as it is updated in each process
  setGrid(Box(grid.getMinCorner().x(), grid.getMinCorner().y(),
              grid.getMinCorner().z(), grid.getMaxCorner().x(),
              grid.getMaxCorner().y(), maxZ + gradation.getPtclMaxRadius()));
}

void
DiscreteElements::findParticleInBox(const Box& container,
                            const DEMParticlePArray& allParticles,
                            DEMParticlePArray& insideParticles)
{
  insideParticles.clear();
  for (const auto& particle : allParticles) {
    // it is critical to use EPS
    if (container.inside(particle->currentPosition(), EPS)) {
      insideParticles.push_back(particle);
    }
  }
}

/*
void
DiscreteElements::findParticleInBox(const Box& container,
                            const DEMParticlePArray& allParticles,
                            DEMParticlePArray& insideParticles)
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
    Vec center = particle->currentPosition();
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
DiscreteElements::writeBoundaryToFile() const
{
  d_writer->writeDomain(&allContainer);
}

void
DiscreteElements::printBoundary() const
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
DiscreteElements::printBdryContact() const
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
DiscreteElements::writeGridToFile() const
{
  d_writer->setMPIComm(s_cartComm);
  d_writer->setMPIProc(s_mpiProcs.x(), s_mpiProcs.y(), s_mpiProcs.z());
  d_writer->writeGrid(&grid);
}

void
DiscreteElements::writeParticlesToFile(int frame) const
{
  d_writer->writeParticles(&allParticleVec, frame);
  d_writer->writeSieves(&gradation);
}

void
DiscreteElements::writeParticlesToFile(DEMParticlePArray& particles, int frame) const
{
  d_writer->writeParticles(&particles, frame);
}

void
DiscreteElements::printParticle(const std::string& fileName, int frame) const
{
  OutputTecplot<DEMParticlePArray> writer(".", 0);
  writer.setParticleFileName(fileName);
  writer.writeParticles(&allParticleVec, frame);
}

void
DiscreteElements::printParticle(const std::string& fileName, 
                                DEMParticlePArray& particles,
                                int frame) const
{
  OutputTecplot<DEMParticlePArray> writer(".", 0);
  writer.setParticleFileName(fileName);
  writer.writeParticles(&particles, frame);
}

void
DiscreteElements::openDepositProg(std::ofstream& ofs, const std::string& str)
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
DiscreteElements::printDepositProg(std::ofstream& ofs)
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
DiscreteElements::tractionErrorTol(REAL sigma, std::string type, REAL sigmaX,
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
DiscreteElements::generateParticle(std::size_t particleLayers,
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
    DEMParticleP newptcl = std::make_shared<DEMParticle>(
      particleNum + 1, 0, Vec(x0, y0, z0), gradation, young, poisson);
    allParticleVec.push_back(newptcl);
    particleNum++;
  } else if (particleLayers == 1) { // a horizontal layer of free particles
    for (x = x1; x - x2 < EPS; x += diameter)
      for (y = y1; y - y2 < EPS; y += diameter) {
        DEMParticleP newptcl = std::make_shared<DEMParticle>(
          particleNum + 1, 0, Vec(x, y, z0), gradation, young, poisson);
        allParticleVec.push_back(newptcl);
        particleNum++;
      }
  } else if (particleLayers == 2) { // multiple layers of free particles
    for (z = z1; z - z2 < EPS; z += diameter) {
      for (x = x1 + offset; x - x2 < EPS; x += diameter)
        for (y = y1 + offset; y - y2 < EPS; y += diameter) {
          DEMParticleP newptcl = std::make_shared<DEMParticle>(
            particleNum + 1, 0, Vec(x, y, z), gradation, young, poisson);
          allParticleVec.push_back(newptcl);
          particleNum++;
        }
      offset *= -1;
    }
  }

  printParticle(genParticle, 0);
}

void
DiscreteElements::trim(bool toRebuild, const std::string& inputParticle,
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
                   [&x1, &y1, &z1, &x2, &y2, &z2, &maxR](DEMParticleP particle) {
                     Vec center = particle->currentPosition();
                     if (center.x() < x1 || center.x() > x2 ||
                         center.y() < y1 || center.y() > y2 ||
                         center.z() < z1 || center.z() + maxR > z2) {
                       return true;
                     }
                     return false;
                   }),
    allParticleVec.end());

  /*
  DEMParticlePArray::iterator itr;
  Vec center;
  for (auto itr = allParticleVec.begin(); itr != allParticleVec.end(); ) {
    center=(*itr)->currentPosition();
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

  printParticle(trmParticle, 0);
}

void
DiscreteElements::removeParticleOutBox()
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
  //std::cout << "MPI Cart rank = " << s_mpiRank << "\n";
  REAL epsilon = EPS;
  particleVec.erase(
    std::remove_if(
      particleVec.begin(), particleVec.end(),
      [&x1, &y1, &z1, &x2, &y2, &z2, &epsilon](DEMParticleP particle) {
        Vec center = particle->currentPosition();
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
  DEMParticlePArray::iterator itr;
  Vec center;
  //std::size_t flag = 0;

  for (itr = particleVec.begin(); itr != particleVec.end(); ) {
    center=(*itr)->currentPosition();
    // it is critical to use EPS
    if ( !(center.x() - x1 >= -EPS && center.x() - x2 < -EPS &&
       center.y() - y1 >= -EPS && center.y() - y2 < -EPS &&
       center.z() - z1 >= -EPS && center.z() - z2 < -EPS) )
  {
    // debugInf << "iter=" << std::setw(8) << iteration << " rank=" <<
  std::setw(2) << s_mpiRank
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
    for (DEMParticlePArray::const_iterator it = particleVec.begin(); it !=
    particleVec.end(); ++it)
    debugInf << std::setw(3) << (*it)->getId();
    debugInf << std::endl;
    }
  */
}


REAL
DiscreteElements::getPtclMaxX(const DEMParticlePArray& inputParticle) const
{
  if (inputParticle.size() == 0)
    return -1 / EPS;

  auto it = inputParticle.cbegin();
  REAL x0 = (*it)->currentPosition().x();
  for (; it != inputParticle.cend(); ++it) {
    if ((*it)->currentPosition().x() > x0)
      x0 = (*it)->currentPosition().x();
  }
  return x0;
}

REAL
DiscreteElements::getPtclMinX(const DEMParticlePArray& inputParticle) const
{
  if (inputParticle.size() == 0)
    return 1 / EPS;

  auto it = inputParticle.cbegin();
  REAL x0 = (*it)->currentPosition().x();
  for (; it != inputParticle.cend(); ++it) {
    if ((*it)->currentPosition().x() < x0)
      x0 = (*it)->currentPosition().x();
  }
  return x0;
}

REAL
DiscreteElements::getPtclMaxY(const DEMParticlePArray& inputParticle) const
{
  if (inputParticle.size() == 0)
    return -1 / EPS;

  auto it = inputParticle.cbegin();
  REAL y0 = (*it)->currentPosition().y();
  for (; it != inputParticle.cend(); ++it) {
    if ((*it)->currentPosition().y() > y0)
      y0 = (*it)->currentPosition().y();
  }
  return y0;
}

REAL
DiscreteElements::getPtclMinY(const DEMParticlePArray& inputParticle) const
{
  if (inputParticle.size() == 0)
    return 1 / EPS;

  auto it = inputParticle.cbegin();
  REAL y0 = (*it)->currentPosition().y();
  for (; it != inputParticle.cend(); ++it) {
    if ((*it)->currentPosition().y() < y0)
      y0 = (*it)->currentPosition().y();
  }
  return y0;
}

REAL
DiscreteElements::getPtclMaxZ(const DEMParticlePArray& inputParticle) const
{
  if (inputParticle.size() == 0)
    return -1 / EPS;

  auto it = inputParticle.cbegin();
  REAL z0 = (*it)->currentPosition().z();
  for (; it != inputParticle.cend(); ++it) {
    if ((*it)->currentPosition().z() > z0)
      z0 = (*it)->currentPosition().z();
  }
  return z0;
}

REAL
DiscreteElements::getPtclMinZ(const DEMParticlePArray& inputParticle) const
{
  if (inputParticle.size() == 0)
    return 1 / EPS;

  auto it = inputParticle.cbegin();
  REAL z0 = (*it)->currentPosition().z();
  for (; it != inputParticle.cend(); ++it) {
    if ((*it)->currentPosition().z() < z0)
      z0 = (*it)->currentPosition().z();
  }
  return z0;
}

void
DiscreteElements::setCommunicator(boost::mpi::communicator& comm)
{
  boostWorld = comm;
  s_mpiWorld = MPI_Comm(comm);
  int mpiProcX = util::getParam<int>("mpiProcX");
  int mpiProcY = util::getParam<int>("mpiProcY");
  int mpiProcZ = util::getParam<int>("mpiProcZ");
  s_mpiProcs = {{ mpiProcX, mpiProcY, mpiProcZ }};

  // create Cartesian virtual topology (unavailable in boost.mpi)
  int ndim = 3;
  int periods[3] = { 0, 0, 0 };
  int reorder = 0; // s_mpiRank not reordered
  MPI_Cart_create(s_mpiWorld, ndim, s_mpiProcs.data(), periods, reorder, &s_cartComm);
  MPI_Comm_rank(s_cartComm, &s_mpiRank);
  MPI_Comm_size(s_cartComm, &s_mpiSize);
  MPI_Cart_coords(s_cartComm, s_mpiRank, ndim, s_mpiCoords.data());
  mpiTag = 0;
  assert(s_mpiRank == boostWorld.rank());
  // debugInf << s_mpiRank << " " << s_mpiCoords[0] << " " << s_mpiCoords[1] << " " <<
  // s_mpiCoords[2] << std::endl;

  for (int iRank = 0; iRank < s_mpiSize; ++iRank) {
    int ndim = 3;
    int coords[3];
    MPI_Cart_coords(s_cartComm, iRank, ndim, coords);
    if (coords[0] == 0 || coords[0] == mpiProcX - 1 || coords[1] == 0 ||
        coords[1] == mpiProcY - 1 || coords[2] == 0 ||
        coords[2] == mpiProcZ - 1)
      bdryProcess.push_back(iRank);
  }
}


bool
DiscreteElements::isBdryProcess()
{
  return (s_mpiCoords.x() == 0 || s_mpiCoords.x() == s_mpiProcs.x() - 1 ||
          s_mpiCoords.y() == 0 || s_mpiCoords.y() == s_mpiProcs.y() - 1 ||
          s_mpiCoords.z() == 0 || s_mpiCoords.z() == s_mpiProcs.z() - 1);
}


void
DiscreteElements::releaseRecvParticle()
{
  // release memory of received particles
  /*
  for (DEMParticlePArray::iterator it = recvParticleVec.begin(); it !=
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
DiscreteElements::migrateParticle()
{
  //std::ostringstream out;
  //out << "Migrate: Rank: " << s_mpiRank << ": in: " << particleVec.size();

  Vec vspan = grid.getMaxCorner() - grid.getMinCorner();
  Vec width = vspan / s_mpiProcs;

  sentParticleVec.clear();
  recvParticleVec.clear();
  d_patchP->sendRecvMigrateXMinus(boostWorld, iteration, width, particleVec);
  d_patchP->sendRecvMigrateXPlus(boostWorld, iteration, width, particleVec);
  d_patchP->waitToFinishX(iteration);
  d_patchP->combineSentParticlesX(iteration, sentParticleVec);
  d_patchP->combineReceivedParticlesX(iteration, recvParticleVec);
  d_patchP->deleteSentParticles<DEMParticleP>(iteration, sentParticleVec, particleVec);
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
  d_patchP->deleteSentParticles<DEMParticleP>(iteration, sentParticleVec, particleVec);
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
  d_patchP->deleteSentParticles<DEMParticleP>(iteration, sentParticleVec, particleVec);
  d_patchP->addReceivedParticles(iteration, recvParticleVec, particleVec);
  //out << " sentZ : " << sentParticleVec.size()
  //    << " recvZ : " << recvParticleVec.size();

  // delete outgoing particles
  d_patchP->removeParticlesOutsidePatch<DEMParticleP>(particleVec);
  //out << " outside: " << particleVec.size();

  //d_patchP->removeDuplicates<DEMParticleP>(particleVec);
  //out <<  ": dup out: " << particleVec.size() << "\n";
  //std::cout << out.str();

}

/*
void
DiscreteElements::migrateParticle()
{
*/
/*
  Vec vspan = grid.getMaxCorner() - grid.getMinCorner();
  Vec seg = vspan / s_mpiProcs;
  REAL segX = seg.x();
  REAL segY = seg.y();
  REAL segZ = seg.z();
  Vec v1 = container.getMinCorner(); // v1, v2 in terms of process
  Vec v2 = container.getMaxCorner();

  // if a neighbor exists, transfer particles crossing the boundary in between.
  DEMParticlePArray particleX1, particleX2;
  DEMParticlePArray particleY1, particleY2;
  DEMParticlePArray particleZ1, particleZ2;
  DEMParticlePArray particleX1Y1, particleX1Y2, particleX1Z1, particleX1Z2;
  DEMParticlePArray particleX2Y1, particleX2Y2, particleX2Z1, particleX2Z2;
  DEMParticlePArray particleY1Z1, particleY1Z2, particleY2Z1, particleY2Z2;
  DEMParticlePArray particleX1Y1Z1, particleX1Y1Z2, particleX1Y2Z1, particleX1Y2Z2;
  DEMParticlePArray particleX2Y1Z1, particleX2Y1Z2, particleX2Y2Z1, particleX2Y2Z2;
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
    << s_mpiRank
    << "   added=";
    for (DEMParticlePArray::const_iterator it = recvParticleVec.begin(); it !=
    recvParticleVec.end(); ++it)
    debugInf << std::setw(3) << (*it)->getId();
    debugInf << " now " << particleVec.size() << ": ";
    for (DEMParticlePArray::const_iterator it = particleVec.begin(); it !=
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
DiscreteElements::gatherParticle()
{
  // update allParticleVec: process 0 collects all updated particles from each
  // other process
  if (s_mpiRank != 0) { // each process except 0
    boostWorld.send(0, mpiTag, particleVec);
  } else { // process 0
    // allParticleVec is cleared before filling with new data
    releaseGatheredParticle();

    // duplicate particleVec so that it is not destroyed by allParticleVec in
    // next iteration,
    // otherwise it causes memory error.
    DEMParticlePArray dupParticleVec(particleVec.size());
    for (std::size_t i = 0; i < dupParticleVec.size(); ++i)
      dupParticleVec[i] = std::make_shared<DEMParticle>(*particleVec[i]);

    // fill allParticleVec with dupParticleVec and received particles
    allParticleVec.insert(allParticleVec.end(), dupParticleVec.begin(),
                          dupParticleVec.end());

    DEMParticlePArray tmpParticleVec;
    long gatherRam = 0;
    for (int iRank = 1; iRank < s_mpiSize; ++iRank) {
      tmpParticleVec.clear(); // do not destroy particles!
      boostWorld.recv(iRank, mpiTag, tmpParticleVec);
      allParticleVec.insert(allParticleVec.end(), tmpParticleVec.begin(),
                            tmpParticleVec.end());
      gatherRam += tmpParticleVec.size();
    }
    // debugInf << "gather: particleNum = " << gatherRam <<  " particleRam = "
    // << gatherRam * sizeof(DEMParticle) << std::endl;
  }
}

void
DiscreteElements::releaseGatheredParticle()
{
  // clear allParticleVec, avoid long time memory footprint.
  /*
  for (DEMParticlePArray::iterator it = allParticleVec.begin(); it !=
  allParticleVec.end(); ++it)
    delete (*it);
  */
  allParticleVec.clear();
  DEMParticlePArray().swap(allParticleVec); // actual memory release
}


void
DiscreteElements::gatherBdryContact()
{
  if (isBdryProcess()) {
    if (s_mpiRank != 0)
      boostWorld.send(0, mpiTag, boundaryVec);
  }

  if (s_mpiRank == 0) {
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
DiscreteElements::gatherEnergy()
{
  calcTransEnergy();
  calcRotatEnergy();
  calcKinetEnergy();
  calcGraviEnergy(allContainer.getMinCorner().z());
  calcMechaEnergy();
}

void
DiscreteElements::closeProg(std::ofstream& ofs)
{
  ofs.close();
}

void
DiscreteElements::getStartDimension(REAL& distX, REAL& distY, REAL& distZ)
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
DiscreteElements::openCompressProg(std::ofstream& ofs, const std::string& str)
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
DiscreteElements::printCompressProg(std::ofstream& ofs, REAL distX, REAL distY,
                            REAL distZ)
{
  REAL x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0, z1 = 0.0, z2 = 0.0;
  for (const auto& boundary : mergeBoundaryVec) {
    switch (boundary->getId()) {
      case 1:
        x1 = boundary->getPoint().x();
        break;
      case 2:
        x2 = boundary->getPoint().x();
        break;
      case 3:
        y1 = boundary->getPoint().y();
        break;
      case 4:
        y2 = boundary->getPoint().y();
        break;
      case 5:
        z1 = boundary->getPoint().z();
        break;
      case 6:
        z2 = boundary->getPoint().z();
        break;
    }
  }
  REAL areaX = (y2 - y1) * (z2 - z1);
  REAL areaY = (z2 - z1) * (x2 - x1);
  REAL areaZ = (x2 - x1) * (y2 - y1);
  REAL bulkVolume = (x2 - x1) * (y2 - y1) * (z2 - z1);
  REAL voidRatio = bulkVolume / getVolume() - 1;

  proc0cout << "Boundary: size = " << mergeBoundaryVec.size()
            << " areas = [" << areaX << "," << areaY << "," << areaZ << "]"
            << " vol = " << bulkVolume << " void ratio = " << voidRatio << "\n";

  REAL var[6], vel[6];
  // normalForce
  for (std::size_t i = 0; i < 6; ++i) {
    var[i] = 0;
    vel[i] = 0;
  }
  for (const auto& boundary : mergeBoundaryVec) {
    std::size_t id = boundary->getId();
    Vec normal = boundary->getNormalForce();
    Vec veloc = boundary->getVeloc();
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
  for (double& i : var) {
    i = 0;
  }
  for (const auto& boundary : mergeBoundaryVec) {
    std::size_t id = boundary->getId();
    var[id - 1] = boundary->getContactNum();
  }
  for (double i : var)
    ofs << std::setw(OWID) << static_cast<std::size_t>(i);
  ofs << std::setw(OWID) << allContactNum;

  // avgPenetr
  for (double& i : var)
    i = 0;
  for (const auto& boundary : mergeBoundaryVec) {
    std::size_t id = boundary->getId();
    var[id - 1] = boundary->getAvgPenetr();
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
DiscreteElements::openParticleProg(std::ofstream& ofs, const std::string& str)
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
DiscreteElements::updateBoundary(REAL sigma, std::string type, REAL sigmaX, REAL sigmaY)
{
  if (s_mpiRank == 0) {
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
DiscreteElements::printContact(const std::string& str) const
{
  // There are two implementions of printContact
  // implementation 1: parallel IO, each process prints to a data file using a
  // shared pointer.
  //                   and use post-processing tool to remove redundant info.
  MPI_Status status;
  MPI_File contactFile;
  MPI_File_open(s_mpiWorld, str.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
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
  combine(csuf, ".p", s_mpiRank, 5);
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
DiscreteElements::calcTransEnergy()
{
  REAL pEngy = 0;
  DEMParticlePArray::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it) {
    if ((*it)->getType() == 0)
      pEngy += (*it)->getTransEnergy();
  }
  MPI_Reduce(&pEngy, &transEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, s_mpiWorld);
}

void
DiscreteElements::calcRotatEnergy()
{
  REAL pEngy = 0;
  DEMParticlePArray::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it) {
    if ((*it)->getType() == 0)
      pEngy += (*it)->getRotatEnergy();
  }
  MPI_Reduce(&pEngy, &rotatEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, s_mpiWorld);
}

void
DiscreteElements::calcKinetEnergy()
{
  REAL pEngy = 0;
  DEMParticlePArray::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it) {
    if ((*it)->getType() == 0)
      pEngy += (*it)->getKinetEnergy();
  }
  MPI_Reduce(&pEngy, &kinetEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, s_mpiWorld);
}

void
DiscreteElements::calcGraviEnergy(REAL ref)
{
  REAL pEngy = 0;
  DEMParticlePArray::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it) {
    if ((*it)->getType() == 0)
      pEngy += (*it)->getPotenEnergy(ref);
  }
  MPI_Reduce(&pEngy, &graviEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, s_mpiWorld);
}

void
DiscreteElements::calcMechaEnergy()
{
  mechaEnergy = kinetEnergy + graviEnergy;
}

REAL
DiscreteElements::getMass() const
{
  REAL var = 0;
  for (const auto& it : allParticleVec)
    var += it->getMass();
  return var;
}

REAL
DiscreteElements::getVolume() const
{
  REAL var = 0;
  for (const auto& it : allParticleVec)
    if (it->getType() == 0)
      var += it->getVolume();
  return var;
}

REAL
DiscreteElements::getAvgTransVelocity() const
{
  REAL avgv = 0;
  std::size_t count = 0;
  DEMParticlePArray::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it)
    if ((*it)->getType() == 0) {
      avgv += vnormL2((*it)->currentVel());
      ++count;
    }
  return avgv /= count;
}

REAL
DiscreteElements::getAvgRotatVelocity() const
{
  REAL avgv = 0;
  std::size_t count = 0;
  DEMParticlePArray::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it)
    if ((*it)->getType() == 0) {
      avgv += vnormL2((*it)->currentOmega());
      ++count;
    }
  return avgv /= count;
}

REAL
DiscreteElements::getAvgForce() const
{
  REAL avgv = 0;
  std::size_t count = 0;
  DEMParticlePArray::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it)
    if ((*it)->getType() == 0) {
      avgv += vnormL2((*it)->getForce());
      ++count;
    }
  return avgv / count;
}

REAL
DiscreteElements::getAvgMoment() const
{
  REAL avgv = 0;
  std::size_t count = 0;
  DEMParticlePArray::const_iterator it;
  for (it = particleVec.begin(); it != particleVec.end(); ++it)
    if ((*it)->getType() == 0) {
      avgv += vnormL2((*it)->getMoment());
      ++count;
    }
  return avgv /= count;
}

void
DiscreteElements::buildBoundary(std::size_t boundaryNum,
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

} // namespace dem ends
