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
#include <DiscreteElements/DEMParticleCreator.h>
#include <Core/Parallel/Patch.h>

#include <Boundary/BoundaryReader.h>
#include <Boundary/CylinderBoundary.h>
#include <Boundary/PlaneBoundary.h>
#include <Boundary/BoundaryFileWriter.h>
#include <Core/Const/Constants.h>
#include <Core/Util/Utility.h>
#include <Core/Math/IntVec.h>
#include <InputOutput/DEMParticleFileReader.h>
#include <InputOutput/DEMParticleFileWriter.h>
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
#include <typeinfo>

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


/**
 * 1) Reads the bounding box from an input boundary file
 * 2) Creates internal boundaries form DEM contact calculations
 * 3) Creates particles
 * 4) Writes the particle data and internal boundary data to the
 *    output files
 * 
 * REQUIRES: d_gradation
 */
void 
DiscreteElements::createAndSaveParticlesAndBoundaries(
  const std::string& outputBoundaryFilename,
  const std::string& outputParticleFilename)
{
  // Get the domain information
  std::string inputBoundaryFilename =
    InputParameter::get().datafile["boundaryFilename"];
  readBoundary(inputBoundaryFilename);
  Box domain = getSpatialDomain();

  // Write the domain and inner boundaries
  BoundaryFileWriter boundaryWriter;
  std::string xmlfile = util::replaceExtension(outputBoundaryFilename, "xml");
  boundaryWriter.writeXML(5, xmlfile, domain);
  std::string csvfile = util::replaceExtension(outputBoundaryFilename, "csv");
  boundaryWriter.writeCSV(5, csvfile, domain);

  auto layerFlag = util::getParam<std::size_t>("particleLayers");
  DEMParticle::DEMParticleShape particleType = 
    DEMParticle::DEMParticleShape::ELLIPSOID;

  DEMParticleCreator creator;
  DEMParticlePArray particles = 
    creator.generateDEMParticles(layerFlag, particleType, domain, d_gradation);
  setAllDEMParticleVec(particles);

  DEMParticleFileWriter writer;
  xmlfile = util::replaceExtension(outputParticleFilename, "xml");
  writer.writeXML(particles, d_gradation, xmlfile);
  csvfile = util::replaceExtension(outputParticleFilename, "csv");
  writer.writeCSV(particles, d_gradation, csvfile);
}

void
DiscreteElements::deposit(const std::string& boundaryFilename,
                  const std::string& particleFilename)
{
  // The output folder (default name is .)
  std::string outputFolder(".");

  // Read the input data
  if (s_mpiRank == 0) {
    readBoundary(boundaryFilename);
    readParticles(particleFilename);
    openProgressOutputFile(progressInf, "deposit_progress");
  }

  //proc0cout << "**NOTICE** Before scatterparticle\n";
  scatterParticles(); // scatter particles only once; also updates grid for the
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
    writePatchGridToFile();
    writeParticlesToFile(iterSnap);
    printBoundaryContacts();
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
    communicateGhostParticles(iteration);

    if (toCheckTime)
      time2 = MPI_Wtime();
    commuT = time2 - time0;

    //proc0cout << "**NOTICE** Before calcTimeStep\n";
    calcTimeStep(); // use values from last step, must call before
                    // findContact (which clears data)

    //proc0cout << "**NOTICE** Before findContact\n";

    findContact();
    if (isBoundaryProcess())
      findBoundaryContacts();

    clearContactForce();

    //proc0cout << "**NOTICE** Before internalForce\n";
    internalForce();

    if (isBoundaryProcess()) {
      //proc0cout << "**NOTICE** Before updateParticle\n";
      boundaryForce();
    }

    //proc0cout << "**NOTICE** Before updateParticle\n";
    updateParticles();
    updatePatchBox(); // universal; updatePatchBoxMaxZ() for deposition only

    /**/ timeCount += timeStep;
    /**/ timeAccrued += timeStep;
    /**/ if (timeCount >= timeTotal / netSnap) {
      // if (iteration % (netStep / netSnap) == 0) {
      if (toCheckTime)
        time1 = MPI_Wtime();

      gatherParticles();

      gatherBoundaryContacts();
      gatherEnergy();
      if (toCheckTime)
        time2 = MPI_Wtime();
      gatherT = time2 - time1;

      if (s_mpiRank == 0) {
        updateFileNames(iterSnap);
        writeBoundaryToFile();
        writePatchGridToFile();
        writeParticlesToFile(iterSnap);
        printBoundaryContacts();
        appendToProgressOutputFile(progressInf);
      }
      printContact(combine(outputFolder, "contact_", iterSnap, 5));

      /**/ timeCount = 0;
      ++iterSnap;
    }

    //proc0cout << "**NOTICE** Before releaseReceivedParticles\n";
    releaseReceivedParticles(); // late release because printContact refers to
                           // received particles
    if (toCheckTime)
      time1 = MPI_Wtime();
    //proc0cout << "**NOTICE** Before migrateParticles\n";
    migrateParticles();
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
    closeProgressOutputFile(progressInf);
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
    reader.readXML(fileName, d_spatialDomain, d_demPatchBox, d_boundaries);
  } else if (firstChar == '{') { // JSON
    //std::cout << "Using the JSON reader\n";
    BoundaryReader reader;
    reader.readJSON(fileName, d_spatialDomain, d_demPatchBox, d_boundaries);
  } else {
    //std::cout << "Using the default text reader\n";
    BoundaryReader reader;
    reader.read(fileName, d_spatialDomain, d_demPatchBox, d_boundaries);
  }
  // std::cout << "d_demPatchBox = " << d_demPatchBox << "\n";

  updateBoundaryAreas(d_boundaries);
}

// Compute boundary areas assuming axis-aligned plane boundaries
// *TODO* Generalize.  For now if any one boundary is not a 
// axis-aligned PlaneBoundary we don't do the area calculation
// and update
void 
DiscreteElements::updateBoundaryAreas(BoundaryPArray& boundaries)
{
  REAL xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0;
  for (const auto& boundary : boundaries) {
    Boundary::BoundaryID id = boundary->getId();
    Vec position = boundary->getPosition();
    switch (id) {
      case Boundary::BoundaryID::XMINUS:
        xmin = position.x();
        break;
      case Boundary::BoundaryID::XPLUS:
        xmax = position.x();
        break;
      case Boundary::BoundaryID::YMINUS:
        ymin = position.y();
        break;
      case Boundary::BoundaryID::YPLUS:
        ymax = position.y();
        break;
      case Boundary::BoundaryID::ZMINUS:
        zmin = position.z();
        break;
      case Boundary::BoundaryID::ZPLUS:
        zmax = position.z();
        break;
      default:
        break;
    }
  }
  REAL areaX = (ymax - ymin) * (zmax - zmin);
  REAL areaY = (zmax - zmin) * (xmax - xmin);
  REAL areaZ = (xmax - xmin) * (ymax - ymin);

  for (auto& boundary : boundaries) {
    Boundary::BoundaryID id = boundary->getId();
    switch (id) {
      case Boundary::BoundaryID::XMINUS:
        boundary->setArea(areaX);
        break;
      case Boundary::BoundaryID::XPLUS:
        boundary->setArea(areaX);
        break;
      case Boundary::BoundaryID::YMINUS:
        boundary->setArea(areaY);
        break;
      case Boundary::BoundaryID::YPLUS:
        boundary->setArea(areaY);
        break;
      case Boundary::BoundaryID::ZMINUS:
        boundary->setArea(areaZ);
        break;
      case Boundary::BoundaryID::ZPLUS:
        boundary->setArea(areaZ);
        break;
      default:
        break;
    }
  }
}

void
DiscreteElements::readParticles(const std::string& particleFilename)
{
  REAL young = util::getParam<REAL>("young");
  REAL poisson = util::getParam<REAL>("poisson");
  bool doInitialize = (util::getParam<int>("demToInitParticle") == 1);

  DEMParticleFileReader reader;
  reader.read(particleFilename, young, poisson, doInitialize, d_allDEMParticles,
              d_gradation);
}

void
DiscreteElements::readBoundaryConditions(const std::string& bcFilename)
{
  d_bc.read(bcFilename);
}

void
DiscreteElements::scatterParticles()
{
  // partition particles and send to each process
  if (s_mpiRank == 0) { // process 0
    setPatchBox(Box(d_demPatchBox.minCorner().x(), d_demPatchBox.minCorner().y(),
                d_demPatchBox.minCorner().z(), d_demPatchBox.maxCorner().x(),
                d_demPatchBox.maxCorner().y(),
                getPtclMaxZ(d_allDEMParticles) + d_gradation.getPtclMaxRadius()));

    Vec v1 = d_demPatchBox.minCorner();
    Vec v2 = d_demPatchBox.maxCorner();
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
      Box domain(lower, upper);
      //std::cout << "lower = " << lower 
      //          << " upper = " << upper << "\n";

      findParticleInBox(domain, d_allDEMParticles, tmpParticleVec);

      if (iRank != 0)
        reqs[iRank - 1] = boostWorld.isend(iRank, mpiTag,
                                           tmpParticleVec); // non-blocking send
      if (iRank == 0) {
        d_patchParticles.resize(tmpParticleVec.size());
        for (auto i = 0u; i < d_patchParticles.size(); ++i)
          d_patchParticles[i] = std::make_shared<DEMParticle>(
            *tmpParticleVec[i]); // default synthesized copy constructor
      } // now d_patchParticles do not share memeory with d_allDEMParticles
    }
    boost::mpi::wait_all(reqs, reqs + s_mpiSize - 1); // for non-blocking send
    delete[] reqs;

  } else { // other processes except 0
    boostWorld.recv(0, mpiTag, d_patchParticles);
  }

  // content of d_allDEMParticles may need to be printed, so do not clear it.
  // if (s_mpiRank == 0) releaseGatheredParticle();
  //proc0cout << "DEM::scatter:: num particles = " << d_patchParticles.size() << "\n";

  // broadcast necessary info
  broadcast(boostWorld, d_gradation, 0);
  broadcast(boostWorld, d_boundaries, 0);
  broadcast(boostWorld, d_spatialDomain, 0);
  broadcast(boostWorld, d_demPatchBox, 0);

  // Create patch for the current process
  REAL ghostWidth = d_gradation.getPtclMaxRadius() * 2;
  createPatch(0, ghostWidth);
}

void 
DiscreteElements::createPatch(int iteration, const REAL& ghostWidth) 
{
  // determine domain of each process
  Vec v1 = d_demPatchBox.minCorner();
  Vec v2 = d_demPatchBox.maxCorner();
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
  // determine domain of each process
  Vec v1 = d_demPatchBox.minCorner();
  Vec v2 = d_demPatchBox.maxCorner();
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
DiscreteElements::communicateGhostParticles()
{
  communicateGhostParticles(-1);
}

void
DiscreteElements::communicateGhostParticles(const int& iteration)
{
  // determine domain of each process
  REAL ghostWidth = d_gradation.getPtclMaxRadius() * 2;
  updatePatch(iteration, ghostWidth);

  // duplicate pointers, pointing to the same memory
  d_mergedParticles.clear();
  d_mergedParticles = d_patchParticles; 

  /*
  std::ostringstream out;
  if (s_mpiRank == 0) {
    out << "Rank: " << s_mpiRank << ": in: " << d_patchParticles.size()
        << " merge: " << d_mergedParticles.size();
  }
  */

  d_patchP->sendRecvGhostXMinus(boostWorld, iteration, d_mergedParticles);
  d_patchP->sendRecvGhostXPlus(boostWorld, iteration, d_mergedParticles);
  d_patchP->waitToFinishX(iteration);
  d_patchP->combineReceivedParticlesX(iteration, d_mergedParticles);

  /*
  if (s_mpiRank == 0) {
    out << " -> " << d_mergedParticles.size();
  }
  */

  d_patchP->sendRecvGhostYMinus(boostWorld, iteration, d_mergedParticles);
  d_patchP->sendRecvGhostYPlus(boostWorld, iteration, d_mergedParticles);
  d_patchP->waitToFinishY(iteration);
  d_patchP->combineReceivedParticlesY(iteration, d_mergedParticles);

  /*
  if (s_mpiRank == 0) {
    out << " -> " << d_mergedParticles.size();
  }
  */

  d_patchP->sendRecvGhostZMinus(boostWorld, iteration, d_mergedParticles);
  d_patchP->sendRecvGhostZPlus(boostWorld, iteration, d_mergedParticles);
  d_patchP->waitToFinishZ(iteration);
  d_patchP->combineReceivedParticlesZ(iteration, d_mergedParticles);

  /*
  if (s_mpiRank == 0) {
    out << " -> " << d_mergedParticles.size();
  }
  */

  //d_patchP->removeDuplicates(d_mergedParticles);

  /*
  if (s_mpiRank == 0) {
    out << " -> " << d_mergedParticles.size() << "\n";

    //std::ostringstream out;
    //out << "Rank: " << s_mpiRank << ": in: " << d_patchParticles.size()
    //    << " recv: " << d_receivedParticles.size();
    //out <<  ": out: " << d_patchParticles.size()
    //    << " merge: " << d_mergedParticles.size() << "\n";
    std::cout << out.str();
  }
  */
}

REAL
DiscreteElements::calcTimeStep()
{
  calcVibraTimeStep();
  calcImpactTimeStep();
  calcContactNum();

  REAL CFL = 0.5;
  std::valarray<REAL> dt(3);
  dt[0] = util::getParam<REAL>("timeStep");
  dt[1] = CFL * d_vibrationTimeStep;
  dt[2] = CFL * d_impactTimeStep;

  timeStep = dt.min();
  //if (s_mpiRank == 0) {
    //std::cout << "Timestep = " << timeStep
    //          << " Vibration timestep = " << d_vibrationTimeStep
    //          << " Impact timestep = " << d_impactTimeStep << "\n";
  //}
  return timeStep;
}

void
DiscreteElements::calcVibraTimeStep()
{
  REAL pTimeStep = 1 / EPS;
  if (d_contacts.size() == 0) {
    pTimeStep = 1 / EPS;
  } else {
    auto it = d_contacts.cbegin();
    pTimeStep = it->getVibraTimeStep();
    for (++it; it != d_contacts.cend(); ++it) {
      REAL val = it->getVibraTimeStep();
      pTimeStep = val < pTimeStep ? val : pTimeStep;
    }
  }

  MPI_Allreduce(&pTimeStep, &d_vibrationTimeStep, 1, MPI_DOUBLE, MPI_MIN, s_mpiWorld);
}

void
DiscreteElements::calcImpactTimeStep()
{
  REAL pTimeStep = 1 / EPS;
  if (d_contacts.size() == 0)
    pTimeStep = 1 / EPS;
  else {
    auto it = d_contacts.cbegin();
    pTimeStep = it->getImpactTimeStep();
    for (++it; it != d_contacts.cend(); ++it) {
      REAL val = it->getImpactTimeStep();
      pTimeStep = val < pTimeStep ? val : pTimeStep;
    }
  }

  MPI_Allreduce(&pTimeStep, &d_impactTimeStep, 1, MPI_DOUBLE, MPI_MIN, s_mpiWorld);
}

void
DiscreteElements::calcContactNum()
{
  std::size_t pContactNum = d_contacts.size();
  MPI_Reduce(&pContactNum, &d_allContactNum, 1, MPI_INT, MPI_SUM, 0, s_mpiWorld);
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

  d_contacts.clear();

#ifdef TIME_PROFILE
  Timer::time_point startOuter, endOuter;
  Timer::time_point startInner, endInner, durationInner = 0;
  startOuter = Timer::now();
#endif

  auto num1 = d_patchParticles.size();      // particles inside domain
  auto num2 = d_mergedParticles.size(); // particles inside domain
                                       // (at front) + particles from
                                       // neighboring blocks (at end)
  // NOT (num1 - 1), in parallel situation where one particle
  // could contact received particles!
  for (auto i = 0; i < num1; ++i) {

    auto particle = d_patchParticles[i];
    Vec u = particle->currentPosition();
    auto particleType = particle->getType();
    auto particleRad = particle->radiusA();

    for (auto j = i + 1; j < num2; ++j) {

      auto mergeParticle = d_mergedParticles[j];
      Vec v = mergeParticle->currentPosition();
      auto mergeParticleType = mergeParticle->getType();
      auto mergeParticleRad = mergeParticle->radiusA();

      /*
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
      */

      if ((vnormL2(v - u) < particleRad + mergeParticleRad) &&
          // not both are fixed particles
          (particleType != DEMParticle::DEMParticleType::FIXED || 
           mergeParticleType != DEMParticle::DEMParticleType::FIXED) &&
          // not both are free boundary particles
          (particleType != DEMParticle::DEMParticleType::BOUNDARY_FREE || 
           mergeParticleType != DEMParticle::DEMParticleType::BOUNDARY_FREE) &&
          // not both are ghost particles
          (particleType != DEMParticle::DEMParticleType::GHOST || 
           mergeParticleType != DEMParticle::DEMParticleType::GHOST)) {

        DEMContact tmpContact(particle.get(), mergeParticle.get());

#ifdef TIME_PROFILE
        startInner = Timer::now();
#endif
        if (tmpContact.isOverlapped()) {
          d_contacts.push_back(tmpContact); // domains use value
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

  d_contacts.clear();

#ifdef TIME_PROFILE
  Timer::time_point startOuter, endOuter;
  startOuter = Timer::now();
#endif

  std::size_t i, j;
  DEMParticle::DEMParticleType particleType, mergeParticleType;
  Vec u, v;
  auto num1 = d_patchParticles.size();      // particles inside domain
  auto num2 = d_mergedParticles.size(); // particles inside domain
                                       // (at front) + particles from
                                       // neighboring blocks (at end)
#pragma omp parallel for num_threads(ompThreads) private(                      \
  i, j, u, v, particleType, mergeParticleType)                                 \
    shared(num1, num2) schedule(dynamic)
  for (i = 0; i < num1; ++i) {
    u = d_patchParticles[i]->currentPosition();
    particleType = d_patchParticles[i]->getType();
    for (j = i + 1; j < num2; ++j) {
      Vec v = d_mergedParticles[j]->currentPosition();
      mergeParticleType = d_mergedParticles[j]->getType();
      if ((vnormL2(v - u) <
           d_patchParticles[i]->radiusA() + d_mergedParticles[j]->radiusA()) &&
          // not both are fixed particles
          (particleType != DEMParticle::DEMParticleType::FIXED || 
           mergeParticleType != DEMParticle::DEMParticleType::FIXED) &&
          // not both are free boundary particles
          (particleType != DEMParticle::DEMParticleType::BOUNDARY_FREE || 
           mergeParticleType != DEMParticle::DEMParticleType::BOUNDARY_FREE) &&
          // not both are ghost particles
          (particleType != DEMParticle::DEMParticleType::GHOST || 
           mergeParticleType != DEMParticle::DEMParticleType::GHOST)) {

        DEMContact tmpContact(d_patchParticles[i].get(), d_mergedParticles[j].get());

        if (tmpContact.isOverlapped()) {
#pragma omp critical
          d_contacts.push_back(tmpContact); // domains use value
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
DiscreteElements::findBoundaryContacts()
{
  for (auto& boundary : d_boundaries) {
    boundary->findBoundaryContacts(d_patchParticles);
  }
}

void
DiscreteElements::clearContactForce()
{
  for (auto& particle : d_patchParticles) {
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
  for (auto& contact : d_contacts) {
    contact.checkinPreviousContactTangents(d_contactTangents);
  }

#ifdef TIME_PROFILE
  auto start = Timer::now();
#endif

  // d_contactTangents must be cleared before filling in new values.
  d_contactTangents.clear();

  for (auto& contact : d_contacts) {
    // cannot be parallelized as it may change a
    // particle's force simultaneously.
    contact.computeContactForces();

    // checkout current tangential force and displacment
    contact.checkoutContactTangents(d_contactTangents);

    pAvg[0] += contact.getNormalForceMagnitude();
    pAvg[1] += contact.getTangentForceMagnitude();
    pAvg[2] += contact.getPenetration();
  }

  if (d_contacts.size() != 0) {
    for (double& pVal : pAvg) {
      pVal /= d_contacts.size();
    }
  }

#ifdef TIME_PROFILE
  auto end = Timer::now();
  debugInf << std::setw(OWID) << "internalForce=" << std::setw(OWID)
           << std::chrono::duration_cast<Seconds>(end - start).count()
           << std::endl;
#endif

  MPI_Reduce(pAvg, sumAvg, 3, MPI_DOUBLE, MPI_SUM, 0, s_mpiWorld);
  d_avgNormalForce = sumAvg[0] / s_mpiSize;
  d_avgShearForce = sumAvg[1] / s_mpiSize;
  d_avgPenetration = sumAvg[2] / s_mpiSize;
}

void
DiscreteElements::boundaryForce()
{
  for (auto& boundary : d_boundaries) {
    boundary->boundaryForce(d_boundaryTangentMap);
  }
}

void
DiscreteElements::updateParticles()
{
  //proc0cout << "Num DEM particles = " << d_patchParticles.size() << "\n";
  for (auto& particle : d_patchParticles)
    particle->update();
}

void
DiscreteElements::updatePatchBox()
{
  updatePatchBoxMinX();
  updatePatchBoxMaxX();
  updatePatchBoxMinY();
  updatePatchBoxMaxY();
  updatePatchBoxMinZ();
  updatePatchBoxMaxZ();
}

void
DiscreteElements::updatePatchBoxMinX()
{
  REAL pMinX = getPtclMinX(d_patchParticles);
  REAL minX = 0;
  MPI_Allreduce(&pMinX, &minX, 1, MPI_DOUBLE, MPI_MIN, s_mpiWorld);

  setPatchBox(Box(minX - d_gradation.getPtclMaxRadius(), d_demPatchBox.minCorner().y(),
              d_demPatchBox.minCorner().z(), d_demPatchBox.maxCorner().x(),
              d_demPatchBox.maxCorner().y(), d_demPatchBox.maxCorner().z()));
}

void
DiscreteElements::updatePatchBoxMaxX()
{
  REAL pMaxX = getPtclMaxX(d_patchParticles);
  REAL maxX = 0;
  MPI_Allreduce(&pMaxX, &maxX, 1, MPI_DOUBLE, MPI_MAX, s_mpiWorld);

  setPatchBox(Box(d_demPatchBox.minCorner().x(), d_demPatchBox.minCorner().y(),
              d_demPatchBox.minCorner().z(), maxX + d_gradation.getPtclMaxRadius(),
              d_demPatchBox.maxCorner().y(), d_demPatchBox.maxCorner().z()));
}

void
DiscreteElements::updatePatchBoxMinY()
{
  REAL pMinY = getPtclMinY(d_patchParticles);
  REAL minY = 0;
  MPI_Allreduce(&pMinY, &minY, 1, MPI_DOUBLE, MPI_MIN, s_mpiWorld);

  setPatchBox(Box(d_demPatchBox.minCorner().x(), minY - d_gradation.getPtclMaxRadius(),
              d_demPatchBox.minCorner().z(), d_demPatchBox.maxCorner().x(),
              d_demPatchBox.maxCorner().y(), d_demPatchBox.maxCorner().z()));
}

void
DiscreteElements::updatePatchBoxMaxY()
{
  REAL pMaxY = getPtclMaxY(d_patchParticles);
  REAL maxY = 0;
  MPI_Allreduce(&pMaxY, &maxY, 1, MPI_DOUBLE, MPI_MAX, s_mpiWorld);

  setPatchBox(Box(d_demPatchBox.minCorner().x(), d_demPatchBox.minCorner().y(),
              d_demPatchBox.minCorner().z(), d_demPatchBox.maxCorner().x(),
              maxY + d_gradation.getPtclMaxRadius(), d_demPatchBox.maxCorner().z()));
}

void
DiscreteElements::updatePatchBoxMinZ()
{
  REAL pMinZ = getPtclMinZ(d_patchParticles);
  REAL minZ = 0;
  MPI_Allreduce(&pMinZ, &minZ, 1, MPI_DOUBLE, MPI_MIN, s_mpiWorld);

  setPatchBox(Box(d_demPatchBox.minCorner().x(), d_demPatchBox.minCorner().y(),
              minZ - d_gradation.getPtclMaxRadius(), d_demPatchBox.maxCorner().x(),
              d_demPatchBox.maxCorner().y(), d_demPatchBox.maxCorner().z()));
}

void
DiscreteElements::updatePatchBoxMaxZ()
{
  // update compute grids adaptively due to particle motion
  REAL pMaxZ = getPtclMaxZ(d_patchParticles);
  REAL maxZ = 0;
  MPI_Allreduce(&pMaxZ, &maxZ, 1, MPI_DOUBLE, MPI_MAX, s_mpiWorld);

  // no need to broadcast grid as it is updated in each process
  setPatchBox(Box(d_demPatchBox.minCorner().x(), d_demPatchBox.minCorner().y(),
              d_demPatchBox.minCorner().z(), d_demPatchBox.maxCorner().x(),
              d_demPatchBox.maxCorner().y(), maxZ + d_gradation.getPtclMaxRadius()));
}

void
DiscreteElements::findParticleInBox(const Box& domain,
                            const DEMParticlePArray& allParticles,
                            DEMParticlePArray& insideParticles)
{
  insideParticles.clear();
  for (const auto& particle : allParticles) {
    // it is critical to use EPS
    if (domain.inside(particle->currentPosition(), EPS)) {
      insideParticles.push_back(particle);
    }
  }
}

void
DiscreteElements::writeBoundaryToFile() const
{
  d_writer->writeDomain(&d_spatialDomain);
}

void
DiscreteElements::printBoundary() const
{
  std::ofstream ofs(d_writer->getBoundaryFilename());
  if (!ofs) {
    debugInf << "stream error: printBoundary" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);

  Vec v1 = d_spatialDomain.minCorner();
  Vec v2 = d_spatialDomain.maxCorner();
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
      << std::setw(OWID) << d_boundaries.size() << std::endl;

  for (const auto& it : d_boundaries) {
    it->print(ofs);
  }

  ofs.close();
}

void
DiscreteElements::printBoundaryContacts() const
{
  std::ofstream ofs(d_writer->getBdryContactFilename());
  if (!ofs) {
    debugInf << "stream error: printBoundaryContacts" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  for (const auto& it : d_mergedBoundaries) {
    it->printContactInfo(ofs);
  }

  ofs.close();
}

void
DiscreteElements::writePatchGridToFile() const
{
  d_writer->setMPIComm(s_cartComm);
  d_writer->setMPIProc(s_mpiProcs.x(), s_mpiProcs.y(), s_mpiProcs.z());
  //std::cout << "d_demPatchBox = " << d_demPatchBox << "\n";
  d_writer->writePatchBoxGrid(&d_demPatchBox);
}

void
DiscreteElements::writeParticlesToFile(int frame) const
{
  d_writer->writeParticles(&d_allDEMParticles, frame);
  d_writer->writeSieves(&d_gradation);
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
  writer.setParticleFilename(fileName);
  std::cout << "d_allDEMParticles = " << d_allDEMParticles.size() << "\n";
  writer.writeParticles(&d_allDEMParticles, frame);
}

void
DiscreteElements::printParticle(const std::string& fileName, 
                                DEMParticlePArray& particles,
                                int frame) const
{
  OutputTecplot<DEMParticlePArray> writer(".", 0);
  writer.setParticleFilename(fileName);
  writer.writeParticles(&particles, frame);
}

bool
DiscreteElements::areBoundaryTractionsEquilibrated(REAL time) 
{
  REAL tol = util::getParam<REAL>("tractionErrorTol");
  for (const auto& boundary : d_boundaries) {
    double area = boundary->getArea();
    Vec normalForce = boundary->getNormalForce();
    Vec appliedForce = boundary->getAppliedTraction(time)*area;
    REAL resultant = (normalForce - appliedForce).length();
    if (resultant > tol) return false;
  }
  return true;
}

bool
DiscreteElements::areBoundaryTractionsEquilibrated(REAL sigma, std::string type, REAL sigmaX,
                           REAL sigmaY)
{
  // sigma implies sigmaZ
  REAL tol = util::getParam<REAL>("tractionErrorTol");

  std::map<std::string, REAL> normalTraction;
  REAL x1, x2, y1, y2, z1, z2;

  // do not use d_mergedBoundaries because each process calls this function.
  for (const auto& boundary : d_boundaries) {
    Boundary::BoundaryID id = boundary->getId();
    Vec force = boundary->getNormalForce();
    Vec point = boundary->getPosition();
    switch (id) {
      case Boundary::BoundaryID::XMINUS:
        normalTraction["x-"] = fabs(force.x());
        x1 = point.x();
        break;
      case Boundary::BoundaryID::XPLUS:
        normalTraction["x+"] = force.x();
        x2 = point.x();
        break;
      case Boundary::BoundaryID::YMINUS:
        normalTraction["y-"] = fabs(force.y());
        y1 = point.y();
        break;
      case Boundary::BoundaryID::YPLUS:
        normalTraction["y+"] = force.y();
        y2 = point.y();
        break;
      case Boundary::BoundaryID::ZMINUS:
        normalTraction["z-"] = fabs(force.z());
        z1 = point.z();
        break;
      case Boundary::BoundaryID::ZPLUS:
        normalTraction["z+"] = force.z();
        z2 = point.z();
        break;
      default:
        break;
    }
  }
  REAL areaX = (y2 - y1) * (z2 - z1);
  REAL areaY = (z2 - z1) * (x2 - x1);
  REAL areaZ = (x2 - x1) * (y2 - y1);
  normalTraction["x-"] /= areaX;
  normalTraction["x+"] /= areaX;
  normalTraction["y-"] /= areaY;
  normalTraction["y+"] /= areaY;
  normalTraction["z-"] /= areaZ;
  normalTraction["z+"] /= areaZ;

  if (type.compare("isotropic") == 0)
    return (fabs(normalTraction["x-"]/sigma - 1) <= tol &&
            fabs(normalTraction["x+"]/sigma - 1) <= tol &&
            fabs(normalTraction["y-"]/sigma - 1) <= tol &&
            fabs(normalTraction["y+"]/sigma - 1) <= tol &&
            fabs(normalTraction["z-"]/sigma - 1) <= tol &&
            fabs(normalTraction["z+"]/sigma - 1) <= tol);

  else if (type.compare("odometer") == 0)
    return (fabs(normalTraction["z-"]/sigma - 1) <= tol &&
            fabs(normalTraction["z+"]/sigma - 1) <= tol);

  else if (type.compare("triaxial") == 0)
    return true; // always near equilibrium

  else if (type.compare("trueTriaxial") == 0)
    return (fabs(normalTraction["x-"]/sigmaX - 1) <= tol &&
            fabs(normalTraction["x+"]/sigmaX - 1) <= tol &&
            fabs(normalTraction["y-"]/sigmaY - 1) <= tol &&
            fabs(normalTraction["y+"]/sigmaY - 1) <= tol &&
            fabs(normalTraction["z-"]/sigma - 1) <= tol &&
            fabs(normalTraction["z+"]/sigma - 1) <= tol);

  return false;
}

void
DiscreteElements::trim(bool toRebuild, const std::string& inputParticle,
               const std::string& trmParticle)
{
  if (toRebuild) {
    readParticles(inputParticle);
  }

  d_trimHistoryNum = d_allDEMParticles.size();

  Vec v1 = d_spatialDomain.minCorner();
  Vec v2 = d_spatialDomain.maxCorner();
  REAL x1 = v1.x();
  REAL y1 = v1.y();
  REAL z1 = v1.z();
  REAL x2 = v2.x();
  REAL y2 = v2.y();
  REAL z2 = v2.z();
  REAL maxR = d_gradation.getPtclMaxRadius();

  // BB: Feb 2, 2017:
  // Not an efficient operation
  // Better approach may be to use a list if random access of vector
  // members is not needed
  d_allDEMParticles.erase(
    std::remove_if(d_allDEMParticles.begin(), d_allDEMParticles.end(),
                   [&x1, &y1, &z1, &x2, &y2, &z2, &maxR](DEMParticleP particle) {
                     Vec center = particle->currentPosition();
                     if (center.x() < x1 || center.x() > x2 ||
                         center.y() < y1 || center.y() > y2 ||
                         center.z() < z1 || center.z() + maxR > z2) {
                       return true;
                     }
                     return false;
                   }),
    d_allDEMParticles.end());

  DEMParticleFileWriter writer;
  writer.writeCSV(d_allDEMParticles, d_gradation, trmParticle+".csv");
  writer.writeXML(d_allDEMParticles, d_gradation, trmParticle+".xml");

  // Also create a VTK output for viz
  // Create the output writer in the master process
  auto folderName =  dem::InputParameter::get().datafile["outputFolder"]+"_trim";
  std::string outputFolder = util::createOutputFolder(folderName);
  //std::cout << "Output folder = " << outputFolder << "\n";
  createOutputWriter(outputFolder, 0);

  writeBoundaryToFile();
  writeParticlesToFile(0);
}

void
DiscreteElements::removeParticleOutBox()
{
  Vec v1 = d_patchDomain.minCorner();
  Vec v2 = d_patchDomain.maxCorner();
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
  d_patchParticles.erase(
    std::remove_if(
      d_patchParticles.begin(), d_patchParticles.end(),
      [&x1, &y1, &z1, &x2, &y2, &z2, &epsilon](DEMParticleP particle) {
        Vec center = particle->currentPosition();
        if (!(center.x() - x1 >= -epsilon && center.x() - x2 < -epsilon &&
              center.y() - y1 >= -epsilon && center.y() - y2 < -epsilon &&
              center.z() - z1 >= -epsilon && center.z() - z2 < -epsilon)) {
          /*
          if (particle->getId() == 2 || particle->getId() == 94) {
            //std::cout << "**WARNING** Removing particle " << particle->getId()
                      << " from domain with \n "
                      << " x : " << center.x() << " not in [" << x1 << ","  << x2 << "]\n"
                      << " y : " << center.y() << " not in [" << y1 << ","  << y2 << "]\n"
                      << " z : " << center.z() << " not in [" << z1 << ","  << z2 << "]\n";
          }
          */
          return true;
        }
        return false;
      }),
    d_patchParticles.end());

  /*
  DEMParticlePArray::iterator itr;
  Vec center;
  //std::size_t flag = 0;

  for (itr = d_patchParticles.begin(); itr != d_patchParticles.end(); ) {
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
    itr = d_patchParticles.erase(itr);
  }
    else
  ++itr;
  }
  */

  /*
    if (flag == 1) {
    debugInf << " now " << d_patchParticles.size() << ": ";
    for (DEMParticlePArray::const_iterator it = d_patchParticles.begin(); it !=
    d_patchParticles.end(); ++it)
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
DiscreteElements::isBoundaryProcess()
{
  return (s_mpiCoords.x() == 0 || s_mpiCoords.x() == s_mpiProcs.x() - 1 ||
          s_mpiCoords.y() == 0 || s_mpiCoords.y() == s_mpiProcs.y() - 1 ||
          s_mpiCoords.z() == 0 || s_mpiCoords.z() == s_mpiProcs.z() - 1);
}


void
DiscreteElements::releaseReceivedParticles()
{
  // release memory of received particles
  d_receivedParticles.clear();
}

void
DiscreteElements::migrateParticles()
{
  //std::ostringstream out;
  //out << "Migrate: Rank: " << s_mpiRank << ": in: " << d_patchParticles.size();

  Vec vspan = d_demPatchBox.maxCorner() - d_demPatchBox.minCorner();
  Vec width = vspan / s_mpiProcs;

  d_sentParticles.clear();
  d_receivedParticles.clear();
  d_patchP->sendRecvMigrateXMinus(boostWorld, iteration, width, d_patchParticles);
  d_patchP->sendRecvMigrateXPlus(boostWorld, iteration, width, d_patchParticles);
  d_patchP->waitToFinishX(iteration);
  d_patchP->combineSentParticlesX(iteration, d_sentParticles);
  d_patchP->combineReceivedParticlesX(iteration, d_receivedParticles);
  d_patchP->deleteSentParticles<DEMParticleP>(iteration, d_sentParticles, d_patchParticles);
  d_patchP->addReceivedParticles(iteration, d_receivedParticles, d_patchParticles);
  //out << " sentX : " << d_sentParticles.size()
  //    << " recvX : " << d_receivedParticles.size();

  d_sentParticles.clear();
  d_receivedParticles.clear();
  d_patchP->sendRecvMigrateYMinus(boostWorld, iteration, width, d_patchParticles);
  d_patchP->sendRecvMigrateYPlus(boostWorld, iteration, width, d_patchParticles);
  d_patchP->waitToFinishY(iteration);
  d_patchP->combineSentParticlesY(iteration, d_sentParticles);
  d_patchP->combineReceivedParticlesY(iteration, d_receivedParticles);
  d_patchP->deleteSentParticles<DEMParticleP>(iteration, d_sentParticles, d_patchParticles);
  d_patchP->addReceivedParticles(iteration, d_receivedParticles, d_patchParticles);
  //out << " sentY : " << d_sentParticles.size()
  //    << " recvY : " << d_receivedParticles.size();

  d_sentParticles.clear();
  d_receivedParticles.clear();
  d_patchP->sendRecvMigrateZMinus(boostWorld, iteration, width, d_patchParticles);
  d_patchP->sendRecvMigrateZPlus(boostWorld, iteration, width, d_patchParticles);
  d_patchP->waitToFinishZ(iteration);
  d_patchP->combineSentParticlesZ(iteration, d_sentParticles);
  d_patchP->combineReceivedParticlesZ(iteration, d_receivedParticles);
  d_patchP->deleteSentParticles<DEMParticleP>(iteration, d_sentParticles, d_patchParticles);
  d_patchP->addReceivedParticles(iteration, d_receivedParticles, d_patchParticles);
  //out << " sentZ : " << d_sentParticles.size()
  //    << " recvZ : " << d_receivedParticles.size();

  // delete outgoing particles
  d_patchP->removeParticlesOutsidePatch<DEMParticleP>(d_patchParticles);
  //out << " outside: " << d_patchParticles.size();

  //d_patchP->removeDuplicates<DEMParticleP>(d_patchParticles);
  //out <<  ": dup out: " << d_patchParticles.size() << "\n";
  //std::cout << out.str();

}


void
DiscreteElements::gatherParticles()
{
  // update d_allDEMParticles: process 0 collects all updated particles from each
  // other process
  if (s_mpiRank != 0) { // each process except 0
    boostWorld.send(0, mpiTag, d_patchParticles);
  } else { // process 0
    // d_allDEMParticles is cleared before filling with new data
    releaseGatheredParticle();

    // duplicate d_patchParticles so that it is not destroyed by d_allDEMParticles in
    // next iteration,
    // otherwise it causes memory error.
    DEMParticlePArray dupParticleVec(d_patchParticles.size());
    for (std::size_t i = 0; i < dupParticleVec.size(); ++i)
      dupParticleVec[i] = std::make_shared<DEMParticle>(*d_patchParticles[i]);

    // fill d_allDEMParticles with dupParticleVec and received particles
    d_allDEMParticles.insert(d_allDEMParticles.end(), dupParticleVec.begin(),
                          dupParticleVec.end());

    DEMParticlePArray tmpParticleVec;
    long gatherRam = 0;
    for (int iRank = 1; iRank < s_mpiSize; ++iRank) {
      tmpParticleVec.clear(); // do not destroy particles!
      boostWorld.recv(iRank, mpiTag, tmpParticleVec);
      d_allDEMParticles.insert(d_allDEMParticles.end(), tmpParticleVec.begin(),
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
  // clear d_allDEMParticles, avoid long time memory footprint.
  /*
  for (DEMParticlePArray::iterator it = d_allDEMParticles.begin(); it !=
  d_allDEMParticles.end(); ++it)
    delete (*it);
  */
  d_allDEMParticles.clear();
  DEMParticlePArray().swap(d_allDEMParticles); // actual memory release
}


void
DiscreteElements::gatherBoundaryContacts()
{
  if (isBoundaryProcess()) {
    if (s_mpiRank != 0)
      boostWorld.send(0, mpiTag, d_boundaries);
  }

  if (s_mpiRank == 0) {
    d_mergedBoundaries.clear();
    BoundaryPArray().swap(d_mergedBoundaries); // actual memory release
    d_mergedBoundaries = d_boundaries;

    BoundaryPArray tmpBoundaryVec;
    for (unsigned long bdryProces : bdryProcess) {
      if (bdryProces != 0) {    // not root process
        tmpBoundaryVec.clear(); // do not destroy particles!
        boostWorld.recv(bdryProces, mpiTag, tmpBoundaryVec);
        // merge tmpBoundaryVec into d_mergedBoundaries
        assert(tmpBoundaryVec.size() == d_mergedBoundaries.size());
        for (std::size_t jt = 0; jt < tmpBoundaryVec.size(); ++jt)
          d_mergedBoundaries[jt]->getBoundaryContacts().insert(
            d_mergedBoundaries[jt]->getBoundaryContacts().end(),
            tmpBoundaryVec[jt]->getBoundaryContacts().begin(),
            tmpBoundaryVec[jt]->getBoundaryContacts().end());
      }
    }

    // must update after collecting all boundary contact info
    for (auto& it : d_mergedBoundaries)
      it->updateStatForce();
  }
}

void
DiscreteElements::gatherEnergy()
{
  calcTranslationalEnergy();
  calcRotationalEnergy();
  calcKineticEnergy();
  calcGravitationalEnergy(d_spatialDomain.minCorner().z());
  calcMechanicalEnergy();
}

void
DiscreteElements::getStartDimension(REAL& distX, REAL& distY, REAL& distZ)
{
  REAL x1 = 0, x2 = 0, y1 = 0, y2 = 0, z1 = 0, z2 = 0;
  // use d_boundaries
  for (const auto& boundary : d_boundaries) {
    switch (boundary->getId()) {
      case Boundary::BoundaryID::XMINUS:
        x1 = boundary->getPosition().x();
        break;
      case Boundary::BoundaryID::XPLUS:
        x2 = boundary->getPosition().x();
        break;
      case Boundary::BoundaryID::YMINUS:
        y1 = boundary->getPosition().y();
        break;
      case Boundary::BoundaryID::YPLUS:
        y2 = boundary->getPosition().y();
        break;
      case Boundary::BoundaryID::ZMINUS:
        z1 = boundary->getPosition().z();
        break;
      case Boundary::BoundaryID::ZPLUS:
        z2 = boundary->getPosition().z();
        break;
      default:
        break;
    }
  }
  distX = x2 - x1;
  distY = y2 - y1;
  distZ = z2 - z1;
}

void
DiscreteElements::openProgressOutputFile(std::ofstream& ofs, const std::string& str)
{
  ofs.open(str);
  if (!ofs) {
    debugInf << "stream error: openProgressOutputFile" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  ofs << std::setw(OWID) << "iteration" << std::setw(OWID) << "normal_x1"
      << std::setw(OWID) << "normal_x2" << std::setw(OWID) << "normal_y1"
      << std::setw(OWID) << "normal_y2" << std::setw(OWID) << "normal_z1"
      << std::setw(OWID) << "normal_z2"
  
      << std::setw(OWID) << "traction_x1"
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

      << std::setw(OWID) << "transEnergy" << std::setw(OWID) << "rotEnergy"
      << std::setw(OWID) << "kinEnergy" << std::setw(OWID) << "gravEnergy"
      << std::setw(OWID) << "mechEnergy"

      << std::setw(OWID) << "vibra_est_dt" << std::setw(OWID) << "impact_est_dt"
      << std::setw(OWID) << "actual_dt"

      << std::endl;
}

void
DiscreteElements::appendToProgressOutputFile(std::ofstream& ofs, REAL distX, REAL distY,
                            REAL distZ)
{
  REAL x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0, z1 = 0.0, z2 = 0.0;
  std::vector<REAL> normalForce = {{0, 0, 0, 0, 0, 0}};
  std::vector<REAL> normalVelocity = {{0, 0, 0, 0, 0, 0}};

  for (const auto& boundary : d_mergedBoundaries) {
    Boundary::BoundaryID id = boundary->getId();
    switch (id) {
      case Boundary::BoundaryID::XMINUS:
        x1 = boundary->getPosition().x();
        normalForce[0] = fabs(boundary->getNormalForce().x());
        normalVelocity[0] = boundary->getVelocity().x();
        break;
      case Boundary::BoundaryID::XPLUS:
        x2 = boundary->getPosition().x();
        normalForce[1] = boundary->getNormalForce().x();
        normalVelocity[1] = boundary->getVelocity().x();
        break;
      case Boundary::BoundaryID::YMINUS:
        y1 = boundary->getPosition().y();
        normalForce[2] = fabs(boundary->getNormalForce().y());
        normalVelocity[2] = boundary->getVelocity().y();
        break;
      case Boundary::BoundaryID::YPLUS:
        y2 = boundary->getPosition().y();
        normalForce[3] = boundary->getNormalForce().y();
        normalVelocity[3] = boundary->getVelocity().y();
        break;
      case Boundary::BoundaryID::ZMINUS:
        z1 = boundary->getPosition().z();
        normalForce[4] = fabs(boundary->getNormalForce().z());
        normalVelocity[4] = boundary->getVelocity().z();
        break;
      case Boundary::BoundaryID::ZPLUS:
        z2 = boundary->getPosition().z();
        normalForce[5] = boundary->getNormalForce().z();
        normalVelocity[5] = boundary->getVelocity().z();
        break;
      default:
        break;
    }
  }

  REAL areaX = (y2 - y1) * (z2 - z1);
  REAL areaY = (z2 - z1) * (x2 - x1);
  REAL areaZ = (x2 - x1) * (y2 - y1);
  REAL bulkVolume = (x2 - x1) * (y2 - y1) * (z2 - z1);
  REAL voidRatio = bulkVolume / volume() - 1;

  // normal traction
  std::vector<REAL> normalTraction = {{0, 0, 0, 0, 0, 0}};
  normalTraction[0] = normalForce[0]/areaX;
  normalTraction[1] = normalForce[1]/areaX;
  normalTraction[2] = normalForce[2]/areaY;
  normalTraction[3] = normalForce[3]/areaY;
  normalTraction[4] = normalForce[4]/areaZ;
  normalTraction[5] = normalForce[5]/areaZ;

  // Write normal force
  ofs << std::setw(OWID) << iteration;
  for (auto force : normalForce)
    ofs << std::setw(OWID) << force;

  // Write normal traction
  REAL avgTraction = 0;
  for (double traction : normalTraction) {
    ofs << std::setw(OWID) << traction;
    avgTraction += traction;
  }
  ofs << std::setw(OWID) << avgTraction / 6;

  // volume
  ofs << std::setw(OWID) << bulkVolume << std::setw(OWID)
      << mass() / bulkVolume << std::setw(OWID) << 1 - (x2 - x1) / distX
      << std::setw(OWID) << 1 - (y2 - y1) / distY << std::setw(OWID)
      << 1 - (z2 - z1) / distZ << std::setw(OWID)
      << 3 - (x2 - x1) / distX - (y2 - y1) / distY - (z2 - z1) / distZ
      << std::setw(OWID) << voidRatio << std::setw(OWID)
      << voidRatio / (1 + voidRatio);

  // velocity
  for (double velocity : normalVelocity)
    ofs << std::setw(OWID) << velocity;

  // b_numContacts
  std::vector<REAL> data = {{0, 0, 0, 0, 0, 0}};
  for (const auto& boundary : d_mergedBoundaries) {
    Boundary::BoundaryID id = boundary->getId();
    data[static_cast<size_t>(id) - 1] = boundary->getNumBoundaryContacts();
  }
  for (double datum : data)
    ofs << std::setw(OWID) << static_cast<std::size_t>(datum);
  ofs << std::setw(OWID) << d_allContactNum;

  // avgPenetr
  data.clear();
  data.reserve(6);
  for (const auto& boundary : d_mergedBoundaries) {
    auto id = boundary->getId();
    data[static_cast<size_t>(id) - 1] = boundary->getAvgPenetration();
  }
  for (double datum : data)
    ofs << std::setw(OWID) << datum;

  // average data
  ofs << std::setw(OWID) << d_avgNormalForce << std::setw(OWID) << d_avgShearForce
      << std::setw(OWID) << d_avgPenetration;

  // energy
  ofs << std::setw(OWID) << d_translationalEnergy << std::setw(OWID) << d_rotationalEnergy
      << std::setw(OWID) << d_kineticEnergy << std::setw(OWID) << d_gravitationalEnergy
      << std::setw(OWID) << d_mechanicalEnergy;

  // time
  ofs << std::setw(OWID) << d_vibrationTimeStep << std::setw(OWID) << d_impactTimeStep
      << std::setw(OWID) << timeStep;

  ofs << std::endl;
}

void
DiscreteElements::closeProgressOutputFile(std::ofstream& ofs)
{
  ofs.close();
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
DiscreteElements::updateBoundary(REAL time, REAL delT, REAL mass)
{
  if (s_mpiRank == 0) {
    for (auto& boundary : d_mergedBoundaries) {
      double area = boundary->getArea();
      boundary->updatePositionAndVelocity(time, delT, area, mass);
    }
    updateBoundaryAreas(d_mergedBoundaries);

    // update d_boundaries from d_mergedBoundaries and remove b_contacts to reduce
    // MPI transmission
    d_boundaries = d_mergedBoundaries;
    for (auto& it : d_boundaries)
      it->clearBoundaryContacts();

    // update d_spatialDomain
    double xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0;
    for (const auto& boundary : d_boundaries) {
      switch (boundary->getId()) {
        case Boundary::BoundaryID::XMINUS:
          xmin = boundary->getPosition().x();
          break;
        case Boundary::BoundaryID::XPLUS:
          xmax = boundary->getPosition().x();
          break;
        case Boundary::BoundaryID::YMINUS:
          ymin = boundary->getPosition().y();
          break;
        case Boundary::BoundaryID::YPLUS:
          ymax = boundary->getPosition().y();
          break;
        case Boundary::BoundaryID::ZMINUS:
          zmin = boundary->getPosition().z();
          break;
        case Boundary::BoundaryID::ZPLUS:
          zmax = boundary->getPosition().z();
          break;
        default:
          break;
      }
    }
    setSpatialDomain(Box(xmin, ymin, zmin, xmax, ymax, zmax));
  }

  broadcast(boostWorld, d_boundaries, 0);
}

void
DiscreteElements::updateBoundary(REAL sigma, std::string type, 
                                 REAL sigmaX, REAL sigmaY)
{
  if (s_mpiRank == 0) {
    REAL x1, x2, y1, y2, z1, z2;
    for (const auto& boundary : d_mergedBoundaries) {
      switch (boundary->getId()) {
        case Boundary::BoundaryID::XMINUS:
          x1 = boundary->getPosition().x();
          break;
        case Boundary::BoundaryID::XPLUS:
          x2 = boundary->getPosition().x();
          break;
        case Boundary::BoundaryID::YMINUS:
          y1 = boundary->getPosition().y();
          break;
        case Boundary::BoundaryID::YPLUS:
          y2 = boundary->getPosition().y();
          break;
        case Boundary::BoundaryID::ZMINUS:
          z1 = boundary->getPosition().z();
          break;
        case Boundary::BoundaryID::ZPLUS:
          z2 = boundary->getPosition().z();
          break;
        default:
          break;
      }
    }
    REAL areaX = (y2 - y1) * (z2 - z1);
    REAL areaY = (z2 - z1) * (x2 - x1);
    REAL areaZ = (x2 - x1) * (y2 - y1);

    if (type.compare("isotropic") == 0) {
      for (auto& it : d_mergedBoundaries)
        it->updateIsotropic(sigma, areaX, areaY, areaZ);
    } else if (type.compare("odometer") == 0) {
      for (auto& it : d_mergedBoundaries)
        it->updateOdometer(sigma, areaX, areaY, areaZ);
    } else if (type.compare("triaxial") == 0) {
      for (auto& it : d_mergedBoundaries)
        it->updateTriaxial(sigma, areaX, areaY, areaZ);
    } else if (type.compare("plnstrn") == 0) {
      for (auto& it : d_mergedBoundaries)
        it->updatePlaneStrain(sigma, areaX, areaY, areaZ);
    } else if (type.compare("trueTriaxial") == 0) {
      for (auto& it : d_mergedBoundaries)
        it->updateTrueTriaxial(sigma, areaX, areaY, areaZ, sigmaX, sigmaY);
    } 

    // update d_boundaries from d_mergedBoundaries and remove b_contacts to reduce
    // MPI transmission
    d_boundaries = d_mergedBoundaries;
    for (auto& it : d_boundaries)
      it->clearBoundaryContacts();

    // update d_spatialDomain
    for (const auto& boundary : d_boundaries) {
      switch (boundary->getId()) {
        case Boundary::BoundaryID::XMINUS:
          x1 = boundary->getPosition().x();
          break;
        case Boundary::BoundaryID::XPLUS:
          x2 = boundary->getPosition().x();
          break;
        case Boundary::BoundaryID::YMINUS:
          y1 = boundary->getPosition().y();
          break;
        case Boundary::BoundaryID::YPLUS:
          y2 = boundary->getPosition().y();
          break;
        case Boundary::BoundaryID::ZMINUS:
          z1 = boundary->getPosition().z();
          break;
        case Boundary::BoundaryID::ZPLUS:
          z2 = boundary->getPosition().z();
          break;
        default:
          break;
      }
    }
    setSpatialDomain(Box(x1, y1, z1, x2, y2, z2));
  }

  broadcast(boostWorld, d_boundaries, 0);
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

  for (const auto& it : d_contacts)
    inf << std::setw(OWID) << it.getP1()->getId() << std::setw(OWID)
        << it.getP2()->getId() << std::setw(OWID) << it.getPoint1().x()
        << std::setw(OWID) << it.getPoint1().y() << std::setw(OWID)
        << it.getPoint1().z() << std::setw(OWID) << it.getPoint2().x()
        << std::setw(OWID) << it.getPoint2().y() << std::setw(OWID)
        << it.getPoint2().z() << std::setw(OWID) << it.radius1()
        << std::setw(OWID) << it.radius2() << std::setw(OWID)
        << it.getPenetration() << std::setw(OWID) << it.getTangentDisplacement()
        << std::setw(OWID) << it.getContactRadius() << std::setw(OWID)
        << it.getR0() << std::setw(OWID) << it.getE0() << std::setw(OWID)
        << it.getNormalForceMagnitude() << std::setw(OWID) << it.getTangentForceMagnitude()
        << std::setw(OWID) << (it.getPoint1().x() + it.getPoint2().x()) / 2
        << std::setw(OWID) << (it.getPoint1().y() + it.getPoint2().y()) / 2
        << std::setw(OWID) << (it.getPoint1().z() + it.getPoint2().z()) / 2
        << std::setw(OWID) << it.getNormalForce().x() << std::setw(OWID)
        << it.getNormalForce().y() << std::setw(OWID) << it.getNormalForce().z()
        << std::setw(OWID) << it.getTangentForce().x() << std::setw(OWID)
        << it.getTangentForce().y() << std::setw(OWID) << it.getTangentForce().z()
        << std::setw(OWID) << it.getVibraTimeStep() << std::setw(OWID)
        << it.getImpactTimeStep() << std::endl;

  int length = (OWID * 28 + 1) * d_contacts.size();
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

  ofs << std::setw(OWID) << d_contacts.size() << std::endl;
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
  for (it = d_contacts.begin(); it != d_contacts.end(); ++it)
    ofs << std::setw(OWID) << it->getP1()->getId()
    << std::setw(OWID) << it->getP2()->getId()
    << std::setw(OWID) << it->getPoint1().x()
    << std::setw(OWID) << it->getPoint1().y()
    << std::setw(OWID) << it->getPoint1().z()
    << std::setw(OWID) << it->getPoint2().x()
    << std::setw(OWID) << it->getPoint2().y()
    << std::setw(OWID) << it->getPoint2().z()
    << std::setw(OWID) << it->radius1()
    << std::setw(OWID) << it->radius2()
    << std::setw(OWID) << it->getPenetration()
    << std::setw(OWID) << it->getTangentDisplacement()
    << std::setw(OWID) << it->getContactRadius()
    << std::setw(OWID) << it->getR0()
    << std::setw(OWID) << it->getE0()
    << std::setw(OWID) << it->getNormalForceMagnitude()
    << std::setw(OWID) << it->getTangentForceMagnitude()
    << std::setw(OWID) << ( it->getPoint1().x() + it->getPoint2().x() )/2
    << std::setw(OWID) << ( it->getPoint1().y() + it->getPoint2().y() )/2
    << std::setw(OWID) << ( it->getPoint1().z() + it->getPoint2().z() )/2
    << std::setw(OWID) << it->getNormalForce().x()
    << std::setw(OWID) << it->getNormalForce().y()
    << std::setw(OWID) << it->getNormalForce().z()
    << std::setw(OWID) << it->getTangentForce().x()
    << std::setw(OWID) << it->getTangentForce().y()
    << std::setw(OWID) << it->getTangentForce().z()
    << std::setw(OWID) << it->getVibraTimeStep()
    << std::setw(OWID) << it->getImpactTimeStep()
    << std::endl;
  ofs.close();
  */
}

void
DiscreteElements::calcTranslationalEnergy()
{
  REAL pEnergy = 0;
  for (const auto& particle : d_patchParticles) {
    if (particle->getType() == DEMParticle::DEMParticleType::FREE) {
      pEnergy += particle->computeTranslationalEnergy();
    }
  }
  MPI_Reduce(&pEnergy, &d_translationalEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, s_mpiWorld);
}

void
DiscreteElements::calcRotationalEnergy()
{
  REAL pEnergy = 0;
  for (const auto& particle : d_patchParticles) {
    if (particle->getType() == DEMParticle::DEMParticleType::FREE) {
      pEnergy += particle->computeRotationalEnergy();
    }
  }
  MPI_Reduce(&pEnergy, &d_rotationalEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, s_mpiWorld);
}

void
DiscreteElements::calcKineticEnergy()
{
  REAL pEnergy = 0;
  for (const auto& particle : d_patchParticles) {
    if (particle->getType() == DEMParticle::DEMParticleType::FREE) {
      pEnergy += particle->computeKineticEnergy();
    }
  }
  MPI_Reduce(&pEnergy, &d_kineticEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, s_mpiWorld);
}

void
DiscreteElements::calcGravitationalEnergy(REAL ref)
{
  REAL pEnergy = 0;
  for (const auto& particle : d_patchParticles) {
    if (particle->getType() == DEMParticle::DEMParticleType::FREE) {
      pEnergy += particle->computePotentialEnergy(ref);
    }
  }
  MPI_Reduce(&pEnergy, &d_gravitationalEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, s_mpiWorld);
}

void
DiscreteElements::calcMechanicalEnergy()
{
  d_mechanicalEnergy = d_kineticEnergy + d_gravitationalEnergy;
}

REAL 
DiscreteElements::getTotalMassFromPatchParticleData() const
{
  return d_patchParticles[0]->getTotalMass();
}

REAL
DiscreteElements::mass() const
{
  REAL var = 0;
  for (const auto& particle : d_allDEMParticles)
    var += particle->mass();
  return var;
}

REAL
DiscreteElements::volume() const
{
  REAL var = 0;
  for (const auto& particle : d_allDEMParticles)
    if (particle->getType() == DEMParticle::DEMParticleType::FREE)
      var += particle->volume();
  return var;
}

REAL
DiscreteElements::getAvgTranslationalVelocity() const
{
  REAL avgv = 0;
  std::size_t count = 0;
  for (const auto& particle : d_patchParticles) {
    if (particle->getType() == DEMParticle::DEMParticleType::FREE) {
      avgv += vnormL2(particle->currentVelocity());
      ++count;
    }
  }
  return avgv /= count;
}

REAL
DiscreteElements::getAvgRotationalVelocity() const
{
  REAL avgv = 0;
  std::size_t count = 0;
  for (const auto& particle : d_patchParticles) {
    if (particle->getType() == DEMParticle::DEMParticleType::FREE) {
      avgv += vnormL2(particle->currentAngularVelocity());
      ++count;
    }
  }
  return avgv /= count;
}

REAL
DiscreteElements::getAvgForce() const
{
  REAL avgv = 0;
  std::size_t count = 0;
  for (const auto& particle : d_patchParticles) {
    if (particle->getType() == DEMParticle::DEMParticleType::FREE) {
      avgv += vnormL2(particle->force());
      ++count;
    }
  }
  return avgv / count;
}

REAL
DiscreteElements::getAvgMoment() const
{
  REAL avgv = 0;
  std::size_t count = 0;
  for (const auto& particle : d_patchParticles) {
    if (particle->getType() == DEMParticle::DEMParticleType::FREE) {
      avgv += vnormL2(particle->moment());
      ++count;
    }
  }
  return avgv /= count;
}

void
DiscreteElements::dragForce() 
{
  for (auto& particle : d_patchParticles) {
    particle->dragForce();
  }
}

} // namespace dem ends
