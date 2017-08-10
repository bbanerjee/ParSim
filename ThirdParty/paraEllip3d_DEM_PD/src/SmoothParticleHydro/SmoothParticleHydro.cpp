#include <DiscreteElements/DEMParticle.h>
#include <SmoothParticleHydro/SmoothParticleHydro.h>

#include <Core/Const/const.h>
#include <Core/Math/IntVec.h>
#include <Core/Math/Matrix.h>
#include <Core/Util/Utility.h>
#include <InputOutput/OutputTecplot.h>
#include <InputOutput/OutputVTK.h>
#include <chrono>

using namespace sph;

using Timer = std::chrono::steady_clock;
using Seconds = std::chrono::seconds;
using IntVec = dem::IntVec;
using Vec = dem::Vec;
using Matrix = dem::Matrix;
using Box = dem::Box;
using DEMParticlePArray = dem::DEMParticlePArray;
using OutputVTK = dem::OutputVTK<SPHParticlePArray>;
using OutputTecplot = dem::OutputTecplot<SPHParticlePArray>;
using communicator = boost::mpi::communicator;

SmoothParticleHydro::SmoothParticleHydro()
{
}

SmoothParticleHydro::~SmoothParticleHydro()
{
  d_allSPHParticleVec.clear();
  d_sphParticleVec.clear();
}

void 
SmoothParticleHydro::setCommunicator(const communicator& boostWorldComm)
{
  d_boostWorld = boostWorldComm;
  d_mpiWorld = MPI_Comm(boostWorldComm);
  int mpiProcX = util::getParam<int>("mpiProcX");
  int mpiProcY = util::getParam<int>("mpiProcY");
  int mpiProcZ = util::getParam<int>("mpiProcZ");
  d_mpiProcs = {{ mpiProcX, mpiProcY, mpiProcZ }};

  // create Cartesian virtual topology (unavailable in boost.mpi)
  int ndim = 3;
  int periods[3] = { 0, 0, 0 };
  int reorder = 0; // d_mpiRank not reordered
  int status = MPI_Cart_create(d_mpiWorld, ndim, d_mpiProcs.data(), 
                               periods, reorder, &d_cartComm);
  if (status != MPI_SUCCESS ) {
    std::cout << "**ERROR** Could not create MPI Cartesian topology.\n";
    exit(-1);
  }
                   
  status = MPI_Comm_rank(d_cartComm, &d_mpiRank);
  if (status != MPI_SUCCESS ) {
    std::cout << "**ERROR** Could not get MPI Cartesian topology rank.\n";
    exit(-1);
  }

  status = MPI_Comm_size(d_cartComm, &d_mpiSize);
  if (status != MPI_SUCCESS ) {
    std::cout << "**ERROR** Could not get MPI Cartesian topology size.\n";
    exit(-1);
  }

  status = MPI_Cart_coords(d_cartComm, d_mpiRank, ndim, d_mpiCoords.data());
  if (status != MPI_SUCCESS ) {
    std::cout << "**ERROR** Could not get MPI Cartesian topology coords.\n";
    exit(-1);
  }

  d_mpiTag = 0;
  assert(d_mpiRank == d_boostWorld.rank());
}

// this is to scatter the dem and sph particle
// two point: (1) the sph ghost particles will be partitioned with dem
// particles together. There is no explicite partition for sph ghost particles.
//        After the partition of the dem particles, the sph ghost particles
// will be paritioned as a class member of dem particle.
//        Before partition, SPHParticle.demParticle points to NULL; after
// receive, SPHParticle.demParticle can point to their dem particle. July 27,
// 2015
//            (2) block partitions for dem domain and sph domain have to be
// the same,
//        otherwise some dem particles and their nearby sph particles will
// belong to different block
//        (grid is expanded to include the boundary sph particles)
//    (3) this partition method here is not very suitable for the free surface
// flow problem, such as bursting dam problem
//        since the total domain is divided, while there are lots of voids in
// the domain. (partition only free sph particles will be better)
//        But our goal is to simulate the porous media in triaxial, insotropic
// or SHPB simulations, particles are filled in the container.
void
SmoothParticleHydro::scatterSPHParticle(const Box& allContainer,
                                        const REAL& ghostWidth,
                                        REAL& bufferLength)
{
  // Comute the buffer size for the computational domain
  // int numLayers = util::getParam<int>("numLayers");
  // REAL spaceInterval = util::getParam<REAL>("spaceInterval");
  // bufferLength = spaceInterval * numLayers;

  // Set d_sphGrid on all procs
  setGrid(Box(allContainer, bufferLength));

  // Create patch for the current process
  int iteration = 0;
  createPatch(iteration, ghostWidth);

  // partition particles and send to each process
  if (d_mpiRank == 0) { // process 0

    /*
    for (auto part : d_allSPHParticleVec) {
      std::cout << "Id:" << part->getId()
                << " Pos:" << part->getInitPosition()
                << " CurPos:" << part->currentPosition() << "\n";
    }
    */

    Vec v1 = d_sphGrid.getMinCorner();
    Vec v2 = d_sphGrid.getMaxCorner();
    Vec vspan = (v2 - v1) / d_mpiProcs;

    boost::mpi::request reqs[d_mpiSize - 1];

    for (int iRank = d_mpiSize - 1; iRank >= 0; --iRank) {

      d_patchSPHParticleVec.clear();

      int ndim = 3;
      IntVec coords;
      MPI_Cart_coords(d_cartComm, iRank, ndim, coords.data());

      Vec lower = v1 + vspan * coords;
      Vec upper = lower + vspan;
      Box container(lower, upper);

      findSPHParticleInBox(container, d_allSPHParticleVec,
                           d_patchSPHParticleVec);

      if (iRank != 0) {

        reqs[iRank - 1] =
          d_boostWorld.isend(iRank, d_mpiTag, d_patchSPHParticleVec);

      } else {

        // Make a deep copy
        d_sphParticleVec.resize(d_patchSPHParticleVec.size());
        for (auto i = 0u; i < d_sphParticleVec.size(); ++i) {
          d_sphParticleVec[i] =
            std::make_shared<SPHParticle>(*d_patchSPHParticleVec[i]);
        }

      } 
    } // end iRank loop

    boost::mpi::wait_all(reqs, reqs + d_mpiSize - 1);

  } else { // other processes except 0

    d_boostWorld.recv(0, d_mpiTag, d_sphParticleVec);
  }

  //std::cout << "d_mpirank = " << d_mpiRank << " Num particles "
  //          << d_sphParticleVec.size() << "\n";

} // scatterDEMSPHParticle

void
SmoothParticleHydro::findSPHParticleInBox(const Box& container,
                                          const SPHParticlePArray& sphParticles,
                                          SPHParticlePArray& insideParticles)
{
  for (const auto& particle : sphParticles) {
    // it is critical to use EPS
    if (container.inside(particle->currentPosition(), dem::EPS)) {
      insideParticles.push_back(particle);
    }
  }
}

void
SmoothParticleHydro::createPatch(int iteration, const REAL& ghostWidth)
{
  // determine container of each process
  Vec v1 = d_sphGrid.getMinCorner();
  Vec v2 = d_sphGrid.getMaxCorner();
  Vec vspan = (v2 - v1) / d_mpiProcs;
  Vec lower = v1 + vspan * d_mpiCoords;
  Vec upper = lower + vspan;
  d_patchP = std::make_unique<SPHPatch>(d_cartComm, d_mpiRank, d_mpiCoords,
                                        lower, upper, ghostWidth, dem::EPS);
}

void
SmoothParticleHydro::updatePatch(int iteration, const REAL& ghostWidth)
{
  // determine container of each process
  Vec v1 = d_sphGrid.getMinCorner();
  Vec v2 = d_sphGrid.getMaxCorner();
  Vec vspan = (v2 - v1) / d_mpiProcs;
  Vec lower = v1 + vspan * d_mpiCoords;
  Vec upper = lower + vspan;
  Box container(lower, upper);
  d_patchP->update(iteration, lower, upper, ghostWidth);
}

/*
void
SmoothParticleHydro::removeSPHParticleOutRectangle()
{
  Vec v1 = container.getMinCorner();
  Vec v2 = container.getMaxCorner();
  REAL x1 = v1.x();
  REAL y1 = v1.y();
  REAL z1 = v1.z();
  REAL x2 = v2.x();
  REAL y2 = v2.y();
  REAL z2 = v2.z();

  std::vector<sph::SPHParticle*>::iterator itr;
  Vec center;
  // std::size_t flag = 0;

  // this for loop may be replaced by for_each, remove, erase and swap (shrink
  // to fit)
  // see item 14 in "effective STL"
  for (itr = d_sphParticleVec.begin(); itr != d_sphParticleVec.end();) {
    center = (*itr)->getCurrPosition();
    // it is critical to use EPS
    if (!(center.x() - x1 >= -EPS && center.x() - x2 < -EPS &&
          center.y() - y1 >= -EPS && center.y() - y2 < -EPS &&
          center.z() - z1 >= -EPS && center.z() - z2 < -EPS)) {
      //  debugInf << "iter=" << std::setw(8) << iteration << " rank=" <<
      // std::setw(2) << d_mpiRank
      //  << " removed=" << std::setw(3) << (*itr)->getId();
      //  flag = 1;
      delete (*itr); // release memory
      itr = d_sphParticleVec.erase(itr);
    } else
      ++itr;
  }
  //  if (flag == 1) {
  //  debugInf << " now " << particleVec.size() << ": ";
  //  for (std::vector<Particle*>::const_iterator it = particleVec.begin(); it
  //!= particleVec.end(); ++it)
  //  debugInf << std::setw(3) << (*it)->getId();
  //  debugInf << std::endl;
  //  }
}
*/

void
SmoothParticleHydro::commuSPHParticle(int iteration,
                                      const REAL& ghostWidth)
{
  updatePatch(iteration, ghostWidth);

  // Initialize merged particle array (particles + ghosts)
  d_mergeSPHParticleVec.clear();
  d_mergeSPHParticleVec = d_sphParticleVec;

  // Create temporary for checking received particles
  SPHParticlePArray recvParticles;

  // Plimpton scheme: x-ghost exchange
  d_patchP->sendRecvGhostXMinus(d_boostWorld, iteration, d_mergeSPHParticleVec);
  d_patchP->sendRecvGhostXPlus(d_boostWorld, iteration, d_mergeSPHParticleVec);
  d_patchP->waitToFinishX(iteration);
  d_patchP->combineReceivedParticlesX(iteration, d_mergeSPHParticleVec);
  d_patchP->combineReceivedParticlesX(iteration, recvParticles);

  /*
  proc0cout << "iteration = " << iteration 
            << " recv X = " << recvParticles.size()
            << " merge = " << d_mergeSPHParticleVec.size() << "\n";
  */

  // Plimpton scheme: y-ghost exchange
  d_patchP->sendRecvGhostYMinus(d_boostWorld, iteration, d_mergeSPHParticleVec);
  d_patchP->sendRecvGhostYPlus(d_boostWorld, iteration, d_mergeSPHParticleVec);
  d_patchP->waitToFinishY(iteration);
  d_patchP->combineReceivedParticlesY(iteration, d_mergeSPHParticleVec);
  d_patchP->combineReceivedParticlesY(iteration, recvParticles);

  /*
  proc0cout << "iteration = " << iteration 
            << " recv Y = " << recvParticles.size()
            << " merge = " << d_mergeSPHParticleVec.size() << "\n";
  */

  // Plimpton scheme: z-ghost exchange
  d_patchP->sendRecvGhostZMinus(d_boostWorld, iteration, d_mergeSPHParticleVec);
  d_patchP->sendRecvGhostZPlus(d_boostWorld, iteration, d_mergeSPHParticleVec);
  d_patchP->waitToFinishZ(iteration);
  d_patchP->combineReceivedParticlesZ(iteration, d_mergeSPHParticleVec);
  d_patchP->combineReceivedParticlesZ(iteration, recvParticles);

  /*
  proc0cout << "iteration = " << iteration 
            << " recv Z = " << recvParticles.size()
            << " merge = " << d_mergeSPHParticleVec.size() << "\n";
  */
}

/*
void
SmoothParticleHydro::releaseRecvParticle()
{
  // release memory of received particles
  for (std::vector<Particle*>::iterator it = recvParticleVec.begin();
       it != recvParticleVec.end(); ++it) {
    for (std::vector<sph::SPHParticle*>::iterator st =
           (*it)->SPHGhostParticleVec.begin();
         st != (*it)->SPHGhostParticleVec.end(); st++) {
      delete (*st); // release memory of these sph ghost particles
    }
    (*it)->SPHGhostParticleVec.clear();
    delete (*it);
  }
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
SmoothParticleHydro::releaseRecvSPHParticle()
{
  // release memory of received particles
  for (std::vector<sph::SPHParticle*>::iterator it =
         d_recvSPHParticleVec.begin();
       it != d_recvSPHParticleVec.end(); ++it)
    delete (*it);
  d_recvSPHParticleVec.clear();
  // 6 surfaces
  rsphParticleX1.clear();
  rsphParticleX2.clear();
  rsphParticleY1.clear();
  rsphParticleY2.clear();
  rsphParticleZ1.clear();
  rsphParticleZ2.clear();
  // 12 edges
  rsphParticleX1Y1.clear();
  rsphParticleX1Y2.clear();
  rsphParticleX1Z1.clear();
  rsphParticleX1Z2.clear();
  rsphParticleX2Y1.clear();
  rsphParticleX2Y2.clear();
  rsphParticleX2Z1.clear();
  rsphParticleX2Z2.clear();
  rsphParticleY1Z1.clear();
  rsphParticleY1Z2.clear();
  rsphParticleY2Z1.clear();
  rsphParticleY2Z2.clear();
  // 8 vertices
  rsphParticleX1Y1Z1.clear();
  rsphParticleX1Y1Z2.clear();
  rsphParticleX1Y2Z1.clear();
  rsphParticleX1Y2Z2.clear();
  rsphParticleX2Y1Z1.clear();
  rsphParticleX2Y1Z2.clear();
  rsphParticleX2Y2Z1.clear();
  rsphParticleX2Y2Z2.clear();
}

void
SmoothParticleHydro::migrateParticle()
{
  Vec vspan = grid.getMaxCorner() - grid.getMinCorner();
  REAL segX = vspan.x() / mpiProcX;
  REAL segY = vspan.y() / mpiProcY;
  REAL segZ = vspan.z() / mpiProcZ;
  Vec v1 = container.getMinCorner(); // v1, v2 in terms of process
  Vec v2 = container.getMaxCorner();

  // if a neighbor exists, transfer particles crossing the boundary in
  // between.
  std::vector<Particle *> particleX1, particleX2;
  std::vector<Particle *> particleY1, particleY2;
  std::vector<Particle *> particleZ1, particleZ2;
  std::vector<Particle *> particleX1Y1, particleX1Y2, particleX1Z1,
    particleX1Z2;
  std::vector<Particle *> particleX2Y1, particleX2Y2, particleX2Z1,
    particleX2Z2;
  std::vector<Particle *> particleY1Z1, particleY1Z2, particleY2Z1,
    particleY2Z2;
  std::vector<Particle *> particleX1Y1Z1, particleX1Y1Z2, particleX1Y2Z1,
    particleX1Y2Z2;
  std::vector<Particle *> particleX2Y1Z1, particleX2Y1Z2, particleX2Y2Z1,
    particleX2Y2Z2;
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
    Rectangle containerX1(v1.x() - segX, v1.y(), v1.z(), v1.x(), v2.y(),
                          v2.z());
    findParticleInRectangle(containerX1, particleVec, particleX1);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleX1.begin();
         it != particleX1.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle(); // now the SPHGhostParticles in
    // part of particleVec are pointing to NULL, don't needed to be back in
    // migrate
    reqX1[0] = boostWorld.isend(rankX1, mpiTag, particleX1);
    reqX1[1] = boostWorld.irecv(rankX1, mpiTag, rParticleX1);
  }
  if (rankX2 >= 0) { // surface x2
    Rectangle containerX2(v2.x(), v1.y(), v1.z(), v2.x() + segX, v2.y(),
                          v2.z());
    findParticleInRectangle(containerX2, particleVec, particleX2);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleX2.begin();
         it != particleX2.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqX2[0] = boostWorld.isend(rankX2, mpiTag, particleX2);
    reqX2[1] = boostWorld.irecv(rankX2, mpiTag, rParticleX2);
  }
  if (rankY1 >= 0) { // surface y1
    Rectangle containerY1(v1.x(), v1.y() - segY, v1.z(), v2.x(), v1.y(),
                          v2.z());
    findParticleInRectangle(containerY1, particleVec, particleY1);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleY1.begin();
         it != particleY1.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqY1[0] = boostWorld.isend(rankY1, mpiTag, particleY1);
    reqY1[1] = boostWorld.irecv(rankY1, mpiTag, rParticleY1);
  }
  if (rankY2 >= 0) { // surface y2
    Rectangle containerY2(v1.x(), v2.y(), v1.z(), v2.x(), v2.y() + segY,
                          v2.z());
    findParticleInRectangle(containerY2, particleVec, particleY2);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleY2.begin();
         it != particleY2.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqY2[0] = boostWorld.isend(rankY2, mpiTag, particleY2);
    reqY2[1] = boostWorld.irecv(rankY2, mpiTag, rParticleY2);
  }
  if (rankZ1 >= 0) { // surface z1
    Rectangle containerZ1(v1.x(), v1.y(), v1.z() - segZ, v2.x(), v2.y(),
                          v1.z());
    findParticleInRectangle(containerZ1, particleVec, particleZ1);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleZ1.begin();
         it != particleZ1.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqZ1[0] = boostWorld.isend(rankZ1, mpiTag, particleZ1);
    reqZ1[1] = boostWorld.irecv(rankZ1, mpiTag, rParticleZ1);
  }
  if (rankZ2 >= 0) { // surface z2
    Rectangle containerZ2(v1.x(), v1.y(), v2.z(), v2.x(), v2.y(),
                          v2.z() + segZ);
    findParticleInRectangle(containerZ2, particleVec, particleZ2);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleZ2.begin();
         it != particleZ2.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqZ2[0] = boostWorld.isend(rankZ2, mpiTag, particleZ2);
    reqZ2[1] = boostWorld.irecv(rankZ2, mpiTag, rParticleZ2);
  }
  // 12 edges
  if (rankX1Y1 >= 0) { // edge x1y1
    Rectangle containerX1Y1(v1.x() - segX, v1.y() - segY, v1.z(), v1.x(),
                            v1.y(), v2.z());
    findParticleInRectangle(containerX1Y1, particleVec, particleX1Y1);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleX1Y1.begin();
         it != particleX1Y1.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqX1Y1[0] = boostWorld.isend(rankX1Y1, mpiTag, particleX1Y1);
    reqX1Y1[1] = boostWorld.irecv(rankX1Y1, mpiTag, rParticleX1Y1);
  }
  if (rankX1Y2 >= 0) { // edge x1y2
    Rectangle containerX1Y2(v1.x() - segX, v2.y(), v1.z(), v1.x(),
                            v2.y() + segY, v2.z());
    findParticleInRectangle(containerX1Y2, particleVec, particleX1Y2);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleX1Y2.begin();
         it != particleX1Y2.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqX1Y2[0] = boostWorld.isend(rankX1Y2, mpiTag, particleX1Y2);
    reqX1Y2[1] = boostWorld.irecv(rankX1Y2, mpiTag, rParticleX1Y2);
  }
  if (rankX1Z1 >= 0) { // edge x1z1
    Rectangle containerX1Z1(v1.x() - segX, v1.y(), v1.z() - segZ, v1.x(),
                            v2.y(), v1.z());
    findParticleInRectangle(containerX1Z1, particleVec, particleX1Z1);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleX1Z1.begin();
         it != particleX1Z1.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqX1Z1[0] = boostWorld.isend(rankX1Z1, mpiTag, particleX1Z1);
    reqX1Z1[1] = boostWorld.irecv(rankX1Z1, mpiTag, rParticleX1Z1);
  }
  if (rankX1Z2 >= 0) { // edge x1z2
    Rectangle containerX1Z2(v1.x() - segX, v1.y(), v2.z(), v1.x(), v2.y(),
                            v2.z() + segZ);
    findParticleInRectangle(containerX1Z2, particleVec, particleX1Z2);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleX1Z2.begin();
         it != particleX1Z2.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqX1Z2[0] = boostWorld.isend(rankX1Z2, mpiTag, particleX1Z2);
    reqX1Z2[1] = boostWorld.irecv(rankX1Z2, mpiTag, rParticleX1Z2);
  }
  if (rankX2Y1 >= 0) { // edge x2y1
    Rectangle containerX2Y1(v2.x(), v1.y() - segY, v1.z(), v2.x() + segX,
                            v1.y(), v2.z());
    findParticleInRectangle(containerX2Y1, particleVec, particleX2Y1);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleX2Y1.begin();
         it != particleX2Y1.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqX2Y1[0] = boostWorld.isend(rankX2Y1, mpiTag, particleX2Y1);
    reqX2Y1[1] = boostWorld.irecv(rankX2Y1, mpiTag, rParticleX2Y1);
  }
  if (rankX2Y2 >= 0) { // edge x2y2
    Rectangle containerX2Y2(v2.x(), v2.y(), v1.z(), v2.x() + segX,
                            v2.y() + segY, v2.z());
    findParticleInRectangle(containerX2Y2, particleVec, particleX2Y2);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleX2Y2.begin();
         it != particleX2Y2.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqX2Y2[0] = boostWorld.isend(rankX2Y2, mpiTag, particleX2Y2);
    reqX2Y2[1] = boostWorld.irecv(rankX2Y2, mpiTag, rParticleX2Y2);
  }
  if (rankX2Z1 >= 0) { // edge x2z1
    Rectangle containerX2Z1(v2.x(), v1.y(), v1.z() - segZ, v2.x() + segX,
                            v2.y(), v1.z());
    findParticleInRectangle(containerX2Z1, particleVec, particleX2Z1);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleX2Z1.begin();
         it != particleX2Z1.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqX2Z1[0] = boostWorld.isend(rankX2Z1, mpiTag, particleX2Z1);
    reqX2Z1[1] = boostWorld.irecv(rankX2Z1, mpiTag, rParticleX2Z1);
  }
  if (rankX2Z2 >= 0) { // edge x2z2
    Rectangle containerX2Z2(v2.x(), v1.y(), v2.z(), v2.x() + segX, v2.y(),
                            v2.z() + segZ);
    findParticleInRectangle(containerX2Z2, particleVec, particleX2Z2);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleX2Z2.begin();
         it != particleX2Z2.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqX2Z2[0] = boostWorld.isend(rankX2Z2, mpiTag, particleX2Z2);
    reqX2Z2[1] = boostWorld.irecv(rankX2Z2, mpiTag, rParticleX2Z2);
  }
  if (rankY1Z1 >= 0) { // edge y1z1
    Rectangle containerY1Z1(v1.x(), v1.y() - segY, v1.z() - segZ, v2.x(),
                            v1.y(), v1.z());
    findParticleInRectangle(containerY1Z1, particleVec, particleY1Z1);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleY1Z1.begin();
         it != particleY1Z1.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqY1Z1[0] = boostWorld.isend(rankY1Z1, mpiTag, particleY1Z1);
    reqY1Z1[1] = boostWorld.irecv(rankY1Z1, mpiTag, rParticleY1Z1);
  }
  if (rankY1Z2 >= 0) { // edge y1z2
    Rectangle containerY1Z2(v1.x(), v1.y() - segY, v2.z(), v2.x(), v1.y(),
                            v2.z() + segZ);
    findParticleInRectangle(containerY1Z2, particleVec, particleY1Z2);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleY1Z2.begin();
         it != particleY1Z2.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqY1Z2[0] = boostWorld.isend(rankY1Z2, mpiTag, particleY1Z2);
    reqY1Z2[1] = boostWorld.irecv(rankY1Z2, mpiTag, rParticleY1Z2);
  }
  if (rankY2Z1 >= 0) { // edge y2z1
    Rectangle containerY2Z1(v1.x(), v2.y(), v1.z() - segZ, v2.x(),
                            v2.y() + segY, v1.z());
    findParticleInRectangle(containerY2Z1, particleVec, particleY2Z1);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleY2Z1.begin();
         it != particleY2Z1.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqY2Z1[0] = boostWorld.isend(rankY2Z1, mpiTag, particleY2Z1);
    reqY2Z1[1] = boostWorld.irecv(rankY2Z1, mpiTag, rParticleY2Z1);
  }
  if (rankY2Z2 >= 0) { // edge y2z2
    Rectangle containerY2Z2(v1.x(), v2.y(), v2.z(), v2.x(), v2.y() + segY,
                            v2.z() + segZ);
    findParticleInRectangle(containerY2Z2, particleVec, particleY2Z2);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleY2Z2.begin();
         it != particleY2Z2.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqY2Z2[0] = boostWorld.isend(rankY2Z2, mpiTag, particleY2Z2);
    reqY2Z2[1] = boostWorld.irecv(rankY2Z2, mpiTag, rParticleY2Z2);
  }
  // 8 vertices
  if (rankX1Y1Z1 >= 0) { // edge x1y1z1
    Rectangle containerX1Y1Z1(v1.x() - segX, v1.y() - segY, v1.z() - segZ,
                              v1.x(), v1.y(), v1.z());
    findParticleInRectangle(containerX1Y1Z1, particleVec, particleX1Y1Z1);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleX1Y1Z1.begin();
         it != particleX1Y1Z1.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqX1Y1Z1[0] = boostWorld.isend(rankX1Y1Z1, mpiTag, particleX1Y1Z1);
    reqX1Y1Z1[1] = boostWorld.irecv(rankX1Y1Z1, mpiTag, rParticleX1Y1Z1);
  }
  if (rankX1Y1Z2 >= 0) { // edge x1y1z2
    Rectangle containerX1Y1Z2(v1.x() - segX, v1.y() - segY, v2.z(), v1.x(),
                              v1.y(), v2.z() + segZ);
    findParticleInRectangle(containerX1Y1Z2, particleVec, particleX1Y1Z2);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleX1Y1Z2.begin();
         it != particleX1Y1Z2.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqX1Y1Z2[0] = boostWorld.isend(rankX1Y1Z2, mpiTag, particleX1Y1Z2);
    reqX1Y1Z2[1] = boostWorld.irecv(rankX1Y1Z2, mpiTag, rParticleX1Y1Z2);
  }
  if (rankX1Y2Z1 >= 0) { // edge x1y2z1
    Rectangle containerX1Y2Z1(v1.x() - segX, v2.y(), v1.z() - segZ, v1.x(),
                              v2.y() + segY, v1.z());
    findParticleInRectangle(containerX1Y2Z1, particleVec, particleX1Y2Z1);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleX1Y2Z1.begin();
         it != particleX1Y2Z1.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqX1Y2Z1[0] = boostWorld.isend(rankX1Y2Z1, mpiTag, particleX1Y2Z1);
    reqX1Y2Z1[1] = boostWorld.irecv(rankX1Y2Z1, mpiTag, rParticleX1Y2Z1);
  }
  if (rankX1Y2Z2 >= 0) { // edge x1y2z2
    Rectangle containerX1Y2Z2(v1.x() - segX, v2.y(), v2.z(), v1.x(),
                              v2.y() + segY, v2.z() + segZ);
    findParticleInRectangle(containerX1Y2Z2, particleVec, particleX1Y2Z2);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleX1Y2Z2.begin();
         it != particleX1Y2Z2.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqX1Y2Z2[0] = boostWorld.isend(rankX1Y2Z2, mpiTag, particleX1Y2Z2);
    reqX1Y2Z2[1] = boostWorld.irecv(rankX1Y2Z2, mpiTag, rParticleX1Y2Z2);
  }
  if (rankX2Y1Z1 >= 0) { // edge x2y1z1
    Rectangle containerX2Y1Z1(v2.x(), v1.y() - segY, v1.z() - segZ,
                              v2.x() + segX, v1.y(), v1.z());
    findParticleInRectangle(containerX2Y1Z1, particleVec, particleX2Y1Z1);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleX2Y1Z1.begin();
         it != particleX2Y1Z1.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqX2Y1Z1[0] = boostWorld.isend(rankX2Y1Z1, mpiTag, particleX2Y1Z1);
    reqX2Y1Z1[1] = boostWorld.irecv(rankX2Y1Z1, mpiTag, rParticleX2Y1Z1);
  }
  if (rankX2Y1Z2 >= 0) { // edge x2y1z2
    Rectangle containerX2Y1Z2(v2.x(), v1.y() - segY, v2.z(), v2.x() + segX,
                              v1.y(), v2.z() + segZ);
    findParticleInRectangle(containerX2Y1Z2, particleVec, particleX2Y1Z2);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleX2Y1Z2.begin();
         it != particleX2Y1Z2.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqX2Y1Z2[0] = boostWorld.isend(rankX2Y1Z2, mpiTag, particleX2Y1Z2);
    reqX2Y1Z2[1] = boostWorld.irecv(rankX2Y1Z2, mpiTag, rParticleX2Y1Z2);
  }
  if (rankX2Y2Z1 >= 0) { // edge x2y2z1
    Rectangle containerX2Y2Z1(v2.x(), v2.y(), v1.z() - segZ, v2.x() + segX,
                              v2.y() + segY, v1.z());
    findParticleInRectangle(containerX2Y2Z1, particleVec, particleX2Y2Z1);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleX2Y2Z1.begin();
         it != particleX2Y2Z1.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
    reqX2Y2Z1[0] = boostWorld.isend(rankX2Y2Z1, mpiTag, particleX2Y2Z1);
    reqX2Y2Z1[1] = boostWorld.irecv(rankX2Y2Z1, mpiTag, rParticleX2Y2Z1);
  }
  if (rankX2Y2Z2 >= 0) { // edge x2y2z2
    Rectangle containerX2Y2Z2(v2.x(), v2.y(), v2.z(), v2.x() + segX,
                              v2.y() + segY, v2.z() + segZ);
    findParticleInRectangle(containerX2Y2Z2, particleVec, particleX2Y2Z2);
    // before send, SPHParticle.demParticle should be NULL
    for (std::vector<Particle*>::iterator it = particleX2Y2Z2.begin();
         it != particleX2Y2Z2.end(); it++)
      (*it)->setNULLDemParticleInSPHParticle();
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

  // do not need to assign SPHGhostParticle.demParticle back in particleX1...,
  // since these particles will be removed
  // delete outgoing particles
  removeParticleOutRectangle();

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

  // after receive, set SPHParticle.demParticle
  for (std::vector<Particle*>::iterator it = recvParticleVec.begin();
       it != recvParticleVec.end(); it++) {
    (*it)->setDemParticleInSPHParticle();
  }

  particleVec.insert(particleVec.end(), recvParticleVec.begin(),
                     recvParticleVec.end());

#ifdef 0
  if (recvParticleVec.size() > 0) {
    debugInf << "iter=" << std::setw(8) << iteration << " rank=" << std::setw(2)
             << d_mpiRank << "   added=";
    for (std::vector<Particle*>::const_iterator it = recvParticleVec.begin();
         it != recvParticleVec.end(); ++it)
      debugInf << std::setw(3) << (*it)->getId();
    debugInf << " now " << particleVec.size() << ": ";
    for (std::vector<Particle*>::const_iterator it = particleVec.begin();
         it != particleVec.end(); ++it)
      debugInf << std::setw(3) << (*it)->getId();
    debugInf << std::endl;
  }
#endif

  // do not release memory of received particles because they are part of and
  managed by particleVec
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

void
SmoothParticleHydro::migrateSPHParticle()
{
  Vec vspan = grid.getMaxCorner() - grid.getMinCorner();
  REAL segX = vspan.x() / mpiProcX;
  REAL segY = vspan.y() / mpiProcY;
  REAL segZ = vspan.z() / mpiProcZ;
  Vec v1 = container.getMinCorner(); // v1, v2 in terms of process
  Vec v2 = container.getMaxCorner();

  // if a neighbor exists, transfer particles crossing the boundary in
  // between.
  std::vector<sph::SPHParticle *> sphParticleX1, sphParticleX2;
  std::vector<sph::SPHParticle *> sphParticleY1, sphParticleY2;
  std::vector<sph::SPHParticle *> sphParticleZ1, sphParticleZ2;
  std::vector<sph::SPHParticle *> sphParticleX1Y1, sphParticleX1Y2,
    sphParticleX1Z1, sphParticleX1Z2;
  std::vector<sph::SPHParticle *> sphParticleX2Y1, sphParticleX2Y2,
    sphParticleX2Z1, sphParticleX2Z2;
  std::vector<sph::SPHParticle *> sphParticleY1Z1, sphParticleY1Z2,
    sphParticleY2Z1, sphParticleY2Z2;
  std::vector<sph::SPHParticle *> sphParticleX1Y1Z1, sphParticleX1Y1Z2,
    sphParticleX1Y2Z1, sphParticleX1Y2Z2;
  std::vector<sph::SPHParticle *> sphParticleX2Y1Z1, sphParticleX2Y1Z2,
    sphParticleX2Y2Z1, sphParticleX2Y2Z2;
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
    Rectangle containerX1(v1.x() - segX, v1.y(), v1.z(), v1.x(), v2.y(),
                          v2.z());
    findSPHParticleInRectangle(containerX1, d_sphParticleVec, sphParticleX1);
    reqX1[0] = boostWorld.isend(rankX1, mpiTag, sphParticleX1);
    reqX1[1] = boostWorld.irecv(rankX1, mpiTag, rsphParticleX1);
  }
  if (rankX2 >= 0) { // surface x2
    Rectangle containerX2(v2.x(), v1.y(), v1.z(), v2.x() + segX, v2.y(),
                          v2.z());
    findSPHParticleInRectangle(containerX2, d_sphParticleVec, sphParticleX2);
    reqX2[0] = boostWorld.isend(rankX2, mpiTag, sphParticleX2);
    reqX2[1] = boostWorld.irecv(rankX2, mpiTag, rsphParticleX2);
  }
  if (rankY1 >= 0) { // surface y1
    Rectangle containerY1(v1.x(), v1.y() - segY, v1.z(), v2.x(), v1.y(),
                          v2.z());
    findSPHParticleInRectangle(containerY1, d_sphParticleVec, sphParticleY1);
    reqY1[0] = boostWorld.isend(rankY1, mpiTag, sphParticleY1);
    reqY1[1] = boostWorld.irecv(rankY1, mpiTag, rsphParticleY1);
  }
  if (rankY2 >= 0) { // surface y2
    Rectangle containerY2(v1.x(), v2.y(), v1.z(), v2.x(), v2.y() + segY,
                          v2.z());
    findSPHParticleInRectangle(containerY2, d_sphParticleVec, sphParticleY2);
    reqY2[0] = boostWorld.isend(rankY2, mpiTag, sphParticleY2);
    reqY2[1] = boostWorld.irecv(rankY2, mpiTag, rsphParticleY2);
  }
  if (rankZ1 >= 0) { // surface z1
    Rectangle containerZ1(v1.x(), v1.y(), v1.z() - segZ, v2.x(), v2.y(),
                          v1.z());
    findSPHParticleInRectangle(containerZ1, d_sphParticleVec, sphParticleZ1);
    reqZ1[0] = boostWorld.isend(rankZ1, mpiTag, sphParticleZ1);
    reqZ1[1] = boostWorld.irecv(rankZ1, mpiTag, rsphParticleZ1);
  }
  if (rankZ2 >= 0) { // surface z2
    Rectangle containerZ2(v1.x(), v1.y(), v2.z(), v2.x(), v2.y(),
                          v2.z() + segZ);
    findSPHParticleInRectangle(containerZ2, d_sphParticleVec, sphParticleZ2);
    reqZ2[0] = boostWorld.isend(rankZ2, mpiTag, sphParticleZ2);
    reqZ2[1] = boostWorld.irecv(rankZ2, mpiTag, rsphParticleZ2);
  }
  // 12 edges
  if (rankX1Y1 >= 0) { // edge x1y1
    Rectangle containerX1Y1(v1.x() - segX, v1.y() - segY, v1.z(), v1.x(),
                            v1.y(), v2.z());
    findSPHParticleInRectangle(containerX1Y1, d_sphParticleVec,
                               sphParticleX1Y1);
    reqX1Y1[0] = boostWorld.isend(rankX1Y1, mpiTag, sphParticleX1Y1);
    reqX1Y1[1] = boostWorld.irecv(rankX1Y1, mpiTag, rsphParticleX1Y1);
  }
  if (rankX1Y2 >= 0) { // edge x1y2
    Rectangle containerX1Y2(v1.x() - segX, v2.y(), v1.z(), v1.x(),
                            v2.y() + segY, v2.z());
    findSPHParticleInRectangle(containerX1Y2, d_sphParticleVec,
                               sphParticleX1Y2);
    reqX1Y2[0] = boostWorld.isend(rankX1Y2, mpiTag, sphParticleX1Y2);
    reqX1Y2[1] = boostWorld.irecv(rankX1Y2, mpiTag, rsphParticleX1Y2);
  }
  if (rankX1Z1 >= 0) { // edge x1z1
    Rectangle containerX1Z1(v1.x() - segX, v1.y(), v1.z() - segZ, v1.x(),
                            v2.y(), v1.z());
    findSPHParticleInRectangle(containerX1Z1, d_sphParticleVec,
                               sphParticleX1Z1);
    reqX1Z1[0] = boostWorld.isend(rankX1Z1, mpiTag, sphParticleX1Z1);
    reqX1Z1[1] = boostWorld.irecv(rankX1Z1, mpiTag, rsphParticleX1Z1);
  }
  if (rankX1Z2 >= 0) { // edge x1z2
    Rectangle containerX1Z2(v1.x() - segX, v1.y(), v2.z(), v1.x(), v2.y(),
                            v2.z() + segZ);
    findSPHParticleInRectangle(containerX1Z2, d_sphParticleVec,
                               sphParticleX1Z2);
    reqX1Z2[0] = boostWorld.isend(rankX1Z2, mpiTag, sphParticleX1Z2);
    reqX1Z2[1] = boostWorld.irecv(rankX1Z2, mpiTag, rsphParticleX1Z2);
  }
  if (rankX2Y1 >= 0) { // edge x2y1
    Rectangle containerX2Y1(v2.x(), v1.y() - segY, v1.z(), v2.x() + segX,
                            v1.y(), v2.z());
    findSPHParticleInRectangle(containerX2Y1, d_sphParticleVec,
                               sphParticleX2Y1);
    reqX2Y1[0] = boostWorld.isend(rankX2Y1, mpiTag, sphParticleX2Y1);
    reqX2Y1[1] = boostWorld.irecv(rankX2Y1, mpiTag, rsphParticleX2Y1);
  }
  if (rankX2Y2 >= 0) { // edge x2y2
    Rectangle containerX2Y2(v2.x(), v2.y(), v1.z(), v2.x() + segX,
                            v2.y() + segY, v2.z());
    findSPHParticleInRectangle(containerX2Y2, d_sphParticleVec,
                               sphParticleX2Y2);
    reqX2Y2[0] = boostWorld.isend(rankX2Y2, mpiTag, sphParticleX2Y2);
    reqX2Y2[1] = boostWorld.irecv(rankX2Y2, mpiTag, rsphParticleX2Y2);
  }
  if (rankX2Z1 >= 0) { // edge x2z1
    Rectangle containerX2Z1(v2.x(), v1.y(), v1.z() - segZ, v2.x() + segX,
                            v2.y(), v1.z());
    findSPHParticleInRectangle(containerX2Z1, d_sphParticleVec,
                               sphParticleX2Z1);
    reqX2Z1[0] = boostWorld.isend(rankX2Z1, mpiTag, sphParticleX2Z1);
    reqX2Z1[1] = boostWorld.irecv(rankX2Z1, mpiTag, rsphParticleX2Z1);
  }
  if (rankX2Z2 >= 0) { // edge x2z2
    Rectangle containerX2Z2(v2.x(), v1.y(), v2.z(), v2.x() + segX, v2.y(),
                            v2.z() + segZ);
    findSPHParticleInRectangle(containerX2Z2, d_sphParticleVec,
                               sphParticleX2Z2);
    reqX2Z2[0] = boostWorld.isend(rankX2Z2, mpiTag, sphParticleX2Z2);
    reqX2Z2[1] = boostWorld.irecv(rankX2Z2, mpiTag, rsphParticleX2Z2);
  }
  if (rankY1Z1 >= 0) { // edge y1z1
    Rectangle containerY1Z1(v1.x(), v1.y() - segY, v1.z() - segZ, v2.x(),
                            v1.y(), v1.z());
    findSPHParticleInRectangle(containerY1Z1, d_sphParticleVec,
                               sphParticleY1Z1);
    reqY1Z1[0] = boostWorld.isend(rankY1Z1, mpiTag, sphParticleY1Z1);
    reqY1Z1[1] = boostWorld.irecv(rankY1Z1, mpiTag, rsphParticleY1Z1);
  }
  if (rankY1Z2 >= 0) { // edge y1z2
    Rectangle containerY1Z2(v1.x(), v1.y() - segY, v2.z(), v2.x(), v1.y(),
                            v2.z() + segZ);
    findSPHParticleInRectangle(containerY1Z2, d_sphParticleVec,
                               sphParticleY1Z2);
    reqY1Z2[0] = boostWorld.isend(rankY1Z2, mpiTag, sphParticleY1Z2);
    reqY1Z2[1] = boostWorld.irecv(rankY1Z2, mpiTag, rsphParticleY1Z2);
  }
  if (rankY2Z1 >= 0) { // edge y2z1
    Rectangle containerY2Z1(v1.x(), v2.y(), v1.z() - segZ, v2.x(),
                            v2.y() + segY, v1.z());
    findSPHParticleInRectangle(containerY2Z1, d_sphParticleVec,
                               sphParticleY2Z1);
    reqY2Z1[0] = boostWorld.isend(rankY2Z1, mpiTag, sphParticleY2Z1);
    reqY2Z1[1] = boostWorld.irecv(rankY2Z1, mpiTag, rsphParticleY2Z1);
  }
  if (rankY2Z2 >= 0) { // edge y2z2
    Rectangle containerY2Z2(v1.x(), v2.y(), v2.z(), v2.x(), v2.y() + segY,
                            v2.z() + segZ);
    findSPHParticleInRectangle(containerY2Z2, d_sphParticleVec,
                               sphParticleY2Z2);
    reqY2Z2[0] = boostWorld.isend(rankY2Z2, mpiTag, sphParticleY2Z2);
    reqY2Z2[1] = boostWorld.irecv(rankY2Z2, mpiTag, rsphParticleY2Z2);
  }
  // 8 vertices
  if (rankX1Y1Z1 >= 0) { // edge x1y1z1
    Rectangle containerX1Y1Z1(v1.x() - segX, v1.y() - segY, v1.z() - segZ,
                              v1.x(), v1.y(), v1.z());
    findSPHParticleInRectangle(containerX1Y1Z1, d_sphParticleVec,
                               sphParticleX1Y1Z1);
    reqX1Y1Z1[0] = boostWorld.isend(rankX1Y1Z1, mpiTag, sphParticleX1Y1Z1);
    reqX1Y1Z1[1] = boostWorld.irecv(rankX1Y1Z1, mpiTag, rsphParticleX1Y1Z1);
  }
  if (rankX1Y1Z2 >= 0) { // edge x1y1z2
    Rectangle containerX1Y1Z2(v1.x() - segX, v1.y() - segY, v2.z(), v1.x(),
                              v1.y(), v2.z() + segZ);
    findSPHParticleInRectangle(containerX1Y1Z2, d_sphParticleVec,
                               sphParticleX1Y1Z2);
    reqX1Y1Z2[0] = boostWorld.isend(rankX1Y1Z2, mpiTag, sphParticleX1Y1Z2);
    reqX1Y1Z2[1] = boostWorld.irecv(rankX1Y1Z2, mpiTag, rsphParticleX1Y1Z2);
  }
  if (rankX1Y2Z1 >= 0) { // edge x1y2z1
    Rectangle containerX1Y2Z1(v1.x() - segX, v2.y(), v1.z() - segZ, v1.x(),
                              v2.y() + segY, v1.z());
    findSPHParticleInRectangle(containerX1Y2Z1, d_sphParticleVec,
                               sphParticleX1Y2Z1);
    reqX1Y2Z1[0] = boostWorld.isend(rankX1Y2Z1, mpiTag, sphParticleX1Y2Z1);
    reqX1Y2Z1[1] = boostWorld.irecv(rankX1Y2Z1, mpiTag, rsphParticleX1Y2Z1);
  }
  if (rankX1Y2Z2 >= 0) { // edge x1y2z2
    Rectangle containerX1Y2Z2(v1.x() - segX, v2.y(), v2.z(), v1.x(),
                              v2.y() + segY, v2.z() + segZ);
    findSPHParticleInRectangle(containerX1Y2Z2, d_sphParticleVec,
                               sphParticleX1Y2Z2);
    reqX1Y2Z2[0] = boostWorld.isend(rankX1Y2Z2, mpiTag, sphParticleX1Y2Z2);
    reqX1Y2Z2[1] = boostWorld.irecv(rankX1Y2Z2, mpiTag, rsphParticleX1Y2Z2);
  }
  if (rankX2Y1Z1 >= 0) { // edge x2y1z1
    Rectangle containerX2Y1Z1(v2.x(), v1.y() - segY, v1.z() - segZ,
                              v2.x() + segX, v1.y(), v1.z());
    findSPHParticleInRectangle(containerX2Y1Z1, d_sphParticleVec,
                               sphParticleX2Y1Z1);
    reqX2Y1Z1[0] = boostWorld.isend(rankX2Y1Z1, mpiTag, sphParticleX2Y1Z1);
    reqX2Y1Z1[1] = boostWorld.irecv(rankX2Y1Z1, mpiTag, rsphParticleX2Y1Z1);
  }
  if (rankX2Y1Z2 >= 0) { // edge x2y1z2
    Rectangle containerX2Y1Z2(v2.x(), v1.y() - segY, v2.z(), v2.x() + segX,
                              v1.y(), v2.z() + segZ);
    findSPHParticleInRectangle(containerX2Y1Z2, d_sphParticleVec,
                               sphParticleX2Y1Z2);
    reqX2Y1Z2[0] = boostWorld.isend(rankX2Y1Z2, mpiTag, sphParticleX2Y1Z2);
    reqX2Y1Z2[1] = boostWorld.irecv(rankX2Y1Z2, mpiTag, rsphParticleX2Y1Z2);
  }
  if (rankX2Y2Z1 >= 0) { // edge x2y2z1
    Rectangle containerX2Y2Z1(v2.x(), v2.y(), v1.z() - segZ, v2.x() + segX,
                              v2.y() + segY, v1.z());
    findSPHParticleInRectangle(containerX2Y2Z1, d_sphParticleVec,
                               sphParticleX2Y2Z1);
    reqX2Y2Z1[0] = boostWorld.isend(rankX2Y2Z1, mpiTag, sphParticleX2Y2Z1);
    reqX2Y2Z1[1] = boostWorld.irecv(rankX2Y2Z1, mpiTag, rsphParticleX2Y2Z1);
  }
  if (rankX2Y2Z2 >= 0) { // edge x2y2z2
    Rectangle containerX2Y2Z2(v2.x(), v2.y(), v2.z(), v2.x() + segX,
                              v2.y() + segY, v2.z() + segZ);
    findSPHParticleInRectangle(containerX2Y2Z2, d_sphParticleVec,
                               sphParticleX2Y2Z2);
    reqX2Y2Z2[0] = boostWorld.isend(rankX2Y2Z2, mpiTag, sphParticleX2Y2Z2);
    reqX2Y2Z2[1] = boostWorld.irecv(rankX2Y2Z2, mpiTag, rsphParticleX2Y2Z2);
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
  removeSPHParticleOutRectangle();

  // add incoming particles
  d_recvSPHParticleVec.clear(); // new use of recvParticleVec
  // 6 surfaces
  if (rankX1 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleX1.begin(), rsphParticleX1.end());
  if (rankX2 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleX2.begin(), rsphParticleX2.end());
  if (rankY1 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleY1.begin(), rsphParticleY1.end());
  if (rankY2 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleY2.begin(), rsphParticleY2.end());
  if (rankZ1 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleZ1.begin(), rsphParticleZ1.end());
  if (rankZ2 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleZ2.begin(), rsphParticleZ2.end());
  // 12 edges
  if (rankX1Y1 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleX1Y1.begin(),
                                rsphParticleX1Y1.end());
  if (rankX1Y2 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleX1Y2.begin(),
                                rsphParticleX1Y2.end());
  if (rankX1Z1 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleX1Z1.begin(),
                                rsphParticleX1Z1.end());
  if (rankX1Z2 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleX1Z2.begin(),
                                rsphParticleX1Z2.end());
  if (rankX2Y1 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleX2Y1.begin(),
                                rsphParticleX2Y1.end());
  if (rankX2Y2 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleX2Y2.begin(),
                                rsphParticleX2Y2.end());
  if (rankX2Z1 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleX2Z1.begin(),
                                rsphParticleX2Z1.end());
  if (rankX2Z2 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleX2Z2.begin(),
                                rsphParticleX2Z2.end());
  if (rankY1Z1 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleY1Z1.begin(),
                                rsphParticleY1Z1.end());
  if (rankY1Z2 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleY1Z2.begin(),
                                rsphParticleY1Z2.end());
  if (rankY2Z1 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleY2Z1.begin(),
                                rsphParticleY2Z1.end());
  if (rankY2Z2 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleY2Z2.begin(),
                                rsphParticleY2Z2.end());
  // 8 vertices
  if (rankX1Y1Z1 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleX1Y1Z1.begin(),
                                rsphParticleX1Y1Z1.end());
  if (rankX1Y1Z2 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleX1Y1Z2.begin(),
                                rsphParticleX1Y1Z2.end());
  if (rankX1Y2Z1 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleX1Y2Z1.begin(),
                                rsphParticleX1Y2Z1.end());
  if (rankX1Y2Z2 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleX1Y2Z2.begin(),
                                rsphParticleX1Y2Z2.end());
  if (rankX2Y1Z1 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleX2Y1Z1.begin(),
                                rsphParticleX2Y1Z1.end());
  if (rankX2Y1Z2 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleX2Y1Z2.begin(),
                                rsphParticleX2Y1Z2.end());
  if (rankX2Y2Z1 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleX2Y2Z1.begin(),
                                rsphParticleX2Y2Z1.end());
  if (rankX2Y2Z2 >= 0)
    d_recvSPHParticleVec.insert(d_recvSPHParticleVec.end(),
                                rsphParticleX2Y2Z2.begin(),
                                rsphParticleX2Y2Z2.end());

  d_sphParticleVec.insert(d_sphParticleVec.end(), d_recvSPHParticleVec.begin(),
                          d_recvSPHParticleVec.end());

#ifdef 0
  if (recvParticleVec.size() > 0) {
    debugInf << "iter=" << std::setw(8) << iteration << " rank=" << std::setw(2)
             << d_mpiRank << "   added=";
    for (std::vector<Particle*>::const_iterator it = recvParticleVec.begin();
         it != recvParticleVec.end(); ++it)
      debugInf << std::setw(3) << (*it)->getId();
    debugInf << " now " << particleVec.size() << ": ";
    for (std::vector<Particle*>::const_iterator it = particleVec.begin();
         it != particleVec.end(); ++it)
      debugInf << std::setw(3) << (*it)->getId();
    debugInf << std::endl;
  }
#endif

  // do not release memory of received particles because they are part of and
  // managed by particleVec
  // 6 surfaces
  rsphParticleX1.clear();
  rsphParticleX2.clear();
  rsphParticleY1.clear();
  rsphParticleY2.clear();
  rsphParticleZ1.clear();
  rsphParticleZ2.clear();
  // 12 edges
  rsphParticleX1Y1.clear();
  rsphParticleX1Y2.clear();
  rsphParticleX1Z1.clear();
  rsphParticleX1Z2.clear();
  rsphParticleX2Y1.clear();
  rsphParticleX2Y2.clear();
  rsphParticleX2Z1.clear();
  rsphParticleX2Z2.clear();
  rsphParticleY1Z1.clear();
  rsphParticleY1Z2.clear();
  rsphParticleY2Z1.clear();
  rsphParticleY2Z2.clear();
  // 8 vertices
  rsphParticleX1Y1Z1.clear();
  rsphParticleX1Y1Z2.clear();
  rsphParticleX1Y2Z1.clear();
  rsphParticleX1Y2Z2.clear();
  rsphParticleX2Y1Z1.clear();
  rsphParticleX2Y1Z2.clear();
  rsphParticleX2Y2Z1.clear();
  rsphParticleX2Y2Z2.clear();

  d_recvSPHParticleVec.clear();
}

void
SmoothParticleHydro::gatherParticle()
{
  // before send, SPHParticle.demParticle should be NULL
  for (std::vector<Particle*>::iterator it = particleVec.begin();
       it != particleVec.end(); it++)
    (*it)->setNULLDemParticleInSPHParticle(); // at this point,
  // SPHGhostParticle.demParticle is not pointing to particleVec

  // update allParticleVec: process 0 collects all updated particles from each
  // other process
  if (d_mpiRank != 0) { // each process except 0
    boostWorld.send(0, mpiTag, particleVec);
  } else { // process 0
    // allParticleVec is cleared before filling with new data
    releaseGatheredParticle();

    // duplicate particleVec so that it is not destroyed by allParticleVec in
    // next iteration,
    // otherwise it causes memory error.
    std::vector<Particle*> dupParticleVec(particleVec.size());
    for (std::size_t i = 0; i < dupParticleVec.size(); ++i) {
      dupParticleVec[i] = new Particle(*particleVec[i]); // at this point,
      // dupParticleVec and particleVec are pointint to the same
      // SPHGhoastParticle
      dupParticleVec[i]->SPHGhostParticleVec.clear(); // at this point,
      // dupParticleVec is pointing to nothing
      for (std::vector<sph::SPHParticle*>::iterator st =
             particleVec[i]->SPHGhostParticleVec.begin();
           st != particleVec[i]->SPHGhostParticleVec.end(); st++) {
        sph::SPHParticle* tmp_sph = new sph::SPHParticle(**st); // create a new
        // SPHGhost particle, which is the same as the one in particleVec
        dupParticleVec[i]->SPHGhostParticleVec.push_back(tmp_sph); // now
        // dupParticleVec points to the new SPHGhostParticle
      }
    }

    // fill allParticleVec with dupParticleVec and received particles
    allParticleVec.insert(allParticleVec.end(), dupParticleVec.begin(),
                          dupParticleVec.end());

    std::vector<Particle*> tmpParticleVec;
    long gatherRam = 0;
    for (int iRank = 1; iRank < d_mpiSize; ++iRank) {

      tmpParticleVec.clear(); // do not destroy particles!
      boostWorld.recv(iRank, mpiTag, tmpParticleVec);
      allParticleVec.insert(allParticleVec.end(), tmpParticleVec.begin(),
                            tmpParticleVec.end());
      gatherRam += tmpParticleVec.size();
    }
    // debugInf << "gather: particleNum = " << gatherRam <<  " particleRam = "
    //<< gatherRam * sizeof(Particle) << std::endl;
  }
  // after receive, set SPHParticle.demParticle
  for (std::vector<Particle*>::iterator it = particleVec.begin();
       it != particleVec.end(); it++) {
    (*it)->setDemParticleInSPHParticle();
  }
}

void
SmoothParticleHydro::releaseGatheredParticle()
{
  // clear allParticleVec, avoid long time memory footprint.
  for (std::vector<Particle*>::iterator it = allParticleVec.begin();
       it != allParticleVec.end(); ++it) {
    for (std::vector<sph::SPHParticle*>::iterator st =
           (*it)->SPHGhostParticleVec.begin();
         st != (*it)->SPHGhostParticleVec.end(); ++st) {
      delete (*st); // this is important to free the memories of sph ghost
                    // particles
    }
    (*it)->SPHGhostParticleVec.clear();
    std::vector<sph::SPHParticle*>().swap((*it)->SPHGhostParticleVec); //
    // actual memory release
    delete (*it);
  }
  allParticleVec.clear();
  std::vector<Particle*>().swap(allParticleVec); // actual memory release
}

void
SmoothParticleHydro::gatherSPHParticle()
{
  // update d_allSPHParticleVec: process 0 collects all updated particles from
  // each other process
  if (d_mpiRank != 0) { // each process except 0
    boostWorld.send(0, mpiTag, d_sphParticleVec);
  } else { // process 0
    // d_allSPHParticleVec is cleared before filling with new data
    releaseGatheredSPHParticle();

    // duplicate d_sphParticleVec so that it is not destroyed by
    // d_allSPHParticleVec in next iteration,
    // otherwise it causes memory error.
    std::vector<sph::SPHParticle*> dupSPHParticleVec(d_sphParticleVec.size());
    for (std::size_t i = 0; i < dupSPHParticleVec.size(); ++i)
      dupSPHParticleVec[i] = new sph::SPHParticle(*d_sphParticleVec[i]);

    // fill allParticleVec with dupParticleVec and received particles
    d_allSPHParticleVec.insert(d_allSPHParticleVec.end(),
                               dupSPHParticleVec.begin(),
                               dupSPHParticleVec.end());

    std::vector<sph::SPHParticle*> tmpSPHParticleVec;
    long gatherRam = 0;
    for (int iRank = 1; iRank < d_mpiSize; ++iRank) {

      tmpSPHParticleVec.clear(); // do not destroy particles!
      boostWorld.recv(iRank, mpiTag, tmpSPHParticleVec);
      d_allSPHParticleVec.insert(d_allSPHParticleVec.end(),
                                 tmpSPHParticleVec.begin(),
                                 tmpSPHParticleVec.end());
      gatherRam += tmpSPHParticleVec.size();
    }
    // debugInf << "gather: particleNum = " << gatherRam <<  " particleRam = "
    //<< gatherRam * sizeof(Particle) << std::endl;
  }
}

void
SmoothParticleHydro::releaseGatheredSPHParticle()
{
  // clear allParticleVec, avoid long time memory footprint.
  for (std::vector<sph::SPHParticle*>::iterator it =
         d_allSPHParticleVec.begin();
       it != d_allSPHParticleVec.end(); ++it)
    delete (*it);
  d_allSPHParticleVec.clear();
  std::vector<sph::SPHParticle*>().swap(d_allSPHParticleVec); // actual memory
                                                              // release
}

void
SmoothParticleHydro::gatherBdryContact()
{
  if (isBdryProcess()) {
    if (d_mpiRank != 0)
      boostWorld.send(0, mpiTag, boundaryVec);
  }

  if (d_mpiRank == 0) {
    mergeBoundaryVec.clear();
    std::vector<Boundary*>().swap(mergeBoundaryVec); // actual memory release
    mergeBoundaryVec = boundaryVec;

    std::vector<Boundary*> tmpBoundaryVec;
    for (std::size_t it = 0; it < bdryProcess.size(); ++it) {
      if (bdryProcess[it] != 0) { // not root process
        tmpBoundaryVec.clear();   // do not destroy particles!
        boostWorld.recv(bdryProcess[it], mpiTag, tmpBoundaryVec);
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
    for (std::vector<Boundary*>::iterator it = mergeBoundaryVec.begin();
         it != mergeBoundaryVec.end(); ++it)
      (*it)->updateStatForce();
  }
}

void
SmoothParticleHydro::openSPHTecplot(std::ofstream& ofs, const char* str)
{
  ofs.open(str);
  if (!ofs) {
    debugInf << "stream error: openSPHTecplot" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  ofs << "Title = \"SPH Particle Information\"" << std::endl;
  ofs << "VARIABLES = \"x\", \"y\",\"z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\" \"Vy\"
\"Vz\" \"Pressure\" \"a_x\" \"a_y\" \"a_z\" \"density_dot\" \"density\" "
      << std::endl;
}

void
SmoothParticleHydro::printSPHTecplot(std::ofstream& ofs, int iframe)
{
  ofs << "ZONE T =\" " << iframe << "-th Load Step\" " << std::endl;
  // Output the coordinates and the array information
  for (std::vector<sph::SPHParticle*>::iterator pt =
         d_allSPHParticleVec.begin();
       pt != d_allSPHParticleVec.end(); pt++) {

    //// not print the most right layer of SPH free particles, 2015/05/19
    // if((*pt)->getInitPosition().getx()==25){
    // continue;
    //}
    //      if((*pt)->getType()==3) continue;  // not print out boundary
    //      particles

    ofs << std::setw(20) << (*pt)->getCurrPosition().x() << std::setw(20)
        << (*pt)->getCurrPosition().y() << std::setw(20)
        << (*pt)->getCurrPosition().z() << std::setw(20)
        << (*pt)->getDisplacement().x() << std::setw(20)
        << (*pt)->getDisplacement().y() << std::setw(20)
        << (*pt)->getDisplacement().z() << std::setw(20)
        << (*pt)->getVelocity().x() << std::setw(20) << (*pt)->getVelocity().y()
        << std::setw(20) << (*pt)->getVelocity().z();

    (*pt)->calculateParticlePressure();
    ofs << std::setw(20) << (*pt)->getParticlePressure();

    ofs << std::setw(20) << (*pt)->getVelocityDot().x() << std::setw(20)
        << (*pt)->getVelocityDot().y() << std::setw(20)
        << (*pt)->getVelocityDot().z() << std::setw(20)
        << (*pt)->getDensityDot() << std::setw(20)
        << (*pt)->getParticleDensity() << std::endl;
  }
}
void
SmoothParticleHydro::printSPHParticle(const char* str) const
{
  std::ofstream ofs(str);
  if (!ofs) {
    debugInf << "stream error: printSPHParticle" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  ofs << "Title = \"SPH Particle Information\"" << std::endl;
  ofs << "VARIABLES = \"x\", \"y\",\"z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\"
\"Vy\" \"Vz\" \"Pressure\" \"a_x\" \"a_y\" \"a_z\" \"density_dot\" "
  "\"density\" "
      << std::endl;

  // Output the coordinates and the array information
  for (std::vector<sph::SPHParticle*>::const_iterator pt =
         d_allSPHParticleVec.begin();
       pt != d_allSPHParticleVec.end(); pt++) {

    //// not print the most right layer of SPH free particles, 2015/05/19
    // if((*pt)->getInitPosition().getx()==25){
    // continue;
    //}
    if ((*pt)->getType() == 3)
      continue; // not print out boundary particles

    ofs << std::setw(20) << (*pt)->getCurrPosition().x() << std::setw(20)
        << (*pt)->getCurrPosition().y() << std::setw(20)
        << (*pt)->getCurrPosition().z() << std::setw(20)
        << (*pt)->getDisplacement().x() << std::setw(20)
        << (*pt)->getDisplacement().y() << std::setw(20)
        << (*pt)->getDisplacement().z() << std::setw(20)
        << (*pt)->getVelocity().x() << std::setw(20) << (*pt)->getVelocity().y()
        << std::setw(20) << (*pt)->getVelocity().z();

    (*pt)->calculateParticlePressure();
    ofs << std::setw(20) << (*pt)->getParticlePressure();

    ofs << std::setw(20) << (*pt)->getVelocityDot().x() << std::setw(20)
        << (*pt)->getVelocityDot().y() << std::setw(20)
        << (*pt)->getVelocityDot().z() << std::setw(20)
        << (*pt)->getDensityDot() << std::setw(20)
        << (*pt)->getParticleDensity() << std::endl;
  }

  ofs.close();
}

void
SmoothParticleHydro::initialSPHVelocity2D()
{
  commuParticle();
  commuSPHParticle(); // this will update container and mergeSPHParticleVec,
                      // both are needed for divideSPHDomain
  calculateSPHDensityDotVelocityDotLinkedList2D(); // calculate velocityDot and
                                                   // densityDot
  // fix the y for all SPH points
  for (std::vector<sph::SPHParticle*>::iterator pt = d_sphParticleVec.begin();
       pt != d_sphParticleVec.end(); pt++) {
    switch ((*pt)->getType()) {
      case 1: // free sph particle
        (*pt)->fixY();
        break;
      case 2: // ghost sph particle
        break;
      case 3: // boundary sph particle
        (*pt)->fixXYZ();
        break;
      default:
        std::cout << "SPH particle type of pta should be 1, 2 or 3!"
                  << std::endl;
        exit(-1);
    } // switch
  }
  initialSPHLeapFrogVelocity(); // initial velocity only for free SPH particles
                                // based on equation (4.3)

  releaseRecvParticle();
  releaseRecvSPHParticle();
} // initialSPHVelocity()

void
SmoothParticleHydro::initialSPHVelocity3D()
{
  commuParticle();
  commuSPHParticle(); // this will update container and mergeSPHParticleVec,
                      // both are needed for divideSPHDomain
  calculateSPHDensityDotVelocityDotLinkedList3D(); // calculate velocityDot and
                                                   // densityDot
  initialSPHLeapFrogVelocity(); // initial velocity only for free SPH particles
                                // based on equation (4.3)

  releaseRecvParticle();
  releaseRecvSPHParticle();
} // initialSPHVelocity()

void
SmoothParticleHydro::initialSPHVelocityCopyDEM3D()
{
  commuSPHParticle(); // this will update container and mergeSPHParticleVec,
                      // both are needed for divideSPHDomain
  calculateSPHDensityDotVelocityDotLinkedList3D(); // calculate velocityDot and
                                                   // densityDot
  initialSPHLeapFrogVelocity(); // initial velocity only for free SPH particles
                                // based on equation (4.3)

  releaseRecvSPHParticle();
} // initialSPHVelocity()

void
SmoothParticleHydro::initialSPHLeapFrogVelocity()
{ // initial particle
  // velocity only for free SPH particles based on equation (4.3)
  for (std::vector<sph::SPHParticle*>::iterator pt = d_sphParticleVec.begin();
       pt != d_sphParticleVec.end(); pt++) {
    if ((*pt)->getType() == 1)
      (*pt)->initialParticleVelocityLeapFrog();
  }
} // end initialSPHLeapFrogVelocity

void
SmoothParticleHydro::updateSPHLeapFrogPositionDensity()
{ // update
  // particle position and density based on equation (4.1)
  for (std::vector<sph::SPHParticle*>::iterator pt = d_sphParticleVec.begin();
       pt != d_sphParticleVec.end(); pt++) {
    switch ((*pt)->getType()) {
      case 1: // free sph particle
        (*pt)->updateParticlePositionDensityLeapFrog();
        break;
      case 2: // ghost sph particle
        break;
      case 3: // boundary sph particle
        (*pt)->updateParticleDensity();
        break;
      default:
        std::cout << "SPH particle type of pta should be 1, 2 or 3!"
                  << std::endl;
        exit(-1);
    } // switch
  }

  std::vector<sph::SPHParticle*>::iterator gt;
  for (std::vector<Particle*>::iterator demt = particleVec.begin();
       demt != particleVec.end(); demt++) {
    (*demt)->updateSPHGhostParticle(); // update position and density by dem
    // particle and sph particle density, and velocity
  }
} // end updateSPHLeapFrogPositionDensity

void
SmoothParticleHydro::updateSPHLeapFrogVelocity()
{ // update particle
  // velocity only for free SPH particles based on equation (4.2)
  for (std::vector<sph::SPHParticle*>::iterator pt = d_sphParticleVec.begin();
       pt != d_sphParticleVec.end(); pt++) {
    if ((*pt)->getType() == 1)
      (*pt)->updateParticleVelocityLeapFrog();
  }
} // end updateSPHLeapFrogVelocity

//  REAL factor = 1.0/(120.0*dem::PI*h*h*h);  // 3D quintic kernel factor
//  REAL factor = 7.0/(478.0*dem::PI*h*h);    // 2D quintic kernel factor
// kernel function, Vec is the position of a b, h is smoothing length
inline REAL
SmoothParticleHydro::kernelFunction(const dem::Vec& a, const dem::Vec& b)
{
  REAL rab = dem::vfabs(a - b);
  REAL s = rab * one_devide_h;
  REAL item1_5 = (3 - s) * (3 - s) * (3 - s) * (3 - s) * (3 - s); // (3-s)^5
  REAL item2_5 = (2 - s) * (2 - s) * (2 - s) * (2 - s) * (2 - s); // (2-s)^5
  REAL item3_5 = (1 - s) * (1 - s) * (1 - s) * (1 - s) * (1 - s); // (1-s)^5
  if (s < 1) {
    return factor_kernel * (item1_5 - 6 * item2_5 + 15 * item3_5);
  } else if (s < 2) {
    return factor_kernel * (item1_5 - 6 * item2_5);
  } else if (s < 3) {
    return factor_kernel * item1_5;
  } else {
    return 0.0;
  }

} // end kernelFunction

// kernel function, s is rab/h, h is smoothing length
inline REAL
SmoothParticleHydro::kernelFunction(REAL s)
{
  REAL item1_5 = (3 - s) * (3 - s) * (3 - s) * (3 - s) * (3 - s); // (3-s)^5
  REAL item2_5 = (2 - s) * (2 - s) * (2 - s) * (2 - s) * (2 - s); // (2-s)^5
  REAL item3_5 = (1 - s) * (1 - s) * (1 - s) * (1 - s) * (1 - s); // (1-s)^5
  if (s < 1) {
    return factor_kernel * (item1_5 - 6 * item2_5 + 15 * item3_5);
  } else if (s < 2) {
    return factor_kernel * (item1_5 - 6 * item2_5);
  } else if (s < 3) {
    return factor_kernel * item1_5;
  } else {
    return 0.0;
  }

} // end kernelFunction

//  REAL factor = 1.0/(120.0*dem::PI*h*h*h*h);  // 3D quintic kernel factor
//  REAL factor = 1.0/(120.0*h*h);  // 1D quintic kernel factor
//  REAL factor = 7.0/(478.0*dem::PI*h*h*h);  // 2D quintic kernel factor
// to calculate delta_aWab, where a is the position of the first particle
inline dem::Vec
SmoothParticleHydro::gradientKernelFunction(const dem::Vec& a,
                                            const dem::Vec& b)
{
  REAL rab = dem::vfabs(a - b);
  REAL s = rab * one_devide_h;
  REAL item1_4 = (3 - s) * (3 - s) * (3 - s) * (3 - s); // (3-s)^4
  REAL item2_4 = (2 - s) * (2 - s) * (2 - s) * (2 - s); // (2-s)^4
  REAL item3_4 = (1 - s) * (1 - s) * (1 - s) * (1 - s); // (1-s)^4
  if (s < 1) {
    return factor_kernel_gradient *
           (-5 * item1_4 + 30 * item2_4 - 75 * item3_4) * (a - b) /
           rab; // it is strange that s is
    // devided here. compare to Liu's SPH code (National University of
    // Singapore)
  } else if (s < 2) {
    return factor_kernel_gradient * (-5 * item1_4 + 30 * item2_4) * (a - b) /
           rab;
  } else if (s < 3) {
    return factor_kernel_gradient * (-5 * item1_4) * (a - b) / rab;
  } else {
    return dem::Vec(0.0);
  }

} // end gradientKernelFunction

// to calculate partial differential dWab_dra
inline REAL
SmoothParticleHydro::partialKernelFunction(const dem::Vec& a, const dem::Vec& b)
{
  REAL rab = dem::vfabs(a - b);
  REAL s = rab * one_devide_h;
  REAL item1_4 = (3 - s) * (3 - s) * (3 - s) * (3 - s); // (3-s)^4
  REAL item2_4 = (2 - s) * (2 - s) * (2 - s) * (2 - s); // (2-s)^4
  REAL item3_4 = (1 - s) * (1 - s) * (1 - s) * (1 - s); // (1-s)^4
  if (s < 1) {
    return factor_kernel_gradient *
           (-5 * item1_4 + 30 * item2_4 - 75 * item3_4);
  } else if (s < 2) {
    return factor_kernel_gradient * (-5 * item1_4 + 30 * item2_4);
  } else if (s < 3) {
    return factor_kernel_gradient * (-5 * item1_4);
  } else {
    return 0.0;
  }

} // end partialKernelFunction

// divide SPH domain in each cpu into different cells in 2D in xz plane.
// see the notes 5/20/2015 and 5/21/2015 or Simpson's paper "Numerical
// techniques for three-dimensional Smoothed Particle Hydrodynamics"
void
SmoothParticleHydro::divideSPHDomain2D()
{

  // clear the std::vector< std::vector<sph::SPHParticle*> > SPHParticleCellVec
  for (int pvec = 0; pvec < SPHParticleCellVec.size(); pvec++) {
    SPHParticleCellVec[pvec].clear();
  }
  SPHParticleCellVec.clear();

  Vec vmin = container.getMinCorner();
  Vec vmax = container.getMaxCorner();
  REAL small_value = 0.01 * spaceInterval;
  REAL xmin = vmin.x() - sphCellSize - small_value;
  REAL zmin = vmin.z() - sphCellSize -
              small_value; // expand the container by sphCellSize for
                           // cells domain is necessary
  REAL xmax = vmax.x() + sphCellSize + small_value;
  REAL zmax = vmax.z() + sphCellSize +
              small_value;             // since the sph domain that we divide is
                                       // mergeSPHParticleVec
  Nx = (xmax - xmin) / kernelSize + 1; // (xmax-xmin)/(3h)+1
  Nz = (zmax - zmin) / kernelSize + 1; // (zmax-zmin)/(3h)+1
  numCell = Nx * Nz;
  SPHParticleCellVec.resize(numCell); // at this point, the cellVec contains
                                      // numCell vectors,
  // each vector is empty but ready to store the pointer of SPH
  // particles (this process will be in the calculation of SPH forces)

  dem::Vec tmp_xyz;
  int num;
  for (std::vector<sph::SPHParticle*>::iterator pt =
         mergeSPHParticleVec.begin();
       pt != mergeSPHParticleVec.end();
       pt++) { // mergeSPHParticle contains also the sph
               // particles from neighboring cpus
    tmp_xyz = (*pt)->getCurrPosition();
    num = int((tmp_xyz.z() - zmin) / kernelSize) * Nx +
          int((tmp_xyz.x() - xmin) / kernelSize);
    SPHParticleCellVec[num].push_back(*pt);
  }

  REAL maxRadius = gradation.getPtclMaxRadius();
  std::vector<sph::SPHParticle*>::iterator gt;
  for (std::vector<Particle*>::iterator pdt = particleVec.begin();
       pdt != particleVec.end();
       ++pdt) { // the sph ghost particles in the current cpu's
                // dem particles
    // will be definitely inside the [xmin, ymin, zmin,
    // xmax ...]
    tmp_xyz = (*pdt)->currentPosition();
    if (tmp_xyz.x() >= xmax + maxRadius || tmp_xyz.z() >= zmax + maxRadius ||
        tmp_xyz.x() <= xmin - maxRadius ||
        tmp_xyz.z() <= zmin - maxRadius) // dem
      particle is outside of the continue;
    for (gt = (*pdt)->SPHGhostParticleVec.begin();
         gt != (*pdt)->SPHGhostParticleVec.end();
         gt++) { // set all SPH ghost particles
      into their cells tmp_xyz = (*gt)->getCurrPosition();
      if (tmp_xyz.x() >= xmax || tmp_xyz.z() >= zmax || tmp_xyz.x() <= xmin ||
          tmp_xyz.z() <= zmin)
        continue; // this sph ghost particle is outside of the container, go to
      next sph ghost particle num =
        int((tmp_xyz.z() - zmin) / kernelSize) * Nx +
        int((tmp_xyz.x() - xmin) / kernelSize);
      SPHParticleCellVec[num].push_back(*gt);
    }
  }

  for (std::vector<Particle*>::iterator pdt = recvParticleVec.begin();
       pdt != recvParticleVec.end();
       pdt++) { // the this time, recvParticleVec are the
    dem particles that are communicated by the commuParticle tmp_xyz =
      (*pdt)->currentPosition();
    if (tmp_xyz.x() >= xmax + maxRadius || tmp_xyz.z() >= zmax + maxRadius ||
        tmp_xyz.x() <= xmin - maxRadius ||
        tmp_xyz.z() <= zmin - maxRadius)   // dem
      particle is outside of the continue; // expanded domain

    for (gt = (*pdt)->SPHGhostParticleVec.begin();
         gt != (*pdt)->SPHGhostParticleVec.end();
         gt++) { // set part SPH ghost particles
      into their cells tmp_xyz = (*gt)->getCurrPosition();
      // not all sph ghost particles in these received particles are inside the
      container if (tmp_xyz.x() >= xmax || tmp_xyz.z() >= zmax ||
                    tmp_xyz.x() <= xmin ||
                    tmp_xyz.z() <= zmin) continue; // this sph ghost particle is
                                                   // outside of the container,
                                                   // go to
      next sph ghost particle num =
        int((tmp_xyz.z() - zmin) / kernelSize) * Nx +
        int((tmp_xyz.x() - xmin) / kernelSize);
      SPHParticleCellVec[num].push_back(*gt);
    }
  }

} // divideSPHDomain2D

//// the momentum equilibrium equation and state equation are implemented as
// in Monaghan's paper (1994), simulate free surface flow using sph
//// here the neighboring list of SPH particles is searched by the cells,
vector<vector<sph::SPHParticle*>> SPHParticleCellVec;
void
SmoothParticleHydro::calculateSPHDensityDotVelocityDotLinkedList2D()
{

  // divide the SPH domain into different cells, each cell will contain SPH
  // particles within it
  divideSPHDomain2D();

  //  checkDivision();  // pass, May 22, 2015
  //  checkNeighborCells2D();  // pass, May 22, 2015

  // initialize the densityDot and velocityDot of all the SPH particles
  dem::Vec tmp_vec = dem::Vec(0,0,
-(util::getParam<REAL>("gravAccel")*(util::getParam<REAL>("gravScale") );
  dem::Vec zero_vec = dem::Vec(0,0,0);
  // in the calculation of forces between sph particles, we need to use the
//mergeSPHParticleVec,
  // since mergeVec contains the sph particles from neighboring cpus. While we
//can only update and migrate and communicate d_sphParticleVec
  for(std::vector<sph::SPHParticle*>::iterator pt=mergeSPHParticleVec.begin();
pt!=mergeSPHParticleVec.end(); pt++){
    (*pt)->setDensityDotVelocityDotZero();
    (*pt)->calculateParticleViscosity(); // everytime when use
    getParticleViscolisity(),
      make sure that pressure has been calculated !!!!(*pt)
        ->calculateParticlePressure(); // everytime when use
    getParticlePressure(),
      make sure that pressure has been calculated !!!!switch ((*pt)->getType())
    {
      case 1: // free sph particle
        (*pt)->addVelocityDot(tmp_vec);
        break;
      case 2: // ghost sph particle
        break;
      case 3: // boundary sph particle
        (*pt)->addVelocityDot(zero_vec);
        break;
      default:
        std::cout << "SPH particle type of pta should be 1, 2 or 3!"
                  << std::endl;
        exit(-1);
    } // switch
  }

  std::vector<sph::SPHParticle*>::iterator gt;
  for(std::vector<Particle*>::iterator pdt=mergeParticleVec.begin();
pdt!=mergeParticleVec.end(); pdt++){  // all sph ghost particles
    for (gt = (*pdt)->SPHGhostParticleVec.begin();
         gt != (*pdt)->SPHGhostParticleVec.end(); gt++) {
      (*gt)->setDensityDotVelocityDotZero();
      (*gt)->calculateParticleViscosity(); // everytime when use
      // getParticleViscolisity(), make sure that pressure has been
      // calculated!!!!
      (*gt)->calculateParticlePressure(); // everytime when use
      // getParticlePressure(), make sure that pressure has been calculated!!!!
    }
  }

  // temporary variables used in the loop
  dem::Vec pta_position;
  dem::Vec ptb_position;
  dem::Vec delta_aWab;
  dem::Vec delta_bWba;
  dem::Vec vab;
  dem::Vec vba;
  dem::Vec vdem;
  dem::Vec dva_dt;
  dem::Vec dvb_dt;
  dem::Vec delta_a;
  dem::Vec delta_b;
  REAL pa, pb, rhoa, rhob, mua, mub;
  REAL rab;
        REAL dWab_dra;
  REAL dWba_drb;
  REAL Wab, Wba;
  REAL da, dB;
  REAL beta;
  REAL xa, ya, xB, yB, k, sqrt_xaya;  // variables for Morris' method to
//calculate da/dB
  std::vector<sph::SPHParticle*>::iterator ptb;
  REAL ra, rb, rc;  // the geometry of the dem particle
  dem::Vec pta_local, ptb_local;  // local position of pta and ptb in the dem
//particle
  dem::Vec demt_curr;
  REAL Gamma_ab, mu_ab, vr_dot;
  REAL alpha = util::getParam<REAL>("alpha");
    REAL alpha_zero = 0;  // the viscous between free and ghost/boundary
  REAL epsilon = util::getParam<REAL>("epsilon");  // parameter for velocity
//correction
  REAL Wq, Ra, Rb, phi_4, coefficient_a, coefficient_b;
  dem::Particle* demt;

        int pnum;
  std::vector<sph::SPHParticle*> tmp_particleVec;  // sph particles in
//neighboring cells
  for(int pvec=0; pvec<SPHParticleCellVec.size(); ++pvec){

    // store all SPH particles in pvec's neighboring cells
    tmp_particleVec.clear(); // it is the same for the particles in the same
                             // cell
    if (pvec + 1 < numCell) {
      pnum = pvec + 1;
      for (std::vector<sph::SPHParticle*>::iterator pt =
             SPHParticleCellVec[pnum].begin();
           pt != SPHParticleCellVec[pnum].end(); pt++) {
        tmp_particleVec.push_back(*pt);
      }
    }
    if (pvec + Nx - 1 < numCell) {
      pnum = pvec + Nx - 1;
      for (std::vector<sph::SPHParticle*>::iterator pt =
             SPHParticleCellVec[pnum].begin();
           pt != SPHParticleCellVec[pnum].end(); pt++) {
        tmp_particleVec.push_back(*pt);
      }
    }
    if (pvec + Nx < numCell) {
      pnum = pvec + Nx;
      for (std::vector<sph::SPHParticle*>::iterator pt =
             SPHParticleCellVec[pnum].begin();
           pt != SPHParticleCellVec[pnum].end(); pt++) {
        tmp_particleVec.push_back(*pt);
      }
    }
    if (pvec + Nx + 1 < numCell) {
      pnum = pvec + Nx + 1;
      for (std::vector<sph::SPHParticle*>::iterator pt =
             SPHParticleCellVec[pnum].begin();
           pt != SPHParticleCellVec[pnum].end(); pt++) {
        tmp_particleVec.push_back(*pt);
      }
    }

    for (std::vector<sph::SPHParticle*>::iterator pta =
           SPHParticleCellVec[pvec].begin();
         pta != SPHParticleCellVec[pvec].end();
         ++pta) { // SPH particles in cell pvec
      pa = (*pta)->getParticlePressure();
      if (pa >= 0) {
        Ra = 0.006;
      } else {
        Ra = 0.6;
      }
      //        mua = (*pta)->getParticleViscosity();
      rhoa = (*pta)->getParticleDensity();
      pta_position = (*pta)->getCurrPosition();

      //    if((*pta)->getType()!=1){  // pta is not free sph particles, i.e.
      //    pta is
      // ghost or boundary particles
      //        continue;  // do not consider pta, as we treat before
      //    }

      for (ptb = pta + 1; ptb != SPHParticleCellVec[pvec].end();
           ++ptb) { // sum over the
                    // SPH particles in the same cell pvec
        ptb_position = (*ptb)->getCurrPosition();
        rab = dem::vfabs(pta_position - ptb_position);
        if (rab <= kernelSize) { // ptb is in the smooth kernel
          pb = (*ptb)->getParticlePressure();
          //          mub = (*ptb)->getParticleViscosity();
          rhob = (*ptb)->getParticleDensity();

          Wq = kernelFunction(rab / smoothLength);
          phi_4 = pow(Wq / Wqmin, 4);
          if (pb >= 0) {
            Rb = 0.006;
          } else {
            Rb = 0.6;
          }
          coefficient_a = 1 + Ra * phi_4;
          coefficient_b = 1 + Rb * phi_4;

          // we have three types of SPH particles: 1, free particle; 2, ghost
          // particle; 3, boundary particle
          // Then we have 3x3 = 9 different types of interactions with the three
          // types of particles
          // so we cannot judge the interaction type only by the type of
          // particle
          // ptb, we need to consider pta also
          // pta          ptb          need to consider or not
          // 1            1              V      free with free
          // 1            2              V      free with ghost
          // 1            3              V      free with boundary

          // 2            1              V      ghost with free
          // 2            2              X      ghost with ghost
          // 2            3              X      ghost with boundary

          // 3            1              V      boundary with free
          // 3            2              X      boundary with ghost
          // 3            3              X       boundary with boundary

          // add the density dot for pta and ptb
          delta_aWab = gradientKernelFunction(pta_position, ptb_position); //
          // this is to add SPH pta
          delta_bWba = -delta_aWab;
          dWab_dra = partialKernelFunction(pta_position, ptb_position); // this
          // is to add SPH pta
          dWba_drb = dWab_dra;
          Wab = kernelFunction(pta_position, ptb_position); // this is to add
                                                            // SPH pta
          Wba = Wab;

          switch ((*ptb)->getType()) {
            case 1: // ptb is free SPH particle
              switch ((*pta)->getType()) {
                case 1: // free with free
                  vab = (*pta)->getVelocity() - (*ptb)->getVelocity();
                  vba = -vab;

                  (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                        (vab * delta_aWab));
                  (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                        (vba * delta_bWba));

                  vr_dot = vab * (pta_position - ptb_position);
                  if (vr_dot < 0) {
                    mu_ab = smoothLength * vr_dot /
                            (rab * rab + 0.01 * smoothLength * smoothLength);
                Gamma_ab =
(-alpha*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rhoa+rhob)*2;
                  } else {
                    Gamma_ab = 0;
                  }

                  dva_dt = -(*ptb)->getParticleMass() *
                           (pa / (rhoa * rhoa) * coefficient_a +
                            pb / (rhob * rhob) * coefficient_b + Gamma_ab) *
                           delta_aWab;
                  (*pta)->addVelocityDot(dva_dt);
                  dvb_dt = -(*pta)->getParticleMass() *
                           (pb / (rhob * rhob) * coefficient_b +
                            pa / (rhoa * rhoa) * coefficient_a + Gamma_ab) *
                           delta_bWba;
                  (*ptb)->addVelocityDot(dvb_dt);

                  delta_a = epsilon * (*ptb)->getParticleMass() * (-vab) * Wab /
                            (rhoa + rhob) * 2;
                  (*pta)->addVelocityCorrection(delta_a);
                  delta_b = epsilon * (*pta)->getParticleMass() * (vab)*Wba /
                            (rhoa + rhob) * 2;
                  (*ptb)->addVelocityCorrection(delta_b);
                  break;
                case 2: // ghost with free
                  demt = (*pta)->getDemParticle();
                  demt_curr = demt->currentPosition();
                  ptb_local =
                    demt->globalToLocal(ptb_position - demt_curr); // the
                  // local position of sph point ptb
                  ra = demt->getA();
                  rc = demt->getC();
                  k = 1.0 / (sqrt(ptb_local.x() * ptb_local.x() / (ra * ra) +
                                  ptb_local.z() * ptb_local.z() / (rc * rc)));
                  da = dem::vfabs(ptb_local -
                                  k * ptb_local); // the distance is the same
                                                  // in rotated coordinates

                  // calculate Vab as the method shown in Morris's paper, 1996

                  // (1) here I wanna use the distance from the point a/b to the
                  // surface of the ellipsoid to simplify the problem
                  pta_local = (*pta)->getLocalPosition();
                  k = 1.0 / (sqrt(pta_local.x() * pta_local.x() / (ra * ra) +
                                  pta_local.z() * pta_local.z() / (rc * rc)));
                  dB = dem::vfabs(pta_local - k * pta_local);

                  //              // (2) here the Morris's method is used
                  //              xa = pta_position.getx(); ya =
                  //              pta_position.gety();
                  //              xB = ptb_position.getx(); yB =
                  //              ptb_position.gety();
                  //              if(ya==0) {da = xa-radius; dB = radius-xB;}
                  //              else if(xa==0) {da = ya-radius; dB =
                  //              radius-yB;}
                  //              else {
                  //          sqrt_xaya = sqrt(xa*xa+ya*ya);
                  //          k = radius/sqrt_xaya;
                  //          da = radius/k - radius;
                  //          dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
                  //              }

                  beta = 1 + dB / da;
                  if (beta > 2.0 || beta < 0 || isnan(beta)) {
                    beta = 2.0;
                  }

                  vdem = (*pta)->getVelocity();
                  vba = beta * ((*ptb)->getVelocity() - vdem);
                  vab = -vba;

                  (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                        (vab * delta_aWab));
                  (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                        (vba * delta_bWba));

                  vr_dot = vab * (pta_position - ptb_position);
                  if (vr_dot < 0) {
                    mu_ab = smoothLength * vr_dot /
                            (rab * rab + 0.01 * smoothLength * smoothLength);
                    Gamma_ab =
(-alpha*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rhoa+rhob)*2;
                  } else {
                    Gamma_ab = 0;
                  }

                  dva_dt = -(*ptb)->getParticleMass() *
                           (pa / (rhoa * rhoa) * coefficient_a +
                            pb / (rhob * rhob) * coefficient_b + Gamma_ab) *
                           delta_aWab;
                  demt->addForce((*pta)->getParticleMass() * dva_dt);
                  demt->addMoment((pta_position - demt_curr) %
                                  ((*pta)->getParticleMass() * dva_dt));
                  //                (*pta)->addVelocityDot(dva_dt);

                  dvb_dt = -(*pta)->getParticleMass() *
                           (pb / (rhob * rhob) * coefficient_b +
                            pa / (rhoa * rhoa) * coefficient_a + Gamma_ab) *
                           delta_bWba;
                  (*ptb)->addVelocityDot(dvb_dt); // the velocities of the ghost
                  // particles will not envolve as the same way as others

                  //                  delta_a =
                  epsilon*(*ptb)->getParticleMass() * (-vab) * Wab /
                    (rhoa + rhob) * 2;
                  //                  (*pta)->addVelocityCorrection(delta_a);
                  delta_b = epsilon * (*pta)->getParticleMass() * (vab)*Wba /
                            (rhoa + rhob) * 2;
                  (*ptb)->addVelocityCorrection(delta_b);

                  break;

                case 3: // boundary with free

                  // calculate Vab as the method shown in Morris's paper, 1996
                  // interact with boundary particles
                  da = ptb_position.z() - allContainer.getMinCorner().z(); //
                  // assume with the bottom boundary
                  dB = allContainer.getMinCorner().z() -
                       pta_position.x(); // assume with
                                         // the bottom boundary
                  if (pta_position.x() <
                      allContainer.getMinCorner().x()) { // with left
                                                         // boundary
                    da = ptb_position.x() - allContainer.getMinCorner().x();
                    dB = allContainer.getMinCorner().x() - pta_position.x();
                  } else if (pta_position.x() >
                             allContainer.getMaxCorner().x()) { // with right
                                                                // boundary
                    da = allContainer.getMaxCorner().x() - ptb_position.x();
                    dB = pta_position.x() - allContainer.getMaxCorner().x();
                  } else if (pta_position.z() >
                             allContainer.getMaxCorner().z()) { // with top
                                                                // boundary
                    da = allContainer.getMaxCorner().z() - ptb_position.z();
                    dB = pta_position.z() - allContainer.getMaxCorner().z();
                  }

                  beta = 1 + dB / da;
                  if (beta > 2 || beta < 0 || isnan(beta)) {
                    beta = 2;
                  }

                  vba = beta * (*ptb)->getVelocity();
                  vab = -vba;

                  (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                        (vab * delta_aWab));
                  (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                        (vba * delta_bWba));

                  vr_dot = vab * (pta_position - ptb_position);
                  if (vr_dot < 0) {
                    mu_ab = smoothLength * vr_dot /
                            (rab * rab + 0.01 * smoothLength * smoothLength);
                    Gamma_ab =
                      (-alpha * (util::getParam<REAL>("soundSpeed") * mu_ab)) /
                      (rhoa + rhob) * 2;
                  } else {
                    Gamma_ab = 0;
                  }

                  //                dva_dt =
                  //-(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
                  //                (*pta)->addVelocityDot(dva_dt);
                  dvb_dt =
                    -(*pta)->getParticleMass() *
                    (pb / (rhob * rhob) + pa / (rhoa * rhoa) + Gamma_ab) *
                    delta_bWba;
                  (*ptb)->addVelocityDot(dvb_dt); // the velocities of the ghost
                  // particles will not envolve as the same way as others

                  //            delta_a =
                  // epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
                  //                (*pta)->addVelocityCorrection(delta_a);
                  delta_b = epsilon * (*pta)->getParticleMass() * (vab)*Wba /
                            (rhoa + rhob) * 2;
                  (*ptb)->addVelocityCorrection(delta_b);

                  //            // apply the boundary forces by Lennard-Jones
                  //            potential as in
                  // Monaghan's paper(1994)
                  //                if(rab<=spaceInterval){ // ptb is in the
                  //                smooth kernel
                  //                dvb_dt = D*(pow(spaceInterval/rab,
                  //                p1)-pow(spaceInterval/rab,
                  // p2))*(ptb_position-pta_position)/(rab*rab);
                  //                (*ptb)->addVelocityDot(dvb_dt);
                  //                } // end if

                  break;
                default:
                  std::cout << "SPH particle type of pta should be 1, 2 or 3!"
                            << std::endl;
                  exit(-1);
              } // end switch pta

              break;
            case 2: // ptb is ghost particle

              if ((*pta)->getType() !=
                  1) { // pta is not free sph particles, i.e. pta
                       // is ghost or boundary particles
                break; // do not consider pta, as we treat before
              }
              demt = (*ptb)->getDemParticle();
              demt_curr = demt->currentPosition();
              pta_local =
                demt->globalToLocal(pta_position - demt_curr); // the local
              // position of sph point pta
              ra = demt->getA();
              rc = demt->getC();
              k = 1.0 / (sqrt(pta_local.x() * pta_local.x() / (ra * ra) +
                              pta_local.z() * pta_local.z() / (rc * rc)));
              da = dem::vfabs(pta_local -
                              k * pta_local); // the distance is the same
                                              // in rotated coordinates

              // calculate Vab as the method shown in Morris's paper, 1996

              // (1) here I wanna use the distance from the point a/b to the
              // surface of the ellipsoid to simplify the problem
              ptb_local = (*ptb)->getLocalPosition();
              k = 1.0 / (sqrt(ptb_local.x() * ptb_local.x() / (ra * ra) +
                              ptb_local.z() * ptb_local.z() / (rc * rc)));
              dB = dem::vfabs(ptb_local - k * ptb_local);

              //              // (2) here the Morris's method is used
              //              xa = pta_position.getx(); ya =
              //              pta_position.gety();
              //              xB = ptb_position.getx(); yB =
              //              ptb_position.gety();
              //              if(ya==0) {da = xa-radius; dB = radius-xB;}
              //              else if(xa==0) {da = ya-radius; dB = radius-yB;}
              //              else {
              //          sqrt_xaya = sqrt(xa*xa+ya*ya);
              //          k = radius/sqrt_xaya;
              //          da = radius/k - radius;
              //          dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
              //              }

              beta = 1 + dB / da;
              if (beta > 2 || beta < 0 || isnan(beta)) {
                beta = 2;
              }

              vdem = (*ptb)->getVelocity();
              vab = beta * ((*pta)->getVelocity() - vdem);
              vba = -vab;

              (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                    (vab * delta_aWab));
              (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                    (vba * delta_bWba));

              vr_dot = vab * (pta_position - ptb_position);
              if (vr_dot < 0) {
                mu_ab = smoothLength * vr_dot /
                        (rab * rab + 0.01 * smoothLength * smoothLength);
                Gamma_ab =
(-alpha*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rhoa+rhob)*2;
              } else {
                Gamma_ab = 0;
              }

              dva_dt = -(*ptb)->getParticleMass() *
                       (pa / (rhoa * rhoa) * coefficient_a +
                        pb / (rhob * rhob) * coefficient_b + Gamma_ab) *
                       delta_aWab;
              (*pta)->addVelocityDot(dva_dt);
              dvb_dt = -(*pta)->getParticleMass() *
                       (pb / (rhob * rhob) * coefficient_b +
                        pa / (rhoa * rhoa) * coefficient_a + Gamma_ab) *
                       delta_bWba;
              demt->addForce((*ptb)->getParticleMass() * dvb_dt);
              demt->addMoment((ptb_position - demt_curr) %
                              ((*ptb)->getParticleMass() * dvb_dt));
              //              (*ptb)->addVelocityDot(dvb_dt);  // the velocities
              //              of the ghost
              // particles will not envolve as the same way as others

              delta_a = epsilon * (*ptb)->getParticleMass() * (-vab) * Wab /
                        (rhoa + rhob) * 2;
              (*pta)->addVelocityCorrection(delta_a);
              //                delta_b =
              // epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
              //                (*ptb)->addVelocityCorrection(delta_b);

              break;
            case 3: // ptb is boundary particle

              if ((*pta)->getType() !=
                  1) { // pta is not free sph particles, i.e. pta
                       // is ghost or boundary particles
                break; // do not consider pta, as we treat before
              }

              // calculate Vab as the method shown in Morris's paper, 1996
              // interact with boundary particles
              da = pta_position.z() - allContainer.getMinCorner().z(); // assume
              // with the bottom boundary
              dB = allContainer.getMinCorner().z() -
                   ptb_position.z(); // assume with
                                     // the bottom boundary
              if (ptb_position.x() <
                  allContainer.getMinCorner().x()) { // with left
                                                     // boundary
                da = pta_position.x() - allContainer.getMinCorner().x();
                dB = allContainer.getMinCorner().x() - ptb_position.x();
              } else if (ptb_position.x() >
                         allContainer.getMaxCorner().x()) { // with
                                                            // right boundary
                da = allContainer.getMaxCorner().x() - pta_position.x();
                dB = ptb_position.x() - allContainer.getMaxCorner().x();
              } else if (ptb_position.z() >
                         allContainer.getMaxCorner().z()) { // with
                                                            // top boundary
                da = allContainer.getMaxCorner().z() - pta_position.z();
                dB = ptb_position.z() - allContainer.getMaxCorner().z();
              }

              beta = 1 + dB / da;
              if (beta > 2 || beta < 0 || isnan(beta)) {
                beta = 2;
              }

              vab = beta * (*pta)->getVelocity();
              vba = -vab;

              (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                    (vab * delta_aWab));
              (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                    (vba * delta_bWba));

              vr_dot = vab * (pta_position - ptb_position);
              if (vr_dot < 0) {
                mu_ab = smoothLength * vr_dot /
                        (rab * rab + 0.01 * smoothLength * smoothLength);
            Gamma_ab =
(-alpha*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rhoa+rhob)*2;
              } else {
                Gamma_ab = 0;
              }

              dva_dt = -(*ptb)->getParticleMass() *
                       (pa / (rhoa * rhoa) * coefficient_a +
                        pb / (rhob * rhob) * coefficient_b + Gamma_ab) *
                       delta_aWab;
              (*pta)->addVelocityDot(dva_dt);
              //              dvb_dt =
              //-(*pta)->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
              //              (*ptb)->addVelocityDot(dvb_dt);  // the velocities
              //              of the ghost
              particles will not envolve as the same way as others

                delta_a = epsilon * (*ptb)->getParticleMass() * (-vab) * Wab /
                          (rhoa + rhob) * 2;
              (*pta)->addVelocityCorrection(delta_a);
              //              delta_b =
              // epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
              //              (*ptb)->addVelocityCorrection(delta_b);

              //          // apply the boundary forces by Lennard-Jones
              //          potential as in
              // Monaghan's paper(1994)
              //              if(rab<=spaceInterval){ // ptb is in the smooth
              //              kernel
              //            dva_dt = D*(pow(spaceInterval/rab,
              //            p1)-pow(spaceInterval/rab,
              // p2))*(pta_position-ptb_position)/(rab*rab);
              //            (*pta)->addVelocityDot(dva_dt);
              //              } // end if
              break;
            default:
              std::cout << "SPH particle type should be 1, 2 or 3!"
                        << std::endl;
              exit(-1);

          } // end swtich type
        }   // end if 3h
      }     // end for ptb in the same cell

      for (ptb = tmp_particleVec.begin(); ptb != tmp_particleVec.end();
           ptb++) { // all
                    // particles in pvec's neighboring cells
        ptb_position = (*ptb)->getCurrPosition();
        rab = dem::vfabs(pta_position - ptb_position);
        if (rab <= kernelSize) { // ptb is in the smooth kernel
          pb = (*ptb)->getParticlePressure();
          //          mub = (*ptb)->getParticleViscosity();
          rhob = (*ptb)->getParticleDensity();

          Wq = kernelFunction(rab / smoothLength);
          phi_4 = pow(Wq / Wqmin, 4);
          if (pb >= 0) {
            Rb = 0.006;
          } else {
            Rb = 0.6;
          }
          coefficient_a = 1 + Ra * phi_4;
          coefficient_b = 1 + Rb * phi_4;

          // add the density dot for pta and ptb
          delta_aWab = gradientKernelFunction(pta_position, ptb_position); //
          // this is to add SPH pta
          delta_bWba = -delta_aWab;
          dWab_dra = partialKernelFunction(pta_position, ptb_position); // this
          // is to add SPH pta
          dWba_drb = dWab_dra;
          Wab = kernelFunction(pta_position, ptb_position); // this is to add
                                                            // SPH pta
          Wba = Wab;
          switch ((*ptb)->getType()) {
            case 1: // ptb is free SPH particle
              switch ((*pta)->getType()) {
                case 1: // free with free
                  vab = (*pta)->getVelocity() - (*ptb)->getVelocity();
                  vba = -vab;

                  (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                        (vab * delta_aWab));
                  (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                        (vba * delta_bWba));

                  vr_dot = vab * (pta_position - ptb_position);
                  if (vr_dot < 0) {
                    mu_ab = smoothLength * vr_dot /
                            (rab * rab + 0.01 * smoothLength * smoothLength);
                Gamma_ab =
(-alpha*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rhoa+rhob)*2;
                  } else {
                    Gamma_ab = 0;
                  }

                  dva_dt = -(*ptb)->getParticleMass() *
                           (pa / (rhoa * rhoa) * coefficient_a +
                            pb / (rhob * rhob) * coefficient_b + Gamma_ab) *
                           delta_aWab;
                  (*pta)->addVelocityDot(dva_dt);
                  dvb_dt = -(*pta)->getParticleMass() *
                           (pb / (rhob * rhob) * coefficient_b +
                            pa / (rhoa * rhoa) * coefficient_a + Gamma_ab) *
                           delta_bWba;
                  (*ptb)->addVelocityDot(dvb_dt);

                  delta_a = epsilon * (*ptb)->getParticleMass() * (-vab) * Wab /
                            (rhoa + rhob) * 2;
                  (*pta)->addVelocityCorrection(delta_a);
                  delta_b = epsilon * (*pta)->getParticleMass() * (vab)*Wba /
                            (rhoa + rhob) * 2;
                  (*ptb)->addVelocityCorrection(delta_b);
                  break;
                case 2: // ghost with free
                  demt = (*pta)->getDemParticle();
                  demt_curr = demt->currentPosition();
                  ptb_local =
                    demt->globalToLocal(ptb_position - demt_curr); // the
                  // local position of sph point ptb
                  ra = demt->getA();
                  rc = demt->getC();
                  k = 1.0 / (sqrt(ptb_local.x() * ptb_local.x() / (ra * ra) +
                                  ptb_local.z() * ptb_local.z() / (rc * rc)));
                  da = dem::vfabs(ptb_local -
                                  k * ptb_local); // the distance is the same
                                                  // in rotated coordinates

                  // calculate Vab as the method shown in Morris's paper, 1996

                  // (1) here I wanna use the distance from the point a/b to the
                  // surface of the ellipsoid to simplify the problem
                  pta_local = (*pta)->getLocalPosition();
                  k = 1.0 / (sqrt(pta_local.x() * pta_local.x() / (ra * ra) +
                                  pta_local.z() * pta_local.z() / (rc * rc)));
                  dB = dem::vfabs(pta_local - k * pta_local);

                  //              // (2) here the Morris's method is used
                  //              xa = pta_position.getx(); ya =
                  //              pta_position.gety();
                  //              xB = ptb_position.getx(); yB =
                  //              ptb_position.gety();
                  //              if(ya==0) {da = xa-radius; dB = radius-xB;}
                  //              else if(xa==0) {da = ya-radius; dB =
                  //              radius-yB;}
                  //              else {
                  //          sqrt_xaya = sqrt(xa*xa+ya*ya);
                  //          k = radius/sqrt_xaya;
                  //          da = radius/k - radius;
                  //          dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
                  //              }

                  beta = 1 + dB / da;
                  if (beta > 2 || beta < 0 || isnan(beta)) {
                    beta = 2;
                  }

                  vdem = (*pta)->getVelocity();
                  vba = beta * ((*ptb)->getVelocity() - vdem);
                  vab = -vba;

                  (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                        (vab * delta_aWab));
                  (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                        (vba * delta_bWba));

                  vr_dot = vab * (pta_position - ptb_position);
                  if (vr_dot < 0) {
                    mu_ab = smoothLength * vr_dot /
                            (rab * rab + 0.01 * smoothLength * smoothLength);
                    Gamma_ab =
(-alpha*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rhoa+rhob)*2;
                  } else {
                    Gamma_ab = 0;
                  }

                  dva_dt = -(*ptb)->getParticleMass() *
                           (pa / (rhoa * rhoa) * coefficient_a +
                            pb / (rhob * rhob) * coefficient_b + Gamma_ab) *
                           delta_aWab;
                  demt->addForce((*pta)->getParticleMass() * dva_dt);
                  demt->addMoment((pta_position - demt_curr) %
                                  ((*pta)->getParticleMass() * dva_dt));
                  //                (*pta)->addVelocityDot(dva_dt);

                  dvb_dt = -(*pta)->getParticleMass() *
                           (pb / (rhob * rhob) * coefficient_b +
                            pa / (rhoa * rhoa) * coefficient_a + Gamma_ab) *
                           delta_bWba;
                  (*ptb)->addVelocityDot(dvb_dt); // the velocities of the ghost
                  particles will not envolve as the same way as others

                    //                  delta_a =
                    epsilon *
                    (*ptb)->getParticleMass() * (-vab) * Wab / (rhoa + rhob) *
                    2;
                  //                  (*pta)->addVelocityCorrection(delta_a);
                  delta_b = epsilon * (*pta)->getParticleMass() * (vab)*Wba /
                            (rhoa + rhob) * 2;
                  (*ptb)->addVelocityCorrection(delta_b);

                  break;

                case 3: // boundary with free

                  // calculate Vab as the method shown in Morris's paper, 1996
                  // interact with boundary particles
                  da = ptb_position.z() - allContainer.getMinCorner().z(); //
                  // assume with the bottom boundary
                  dB = allContainer.getMinCorner().z() -
                       pta_position.x(); // assume with
                                         // the bottom boundary
                  if (pta_position.x() <
                      allContainer.getMinCorner().x()) { // with left
                                                         // boundary
                    da = ptb_position.x() - allContainer.getMinCorner().x();
                    dB = allContainer.getMinCorner().x() - pta_position.x();
                  } else if (pta_position.x() >
                             allContainer.getMaxCorner().x()) { // with right
                                                                // boundary
                    da = allContainer.getMaxCorner().x() - ptb_position.x();
                    dB = pta_position.x() - allContainer.getMaxCorner().x();
                  } else if (pta_position.z() >
                             allContainer.getMaxCorner().z()) { // with top
                                                                // boundary
                    da = allContainer.getMaxCorner().z() - ptb_position.z();
                    dB = pta_position.z() - allContainer.getMaxCorner().z();
                  }

                  beta = 1 + dB / da;
                  if (beta > 2 || beta < 0 || isnan(beta)) {
                    beta = 2;
                  }

                  vba = beta * (*ptb)->getVelocity();
                  vab = -vba;

                  (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                        (vab * delta_aWab));
                  (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                        (vba * delta_bWba));

                  vr_dot = vab * (pta_position - ptb_position);
                  if (vr_dot < 0) {
                    mu_ab = smoothLength * vr_dot /
                            (rab * rab + 0.01 * smoothLength * smoothLength);
                Gamma_ab =
(-alpha*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rhoa+rhob)*2;
                  } else {
                    Gamma_ab = 0;
                  }

                  //                dva_dt =
                  //-(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
                  //                (*pta)->addVelocityDot(dva_dt);
                  dvb_dt =
                    -(*pta)->getParticleMass() *
                    (pb / (rhob * rhob) + pa / (rhoa * rhoa) + Gamma_ab) *
                    delta_bWba;
                  (*ptb)->addVelocityDot(dvb_dt); // the velocities of the ghost
                  particles will not envolve as the same way as others

                    //            delta_a =
                    // epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
                    //                (*pta)->addVelocityCorrection(delta_a);
                    delta_b = epsilon * (*pta)->getParticleMass() * (vab)*Wba /
                              (rhoa + rhob) * 2;
                  (*ptb)->addVelocityCorrection(delta_b);

                  // apply the boundary forces by Lennard-Jones potential as in
                  // Monaghan's paper(1994)
                  if (rab <= spaceInterval) { // ptb is in the smooth kernel
                    dvb_dt = D * (pow(spaceInterval / rab, p1) -
                                  pow(spaceInterval / rab, p2)) *
                             (ptb_position - pta_position) / (rab * rab);
                    (*ptb)->addVelocityDot(dvb_dt);
                  } // end if

                  break;
                default:
                  std::cout << "SPH particle type of pta should be 1, 2 or 3!"
                            << std::endl;
                  exit(-1);
              } // end switch pta

              break;
            case 2: // ptb is ghost particle

              if ((*pta)->getType() !=
                  1) { // pta is not free sph particles, i.e. pta
                       // is ghost or boundary particles
                break; // do not consider pta, as we treat before
              }
              demt = (*ptb)->getDemParticle();
              demt_curr = demt->currentPosition();
              pta_local =
                demt->globalToLocal(pta_position - demt_curr); // the local
              // position of sph point pta
              ra = demt->getA();
              rc = demt->getC();
              k = 1.0 / (sqrt(pta_local.x() * pta_local.x() / (ra * ra) +
                              pta_local.z() * pta_local.z() / (rc * rc)));
              da = dem::vfabs(pta_local -
                              k * pta_local); // the distance is the same
                                              // in rotated coordinates

              // calculate Vab as the method shown in Morris's paper, 1996

              // (1) here I wanna use the distance from the point a/b to the
              // surface of the ellipsoid to simplify the problem
              ptb_local = (*ptb)->getLocalPosition();
              k = 1.0 / (sqrt(ptb_local.x() * ptb_local.x() / (ra * ra) +
                              ptb_local.z() * ptb_local.z() / (rc * rc)));
              dB = dem::vfabs(ptb_local - k * ptb_local);

              //              // (2) here the Morris's method is used
              //              xa = pta_position.getx(); ya =
              //              pta_position.gety();
              //              xB = ptb_position.getx(); yB =
              //              ptb_position.gety();
              //              if(ya==0) {da = xa-radius; dB = radius-xB;}
              //              else if(xa==0) {da = ya-radius; dB = radius-yB;}
              //              else {
              //          sqrt_xaya = sqrt(xa*xa+ya*ya);
              //          k = radius/sqrt_xaya;
              //          da = radius/k - radius;
              //          dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
              //              }

              beta = 1 + dB / da;
              if (beta > 2 || beta < 0 || isnan(beta)) {
                beta = 2;
              }

              vdem = (*ptb)->getVelocity();
              vab = beta * ((*pta)->getVelocity() - vdem);
              vba = -vab;

              (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                    (vab * delta_aWab));
              (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                    (vba * delta_bWba));

              vr_dot = vab * (pta_position - ptb_position);
              if (vr_dot < 0) {
                mu_ab = smoothLength * vr_dot /
                        (rab * rab + 0.01 * smoothLength * smoothLength);
                Gamma_ab =
(-alpha*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rhoa+rhob)*2;
              } else {
                Gamma_ab = 0;
              }

              dva_dt = -(*ptb)->getParticleMass() *
                       (pa / (rhoa * rhoa) * coefficient_a +
                        pb / (rhob * rhob) * coefficient_b + Gamma_ab) *
                       delta_aWab;
              (*pta)->addVelocityDot(dva_dt);
              dvb_dt = -(*pta)->getParticleMass() *
                       (pb / (rhob * rhob) * coefficient_b +
                        pa / (rhoa * rhoa) * coefficient_a + Gamma_ab) *
                       delta_bWba;
              demt->addForce((*ptb)->getParticleMass() * dvb_dt);
              demt->addMoment((ptb_position - demt_curr) %
                              ((*ptb)->getParticleMass() * dvb_dt));
              //              (*ptb)->addVelocityDot(dvb_dt);  // the velocities
              //              of the ghost
              // particles will not envolve as the same way as others

              delta_a = epsilon * (*ptb)->getParticleMass() * (-vab) * Wab /
                        (rhoa + rhob) * 2;
              (*pta)->addVelocityCorrection(delta_a);
              //                delta_b =
              // epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
              //                (*ptb)->addVelocityCorrection(delta_b);

              break;

            case 3: // ptb is boundary particle

              if ((*pta)->getType() !=
                  1) { // pta is not free sph particles, i.e. pta
                       // is ghost or boundary particles
                break; // do not consider pta, as we treat before
              }

              // calculate Vab as the method shown in Morris's paper, 1996
              // interact with boundary particles
              da = pta_position.z() - allContainer.getMinCorner().z(); // assume
              // with the bottom boundary
              dB = allContainer.getMinCorner().z() -
                   ptb_position.z(); // assume with
                                     // the bottom boundary
              if (ptb_position.x() <
                  allContainer.getMinCorner().x()) { // with left
                                                     // boundary
                da = pta_position.x() - allContainer.getMinCorner().x();
                dB = allContainer.getMinCorner().x() - ptb_position.x();
              } else if (ptb_position.x() >
                         allContainer.getMaxCorner().x()) { // with
                                                            // right boundary
                da = allContainer.getMaxCorner().x() - pta_position.x();
                dB = ptb_position.x() - allContainer.getMaxCorner().x();
              } else if (ptb_position.z() >
                         allContainer.getMaxCorner().z()) { // with
                                                            // top boundary
                da = allContainer.getMaxCorner().z() - pta_position.z();
                dB = ptb_position.z() - allContainer.getMaxCorner().z();
              }

              beta = 1 + dB / da;
              if (beta > 2 || beta < 0 || isnan(beta)) {
                beta = 2;
              }

              vab = beta * (*pta)->getVelocity();
              vba = -vab;

              (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                    (vab * delta_aWab));
              (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                    (vba * delta_bWba));

              vr_dot = vab * (pta_position - ptb_position);
              if (vr_dot < 0) {
                mu_ab = smoothLength * vr_dot /
                        (rab * rab + 0.01 * smoothLength * smoothLength);
            Gamma_ab =
(-alpha*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rhoa+rhob)*2;
              } else {
                Gamma_ab = 0;
              }

              dva_dt = -(*ptb)->getParticleMass() *
                       (pa / (rhoa * rhoa) * coefficient_a +
                        pb / (rhob * rhob) * coefficient_b + Gamma_ab) *
                       delta_aWab;
              (*pta)->addVelocityDot(dva_dt);
              //              dvb_dt =
              //-(*pta)->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
              //              (*ptb)->addVelocityDot(dvb_dt);  // the velocities
              //              of the ghost
              particles will not envolve as the same way as others

                delta_a = epsilon * (*ptb)->getParticleMass() * (-vab) * Wab /
                          (rhoa + rhob) * 2;
              (*pta)->addVelocityCorrection(delta_a);
              //              delta_b =
              epsilon*(*pta)->getParticleMass() * (vab)*Wba / (rhoa + rhob) * 2;
              //              (*ptb)->addVelocityCorrection(delta_b);

              // apply the boundary forces by Lennard-Jones potential as in
              // Monaghan's paper(1994)
              if (rab <= spaceInterval) { // ptb is in the smooth kernel
                dva_dt = D * (pow(spaceInterval / rab, p1) -
                              pow(spaceInterval / rab, p2)) *
                         (pta_position - ptb_position) / (rab * rab);
                (*pta)->addVelocityDot(dva_dt);
              } // end if

              break;
            default:
              std::cout << "SPH particle type should be 1, 2 or 3!"
                        << std::endl;
              exit(-1);

          } // end swtich type
        }   // end if 3h
      }     // end for ptb in neighbor cells
    }       // end for pta

    tmp_particleVec.clear(); // clear elements in tmp-vector for particles
                             // neighboring cells, it is important

  } // end for pvec, different cells

  //      // apply the boundary forces by Lennard-Jones potential as in
  //      Monaghan's
  // paper(1994)
  //      for(ptb=SPHBoundaryParticleVec.begin();
  // ptb!=SPHBoundaryParticleVec.end(); ptb++){
  //    ptb_position = (*ptb)->getCurrPosition();
  //    rab = dem::vfabs(pta_position-ptb_position);
  //    if(rab<=spaceInterval){ // ptb is in the smooth kernel
  //        dva_dt = D*(pow(spaceInterval/rab, p1)-pow(spaceInterval/rab,
  // p2))*(pta_position-ptb_position)/(rab*rab);
  //        (*pta)->addVelocityDot(dva_dt);
  //    } // end if
  //      } // end ptb

} // end calculateSPHDensityDotVelocityDotLinkedList2D()

// divide SPH domain in each cpu into different cells in 2D in xz plane.
// see the notes 5/20/2015 and 5/21/2015 or Simpson's paper "Numerical
// techniques for three-dimensional Smoothed Particle Hydrodynamics"
void
SmoothParticleHydro::divideSPHDomain3D()
{

  // clear the std::vector< std::vector<sph::SPHParticle*> > SPHParticleCellVec
  for (int pvec = 0; pvec < SPHParticleCellVec.size(); pvec++) {
    SPHParticleCellVec[pvec].clear();
  }
  SPHParticleCellVec.clear();

  Vec vmin = container.getMinCorner();
  Vec vmax = container.getMaxCorner();
  REAL small_value = 0.01 * spaceInterval;
  REAL xmin = vmin.x() - sphCellSize - small_value;
  REAL ymin = vmin.y() - sphCellSize - small_value;
  REAL zmin = vmin.z() - sphCellSize - small_value;
  // expand the container by sphCellSize for cells domain is necessary
  REAL xmax = vmax.x() + sphCellSize + small_value;
  REAL ymax = vmax.y() + sphCellSize + small_value;
  REAL zmax = vmax.z() + sphCellSize + small_value;
  // since the sph domain that we divide is mergeSPHParticleVec
  Nx = (xmax - xmin) / kernelSize + 1; // (xmax-xmin)/(3h)+1
  Nz = (zmax - zmin) / kernelSize + 1; // (zmax-zmin)/(3h)+1
  Ny = (ymax - ymin) / kernelSize + 1; // (ymax-ymin)/(3h)+1
  numCell = Nx * Nz * Ny;
  SPHParticleCellVec.resize(numCell); // at this point, the cellVec contains
                                      // numCell vectors,
  // each vector is empty but ready to store the pointer of SPH
  // particles (this process will be in the calculation of SPH forces)
  pnum_vec.clear();
  pnum_vec.push_back(1);
  pnum_vec.push_back(Nx - 1);
  pnum_vec.push_back(Nx);
  pnum_vec.push_back(Nx + 1);
  pnum_vec.push_back(Nx * Nz);
  pnum_vec.push_back(Nx * Nz - Nx);
  pnum_vec.push_back(Nx * Nz + Nx);
  pnum_vec.push_back(Nx * Nz - 1);
  pnum_vec.push_back(Nx * Nz - 1 - Nx);
  pnum_vec.push_back(Nx * Nz - 1 + Nx);
  pnum_vec.push_back(Nx * Nz + 1);
  pnum_vec.push_back(Nx * Nz + 1 - Nx);
  pnum_vec.push_back(Nx * Nz + 1 + Nx);

  dem::Vec tmp_xyz;
  int num;
  for (std::vector<sph::SPHParticle*>::iterator pt =
         mergeSPHParticleVec.begin();
       pt != mergeSPHParticleVec.end();
       pt++) { // mergeSPHParticle contains also the sph
               // particles from neighboring cpus
    tmp_xyz = (*pt)->getCurrPosition();
    num = int((tmp_xyz.y() - ymin) / kernelSize) * Nx * Nz +
          int((tmp_xyz.z() - zmin) / kernelSize) * Nx +
          int((tmp_xyz.x() - xmin) / kernelSize);
    SPHParticleCellVec[num].push_back(*pt);
  }

  REAL maxRadius = gradation.getPtclMaxRadius();
  std::vector<sph::SPHParticle*>::iterator gt;
  for (std::vector<Particle*>::iterator pdt = particleVec.begin();
       pdt != particleVec.end();
       pdt++) { // the sph ghost particles in the current cpu's
                // dem particles
    // will be definitely inside the [xmin, ymin, zmin,
    // xmax ...]
    tmp_xyz = (*pdt)->currentPosition();
    if (tmp_xyz.x() >= xmax + maxRadius || tmp_xyz.y() >= ymax + maxRadius ||
        tmp_xyz.z() >= zmax + maxRadius || tmp_xyz.x() <= xmin - maxRadius ||
        tmp_xyz.y() <= ymin - maxRadius ||
        tmp_xyz.z() <= zmin - maxRadius) // dem particle is outside of the
      continue;
    for (gt = (*pdt)->SPHGhostParticleVec.begin();
         gt != (*pdt)->SPHGhostParticleVec.end();
         gt++) { // set all SPH ghost particles
                 // into their cells
      tmp_xyz = (*gt)->getCurrPosition();
      if (tmp_xyz.x() >= xmax || tmp_xyz.y() >= ymax || tmp_xyz.z() >= zmax ||
          tmp_xyz.x() <= xmin || tmp_xyz.y() <= ymin || tmp_xyz.z() <= zmin)
        continue; // this sph ghost particle is outside of the container, go to
                  // next sph ghost particle
      num = int((tmp_xyz.y() - ymin) / kernelSize) * Nx * Nz +
            int((tmp_xyz.z() - zmin) / kernelSize) * Nx +
            int((tmp_xyz.x() - xmin) / kernelSize);
      SPHParticleCellVec[num].push_back(*gt);
    }
  }

  for (std::vector<Particle*>::iterator pdt = recvParticleVec.begin();
       pdt != recvParticleVec.end();
       pdt++) { // the this time, recvParticleVec are the
                // dem particles that are communicated by the commuParticle
    tmp_xyz = (*pdt)->currentPosition();
    if (tmp_xyz.x() >= xmax + maxRadius || tmp_xyz.y() >= ymax + maxRadius ||
        tmp_xyz.z() >= zmax + maxRadius || tmp_xyz.x() <= xmin - maxRadius ||
        tmp_xyz.y() <= ymin - maxRadius ||
        tmp_xyz.z() <= zmin - maxRadius) // dem particle is outside of the
      continue;                          // expanded domain

    for (gt = (*pdt)->SPHGhostParticleVec.begin();
         gt != (*pdt)->SPHGhostParticleVec.end();
         gt++) { // set part SPH ghost particles
                 // into their cells
      tmp_xyz = (*gt)->getCurrPosition();
      // not all sph ghost particles in these received particles are inside the
      // container
      if (tmp_xyz.x() >= xmax || tmp_xyz.y() >= ymax || tmp_xyz.z() >= zmax ||
          tmp_xyz.x() <= xmin || tmp_xyz.y() <= ymin || tmp_xyz.z() <= zmin)
        continue; // this sph ghost particle is outside of the container, go to
      next sph ghost particle num =
        int((tmp_xyz.y() - ymin) / kernelSize) * Nx * Nz +
        int((tmp_xyz.z() - zmin) / kernelSize) * Nx +
        int((tmp_xyz.x() - xmin) / kernelSize);
      SPHParticleCellVec[num].push_back(*gt);
    }
  }

} // divideSPHDomain3D

//// the momentum equilibrium equation and state equation are implemented as
// in Monaghan's paper (1994), simulate free surface flow using sph
//// here the neighboring list of SPH particles is searched by the cells,
vector<vector<sph::SPHParticle*>> SPHParticleCellVec;
void
SmoothParticleHydro::calculateSPHDensityDotVelocityDotLinkedList3D()
{

  // divide the SPH domain into different cells, each cell will contain SPH
  // particles within it
  divideSPHDomain3D();

  //  checkDivision();  // pass, May 22, 2015
  //  checkNeighborCells3D();  // pass, May 22, 2015

  // initialize the densityDot and velocityDot of all the SPH particles
  dem::Vec tmp_vec = dem::Vec(0,0,
-(util::getParam<REAL>("gravAccel")*(util::getParam<REAL>("gravScale") );
  dem::Vec zero_vec = dem::Vec(0,0,0);
  // in the calculation of forces between sph particles, we need to use the
//mergeSPHParticleVec,
  // since mergeVec contains the sph particles from neighboring cpus. While we
//can only update and migrate and communicate d_sphParticleVec
  for(std::vector<sph::SPHParticle*>::iterator pt=mergeSPHParticleVec.begin();
pt!=mergeSPHParticleVec.end(); pt++){
    (*pt)->setDensityDotVelocityDotZero();
    (*pt)->calculateParticleViscosity(); // everytime when use
    // getParticleViscolisity(), make sure that pressure has been calculated!!!!
    (*pt)->calculateParticlePressure(); // everytime when use
    // getParticlePressure(), make sure that pressure has been calculated!!!!
    switch ((*pt)->getType()) {
      case 1: // free sph particle
        (*pt)->addVelocityDot(tmp_vec);
        break;
      case 2: // ghost sph particle
        break;
      case 3: // boundary sph particle
        (*pt)->addVelocityDot(zero_vec);
        break;
      default:
        std::cout << "SPH particle type of pta should be 1, 2 or 3!"
                  << std::endl;
        exit(-1);
    } // switch
  }

  std::vector<sph::SPHParticle*>::iterator gt;
  for(std::vector<Particle*>::iterator pdt=mergeParticleVec.begin();
pdt!=mergeParticleVec.end(); pdt++){  // all sph ghost particles
    for (gt = (*pdt)->SPHGhostParticleVec.begin();
         gt != (*pdt)->SPHGhostParticleVec.end(); gt++) {
      (*gt)->setDensityDotVelocityDotZero();
      (*gt)->calculateParticleViscosity(); // everytime when use
      // getParticleViscolisity(), make sure that pressure has been
      // calculated!!!!
      (*gt)->calculateParticlePressure(); // everytime when use
      // getParticlePressure(), make sure that pressure has been calculated!!!!
    }
  }

  // temporary variables used in the loop
  dem::Vec pta_position;
  dem::Vec ptb_position;
  dem::Vec delta_aWab;
  dem::Vec delta_bWba;
  dem::Vec vab;
  dem::Vec vba;
  dem::Vec vdem;
  dem::Vec dva_dt;
  dem::Vec dvb_dt;
  dem::Vec delta_a;
  dem::Vec delta_b;
  REAL pa, pb, rhoa, rhob, mua, mub;
  REAL rab;
        REAL dWab_dra;
  REAL dWba_drb;
  REAL Wab, Wba;
  REAL da, dB;
  REAL beta;
  REAL xa, ya, xB, yB, k, sqrt_xaya;  // variables for Morris' method to
//calculate da/dB
  std::vector<sph::SPHParticle*>::iterator ptb;
  REAL ra, rb, rc;  // the geometry of the dem particle
  dem::Vec pta_local, ptb_local;  // local position of pta and ptb in the dem
//particle
  dem::Vec demt_curr;
  REAL Gamma_ab, mu_ab, vr_dot;
  REAL alpha = util::getParam<REAL>("alpha");
     REAL alpha_zero = 0;  // the viscous between free and ghost/boundary
  REAL epsilon = util::getParam<REAL>("epsilon");  // parameter for velocity
//correction
  REAL Wq, Ra, Rb, phi_4, coefficient_a, coefficient_b;
  dem::Particle* demt;

  int pnum;
  std::vector<sph::SPHParticle*> tmp_particleVec;  // sph particles in
//neighboring cells
  for(int pvec=0; pvec<SPHParticleCellVec.size(); pvec++){

    // store all SPH particles in pvec's neighboring cells
    tmp_particleVec.clear(); // it is the same for the particles in the same
                             // cell
    for (std::vector<int>::const_iterator pint = pnum_vec.begin();
         pint != pnum_vec.end(); pint++) {
      pnum = pvec + (*pint);
      if (pnum < numCell) {
        for (std::vector<sph::SPHParticle*>::iterator pt =
               SPHParticleCellVec[pnum].begin();
             pt != SPHParticleCellVec[pnum].end(); pt++) {
          tmp_particleVec.push_back(*pt);
        }
      }
    }

    for (std::vector<sph::SPHParticle*>::iterator pta =
           SPHParticleCellVec[pvec].begin();
         pta != SPHParticleCellVec[pvec].end();
         pta++) { // SPH particles in cell pvec
      pa = (*pta)->getParticlePressure();
      if (pa >= 0) {
        Ra = 0.006;
      } else {
        Ra = 0.6;
      }
      //        mua = (*pta)->getParticleViscosity();
      rhoa = (*pta)->getParticleDensity();
      pta_position = (*pta)->getCurrPosition();

      //    if((*pta)->getType()!=1){  // pta is not free sph particles, i.e.
      //    pta is
      // ghost or boundary particles
      //        continue;  // do not consider pta, as we treat before
      //    }

      for (ptb = pta + 1; ptb != SPHParticleCellVec[pvec].end();
           ptb++) { // sum over the
                    // SPH particles in the same cell pvec
        ptb_position = (*ptb)->getCurrPosition();
        rab = dem::vfabs(pta_position - ptb_position);
        if (rab <= kernelSize) { // ptb is in the smooth kernel
          pb = (*ptb)->getParticlePressure();
          //          mub = (*ptb)->getParticleViscosity();
          rhob = (*ptb)->getParticleDensity();

          Wq = kernelFunction(rab / smoothLength);
          phi_4 = pow(Wq / Wqmin, 4);
          if (pb >= 0) {
            Rb = 0.006;
          } else {
            Rb = 0.6;
          }
          coefficient_a = 1 + Ra * phi_4;
          coefficient_b = 1 + Rb * phi_4;

          // we have three types of SPH particles: 1, free particle; 2, ghost
          // particle; 3, boundary particle
          // Then we have 3x3 = 9 different types of interactions with the three
          // types of particles
          // so we cannot judge the interaction type only by the type of
          // particle
          // ptb, we need to consider pta also
          // pta          ptb          need to consider or not
          // 1            1              V      free with free
          // 1            2              V      free with ghost
          // 1            3              V      free with boundary

          // 2            1              V      ghost with free
          // 2            2              X      ghost with ghost
          // 2            3              X      ghost with boundary

          // 3            1              V      boundary with free
          // 3            2              X      boundary with ghost
          // 3            3              X       boundary with boundary

          // add the density dot for pta and ptb
          delta_aWab = gradientKernelFunction(pta_position, ptb_position); //
          // this is to add SPH pta
          delta_bWba = -delta_aWab;
          dWab_dra = partialKernelFunction(pta_position, ptb_position); // this
          // is to add SPH pta
          dWba_drb = dWab_dra;
          Wab = kernelFunction(pta_position, ptb_position); // this is to add
                                                            // SPH pta
          Wba = Wab;

          switch ((*ptb)->getType()) {
            case 1: // ptb is free SPH particle
              switch ((*pta)->getType()) {
                case 1: // free with free
                  vab = (*pta)->getVelocity() - (*ptb)->getVelocity();
                  vba = -vab;

                  (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                        (vab * delta_aWab));
                  (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                        (vba * delta_bWba));

                  vr_dot = vab * (pta_position - ptb_position);
                  if (vr_dot < 0) {
                    mu_ab = smoothLength * vr_dot /
                            (rab * rab + 0.01 * smoothLength * smoothLength);
                Gamma_ab =
(-alpha*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rhoa+rhob)*2;
                  } else {
                    Gamma_ab = 0;
                  }

                  dva_dt = -(*ptb)->getParticleMass() *
                           (pa / (rhoa * rhoa) * coefficient_a +
                            pb / (rhob * rhob) * coefficient_b + Gamma_ab) *
                           delta_aWab;
                  (*pta)->addVelocityDot(dva_dt);
                  dvb_dt = -(*pta)->getParticleMass() *
                           (pb / (rhob * rhob) * coefficient_b +
                            pa / (rhoa * rhoa) * coefficient_a + Gamma_ab) *
                           delta_bWba;
                  (*ptb)->addVelocityDot(dvb_dt);

                  delta_a = epsilon * (*ptb)->getParticleMass() * (-vab) * Wab /
                            (rhoa + rhob) * 2;
                  (*pta)->addVelocityCorrection(delta_a);
                  delta_b = epsilon * (*pta)->getParticleMass() * (vab)*Wba /
                            (rhoa + rhob) * 2;
                  (*ptb)->addVelocityCorrection(delta_b);
                  break;
                case 2: // ghost with free
                  demt = (*pta)->getDemParticle();
                  demt_curr = demt->currentPosition();
                  ptb_local =
                    demt->globalToLocal(ptb_position - demt_curr); // the
                  // local position of sph point ptb
                  ra = demt->getA();
                  rb = demt->getB();
                  rc = demt->getC();
                  k = 1.0 / (sqrt(ptb_local.x() * ptb_local.x() / (ra * ra) +
                                  ptb_local.y() * ptb_local.y() / (rb * rb) +
                                  ptb_local.z() * ptb_local.z() / (rc * rc)));
                  da = dem::vfabs(ptb_local -
                                  k * ptb_local); // the distance is the same
                                                  // in rotated coordinates

                  // calculate Vab as the method shown in Morris's paper, 1996

                  // (1) here I wanna use the distance from the point a/b to the
                  // surface of the ellipsoid to simplify the problem
                  pta_local = (*pta)->getLocalPosition();
                  k = 1.0 / (sqrt(pta_local.x() * pta_local.x() / (ra * ra) +
                                  pta_local.y() * pta_local.y() / (rb * rb) +
                                  pta_local.z() * pta_local.z() / (rc * rc)));
                  dB = dem::vfabs(pta_local - k * pta_local);

                  //              // (2) here the Morris's method is used
                  //              xa = pta_position.x(); ya = pta_position.y();
                  //              xB = ptb_position.x(); yB = ptb_position.y();
                  //              if(ya==0) {da = xa-radius; dB = radius-xB;}
                  //              else if(xa==0) {da = ya-radius; dB =
                  //              radius-yB;}
                  //              else {
                  //          sqrt_xaya = sqrt(xa*xa+ya*ya);
                  //          k = radius/sqrt_xaya;
                  //          da = radius/k - radius;
                  //          dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
                  //              }

                  beta = 1 + dB / da;
                  if (beta > 2 || beta < 0 || isnan(beta)) {
                    beta = 2;
                  }

                  vdem = (*pta)->getVelocity();
                  vba = beta * ((*ptb)->getVelocity() - vdem);
                  vab = -vba;

                  (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                        (vab * delta_aWab));
                  (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                        (vba * delta_bWba));

                  vr_dot = vab * (pta_position - ptb_position);
                  if (vr_dot < 0) {
                    mu_ab = smoothLength * vr_dot /
                            (rab * rab + 0.01 * smoothLength * smoothLength);
                    Gamma_ab =
(-alpha_zero*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rhoa+rhob)*2;
                  } else {
                    Gamma_ab = 0;
                  }

                  dva_dt = -(*ptb)->getParticleMass() *
                           (pa / (rhoa * rhoa) * coefficient_a +
                            pb / (rhob * rhob) * coefficient_b + Gamma_ab) *
                           delta_aWab;
                  demt->addForce((*pta)->getParticleMass() * dva_dt);
                  demt->addMoment((pta_position - demt_curr) %
                                  ((*pta)->getParticleMass() * dva_dt));
                  //                (*pta)->addVelocityDot(dva_dt);

                  dvb_dt = -(*pta)->getParticleMass() *
                           (pb / (rhob * rhob) * coefficient_b +
                            pa / (rhoa * rhoa) * coefficient_a + Gamma_ab) *
                           delta_bWba;
                  (*ptb)->addVelocityDot(dvb_dt); // the velocities of the ghost
                  particles will not envolve as the same way as others

                    //                  delta_a =
                    // epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
                    //                  (*pta)->addVelocityCorrection(delta_a);
                    delta_b = epsilon * (*pta)->getParticleMass() * (vab)*Wba /
                              (rhoa + rhob) * 2;
                  (*ptb)->addVelocityCorrection(delta_b);
                  break;

                case 3: // boundary with free

                  // calculate Vab as the method shown in Morris's paper, 1996
                  // interact with boundary particles
                  da = ptb_position.z() - allContainer.getMinCorner().z(); //
                  // assume with the bottom boundary
                  dB = allContainer.getMinCorner().z() -
                       pta_position.x(); // assume with
                                         // the bottom boundary
                  if (pta_position.x() <
                      allContainer.getMinCorner().x()) { // with left
                                                         // boundary
                    da = ptb_position.x() - allContainer.getMinCorner().x();
                    dB = allContainer.getMinCorner().x() - pta_position.x();
                  } else if (pta_position.x() >
                             allContainer.getMaxCorner().x()) { // with right
                                                                // boundary
                    da = allContainer.getMaxCorner().x() - ptb_position.x();
                    dB = pta_position.x() - allContainer.getMaxCorner().x();
                  } else if (pta_position.z() >
                             allContainer.getMaxCorner().z()) { // with top
                                                                // boundary
                    da = allContainer.getMaxCorner().z() - ptb_position.z();
                    dB = pta_position.z() - allContainer.getMaxCorner().z();
                  }

                  beta = 1 + dB / da;
                  if (beta > 2 || beta < 0 || isnan(beta)) {
                    beta = 2;
                  }

                  vba = beta * (*ptb)->getVelocity();
                  vab = -vba;

                  (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                        (vab * delta_aWab));
                  (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                        (vba * delta_bWba));

                  vr_dot = vab * (pta_position - ptb_position);
                  if (vr_dot < 0) {
                    mu_ab = smoothLength * vr_dot /
                            (rab * rab + 0.01 * smoothLength * smoothLength);
                Gamma_ab =
(-alpha_zero*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rhoa+rhob)*2;
                  } else {
                    Gamma_ab = 0;
                  }

                  //                dva_dt =
                  //-(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
                  //                (*pta)->addVelocityDot(dva_dt);
                  dvb_dt =
                    -(*pta)->getParticleMass() *
                    (pb / (rhob * rhob) + pa / (rhoa * rhoa) + Gamma_ab) *
                    delta_bWba;
                  (*ptb)->addVelocityDot(dvb_dt); // the velocities of the ghost
                  particles will not envolve as the same way as others

                    //            delta_a =
                    // epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
                    //                (*pta)->addVelocityCorrection(delta_a);
                    delta_b = epsilon * (*pta)->getParticleMass() * (vab)*Wba /
                              (rhoa + rhob) * 2;
                  (*ptb)->addVelocityCorrection(delta_b);

                  //            // apply the boundary forces by Lennard-Jones
                  //            potential as in
                  // Monaghan's paper(1994)
                  //                if(rab<=spaceInterval){ // ptb is in the
                  //                smooth kernel
                  //                dvb_dt = D*(pow(spaceInterval/rab,
                  //                p1)-pow(spaceInterval/rab,
                  // p2))*(ptb_position-pta_position)/(rab*rab);
                  //                (*ptb)->addVelocityDot(dvb_dt);
                  //                } // end if

                  break;
                default:
                  std::cout << "SPH particle type of pta should be 1, 2 or 3!"
                            << std::endl;
                  exit(-1);
              } // end switch pta

              break;
            case 2: // ptb is ghost particle

              if ((*pta)->getType() !=
                  1) { // pta is not free sph particles, i.e. pta
                       // is ghost or boundary particles
                break; // do not consider pta, as we treat before
              }
              demt = (*ptb)->getDemParticle();
              demt_curr = demt->currentPosition();
              pta_local =
                demt->globalToLocal(pta_position - demt_curr); // the local
              // position of sph point pta
              ra = demt->getA();
              rb = demt->getB();
              rc = demt->getC();
              k = 1.0 / (sqrt(pta_local.x() * pta_local.x() / (ra * ra) +
                              pta_local.y() * pta_local.y() / (rb * rb) +
                              pta_local.z() * pta_local.z() / (rc * rc)));
              da = dem::vfabs(pta_local -
                              k * pta_local); // the distance is the same
                                              // in rotated coordinates

              // calculate Vab as the method shown in Morris's paper, 1996

              // (1) here I wanna use the distance from the point a/b to the
              // surface of the ellipsoid to simplify the problem
              ptb_local = (*ptb)->getLocalPosition();
              k = 1.0 / (sqrt(ptb_local.x() * ptb_local.x() / (ra * ra) +
                              ptb_local.y() * ptb_local.y() / (rb * rb) +
                              ptb_local.z() * ptb_local.z() / (rc * rc)));
              dB = dem::vfabs(ptb_local - k * ptb_local);

              //              // (2) here the Morris's method is used
              //              xa = pta_position.x(); ya = pta_position.gety();
              //              xB = ptb_position.x(); yB = ptb_position.gety();
              //              if(ya==0) {da = xa-radius; dB = radius-xB;}
              //              else if(xa==0) {da = ya-radius; dB = radius-yB;}
              //              else {
              //          sqrt_xaya = sqrt(xa*xa+ya*ya);
              //          k = radius/sqrt_xaya;
              //          da = radius/k - radius;
              //          dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
              //              }

              beta = 1 + dB / da;
              if (beta > 2 || beta < 0 || isnan(beta)) {
                beta = 2;
              }

              vdem = (*ptb)->getVelocity();
              vab = beta * ((*pta)->getVelocity() - vdem);
              vba = -vab;

              (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                    (vab * delta_aWab));
              (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                    (vba * delta_bWba));

              vr_dot = vab * (pta_position - ptb_position);
              if (vr_dot < 0) {
                mu_ab = smoothLength * vr_dot /
                        (rab * rab + 0.01 * smoothLength * smoothLength);
                Gamma_ab =
(-alpha_zero*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rhoa+rhob)*2;
              } else {
                Gamma_ab = 0;
              }

              dva_dt = -(*ptb)->getParticleMass() *
                       (pa / (rhoa * rhoa) * coefficient_a +
                        pb / (rhob * rhob) * coefficient_b + Gamma_ab) *
                       delta_aWab;
              (*pta)->addVelocityDot(dva_dt);
              dvb_dt = -(*pta)->getParticleMass() *
                       (pb / (rhob * rhob) * coefficient_b +
                        pa / (rhoa * rhoa) * coefficient_a + Gamma_ab) *
                       delta_bWba;
              demt->addForce((*ptb)->getParticleMass() * dvb_dt);
              demt->addMoment((ptb_position - demt_curr) %
                              ((*ptb)->getParticleMass() * dvb_dt));
              //              (*ptb)->addVelocityDot(dvb_dt);  // the velocities
              //              of the ghost
              // particles will not envolve as the same way as others

              delta_a = epsilon * (*ptb)->getParticleMass() * (-vab) * Wab /
                        (rhoa + rhob) * 2;
              (*pta)->addVelocityCorrection(delta_a);
              //                delta_b =
              // epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
              //                (*ptb)->addVelocityCorrection(delta_b);
              break;
            case 3: // ptb is boundary particle

              if ((*pta)->getType() !=
                  1) { // pta is not free sph particles, i.e. pta
                       // is ghost or boundary particles
                break; // do not consider pta, as we treat before
              }

              // calculate Vab as the method shown in Morris's paper, 1996
              // interact with boundary particles
              da = pta_position.z() - allContainer.getMinCorner().z(); // assume
              // with the bottom boundary
              dB = allContainer.getMinCorner().z() -
                   ptb_position.z(); // assume with
                                     // the bottom boundary
              if (ptb_position.x() <
                  allContainer.getMinCorner().x()) { // with left
                                                     // boundary
                da = pta_position.x() - allContainer.getMinCorner().x();
                dB = allContainer.getMinCorner().x() - ptb_position.x();
              } else if (ptb_position.x() >
                         allContainer.getMaxCorner().x()) { // with
                                                            // right boundary
                da = allContainer.getMaxCorner().x() - pta_position.x();
                dB = ptb_position.x() - allContainer.getMaxCorner().x();
              } else if (ptb_position.z() >
                         allContainer.getMaxCorner().z()) { // with
                                                            // top boundary
                da = allContainer.getMaxCorner().z() - pta_position.z();
                dB = ptb_position.z() - allContainer.getMaxCorner().z();
              }

              beta = 1 + dB / da;
              if (beta > 2 || beta < 0 || isnan(beta)) {
                beta = 2;
              }

              vab = beta * (*pta)->getVelocity();
              vba = -vab;

              (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                    (vab * delta_aWab));
              (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                    (vba * delta_bWba));

              vr_dot = vab * (pta_position - ptb_position);
              if (vr_dot < 0) {
                mu_ab = smoothLength * vr_dot /
                        (rab * rab + 0.01 * smoothLength * smoothLength);
            Gamma_ab =
(-alpha_zero*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rhoa+rhob)*2;
              } else {
                Gamma_ab = 0;
              }

              dva_dt = -(*ptb)->getParticleMass() *
                       (pa / (rhoa * rhoa) * coefficient_a +
                        pb / (rhob * rhob) * coefficient_b + Gamma_ab) *
                       delta_aWab;
              (*pta)->addVelocityDot(dva_dt);
              //              dvb_dt =
              //-(*pta)->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
              //              (*ptb)->addVelocityDot(dvb_dt);  // the velocities
              //              of the ghost
              // particles will not envolve as the same way as others

              delta_a = epsilon * (*ptb)->getParticleMass() * (-vab) * Wab /
                        (rhoa + rhob) * 2;
              (*pta)->addVelocityCorrection(delta_a);
              //              delta_b =
              // epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
              //              (*ptb)->addVelocityCorrection(delta_b);

              //          // apply the boundary forces by Lennard-Jones
              //          potential as in
              // Monaghan's paper(1994)
              //              if(rab<=spaceInterval){ // ptb is in the smooth
              //              kernel
              //            dva_dt = D*(pow(spaceInterval/rab,
              //            p1)-pow(spaceInterval/rab,
              // p2))*(pta_position-ptb_position)/(rab*rab);
              //            (*pta)->addVelocityDot(dva_dt);
              //              } // end if
              break;
            default:
              std::cout << "SPH particle type should be 1, 2 or 3!"
                        << std::endl;
              exit(-1);

          } // end swtich type
        }   // end if 3h
      }     // end for ptb in the same cell

      for (ptb = tmp_particleVec.begin(); ptb != tmp_particleVec.end();
           ptb++) { // all
                    // particles in pvec's neighboring cells
        ptb_position = (*ptb)->getCurrPosition();
        rab = dem::vfabs(pta_position - ptb_position);
        if (rab <= kernelSize) { // ptb is in the smooth kernel
          pb = (*ptb)->getParticlePressure();
          //          mub = (*ptb)->getParticleViscosity();
          rhob = (*ptb)->getParticleDensity();

          Wq = kernelFunction(rab / smoothLength);
          phi_4 = pow(Wq / Wqmin, 4);
          if (pb >= 0) {
            Rb = 0.006;
          } else {
            Rb = 0.6;
          }
          coefficient_a = 1 + Ra * phi_4;
          coefficient_b = 1 + Rb * phi_4;

          // add the density dot for pta and ptb
          delta_aWab = gradientKernelFunction(pta_position, ptb_position); //
          // this is to add SPH pta
          delta_bWba = -delta_aWab;
          dWab_dra = partialKernelFunction(pta_position, ptb_position); // this
          // is to add SPH pta
          dWba_drb = dWab_dra;
          Wab = kernelFunction(pta_position, ptb_position); // this is to add
                                                            // SPH pta
          Wba = Wab;
          switch ((*ptb)->getType()) {
            case 1: // ptb is free SPH particle
              switch ((*pta)->getType()) {
                case 1: // free with free
                  vab = (*pta)->getVelocity() - (*ptb)->getVelocity();
                  vba = -vab;

                  (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                        (vab * delta_aWab));
                  (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                        (vba * delta_bWba));

                  vr_dot = vab * (pta_position - ptb_position);
                  if (vr_dot < 0) {
                    mu_ab = smoothLength * vr_dot /
                            (rab * rab + 0.01 * smoothLength * smoothLength);
                Gamma_ab =
(-alpha*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rhoa+rhob)*2;
                  } else {
                    Gamma_ab = 0;
                  }

                  dva_dt = -(*ptb)->getParticleMass() *
                           (pa / (rhoa * rhoa) * coefficient_a +
                            pb / (rhob * rhob) * coefficient_b + Gamma_ab) *
                           delta_aWab;
                  (*pta)->addVelocityDot(dva_dt);
                  dvb_dt = -(*pta)->getParticleMass() *
                           (pb / (rhob * rhob) * coefficient_b +
                            pa / (rhoa * rhoa) * coefficient_a + Gamma_ab) *
                           delta_bWba;
                  (*ptb)->addVelocityDot(dvb_dt);

                  delta_a = epsilon * (*ptb)->getParticleMass() * (-vab) * Wab /
                            (rhoa + rhob) * 2;
                  (*pta)->addVelocityCorrection(delta_a);
                  delta_b = epsilon * (*pta)->getParticleMass() * (vab)*Wba /
                            (rhoa + rhob) * 2;
                  (*ptb)->addVelocityCorrection(delta_b);
                  break;
                case 2: // ghost with free

                  demt = (*pta)->getDemParticle();
                  demt_curr = demt->currentPosition();
                  ptb_local =
                    demt->globalToLocal(ptb_position - demt_curr); // the
                  // local position of sph point ptb
                  ra = demt->getA();
                  rb = demt->getB();
                  rc = demt->getC();
                  k = 1.0 / (sqrt(ptb_local.x() * ptb_local.x() / (ra * ra) +
                                  ptb_local.y() * ptb_local.y() / (rb * rb) +
                                  ptb_local.z() * ptb_local.z() / (rc * rc)));
                  da = dem::vfabs(ptb_local -
                                  k * ptb_local); // the distance is the same
                                                  // in rotated coordinates

                  // calculate Vab as the method shown in Morris's paper, 1996

                  // (1) here I wanna use the distance from the point a/b to the
                  // surface of the ellipsoid to simplify the problem
                  pta_local = (*pta)->getLocalPosition();
                  k = 1.0 / (sqrt(pta_local.x() * pta_local.x() / (ra * ra) +
                                  pta_local.y() * pta_local.y() / (rb * rb) +
                                  pta_local.z() * pta_local.z() / (rc * rc)));
                  dB = dem::vfabs(pta_local - k * pta_local);

                  //              // (2) here the Morris's method is used
                  //              xa = pta_position.x(); ya = pta_position.y();
                  //              xB = ptb_position.x(); yB = ptb_position.y();
                  //              if(ya==0) {da = xa-radius; dB = radius-xB;}
                  //              else if(xa==0) {da = ya-radius; dB =
                  //              radius-yB;}
                  //              else {
                  //          sqrt_xaya = sqrt(xa*xa+ya*ya);
                  //          k = radius/sqrt_xaya;
                  //          da = radius/k - radius;
                  //          dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
                  //              }

                  beta = 1 + dB / da;
                  if (beta > 2 || beta < 0 || isnan(beta)) {
                    beta = 2;
                  }

                  vdem = (*pta)->getVelocity();
                  vba = beta * ((*ptb)->getVelocity() - vdem);
                  vab = -vba;

                  (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                        (vab * delta_aWab));
                  (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                        (vba * delta_bWba));

                  vr_dot = vab * (pta_position - ptb_position);
                  if (vr_dot < 0) {
                    mu_ab = smoothLength * vr_dot /
                            (rab * rab + 0.01 * smoothLength * smoothLength);
                    Gamma_ab =
(-alpha_zero*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rhoa+rhob)*2;
                  } else {
                    Gamma_ab = 0;
                  }

                  dva_dt = -(*ptb)->getParticleMass() *
                           (pa / (rhoa * rhoa) * coefficient_a +
                            pb / (rhob * rhob) * coefficient_b + Gamma_ab) *
                           delta_aWab;
                  demt->addForce((*pta)->getParticleMass() * dva_dt);
                  demt->addMoment((pta_position - demt_curr) %
                                  ((*pta)->getParticleMass() * dva_dt));
                  //                (*pta)->addVelocityDot(dva_dt);

                  dvb_dt = -(*pta)->getParticleMass() *
                           (pb / (rhob * rhob) * coefficient_b +
                            pa / (rhoa * rhoa) * coefficient_a + Gamma_ab) *
                           delta_bWba;
                  (*ptb)->addVelocityDot(dvb_dt); // the velocities of the ghost
                  // particles will not envolve as the same way as others

                  //                  delta_a =
                  // epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
                  //                  (*pta)->addVelocityCorrection(delta_a);
                  delta_b = epsilon * (*pta)->getParticleMass() * (vab)*Wba /
                            (rhoa + rhob) * 2;
                  (*ptb)->addVelocityCorrection(delta_b);
                  break;
                case 3: // boundary with free

                  // calculate Vab as the method shown in Morris's paper, 1996
                  // interact with boundary particles
                  da = ptb_position.z() - allContainer.getMinCorner().z(); //
                  // assume with the bottom boundary
                  dB = allContainer.getMinCorner().z() -
                       pta_position.x(); // assume with
                                         // the bottom boundary
                  if (pta_position.x() <
                      allContainer.getMinCorner().x()) { // with left
                                                         // boundary
                    da = ptb_position.x() - allContainer.getMinCorner().x();
                    dB = allContainer.getMinCorner().x() - pta_position.x();
                  } else if (pta_position.x() >
                             allContainer.getMaxCorner().x()) { // with right
                                                                // boundary
                    da = allContainer.getMaxCorner().x() - ptb_position.x();
                    dB = pta_position.x() - allContainer.getMaxCorner().x();
                  } else if (pta_position.z() >
                             allContainer.getMaxCorner().z()) { // with top
                                                                // boundary
                    da = allContainer.getMaxCorner().z() - ptb_position.z();
                    dB = pta_position.z() - allContainer.getMaxCorner().z();
                  }

                  beta = 1 + dB / da;
                  if (beta > 2 || beta < 0 || isnan(beta)) {
                    beta = 2;
                  }

                  vba = beta * (*ptb)->getVelocity();
                  vab = -vba;

                  (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                        (vab * delta_aWab));
                  (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                        (vba * delta_bWba));

                  vr_dot = vab * (pta_position - ptb_position);
                  if (vr_dot < 0) {
                    mu_ab = smoothLength * vr_dot /
                            (rab * rab + 0.01 * smoothLength * smoothLength);
                Gamma_ab =
(-alpha_zero*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rhoa+rhob)*2;
                  } else {
                    Gamma_ab = 0;
                  }

                  //                dva_dt =
                  //-(*ptb)->getParticleMass()*(pa/(rhoa*rhoa)*coefficient_a+pb/(rhob*rhob)*coefficient_b+Gamma_ab)*delta_aWab;
                  //                (*pta)->addVelocityDot(dva_dt);
                  dvb_dt =
                    -(*pta)->getParticleMass() *
                    (pb / (rhob * rhob) + pa / (rhoa * rhoa) + Gamma_ab) *
                    delta_bWba;
                  (*ptb)->addVelocityDot(dvb_dt); // the velocities of the ghost
                  particles will not envolve as the same way as others

                    //            delta_a =
                    // epsilon*(*ptb)->getParticleMass()*(-vab)*Wab/(rhoa+rhob)*2;
                    //                (*pta)->addVelocityCorrection(delta_a);
                    delta_b = epsilon * (*pta)->getParticleMass() * (vab)*Wba /
                              (rhoa + rhob) * 2;
                  (*ptb)->addVelocityCorrection(delta_b);

                  // apply the boundary forces by Lennard-Jones potential as in
                  // Monaghan's paper(1994)
                  if (rab <= spaceInterval) { // ptb is in the smooth kernel
                    dvb_dt = D * (pow(spaceInterval / rab, p1) -
                                  pow(spaceInterval / rab, p2)) *
                             (ptb_position - pta_position) / (rab * rab);
                    (*ptb)->addVelocityDot(dvb_dt);
                  } // end if

                  break;
                default:
                  std::cout << "SPH particle type of pta should be 1, 2 or 3!"
                            << std::endl;
                  exit(-1);
              } // end switch pta

              break;
            case 2: // ptb is ghost particle

              if ((*pta)->getType() !=
                  1) { // pta is not free sph particles, i.e. pta
                       // is ghost or boundary particles
                break; // do not consider pta, as we treat before
              }

              demt = (*ptb)->getDemParticle();
              demt_curr = demt->currentPosition();
              pta_local =
                demt->globalToLocal(pta_position - demt_curr); // the local
              // position of sph point pta
              ra = demt->getA();
              rb = demt->getB();
              rc = demt->getC();
              k = 1.0 / (sqrt(pta_local.x() * pta_local.x() / (ra * ra) +
                              pta_local.y() * pta_local.y() / (rb * rb) +
                              pta_local.z() * pta_local.z() / (rc * rc)));
              da = dem::vfabs(pta_local -
                              k * pta_local); // the distance is the same
                                              // in rotated coordinates

              // calculate Vab as the method shown in Morris's paper, 1996

              // (1) here I wanna use the distance from the point a/b to the
              // surface of the ellipsoid to simplify the problem
              ptb_local = (*ptb)->getLocalPosition();
              k = 1.0 / (sqrt(ptb_local.x() * ptb_local.x() / (ra * ra) +
                              ptb_local.y() * ptb_local.y() / (rb * rb) +
                              ptb_local.z() * ptb_local.z() / (rc * rc)));
              dB = dem::vfabs(ptb_local - k * ptb_local);

              //              // (2) here the Morris's method is used
              //              xa = pta_position.x(); ya = pta_position.y();
              //              xB = ptb_position.x(); yB = ptb_position.y();
              //              if(ya==0) {da = xa-radius; dB = radius-xB;}
              //              else if(xa==0) {da = ya-radius; dB = radius-yB;}
              //              else {
              //          sqrt_xaya = sqrt(xa*xa+ya*ya);
              //          k = radius/sqrt_xaya;
              //          da = radius/k - radius;
              //          dB = fabs(xa*xB+ya*yB-k*radius*radius)/sqrt_xaya;
              //              }

              beta = 1 + dB / da;
              if (beta > 2 || beta < 0 || isnan(beta)) {
                beta = 2;
              }

              vdem = (*ptb)->getVelocity();
              vab = beta * ((*pta)->getVelocity() - vdem);
              vba = -vab;

              (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                    (vab * delta_aWab));
              (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                    (vba * delta_bWba));

              vr_dot = vab * (pta_position - ptb_position);
              if (vr_dot < 0) {
                mu_ab = smoothLength * vr_dot /
                        (rab * rab + 0.01 * smoothLength * smoothLength);
                Gamma_ab =
(-alpha_zero*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rhoa+rhob)*2;
              } else {
                Gamma_ab = 0;
              }

              dva_dt = -(*ptb)->getParticleMass() *
                       (pa / (rhoa * rhoa) * coefficient_a +
                        pb / (rhob * rhob) * coefficient_b + Gamma_ab) *
                       delta_aWab;
              (*pta)->addVelocityDot(dva_dt);
              dvb_dt = -(*pta)->getParticleMass() *
                       (pb / (rhob * rhob) * coefficient_b +
                        pa / (rhoa * rhoa) * coefficient_a + Gamma_ab) *
                       delta_bWba;
              demt->addForce((*ptb)->getParticleMass() * dvb_dt);
              demt->addMoment((ptb_position - demt_curr) %
                              ((*ptb)->getParticleMass() * dvb_dt));
              //              (*ptb)->addVelocityDot(dvb_dt);  // the velocities
              //              of the ghost
              // particles will not envolve as the same way as others

              delta_a = epsilon * (*ptb)->getParticleMass() * (-vab) * Wab /
                        (rhoa + rhob) * 2;
              (*pta)->addVelocityCorrection(delta_a);
              //                delta_b =
              // epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
              //                (*ptb)->addVelocityCorrection(delta_b);
              break;
            case 3: // ptb is boundary particle

              if ((*pta)->getType() !=
                  1) { // pta is not free sph particles, i.e. pta
                       // is ghost or boundary particles
                break; // do not consider pta, as we treat before
              }

              // calculate Vab as the method shown in Morris's paper, 1996
              // interact with boundary particles
              da = pta_position.z() - allContainer.getMinCorner().z(); // assume
              // with the bottom boundary
              dB = allContainer.getMinCorner().z() -
                   ptb_position.z(); // assume with
                                     // the bottom boundary
              if (ptb_position.x() <
                  allContainer.getMinCorner().x()) { // with left
                                                     // boundary
                da = pta_position.x() - allContainer.getMinCorner().x();
                dB = allContainer.getMinCorner().x() - ptb_position.x();
              } else if (ptb_position.x() >
                         allContainer.getMaxCorner().x()) { // with
                                                            // right boundary
                da = allContainer.getMaxCorner().x() - pta_position.x();
                dB = ptb_position.x() - allContainer.getMaxCorner().x();
              } else if (ptb_position.z() >
                         allContainer.getMaxCorner().z()) { // with
                                                            // top boundary
                da = allContainer.getMaxCorner().z() - pta_position.z();
                dB = ptb_position.z() - allContainer.getMaxCorner().z();
              }

              beta = 1 + dB / da;
              if (beta > 2 || beta < 0 || isnan(beta)) {
                beta = 2;
              }

              vab = beta * (*pta)->getVelocity();
              vba = -vab;

              (*pta)->addDensityDot((*ptb)->getParticleMass() *
                                    (vab * delta_aWab));
              (*ptb)->addDensityDot((*pta)->getParticleMass() *
                                    (vba * delta_bWba));

              vr_dot = vab * (pta_position - ptb_position);
              if (vr_dot < 0) {
                mu_ab = smoothLength * vr_dot /
                        (rab * rab + 0.01 * smoothLength * smoothLength);
            Gamma_ab =
(-alpha_zero*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rhoa+rhob)*2;
              } else {
                Gamma_ab = 0;
              }

              dva_dt = -(*ptb)->getParticleMass() *
                       (pa / (rhoa * rhoa) * coefficient_a +
                        pb / (rhob * rhob) * coefficient_b + Gamma_ab) *
                       delta_aWab;
              (*pta)->addVelocityDot(dva_dt);
              //              dvb_dt =
              //-(*pta)->getParticleMass()*(pb/(rhob*rhob)+pa/(rhoa*rhoa)+Gamma_ab)*delta_bWba;
              //              (*ptb)->addVelocityDot(dvb_dt);  // the velocities
              //              of the ghost
              // particles will not envolve as the same way as others

              delta_a = epsilon * (*ptb)->getParticleMass() * (-vab) * Wab /
                        (rhoa + rhob) * 2;
              (*pta)->addVelocityCorrection(delta_a);
              //              delta_b =
              // epsilon*(*pta)->getParticleMass()*(vab)*Wba/(rhoa+rhob)*2;
              //              (*ptb)->addVelocityCorrection(delta_b);

              // apply the boundary forces by Lennard-Jones potential as in
              // Monaghan's paper(1994)
              if (rab <= spaceInterval) { // ptb is in the smooth kernel
                dva_dt = D * (pow(spaceInterval / rab, p1) -
                              pow(spaceInterval / rab, p2)) *
                         (pta_position - ptb_position) / (rab * rab);
                (*pta)->addVelocityDot(dva_dt);
              } // end if

              break;
            default:
              std::cout << "SPH particle type should be 1, 2 or 3!"
                        << std::endl;
              exit(-1);

          } // end swtich type
        }   // end if 3h
      }     // end for ptb in neighbor cells
    }       // end for pta

    tmp_particleVec.clear(); // clear elements in tmp-vector for particles
                             // neighboring cells, it is important

  } // end for pvec, different cells

  //      // apply the boundary forces by Lennard-Jones potential as in
  //      Monaghan's
  // paper(1994)
  //      for(ptb=SPHBoundaryParticleVec.begin();
  // ptb!=SPHBoundaryParticleVec.end(); ptb++){
  //    ptb_position = (*ptb)->getCurrPosition();
  //    rab = dem::vfabs(pta_position-ptb_position);
  //    if(rab<=spaceInterval){ // ptb is in the smooth kernel
  //        dva_dt = D*(pow(spaceInterval/rab, p1)-pow(spaceInterval/rab,
  // p2))*(pta_position-ptb_position)/(rab*rab);
  //        (*pta)->addVelocityDot(dva_dt);
  //    } // end if
  //      } // end ptb

} // end calculateSPHDensityDotVelocityDotLinkedList3D()
*/
