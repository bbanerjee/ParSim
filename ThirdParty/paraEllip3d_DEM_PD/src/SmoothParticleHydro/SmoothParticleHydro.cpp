#include <DiscreteElements/DEMParticle.h>
#include <SmoothParticleHydro/SmoothParticleHydro.h>
#include <SmoothParticleHydro/SPHInteractions.h>

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
  d_mpiProcs = { { mpiProcX, mpiProcY, mpiProcZ } };

  // create Cartesian virtual topology (unavailable in boost.mpi)
  int ndim = 3;
  int periods[3] = { 0, 0, 0 };
  int reorder = 0; // d_mpiRank not reordered
  int status = MPI_Cart_create(d_mpiWorld, ndim, d_mpiProcs.data(), periods,
                               reorder, &d_cartComm);
  if (status != MPI_SUCCESS) {
    std::cout << "**ERROR** Could not create MPI Cartesian topology.\n";
    exit(-1);
  }

  status = MPI_Comm_rank(d_cartComm, &d_mpiRank);
  if (status != MPI_SUCCESS) {
    std::cout << "**ERROR** Could not get MPI Cartesian topology rank.\n";
    exit(-1);
  }

  status = MPI_Comm_size(d_cartComm, &d_mpiSize);
  if (status != MPI_SUCCESS) {
    std::cout << "**ERROR** Could not get MPI Cartesian topology size.\n";
    exit(-1);
  }

  status = MPI_Cart_coords(d_cartComm, d_mpiRank, ndim, d_mpiCoords.data());
  if (status != MPI_SUCCESS) {
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

  // Set d_sphPatchBox on all procs
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

    Vec v1 = d_sphPatchBox.getMinCorner();
    Vec v2 = d_sphPatchBox.getMaxCorner();
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

  // std::cout << "d_mpirank = " << d_mpiRank << " Num particles "
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
  Vec v1 = d_sphPatchBox.getMinCorner();
  Vec v2 = d_sphPatchBox.getMaxCorner();
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
  Vec v1 = d_sphPatchBox.getMinCorner();
  Vec v2 = d_sphPatchBox.getMaxCorner();
  Vec vspan = (v2 - v1) / d_mpiProcs;
  Vec lower = v1 + vspan * d_mpiCoords;
  Vec upper = lower + vspan;
  Box container(lower, upper);
  d_patchP->update(iteration, lower, upper, ghostWidth);
}

void
SmoothParticleHydro::commuSPHParticle(int iteration, const REAL& ghostWidth)
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

void
SmoothParticleHydro::migrateSPHParticle(int iteration)
{
  Vec domainWidth = d_sphPatchBox.getMaxCorner() - d_sphPatchBox.getMinCorner();
  Vec patchWidth = domainWidth / d_mpiProcs;

  // Migrate particles in the x-direction
  dem::ParticleIDHashMap sentParticles;
  SPHParticlePArray recvParticles;
  d_patchP->sendRecvMigrateXMinus(d_boostWorld, iteration, patchWidth,
                                  d_sphParticleVec);
  d_patchP->sendRecvMigrateXPlus(d_boostWorld, iteration, patchWidth,
                                 d_sphParticleVec);
  d_patchP->waitToFinishX(iteration);
  d_patchP->combineSentParticlesX(iteration, sentParticles);
  d_patchP->combineReceivedParticlesX(iteration, recvParticles);
  d_patchP->deleteSentParticles<SPHParticleP>(iteration, sentParticles,
                                              d_sphParticleVec);
  d_patchP->addReceivedParticles(iteration, recvParticles, d_sphParticleVec);

  /*
  std::ostringstream out;
  out << "Iter = " << iteration << " Rank = " << d_mpiRank
      << " Recv X = " << recvParticles.size()
      << " Sent X = " << sentParticles.size()
      << " All = " << d_sphParticleVec.size() << "\n";
  */

  // Migrate particles in the y-direction
  sentParticles.clear();
  recvParticles.clear();
  d_patchP->sendRecvMigrateYMinus(d_boostWorld, iteration, patchWidth,
                                  d_sphParticleVec);
  d_patchP->sendRecvMigrateYPlus(d_boostWorld, iteration, patchWidth,
                                 d_sphParticleVec);
  d_patchP->waitToFinishY(iteration);
  d_patchP->combineSentParticlesY(iteration, sentParticles);
  d_patchP->combineReceivedParticlesY(iteration, recvParticles);
  d_patchP->deleteSentParticles<SPHParticleP>(iteration, sentParticles,
                                              d_sphParticleVec);
  d_patchP->addReceivedParticles(iteration, recvParticles, d_sphParticleVec);

  /*
  out << "Iter = " << iteration << " Rank = " << d_mpiRank
      << " Recv Y = " << recvParticles.size()
      << " Sent Y = " << sentParticles.size()
      << " All = " << sphParticleVec.size() << "\n";
  */

  // Migrate particles in the z-direction
  sentParticles.clear();
  recvParticles.clear();
  d_patchP->sendRecvMigrateZMinus(d_boostWorld, iteration, patchWidth,
                                  d_sphParticleVec);
  d_patchP->sendRecvMigrateZPlus(d_boostWorld, iteration, patchWidth,
                                 d_sphParticleVec);
  d_patchP->waitToFinishZ(iteration);
  d_patchP->combineSentParticlesZ(iteration, sentParticles);
  d_patchP->combineReceivedParticlesZ(iteration, recvParticles);
  d_patchP->deleteSentParticles<SPHParticleP>(iteration, sentParticles,
                                              d_sphParticleVec);
  d_patchP->addReceivedParticles(iteration, recvParticles, d_sphParticleVec);

  /*
  out << "Iter = " << iteration << " Rank = " << d_mpiRank
      << " Recv Z = " << recvParticles.size()
      << " Sent Z = " << sentParticles.size()
      << " All = " << sphParticleVec.size() << "\n";
  std::cout << out.str();
  */

  // delete outgoing particles
  // d_patchP->removeParticlesOutsidePatch<SPHParticleP>(d_sphParticleVec);

  // out << "After del: " << sphParticleVec.size() << "\n";
  /*
  std::cout << out.str();
  */
}

// update d_allSPHParticleVec: process 0 collects all updated particles from
// each other process
void
SmoothParticleHydro::gatherSPHParticle()
{
  if (d_mpiRank != 0) {

    d_boostWorld.send(0, d_mpiTag, d_sphParticleVec);

  } else { 

    d_allSPHParticleVec.clear();

    // d_allSPHParticleVec is cleared before filling with new data
    //releaseGatheredSPHParticle();

    // duplicate d_sphParticleVec so that it is not destroyed by
    // d_allSPHParticleVec in next iteration,
    // otherwise it causes memory error.
    SPHParticlePArray dupSPHParticleVec(d_sphParticleVec.size());
    std::size_t index = 0;
    for (const auto& particle :  d_sphParticleVec) {
      dupSPHParticleVec[index] = std::make_shared<SPHParticle>(*particle);
      index++;
    }

    // fill allParticleVec with dupParticleVec
    d_allSPHParticleVec.insert(d_allSPHParticleVec.end(),
                               dupSPHParticleVec.begin(),
                               dupSPHParticleVec.end());

    // Add received particles
    //long gatherRam = 0;
    for (int iRank = 1; iRank < d_mpiSize; ++iRank) {

      SPHParticlePArray recvParticles;
      d_boostWorld.recv(iRank, d_mpiTag, recvParticles);
      d_allSPHParticleVec.insert(d_allSPHParticleVec.end(),
                                 recvParticles.begin(),
                                 recvParticles.end());
      //gatherRam += recvParticles.size();
    }
    // debugInf << "gather: particleNum = " << gatherRam <<  " particleRam = "
    //<< gatherRam * sizeof(SPHParticle) << std::endl;
  }
}

// initialize the densityRate and velocityDot of all the SPH particles
void 
SmoothParticleHydro::initializeDensityRateAndAcceleration()
{
  REAL gravAccel = util::getParam<REAL>("gravAccel");
  REAL gravScale = util::getParam<REAL>("gravScale");

  // Assume gravity is in the z-direction
  dem::Vec bodyForceAcc(0,0, -gravAccel*gravScale);
  dem::Vec zero(0,0,0);

  for (auto& particle : d_mergeSPHParticleVec) {
    particle->setDensityRateAccelerationZero();
    particle->calculateViscosity();
    particle->calculatePressure();
    switch (particle->getType()) {
      case SPHParticleType::FREE: 
        particle->incAcceleration(bodyForceAcc);
        break;
      case SPHParticleType::GHOST: 
        particle->incAcceleration(zero);
        break;
      case SPHParticleType::BOUNDARY: 
        particle->incAcceleration(zero);
        break;
      default:
        std::cout << "SPH particle type should be FREE, GHOST or BOUNDARY!"
                  << std::endl;
        exit(-1);
    } // switch
  }
}

// divide SPH domain in each cpu into different cells in 2D in xz plane.
// see the notes 5/20/2015 and 5/21/2015 or Simpson's paper "Numerical
// techniques for three-dimensional Smoothed Particle Hydrodynamics"
//
// **NOTE** bufferWidth = spaceInterval;
template <int dim>
void
SmoothParticleHydro::assignParticlesToPatchGrid(const Box& container,
                                                const REAL& bufferWidth,
                                                const REAL& ghostWidth,
                                                const REAL& kernelSize)
{
  REAL small_value = 0.01 * bufferWidth;
  REAL ghostBuffer = ghostWidth + small_value;

  // expand the container by bufferWidth + small_value
  Box expandedContainer(container, ghostBuffer);

  // Compute the number of cells in each dimension
  int nx = std::round(expandedContainer.getDimx()/kernelSize) + 1;
  int ny = std::round(expandedContainer.getDimy()/kernelSize) + 1;
  int nz = std::round(expandedContainer.getDimz()/kernelSize) + 1;
  d_numGridCells.setX(nx);
  d_numGridCells.setY(ny); 
  d_numGridCells.setZ(nz);

  d_sphPatchGridIndex.resize(nx*ny*nz); 
  int cellID = 0;
  for (int kk = 0; kk < nz; kk++) {
    for (int jj = 0; jj < ny; jj++) {
      for (int ii = 0; ii < nx; ii++) {
        d_sphPatchGridIndex[cellID] = IntVec(ii, jj, kk);
        cellID++;
      }
    }
  }

  d_sphPatchGrid.resize(nx*ny*nz); 
  for (const auto& particle : d_mergeSPHParticleVec) {
    int cellIndex = getCellIndex<dim>(expandedContainer.getMinCorner(),
                                      kernelSize, 
                                      d_numGridCells,
                                      particle->currentPosition());
    d_sphPatchGrid[cellIndex].push_back(particle);
  }

} // assignParticlesToPatchGrid2D

template <>
int 
SmoothParticleHydro::getCellIndex<1>(const dem::Vec& cellMinCorner,
                                     const REAL& cellWidth,
                                     const IntVec& numGridCells,
                                     const Vec& pos)
{
  Vec relPos = (pos - cellMinCorner)/cellWidth;
  int xindex = std::floor(relPos.x());
  int cellIndex = xindex;
  return cellIndex;  
}

template <>
int 
SmoothParticleHydro::getCellIndex<2>(const dem::Vec& cellMinCorner,
                                     const REAL& cellWidth,
                                     const IntVec& numGridCells,
                                     const Vec& pos)
{
  Vec relPos = (pos - cellMinCorner)/cellWidth;
  int xindex = std::floor(relPos.x());
  int zindex = std::floor(relPos.z());
  int cellIndex = zindex * numGridCells.x() + xindex;
  return cellIndex;  
}

template <>
int 
SmoothParticleHydro::getCellIndex<3>(const dem::Vec& cellMinCorner,
                                     const REAL& cellWidth,
                                     const IntVec& numGridCells,
                                     const Vec& pos)
{
  Vec relPos = (pos - cellMinCorner)/cellWidth;
  int xindex = std::floor(relPos.x());
  int yindex = std::floor(relPos.y());
  int zindex = std::floor(relPos.z());
  int cellIndex = zindex * numGridCells.y() * numGridCells.x() +
                  yindex * numGridCells.x() + 
                  xindex;
  return cellIndex;  
}

template <>
int
SmoothParticleHydro::getIndex<1>(const IntVec& cell, const int& nx, 
                                 const int& ny, const int& nz) const
{
  return cell.x();
}

template <>
int
SmoothParticleHydro::getIndex<2>(const IntVec& cell, const int& nx,
                                 const int& ny, const int& nz) const
{
  return cell.y()*nx + cell.x();
}

template <>
int
SmoothParticleHydro::getIndex<3>(const IntVec& cell, const int& nx,
                                 const int& ny, const int& nz) const
{
  return cell.z()*ny*nx + cell.y()*nx + cell.x();
}

template <>
std::vector<int>
SmoothParticleHydro::getAdjacentCellIndices<1>(const IntVec& cell,
                                               const IntVec& numGridCells)
{
  // The adjacent cell indices
  std::vector<int> neighbors = {-1, 1};
  int nx = numGridCells.x();

  std::vector<int> neighborIndices;
  for (auto& xneighbor : neighbors) {
    int xindex = cell.x() + xneighbor; 
    if ((xindex < 0) || (xindex > nx - 1)) continue;
    IntVec cellNeighbor(xindex, 0, 0);
    neighborIndices.push_back(getIndex<1>(cellNeighbor, nx, 0, 0));
  }

  return neighborIndices;
}

template <>
std::vector<int>
SmoothParticleHydro::getAdjacentCellIndices<2>(const IntVec& cell,
                                               const IntVec& numGridCells)
{
  std::vector<int> neighbors = {-1, 1};
  int nx = numGridCells.x();
  int nz = numGridCells.z();

  // The adjacent cell indices
  std::vector<int> neighborIndices;
  for (auto& zneighbor : neighbors) {
    int zindex = cell.z() + zneighbor; 
    if ((zindex < 0) || (zindex > nz - 1)) continue;
    for (auto& xneighbor : neighbors) {
      int xindex = cell.x() + xneighbor; 
      if ((xindex < 0) || (xindex > nx - 1)) continue;
      IntVec cellNeighbor(xindex, 0, zindex);
      neighborIndices.push_back(getIndex<2>(cellNeighbor, nx, 0, nz));
    }
  }

  return neighborIndices;
}

template <>
std::vector<int>
SmoothParticleHydro::getAdjacentCellIndices<3>(const IntVec& cell,
                                               const IntVec& numGridCells)
{
  std::vector<int> neighbors = {-1, 1};
  int nx = numGridCells.x();
  int ny = numGridCells.y();
  int nz = numGridCells.z();

  // The adjacent cell indices
  std::vector<int> neighborIndices;
  for (auto& zneighbor : neighbors) {
    int zindex = cell.z() + zneighbor; 
    if ((zindex < 0) || (zindex > nz - 1)) continue;
    for (auto& yneighbor : neighbors) {
      int yindex = cell.y() + yneighbor; 
      if ((yindex < 0) || (yindex > ny - 1)) continue;
      for (auto& xneighbor : neighbors) {
        int xindex = cell.x() + xneighbor; 
        if ((xindex < 0) || (xindex > nx - 1)) continue;
        IntVec cellNeighbor(xindex, yindex, zindex);
        neighborIndices.push_back(getIndex<3>(cellNeighbor, nx, ny, nz));
      }
    }
  }

  return neighborIndices;

}

template <int dim>
SPHParticlePArray 
SmoothParticleHydro::getParticlesInAdjacentCells(const int& cellIndex,
                                                 const IntVec& numGridCells) const
{
  // Current cell
  IntVec cell = d_sphPatchGridIndex[cellIndex];

  // Find adjacent cells
  auto neighborIndices = getAdjacentCellIndices<dim>(cell, numGridCells);

  // Create an array of particles from adjacent cells
  SPHParticlePArray neighborParticles;
  for (const auto& neighborIndex : neighborIndices) {
    for (const auto& particle : d_sphPatchGrid[neighborIndex]) {
      neighborParticles.push_back(particle);
    }
  }

  return neighborParticles;
}

template <int dim>
void
SmoothParticleHydro::computeParticleInteractions(SPHParticleP& sph_part_a,
                                                 SPHParticleP& sph_part_b,
                                                 const REAL& kernelSize,
                                                 const REAL& smoothLength,
                                                 const Box& allContainer) const
{
  REAL alpha = util::getParam<REAL>("alpha");
  REAL epsilon = util::getParam<REAL>("epsilon");  // parameter for velocity
  REAL soundSpeed = util::getParam<REAL>("soundSpeed");

  SPHInteractions interaction;
  interaction.setDataParticleA(sph_part_a);
  interaction.setDataParticleB(sph_part_b);
  interaction.initializeInteractionData();

  interaction.computeInteractionKernel<dim>(smoothLength);
  interaction.updateParticleCoeffs();

  // we have three types of SPH particles: 1, free particle; 2, ghost
  // particle; 3, boundary particle.  Then we have 3x3 = 9 different 
  // types of interactions with the three types of particles
  // sph_part_a   sph_part_b     need to consider or not
  // 1            1              V      free with free
  // 1            2              V      free with ghost
  // 1            3              V      free with boundary

  // 2            1              V      ghost with free
  // 2            2              X      ghost with ghost
  // 2            3              X      ghost with boundary

  // 3            1              V      boundary with free
  // 3            2              X      boundary with ghost
  // 3            3              X       boundary with boundary
  if (sph_part_b->getType() == SPHParticleType::FREE) {

    switch(sph_part_a->getType()) {

    case SPHParticleType::FREE:

      interaction.updateMomentumExchangeCoeffs(smoothLength, alpha, soundSpeed);
      interaction.updateInteractionFreeFree(sph_part_a, sph_part_b, epsilon);
      break;        

    case SPHParticleType::GHOST:

      interaction.doGhostBoundaryVelocityCorrection(sph_part_a);
      interaction.updateMomentumExchangeCoeffs(smoothLength, alpha, soundSpeed);
      interaction.updateInteractionGhostFree(sph_part_a, sph_part_b, epsilon);
      break;        

    case SPHParticleType::BOUNDARY:

      // calculate Vab as the method shown in Morris's paper, 1996
      // interact with boundary particles
      interaction.doDomainBoundaryVelocityCorrection(allContainer);
      interaction.updateMomentumExchangeCoeffs(smoothLength, alpha, soundSpeed);
      interaction.updateInteractionBoundaryFree(sph_part_a, sph_part_b, epsilon);
      break;        

    default:
      break;
    }

  } else if (sph_part_b->getType() == SPHParticleType::GHOST &&
              sph_part_a->getType() == SPHParticleType::FREE) {

    interaction.swapParticles();
    interaction.doGhostBoundaryVelocityCorrection(sph_part_b);
    interaction.updateMomentumExchangeCoeffs(smoothLength, alpha, soundSpeed);
    interaction.updateInteractionGhostFree(sph_part_b, sph_part_a, epsilon);

  } else if (sph_part_b->getType() == SPHParticleType::BOUNDARY &&
              sph_part_a->getType() == SPHParticleType::FREE) {

    interaction.swapParticles();
    interaction.doDomainBoundaryVelocityCorrection(allContainer);
    interaction.updateMomentumExchangeCoeffs(smoothLength, alpha, soundSpeed);
    interaction.updateInteractionBoundaryFree(sph_part_b, sph_part_a, epsilon);

  }
}

//// the momentum equilibrium equation and state equation are implemented as
// in Monaghan's paper (1994), simulate free surface flow using sph
//// here the neighboring list of SPH particles is searched by the cells,
template <int dim>
void
SmoothParticleHydro::updateParticleInteractions(const REAL& kernelSize,
                                                const REAL& smoothLength,
                                                const Box& allContainer)
{
  // Initialize the rate quantities that are integrated
  initializeDensityRateAndAcceleration();

  // divide the SPH domain into different cells, each cell may contain SPH
  // particles within it
  assignParticlesToPatchGrid<dim>();

  for (int cellIndex=0; cellIndex < d_sphPatchGrid.size(); ++cellIndex) {

    SPHParticlePArray cellParticles = d_sphPatchGrid[cellIndex];
    SPHParticlePArray neighborParticles = 
      getParticlesInAdjacentCells<dim>(cellIndex, d_numGridCells);

    for (auto iter_a = cellParticles.begin(); iter_a != cellParticles.end(); iter_a++) {
      auto sph_part_a = *iter_a;

      for (auto iter_b = iter_a+1; iter_b != cellParticles.end(); iter_b++) {
        auto sph_part_b = *iter_b;

        // Go to next particle if the distance is too large
        if (sph_part_b->isOutsideInfluenceZone(*sph_part_a, kernelSize)) continue;

        computeParticleInteractions<dim>(sph_part_a, sph_part_b,
                                         kernelSize, smoothLength, allContainer);

      } // end for sph_part_b in the same cell

      for (auto& sph_part_b : neighborParticles) {

        // Go to next particle if the distance is too large
        if (sph_part_b->isOutsideInfluenceZone(*sph_part_a, kernelSize)) continue;

        computeParticleInteractions<dim>(sph_part_a, sph_part_b,
                                         kernelSize, smoothLength, allContainer);

      } // end for sph_part_b in neighbor cells
    } // end for sph_part_a
  } // end for cellIndex, different cells
} // end updateParticleInteractions()

/*
void
SmoothParticleHydro::initialSPHVelocity2D()
{
  commuParticle();
  commuSPHParticle(); // this will update container and mergeSPHParticleVec,
                      // both are needed for assignParticlesToPatchGrid
  calculateSPHDensityRateAccelerationLinkedList2D(); // calculate velocityDot and
                                                   // densityRate
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
        std::cout << "SPH particle type of sph_part_a should be 1, 2 or 3!"
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
                      // both are needed for assignParticlesToPatchGrid
  calculateSPHDensityRateAccelerationLinkedList3D(); // calculate velocityDot and
                                                   // densityRate
  initialSPHLeapFrogVelocity(); // initial velocity only for free SPH particles
                                // based on equation (4.3)

  releaseRecvParticle();
  releaseRecvSPHParticle();
} // initialSPHVelocity()

void
SmoothParticleHydro::initialSPHVelocityCopyDEM3D()
{
  commuSPHParticle(); // this will update container and mergeSPHParticleVec,
                      // both are needed for assignParticlesToPatchGrid
  calculateSPHDensityRateAccelerationLinkedList3D(); // calculate velocityDot and
                                                   // densityRate
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
        std::cout << "SPH particle type of sph_part_a should be 1, 2 or 3!"
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
*/

/*
// divide SPH domain in each cpu into different cells in 2D in xz plane.
// see the notes 5/20/2015 and 5/21/2015 or Simpson's paper "Numerical
// techniques for three-dimensional Smoothed Particle Hydrodynamics"
void
SmoothParticleHydro::assignParticlesToPatchGrid3D()
{

  // clear the std::vector< std::vector<sph::SPHParticle*> > d_sphPatchGrid
  for (int cellIndex = 0; cellIndex < d_sphPatchGrid.size(); cellIndex++) {
    d_sphPatchGrid[cellIndex].clear();
  }
  d_sphPatchGrid.clear();

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
  d_sphPatchGrid.resize(numCell); // at this point, the cellVec contains
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
    tmp_xyz = (*pt)->currentPosition();
    num = int((tmp_xyz.y() - ymin) / kernelSize) * Nx * Nz +
          int((tmp_xyz.z() - zmin) / kernelSize) * Nx +
          int((tmp_xyz.x() - xmin) / kernelSize);
    d_sphPatchGrid[num].push_back(*pt);
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
      tmp_xyz = (*gt)->currentPosition();
      if (tmp_xyz.x() >= xmax || tmp_xyz.y() >= ymax || tmp_xyz.z() >= zmax ||
          tmp_xyz.x() <= xmin || tmp_xyz.y() <= ymin || tmp_xyz.z() <= zmin)
        continue; // this sph ghost particle is outside of the container, go to
                  // next sph ghost particle
      num = int((tmp_xyz.y() - ymin) / kernelSize) * Nx * Nz +
            int((tmp_xyz.z() - zmin) / kernelSize) * Nx +
            int((tmp_xyz.x() - xmin) / kernelSize);
      d_sphPatchGrid[num].push_back(*gt);
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
      tmp_xyz = (*gt)->currentPosition();
      // not all sph ghost particles in these received particles are inside the
      // container
      if (tmp_xyz.x() >= xmax || tmp_xyz.y() >= ymax || tmp_xyz.z() >= zmax ||
          tmp_xyz.x() <= xmin || tmp_xyz.y() <= ymin || tmp_xyz.z() <= zmin)
        continue; // this sph ghost particle is outside of the container, go to
      next sph ghost particle num =
        int((tmp_xyz.y() - ymin) / kernelSize) * Nx * Nz +
        int((tmp_xyz.z() - zmin) / kernelSize) * Nx +
        int((tmp_xyz.x() - xmin) / kernelSize);
      d_sphPatchGrid[num].push_back(*gt);
    }
  }

} // assignParticlesToPatchGrid3D
*/

/*
//// the momentum equilibrium equation and state equation are implemented as
// in Monaghan's paper (1994), simulate free surface flow using sph
//// here the neighboring list of SPH particles is searched by the cells,
vector<vector<sph::SPHParticle*>> d_sphPatchGrid;
void
SmoothParticleHydro::calculateSPHDensityRateAccelerationLinkedList3D()
{

  // divide the SPH domain into different cells, each cell will contain SPH
  // particles within it
  assignParticlesToPatchGrid3D();

  //  checkDivision();  // pass, May 22, 2015
  //  checkNeighborCells3D();  // pass, May 22, 2015

  // initialize the densityRate and velocityDot of all the SPH particles
  dem::Vec tmp_vec = dem::Vec(0,0,
-(util::getParam<REAL>("gravAccel")*(util::getParam<REAL>("gravScale") );
  dem::Vec zero_vec = dem::Vec(0,0,0);
  // in the calculation of forces between sph particles, we need to use the
//mergeSPHParticleVec,
  // since mergeVec contains the sph particles from neighboring cpus. While we
//can only update and migrate and communicate d_sphParticleVec
  for(std::vector<sph::SPHParticle*>::iterator pt=mergeSPHParticleVec.begin();
pt!=mergeSPHParticleVec.end(); pt++){
    (*pt)->setDensityRateAccelerationZero();
    (*pt)->calculateParticleViscosity(); // everytime when use
    // getParticleViscolisity(), make sure that pressure has been calculated!!!!
    (*pt)->calculateParticlePressure(); // everytime when use
    // getPressure(), make sure that pressure has been calculated!!!!
    switch ((*pt)->getType()) {
      case 1: // free sph particle
        (*pt)->incAcceleration(tmp_vec);
        break;
      case 2: // ghost sph particle
        break;
      case 3: // boundary sph particle
        (*pt)->incAcceleration(zero_vec);
        break;
      default:
        std::cout << "SPH particle type of sph_part_a should be 1, 2 or 3!"
                  << std::endl;
        exit(-1);
    } // switch
  }

  std::vector<sph::SPHParticle*>::iterator gt;
  for(std::vector<Particle*>::iterator pdt=mergeParticleVec.begin();
pdt!=mergeParticleVec.end(); pdt++){  // all sph ghost particles
    for (gt = (*pdt)->SPHGhostParticleVec.begin();
         gt != (*pdt)->SPHGhostParticleVec.end(); gt++) {
      (*gt)->setDensityRateAccelerationZero();
      (*gt)->calculateParticleViscosity(); // everytime when use
      // getParticleViscolisity(), make sure that pressure has been
      // calculated!!!!
      (*gt)->calculateParticlePressure(); // everytime when use
      // getPressure(), make sure that pressure has been calculated!!!!
    }
  }

  // temporary variables used in the loop
  dem::Vec pos_a;
  dem::Vec pos_b;
  dem::Vec gradWab_a;
  dem::Vec gradWba_b;
  dem::Vec vab;
  dem::Vec vel_ba;
  dem::Vec vdem;
  dem::Vec dva_dt;
  dem::Vec dvb_dt;
  dem::Vec delta_a;
  dem::Vec delta_b;
  REAL press_a, press_b, rho_a, rho_b, mu_a, mu_b;
  REAL rab;
        REAL dWab_dra;
  REAL dWba_drb;
  REAL Wab, Wba;
  REAL da, dB;
  REAL beta;
  REAL xa, ya, xB, yB, k, sqrt_xaya;  // variables for Morris' method to
//calculate da/dB
  std::vector<sph::SPHParticle*>::iterator sph_part_b;
  REAL ra, rb, rc;  // the geometry of the dem particle
  dem::Vec local_a, local_b;  // local position of sph_part_a and sph_part_b in the dem
//particle
  dem::Vec dem_pos;
  REAL Gamma_ab, mu_ab, vel_radial;
  REAL alpha = util::getParam<REAL>("alpha");
     REAL alpha_zero = 0;  // the viscous between free and ghost/boundary
  REAL epsilon = util::getParam<REAL>("epsilon");  // parameter for velocity
//correction
  REAL Wq, Ra, Rb, phi_4, coeff_a, coeff_b;
  dem::Particle* demt;

  int pnum;
  std::vector<sph::SPHParticle*> tmp_particleVec;  // sph particles in
//neighboring cells
  for(int cellIndex=0; cellIndex<d_sphPatchGrid.size(); cellIndex++){

    // store all SPH particles in cellIndex's neighboring cells
    tmp_particleVec.clear(); // it is the same for the particles in the same
                             // cell
    for (std::vector<int>::const_iterator pint = pnum_vec.begin();
         pint != pnum_vec.end(); pint++) {
      pnum = cellIndex + (*pint);
      if (pnum < numCell) {
        for (std::vector<sph::SPHParticle*>::iterator pt =
               d_sphPatchGrid[pnum].begin();
             pt != d_sphPatchGrid[pnum].end(); pt++) {
          tmp_particleVec.push_back(*pt);
        }
      }
    }

    for (std::vector<sph::SPHParticle*>::iterator sph_part_a =
           d_sphPatchGrid[cellIndex].begin();
         sph_part_a != d_sphPatchGrid[cellIndex].end();
         sph_part_a++) { // SPH particles in cell cellIndex
      press_a = sph_part_a->getPressure();
      if (press_a >= 0) {
        Ra = 0.006;
      } else {
        Ra = 0.6;
      }
      //        mu_a = sph_part_a->getViscosity();
      rho_a = sph_part_a->getDensity();
      pos_a = sph_part_a->currentPosition();

      //    if(sph_part_a->getType()!=1){  // sph_part_a is not free sph particles, i.e.
      //    sph_part_a is
      // ghost or boundary particles
      //        continue;  // do not consider sph_part_a, as we treat before
      //    }

      for (sph_part_b = sph_part_a + 1; sph_part_b != d_sphPatchGrid[cellIndex].end();
           sph_part_b++) { // sum over the
                    // SPH particles in the same cell cellIndex
        pos_b = sph_part_b->currentPosition();
        rab = dem::vfabs(pos_a - pos_b);
        if (rab <= kernelSize) { // sph_part_b is in the smooth kernel
          press_b = sph_part_b->getPressure();
          //          mu_b = sph_part_b->getViscosity();
          rho_b = sph_part_b->getDensity();

          Wq = kernelFunction(rab / smoothLength);
          phi_4 = pow(Wq / Wqmin, 4);
          if (press_b >= 0) {
            Rb = 0.006;
          } else {
            Rb = 0.6;
          }
          coeff_a = 1 + Ra * phi_4;
          coeff_b = 1 + Rb * phi_4;

          // we have three types of SPH particles: 1, free particle; 2, ghost
          // particle; 3, boundary particle
          // Then we have 3x3 = 9 different types of interactions with the three
          // types of particles
          // so we cannot judge the interaction type only by the type of
          // particle
          // sph_part_b, we need to consider sph_part_a also
          // sph_part_a          sph_part_b          need to consider or not
          // 1            1              V      free with free
          // 1            2              V      free with ghost
          // 1            3              V      free with boundary

          // 2            1              V      ghost with free
          // 2            2              X      ghost with ghost
          // 2            3              X      ghost with boundary

          // 3            1              V      boundary with free
          // 3            2              X      boundary with ghost
          // 3            3              X       boundary with boundary

          // add the density dot for sph_part_a and sph_part_b
          gradWab_a = gradientKernelFunction(pos_a, pos_b); //
          // this is to add SPH sph_part_a
          gradWba_b = -gradWab_a;
          dWab_dra = partialKernelFunction(pos_a, pos_b); // this
          // is to add SPH sph_part_a
          dWba_drb = dWab_dra;
          Wab = kernelFunction(pos_a, pos_b); // this is to add
                                                            // SPH sph_part_a
          Wba = Wab;

          switch (sph_part_b->getType()) {
            case 1: // sph_part_b is free SPH particle
              switch (sph_part_a->getType()) {
                case 1: // free with free
                  vab = sph_part_a->getVelocity() - sph_part_b->getVelocity();
                  vel_ba = -vab;

                  sph_part_a->incDensityRate(mass_b *
                                        (vab * gradWab_a));
                  sph_part_b->incDensityRate(mass_a *
                                        (vel_ba * gradWba_b));

                  vel_radial = vab * (pos_a - pos_b);
                  if (vel_radial < 0) {
                    mu_ab = smoothLength * vel_radial /
                            (rab * rab + 0.01 * smoothLength * smoothLength);
                Gamma_ab =
(-alpha*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rho_a+rho_b)*2;
                  } else {
                    Gamma_ab = 0;
                  }

                  dva_dt = -mass_b *
                           (p_rhoSq_a * coeff_a +
                            p_rhoSq_b * coeff_b + Gamma_ab) *
                           gradWab_a;
                  sph_part_a->incAcceleration(dva_dt);
                  dvb_dt = -mass_a *
                           (p_rhoSq_b * coeff_b +
                            p_rhoSq_a * coeff_a + Gamma_ab) *
                           gradWba_b;
                  sph_part_b->incAcceleration(dvb_dt);

                  delta_a = epsilon * mass_b * (-vab) * Wab /
                            (rho_a + rho_b) * 2;
                  sph_part_a->incVelocityCorrection(delta_a);
                  delta_b = epsilon * mass_a * (vab)*Wba /
                            (rho_a + rho_b) * 2;
                  sph_part_b->incVelocityCorrection(delta_b);
                  break;
                case 2: // ghost with free
                  demt = sph_part_a->getDemParticle();
                  dem_pos = demt->currentPosition();
                  local_b =
                    demt->globalToLocal(pos_b - dem_pos); // the
                  // local position of sph point sph_part_b
                  ra = demt->getA();
                  rb = demt->getB();
                  rc = demt->getC();
                  k = 1.0 / (sqrt(local_b.x() * local_b.x() / (ra * ra) +
                                  local_b.y() * local_b.y() / (rb * rb) +
                                  local_b.z() * local_b.z() / (rc * rc)));
                  da = dem::vfabs(local_b -
                                  k * local_b); // the distance is the same
                                                  // in rotated coordinates

                  // calculate Vab as the method shown in Morris's paper, 1996

                  // (1) here I wanna use the distance from the point a/b to the
                  // surface of the ellipsoid to simplify the problem
                  local_a = sph_part_a->getLocalPosition();
                  k = 1.0 / (sqrt(local_a.x() * local_a.x() / (ra * ra) +
                                  local_a.y() * local_a.y() / (rb * rb) +
                                  local_a.z() * local_a.z() / (rc * rc)));
                  dB = dem::vfabs(local_a - k * local_a);

                  //              // (2) here the Morris's method is used
                  //              xa = pos_a.x(); ya = pos_a.y();
                  //              xB = pos_b.x(); yB = pos_b.y();
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

                  vdem = sph_part_a->getVelocity();
                  vel_ba = beta * (sph_part_b->getVelocity() - vdem);
                  vab = -vel_ba;

                  sph_part_a->incDensityRate(mass_b *
                                        (vab * gradWab_a));
                  sph_part_b->incDensityRate(mass_a *
                                        (vel_ba * gradWba_b));

                  vel_radial = vab * (pos_a - pos_b);
                  if (vel_radial < 0) {
                    mu_ab = smoothLength * vel_radial /
                            (rab * rab + 0.01 * smoothLength * smoothLength);
                    Gamma_ab =
(-alpha_zero*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rho_a+rho_b)*2;
                  } else {
                    Gamma_ab = 0;
                  }

                  dva_dt = -mass_b *
                           (p_rhoSq_a * coeff_a +
                            p_rhoSq_b * coeff_b + Gamma_ab) *
                           gradWab_a;
                  demt->addForce(mass_a * dva_dt);
                  demt->addMoment((pos_a - dem_pos) %
                                  (mass_a * dva_dt));
                  //                sph_part_a->incAcceleration(dva_dt);

                  dvb_dt = -mass_a *
                           (p_rhoSq_b * coeff_b +
                            p_rhoSq_a * coeff_a + Gamma_ab) *
                           gradWba_b;
                  sph_part_b->incAcceleration(dvb_dt); // the velocities of the ghost
                  particles will not envolve as the same way as others

                    //                  delta_a =
                    // epsilon*mass_b*(-vab)*Wab/(rho_a+rho_b)*2;
                    //                  sph_part_a->incVelocityCorrection(delta_a);
                    delta_b = epsilon * mass_a * (vab)*Wba /
                              (rho_a + rho_b) * 2;
                  sph_part_b->incVelocityCorrection(delta_b);
                  break;

                case 3: // boundary with free

                  // calculate Vab as the method shown in Morris's paper, 1996
                  // interact with boundary particles
                  da = pos_b.z() - allContainer.getMinCorner().z(); //
                  // assume with the bottom boundary
                  dB = allContainer.getMinCorner().z() -
                       pos_a.x(); // assume with
                                         // the bottom boundary
                  if (pos_a.x() <
                      allContainer.getMinCorner().x()) { // with left
                                                         // boundary
                    da = pos_b.x() - allContainer.getMinCorner().x();
                    dB = allContainer.getMinCorner().x() - pos_a.x();
                  } else if (pos_a.x() >
                             allContainer.getMaxCorner().x()) { // with right
                                                                // boundary
                    da = allContainer.getMaxCorner().x() - pos_b.x();
                    dB = pos_a.x() - allContainer.getMaxCorner().x();
                  } else if (pos_a.z() >
                             allContainer.getMaxCorner().z()) { // with top
                                                                // boundary
                    da = allContainer.getMaxCorner().z() - pos_b.z();
                    dB = pos_a.z() - allContainer.getMaxCorner().z();
                  }

                  beta = 1 + dB / da;
                  if (beta > 2 || beta < 0 || isnan(beta)) {
                    beta = 2;
                  }

                  vel_ba = beta * sph_part_b->getVelocity();
                  vab = -vel_ba;

                  sph_part_a->incDensityRate(mass_b *
                                        (vab * gradWab_a));
                  sph_part_b->incDensityRate(mass_a *
                                        (vel_ba * gradWba_b));

                  vel_radial = vab * (pos_a - pos_b);
                  if (vel_radial < 0) {
                    mu_ab = smoothLength * vel_radial /
                            (rab * rab + 0.01 * smoothLength * smoothLength);
                Gamma_ab =
(-alpha_zero*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rho_a+rho_b)*2;
                  } else {
                    Gamma_ab = 0;
                  }

                  //                dva_dt =
                  //-mass_b*(press_a/(rho_a*rho_a)*coeff_a+press_b/(rho_b*rho_b)*coeff_b+Gamma_ab)*gradWab_a;
                  //                sph_part_a->incAcceleration(dva_dt);
                  dvb_dt =
                    -mass_a *
                    (p_rhoSq_b + p_rhoSq_a + Gamma_ab) *
                    gradWba_b;
                  sph_part_b->incAcceleration(dvb_dt); // the velocities of the ghost
                  particles will not envolve as the same way as others

                    //            delta_a =
                    // epsilon*mass_b*(-vab)*Wab/(rho_a+rho_b)*2;
                    //                sph_part_a->incVelocityCorrection(delta_a);
                    delta_b = epsilon * mass_a * (vab)*Wba /
                              (rho_a + rho_b) * 2;
                  sph_part_b->incVelocityCorrection(delta_b);

                  //            // apply the boundary forces by Lennard-Jones
                  //            potential as in
                  // Monaghan's paper(1994)
                  //                if(rab<=spaceInterval){ // sph_part_b is in the
                  //                smooth kernel
                  //                dvb_dt = D*(pow(spaceInterval/rab,
                  //                p1)-pow(spaceInterval/rab,
                  // p2))*(pos_b-pos_a)/(rab*rab);
                  //                sph_part_b->incAcceleration(dvb_dt);
                  //                } // end if

                  break;
                default:
                  std::cout << "SPH particle type of sph_part_a should be 1, 2 or 3!"
                            << std::endl;
                  exit(-1);
              } // end switch sph_part_a

              break;
            case 2: // sph_part_b is ghost particle

              if (sph_part_a->getType() !=
                  1) { // sph_part_a is not free sph particles, i.e. sph_part_a
                       // is ghost or boundary particles
                break; // do not consider sph_part_a, as we treat before
              }
              demt = sph_part_b->getDemParticle();
              dem_pos = demt->currentPosition();
              local_a =
                demt->globalToLocal(pos_a - dem_pos); // the local
              // position of sph point sph_part_a
              ra = demt->getA();
              rb = demt->getB();
              rc = demt->getC();
              k = 1.0 / (sqrt(local_a.x() * local_a.x() / (ra * ra) +
                              local_a.y() * local_a.y() / (rb * rb) +
                              local_a.z() * local_a.z() / (rc * rc)));
              da = dem::vfabs(local_a -
                              k * local_a); // the distance is the same
                                              // in rotated coordinates

              // calculate Vab as the method shown in Morris's paper, 1996

              // (1) here I wanna use the distance from the point a/b to the
              // surface of the ellipsoid to simplify the problem
              local_b = sph_part_b->getLocalPosition();
              k = 1.0 / (sqrt(local_b.x() * local_b.x() / (ra * ra) +
                              local_b.y() * local_b.y() / (rb * rb) +
                              local_b.z() * local_b.z() / (rc * rc)));
              dB = dem::vfabs(local_b - k * local_b);

              //              // (2) here the Morris's method is used
              //              xa = pos_a.x(); ya = pos_a.gety();
              //              xB = pos_b.x(); yB = pos_b.gety();
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

              vdem = sph_part_b->getVelocity();
              vab = beta * (sph_part_a->getVelocity() - vdem);
              vel_ba = -vab;

              sph_part_a->incDensityRate(mass_b *
                                    (vab * gradWab_a));
              sph_part_b->incDensityRate(mass_a *
                                    (vel_ba * gradWba_b));

              vel_radial = vab * (pos_a - pos_b);
              if (vel_radial < 0) {
                mu_ab = smoothLength * vel_radial /
                        (rab * rab + 0.01 * smoothLength * smoothLength);
                Gamma_ab =
(-alpha_zero*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rho_a+rho_b)*2;
              } else {
                Gamma_ab = 0;
              }

              dva_dt = -mass_b *
                       (p_rhoSq_a * coeff_a +
                        p_rhoSq_b * coeff_b + Gamma_ab) *
                       gradWab_a;
              sph_part_a->incAcceleration(dva_dt);
              dvb_dt = -mass_a *
                       (p_rhoSq_b * coeff_b +
                        p_rhoSq_a * coeff_a + Gamma_ab) *
                       gradWba_b;
              demt->addForce(mass_b * dvb_dt);
              demt->addMoment((pos_b - dem_pos) %
                              (mass_b * dvb_dt));
              //              sph_part_b->incAcceleration(dvb_dt);  // the velocities
              //              of the ghost
              // particles will not envolve as the same way as others

              delta_a = epsilon * mass_b * (-vab) * Wab /
                        (rho_a + rho_b) * 2;
              sph_part_a->incVelocityCorrection(delta_a);
              //                delta_b =
              // epsilon*mass_a*(vab)*Wba/(rho_a+rho_b)*2;
              //                sph_part_b->incVelocityCorrection(delta_b);
              break;
            case 3: // sph_part_b is boundary particle

              if (sph_part_a->getType() !=
                  1) { // sph_part_a is not free sph particles, i.e. sph_part_a
                       // is ghost or boundary particles
                break; // do not consider sph_part_a, as we treat before
              }

              // calculate Vab as the method shown in Morris's paper, 1996
              // interact with boundary particles
              da = pos_a.z() - allContainer.getMinCorner().z(); // assume
              // with the bottom boundary
              dB = allContainer.getMinCorner().z() -
                   pos_b.z(); // assume with
                                     // the bottom boundary
              if (pos_b.x() <
                  allContainer.getMinCorner().x()) { // with left
                                                     // boundary
                da = pos_a.x() - allContainer.getMinCorner().x();
                dB = allContainer.getMinCorner().x() - pos_b.x();
              } else if (pos_b.x() >
                         allContainer.getMaxCorner().x()) { // with
                                                            // right boundary
                da = allContainer.getMaxCorner().x() - pos_a.x();
                dB = pos_b.x() - allContainer.getMaxCorner().x();
              } else if (pos_b.z() >
                         allContainer.getMaxCorner().z()) { // with
                                                            // top boundary
                da = allContainer.getMaxCorner().z() - pos_a.z();
                dB = pos_b.z() - allContainer.getMaxCorner().z();
              }

              beta = 1 + dB / da;
              if (beta > 2 || beta < 0 || isnan(beta)) {
                beta = 2;
              }

              vab = beta * sph_part_a->getVelocity();
              vel_ba = -vab;

              sph_part_a->incDensityRate(mass_b *
                                    (vab * gradWab_a));
              sph_part_b->incDensityRate(mass_a *
                                    (vel_ba * gradWba_b));

              vel_radial = vab * (pos_a - pos_b);
              if (vel_radial < 0) {
                mu_ab = smoothLength * vel_radial /
                        (rab * rab + 0.01 * smoothLength * smoothLength);
            Gamma_ab =
(-alpha_zero*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rho_a+rho_b)*2;
              } else {
                Gamma_ab = 0;
              }

              dva_dt = -mass_b *
                       (p_rhoSq_a * coeff_a +
                        p_rhoSq_b * coeff_b + Gamma_ab) *
                       gradWab_a;
              sph_part_a->incAcceleration(dva_dt);
              //              dvb_dt =
              //-mass_a*(press_b/(rho_b*rho_b)+press_a/(rho_a*rho_a)+Gamma_ab)*gradWba_b;
              //              sph_part_b->incAcceleration(dvb_dt);  // the velocities
              //              of the ghost
              // particles will not envolve as the same way as others

              delta_a = epsilon * mass_b * (-vab) * Wab /
                        (rho_a + rho_b) * 2;
              sph_part_a->incVelocityCorrection(delta_a);
              //              delta_b =
              // epsilon*mass_a*(vab)*Wba/(rho_a+rho_b)*2;
              //              sph_part_b->incVelocityCorrection(delta_b);

              //          // apply the boundary forces by Lennard-Jones
              //          potential as in
              // Monaghan's paper(1994)
              //              if(rab<=spaceInterval){ // sph_part_b is in the smooth
              //              kernel
              //            dva_dt = D*(pow(spaceInterval/rab,
              //            p1)-pow(spaceInterval/rab,
              // p2))*(pos_a-pos_b)/(rab*rab);
              //            sph_part_a->incAcceleration(dva_dt);
              //              } // end if
              break;
            default:
              std::cout << "SPH particle type should be 1, 2 or 3!"
                        << std::endl;
              exit(-1);

          } // end swtich type
        }   // end if 3h
      }     // end for sph_part_b in the same cell

      for (sph_part_b = tmp_particleVec.begin(); sph_part_b != tmp_particleVec.end();
           sph_part_b++) { // all
                    // particles in cellIndex's neighboring cells
        pos_b = sph_part_b->currentPosition();
        rab = dem::vfabs(pos_a - pos_b);
        if (rab <= kernelSize) { // sph_part_b is in the smooth kernel
          press_b = sph_part_b->getPressure();
          //          mu_b = sph_part_b->getViscosity();
          rho_b = sph_part_b->getDensity();

          Wq = kernelFunction(rab / smoothLength);
          phi_4 = pow(Wq / Wqmin, 4);
          if (press_b >= 0) {
            Rb = 0.006;
          } else {
            Rb = 0.6;
          }
          coeff_a = 1 + Ra * phi_4;
          coeff_b = 1 + Rb * phi_4;

          // add the density dot for sph_part_a and sph_part_b
          gradWab_a = gradientKernelFunction(pos_a, pos_b); //
          // this is to add SPH sph_part_a
          gradWba_b = -gradWab_a;
          dWab_dra = partialKernelFunction(pos_a, pos_b); // this
          // is to add SPH sph_part_a
          dWba_drb = dWab_dra;
          Wab = kernelFunction(pos_a, pos_b); // this is to add
                                                            // SPH sph_part_a
          Wba = Wab;
          switch (sph_part_b->getType()) {
            case 1: // sph_part_b is free SPH particle
              switch (sph_part_a->getType()) {
                case 1: // free with free
                  vab = sph_part_a->getVelocity() - sph_part_b->getVelocity();
                  vel_ba = -vab;

                  sph_part_a->incDensityRate(mass_b *
                                        (vab * gradWab_a));
                  sph_part_b->incDensityRate(mass_a *
                                        (vel_ba * gradWba_b));

                  vel_radial = vab * (pos_a - pos_b);
                  if (vel_radial < 0) {
                    mu_ab = smoothLength * vel_radial /
                            (rab * rab + 0.01 * smoothLength * smoothLength);
                Gamma_ab =
(-alpha*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rho_a+rho_b)*2;
                  } else {
                    Gamma_ab = 0;
                  }

                  dva_dt = -mass_b *
                           (p_rhoSq_a * coeff_a +
                            p_rhoSq_b * coeff_b + Gamma_ab) *
                           gradWab_a;
                  sph_part_a->incAcceleration(dva_dt);
                  dvb_dt = -mass_a *
                           (p_rhoSq_b * coeff_b +
                            p_rhoSq_a * coeff_a + Gamma_ab) *
                           gradWba_b;
                  sph_part_b->incAcceleration(dvb_dt);

                  delta_a = epsilon * mass_b * (-vab) * Wab /
                            (rho_a + rho_b) * 2;
                  sph_part_a->incVelocityCorrection(delta_a);
                  delta_b = epsilon * mass_a * (vab)*Wba /
                            (rho_a + rho_b) * 2;
                  sph_part_b->incVelocityCorrection(delta_b);
                  break;
                case 2: // ghost with free

                  demt = sph_part_a->getDemParticle();
                  dem_pos = demt->currentPosition();
                  local_b =
                    demt->globalToLocal(pos_b - dem_pos); // the
                  // local position of sph point sph_part_b
                  ra = demt->getA();
                  rb = demt->getB();
                  rc = demt->getC();
                  k = 1.0 / (sqrt(local_b.x() * local_b.x() / (ra * ra) +
                                  local_b.y() * local_b.y() / (rb * rb) +
                                  local_b.z() * local_b.z() / (rc * rc)));
                  da = dem::vfabs(local_b -
                                  k * local_b); // the distance is the same
                                                  // in rotated coordinates

                  // calculate Vab as the method shown in Morris's paper, 1996

                  // (1) here I wanna use the distance from the point a/b to the
                  // surface of the ellipsoid to simplify the problem
                  local_a = sph_part_a->getLocalPosition();
                  k = 1.0 / (sqrt(local_a.x() * local_a.x() / (ra * ra) +
                                  local_a.y() * local_a.y() / (rb * rb) +
                                  local_a.z() * local_a.z() / (rc * rc)));
                  dB = dem::vfabs(local_a - k * local_a);

                  //              // (2) here the Morris's method is used
                  //              xa = pos_a.x(); ya = pos_a.y();
                  //              xB = pos_b.x(); yB = pos_b.y();
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

                  vdem = sph_part_a->getVelocity();
                  vel_ba = beta * (sph_part_b->getVelocity() - vdem);
                  vab = -vel_ba;

                  sph_part_a->incDensityRate(mass_b *
                                        (vab * gradWab_a));
                  sph_part_b->incDensityRate(mass_a *
                                        (vel_ba * gradWba_b));

                  vel_radial = vab * (pos_a - pos_b);
                  if (vel_radial < 0) {
                    mu_ab = smoothLength * vel_radial /
                            (rab * rab + 0.01 * smoothLength * smoothLength);
                    Gamma_ab =
(-alpha_zero*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rho_a+rho_b)*2;
                  } else {
                    Gamma_ab = 0;
                  }

                  dva_dt = -mass_b *
                           (p_rhoSq_a * coeff_a +
                            p_rhoSq_b * coeff_b + Gamma_ab) *
                           gradWab_a;
                  demt->addForce(mass_a * dva_dt);
                  demt->addMoment((pos_a - dem_pos) %
                                  (mass_a * dva_dt));
                  //                sph_part_a->incAcceleration(dva_dt);

                  dvb_dt = -mass_a *
                           (p_rhoSq_b * coeff_b +
                            p_rhoSq_a * coeff_a + Gamma_ab) *
                           gradWba_b;
                  sph_part_b->incAcceleration(dvb_dt); // the velocities of the ghost
                  // particles will not envolve as the same way as others

                  //                  delta_a =
                  // epsilon*mass_b*(-vab)*Wab/(rho_a+rho_b)*2;
                  //                  sph_part_a->incVelocityCorrection(delta_a);
                  delta_b = epsilon * mass_a * (vab)*Wba /
                            (rho_a + rho_b) * 2;
                  sph_part_b->incVelocityCorrection(delta_b);
                  break;
                case 3: // boundary with free

                  // calculate Vab as the method shown in Morris's paper, 1996
                  // interact with boundary particles
                  da = pos_b.z() - allContainer.getMinCorner().z(); //
                  // assume with the bottom boundary
                  dB = allContainer.getMinCorner().z() -
                       pos_a.x(); // assume with
                                         // the bottom boundary
                  if (pos_a.x() <
                      allContainer.getMinCorner().x()) { // with left
                                                         // boundary
                    da = pos_b.x() - allContainer.getMinCorner().x();
                    dB = allContainer.getMinCorner().x() - pos_a.x();
                  } else if (pos_a.x() >
                             allContainer.getMaxCorner().x()) { // with right
                                                                // boundary
                    da = allContainer.getMaxCorner().x() - pos_b.x();
                    dB = pos_a.x() - allContainer.getMaxCorner().x();
                  } else if (pos_a.z() >
                             allContainer.getMaxCorner().z()) { // with top
                                                                // boundary
                    da = allContainer.getMaxCorner().z() - pos_b.z();
                    dB = pos_a.z() - allContainer.getMaxCorner().z();
                  }

                  beta = 1 + dB / da;
                  if (beta > 2 || beta < 0 || isnan(beta)) {
                    beta = 2;
                  }

                  vel_ba = beta * sph_part_b->getVelocity();
                  vab = -vel_ba;

                  sph_part_a->incDensityRate(mass_b *
                                        (vab * gradWab_a));
                  sph_part_b->incDensityRate(mass_a *
                                        (vel_ba * gradWba_b));

                  vel_radial = vab * (pos_a - pos_b);
                  if (vel_radial < 0) {
                    mu_ab = smoothLength * vel_radial /
                            (rab * rab + 0.01 * smoothLength * smoothLength);
                Gamma_ab =
(-alpha_zero*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rho_a+rho_b)*2;
                  } else {
                    Gamma_ab = 0;
                  }

                  //                dva_dt =
                  //-mass_b*(press_a/(rho_a*rho_a)*coeff_a+press_b/(rho_b*rho_b)*coeff_b+Gamma_ab)*gradWab_a;
                  //                sph_part_a->incAcceleration(dva_dt);
                  dvb_dt =
                    -mass_a *
                    (p_rhoSq_b + p_rhoSq_a + Gamma_ab) *
                    gradWba_b;
                  sph_part_b->incAcceleration(dvb_dt); // the velocities of the ghost
                  particles will not envolve as the same way as others

                    //            delta_a =
                    // epsilon*mass_b*(-vab)*Wab/(rho_a+rho_b)*2;
                    //                sph_part_a->incVelocityCorrection(delta_a);
                    delta_b = epsilon * mass_a * (vab)*Wba /
                              (rho_a + rho_b) * 2;
                  sph_part_b->incVelocityCorrection(delta_b);

                  // apply the boundary forces by Lennard-Jones potential as in
                  // Monaghan's paper(1994)
                  if (rab <= spaceInterval) { // sph_part_b is in the smooth kernel
                    dvb_dt = D * (pow(spaceInterval / rab, p1) -
                                  pow(spaceInterval / rab, p2)) *
                             (pos_b - pos_a) / (rab * rab);
                    sph_part_b->incAcceleration(dvb_dt);
                  } // end if

                  break;
                default:
                  std::cout << "SPH particle type of sph_part_a should be 1, 2 or 3!"
                            << std::endl;
                  exit(-1);
              } // end switch sph_part_a

              break;
            case 2: // sph_part_b is ghost particle

              if (sph_part_a->getType() !=
                  1) { // sph_part_a is not free sph particles, i.e. sph_part_a
                       // is ghost or boundary particles
                break; // do not consider sph_part_a, as we treat before
              }

              demt = sph_part_b->getDemParticle();
              dem_pos = demt->currentPosition();
              local_a =
                demt->globalToLocal(pos_a - dem_pos); // the local
              // position of sph point sph_part_a
              ra = demt->getA();
              rb = demt->getB();
              rc = demt->getC();
              k = 1.0 / (sqrt(local_a.x() * local_a.x() / (ra * ra) +
                              local_a.y() * local_a.y() / (rb * rb) +
                              local_a.z() * local_a.z() / (rc * rc)));
              da = dem::vfabs(local_a -
                              k * local_a); // the distance is the same
                                              // in rotated coordinates

              // calculate Vab as the method shown in Morris's paper, 1996

              // (1) here I wanna use the distance from the point a/b to the
              // surface of the ellipsoid to simplify the problem
              local_b = sph_part_b->getLocalPosition();
              k = 1.0 / (sqrt(local_b.x() * local_b.x() / (ra * ra) +
                              local_b.y() * local_b.y() / (rb * rb) +
                              local_b.z() * local_b.z() / (rc * rc)));
              dB = dem::vfabs(local_b - k * local_b);

              //              // (2) here the Morris's method is used
              //              xa = pos_a.x(); ya = pos_a.y();
              //              xB = pos_b.x(); yB = pos_b.y();
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

              vdem = sph_part_b->getVelocity();
              vab = beta * (sph_part_a->getVelocity() - vdem);
              vel_ba = -vab;

              sph_part_a->incDensityRate(mass_b *
                                    (vab * gradWab_a));
              sph_part_b->incDensityRate(mass_a *
                                    (vel_ba * gradWba_b));

              vel_radial = vab * (pos_a - pos_b);
              if (vel_radial < 0) {
                mu_ab = smoothLength * vel_radial /
                        (rab * rab + 0.01 * smoothLength * smoothLength);
                Gamma_ab =
(-alpha_zero*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rho_a+rho_b)*2;
              } else {
                Gamma_ab = 0;
              }

              dva_dt = -mass_b *
                       (p_rhoSq_a * coeff_a +
                        p_rhoSq_b * coeff_b + Gamma_ab) *
                       gradWab_a;
              sph_part_a->incAcceleration(dva_dt);
              dvb_dt = -mass_a *
                       (p_rhoSq_b * coeff_b +
                        p_rhoSq_a * coeff_a + Gamma_ab) *
                       gradWba_b;
              demt->addForce(mass_b * dvb_dt);
              demt->addMoment((pos_b - dem_pos) %
                              (mass_b * dvb_dt));
              //              sph_part_b->incAcceleration(dvb_dt);  // the velocities
              //              of the ghost
              // particles will not envolve as the same way as others

              delta_a = epsilon * mass_b * (-vab) * Wab /
                        (rho_a + rho_b) * 2;
              sph_part_a->incVelocityCorrection(delta_a);
              //                delta_b =
              // epsilon*mass_a*(vab)*Wba/(rho_a+rho_b)*2;
              //                sph_part_b->incVelocityCorrection(delta_b);
              break;
            case 3: // sph_part_b is boundary particle

              if (sph_part_a->getType() !=
                  1) { // sph_part_a is not free sph particles, i.e. sph_part_a
                       // is ghost or boundary particles
                break; // do not consider sph_part_a, as we treat before
              }

              // calculate Vab as the method shown in Morris's paper, 1996
              // interact with boundary particles
              da = pos_a.z() - allContainer.getMinCorner().z(); // assume
              // with the bottom boundary
              dB = allContainer.getMinCorner().z() -
                   pos_b.z(); // assume with
                                     // the bottom boundary
              if (pos_b.x() <
                  allContainer.getMinCorner().x()) { // with left
                                                     // boundary
                da = pos_a.x() - allContainer.getMinCorner().x();
                dB = allContainer.getMinCorner().x() - pos_b.x();
              } else if (pos_b.x() >
                         allContainer.getMaxCorner().x()) { // with
                                                            // right boundary
                da = allContainer.getMaxCorner().x() - pos_a.x();
                dB = pos_b.x() - allContainer.getMaxCorner().x();
              } else if (pos_b.z() >
                         allContainer.getMaxCorner().z()) { // with
                                                            // top boundary
                da = allContainer.getMaxCorner().z() - pos_a.z();
                dB = pos_b.z() - allContainer.getMaxCorner().z();
              }

              beta = 1 + dB / da;
              if (beta > 2 || beta < 0 || isnan(beta)) {
                beta = 2;
              }

              vab = beta * sph_part_a->getVelocity();
              vel_ba = -vab;

              sph_part_a->incDensityRate(mass_b *
                                    (vab * gradWab_a));
              sph_part_b->incDensityRate(mass_a *
                                    (vel_ba * gradWba_b));

              vel_radial = vab * (pos_a - pos_b);
              if (vel_radial < 0) {
                mu_ab = smoothLength * vel_radial /
                        (rab * rab + 0.01 * smoothLength * smoothLength);
            Gamma_ab =
(-alpha_zero*(util::getParam<REAL>("soundSpeed")*mu_ab)/(rho_a+rho_b)*2;
              } else {
                Gamma_ab = 0;
              }

              dva_dt = -mass_b *
                       (p_rhoSq_a * coeff_a +
                        p_rhoSq_b * coeff_b + Gamma_ab) *
                       gradWab_a;
              sph_part_a->incAcceleration(dva_dt);
              //              dvb_dt =
              //-mass_a*(press_b/(rho_b*rho_b)+press_a/(rho_a*rho_a)+Gamma_ab)*gradWba_b;
              //              sph_part_b->incAcceleration(dvb_dt);  // the velocities
              //              of the ghost
              // particles will not envolve as the same way as others

              delta_a = epsilon * mass_b * (-vab) * Wab /
                        (rho_a + rho_b) * 2;
              sph_part_a->incVelocityCorrection(delta_a);
              //              delta_b =
              // epsilon*mass_a*(vab)*Wba/(rho_a+rho_b)*2;
              //              sph_part_b->incVelocityCorrection(delta_b);

              // apply the boundary forces by Lennard-Jones potential as in
              // Monaghan's paper(1994)
              if (rab <= spaceInterval) { // sph_part_b is in the smooth kernel
                dva_dt = D * (pow(spaceInterval / rab, p1) -
                              pow(spaceInterval / rab, p2)) *
                         (pos_a - pos_b) / (rab * rab);
                sph_part_a->incAcceleration(dva_dt);
              } // end if

              break;
            default:
              std::cout << "SPH particle type should be 1, 2 or 3!"
                        << std::endl;
              exit(-1);

          } // end swtich type
        }   // end if 3h
      }     // end for sph_part_b in neighbor cells
    }       // end for sph_part_a

    tmp_particleVec.clear(); // clear elements in tmp-vector for particles
                             // neighboring cells, it is important

  } // end for cellIndex, different cells

  //      // apply the boundary forces by Lennard-Jones potential as in
  //      Monaghan's
  // paper(1994)
  //      for(sph_part_b=SPHBoundaryParticleVec.begin();
  // sph_part_b!=SPHBoundaryParticleVec.end(); sph_part_b++){
  //    pos_b = sph_part_b->currentPosition();
  //    rab = dem::vfabs(pos_a-pos_b);
  //    if(rab<=spaceInterval){ // sph_part_b is in the smooth kernel
  //        dva_dt = D*(pow(spaceInterval/rab, p1)-pow(spaceInterval/rab,
  // p2))*(pos_a-pos_b)/(rab*rab);
  //        sph_part_a->incAcceleration(dva_dt);
  //    } // end if
  //      } // end sph_part_b

} // end calculateSPHDensityRateAccelerationLinkedList3D()
*/

/*
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

    ofs << std::setw(20) << (*pt)->currentPosition().x() << std::setw(20)
        << (*pt)->currentPosition().y() << std::setw(20)
        << (*pt)->currentPosition().z() << std::setw(20)
        << (*pt)->getDisplacement().x() << std::setw(20)
        << (*pt)->getDisplacement().y() << std::setw(20)
        << (*pt)->getDisplacement().z() << std::setw(20)
        << (*pt)->getVelocity().x() << std::setw(20) << (*pt)->getVelocity().y()
        << std::setw(20) << (*pt)->getVelocity().z();

    (*pt)->calculateParticlePressure();
    ofs << std::setw(20) << (*pt)->getPressure();

    ofs << std::setw(20) << (*pt)->getAcceleration().x() << std::setw(20)
        << (*pt)->getAcceleration().y() << std::setw(20)
        << (*pt)->getAcceleration().z() << std::setw(20)
        << (*pt)->getDensityRate() << std::setw(20)
        << (*pt)->getDensity() << std::endl;
  }

  ofs.close();
}
*/

namespace sph {
template void
SmoothParticleHydro::assignParticlesToPatchGrid<2>(const Box& container,
  const REAL& bufferWidth, const REAL& ghostWidth, const REAL& kernelSize);
template void
SmoothParticleHydro::assignParticlesToPatchGrid<3>(const Box& container,
  const REAL& bufferWidth, const REAL& ghostWidth, const REAL& kernelSize);
}