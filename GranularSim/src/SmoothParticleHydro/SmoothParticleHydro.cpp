#include <DiscreteElements/DEMParticle.h>
#include <SmoothParticleHydro/SmoothParticleHydro.h>
#include <SmoothParticleHydro/SPHInteractions.h>

#include <Core/Const/Constants.h>
#include <Core/Math/IntVec.h>
#include <Core/Math/Matrix.h>
#include <Core/Util/Utility.h>
#include <InputOutput/OutputTecplot.h>
#include <InputOutput/OutputVTK.h>
#include <chrono>


using Timer = std::chrono::steady_clock;
using Seconds = std::chrono::seconds;
using IntVec = dem::IntVec;
using Vec = dem::Vec;
using Matrix = dem::Matrix;
using Box = dem::Box;
using DEMParticlePArray = dem::DEMParticlePArray;
using communicator = boost::mpi::communicator;

namespace sph {

using OutputVTK = dem::OutputVTK<SPHParticlePArray>;
using OutputTecplot = dem::OutputTecplot<SPHParticlePArray>;

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
//        (patchGrid is expanded to include the boundary sph particles)
//    (3) this partition method here is not very suitable for the free surface
// flow problem, such as bursting dam problem
//        since the total domain is divided, while there are lots of voids in
// the domain. (partition only free sph particles will be better)
//        But our goal is to simulate the porous media in triaxial, insotropic
// or SHPB simulations, particles are filled in the domain.
void
SmoothParticleHydro::scatterSPHParticle(const Box& spatialDomain,
                                        const REAL& ghostWidth,
                                        REAL& bufferLength)
{
  // Comute the buffer size for the computational domain
  // int numLayers = util::getParam<int>("numLayers");
  // REAL spaceInterval = util::getParam<REAL>("spaceInterval");
  // bufferLength = spaceInterval * numLayers;

  // Set d_sphPatchBox on all procs
  setPatchBox(Box(spatialDomain, bufferLength));

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

    Vec v1 = d_sphPatchBox.minCorner();
    Vec v2 = d_sphPatchBox.maxCorner();
    Vec vspan = (v2 - v1) / d_mpiProcs;

    auto reqs = new boost::mpi::request[d_mpiSize - 1];

    for (int iRank = d_mpiSize - 1; iRank >= 0; --iRank) {

      d_patchSPHParticleVec.clear();

      int ndim = 3;
      IntVec coords;
      MPI_Cart_coords(d_cartComm, iRank, ndim, coords.data());

      Vec lower = v1 + vspan * coords;
      Vec upper = lower + vspan;
      Box domain(lower, upper);

      findSPHParticleInBox(domain, d_allSPHParticleVec,
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
SmoothParticleHydro::findSPHParticleInBox(const Box& domain,
                                          const SPHParticlePArray& sphParticles,
                                          SPHParticlePArray& insideParticles)
{
  for (const auto& particle : sphParticles) {
    // it is critical to use EPS
    if (domain.inside(particle->currentPosition(), dem::EPS)) {
      insideParticles.push_back(particle);
    }
  }
}

void
SmoothParticleHydro::createPatch(int iteration, const REAL& ghostWidth)
{
  // determine domain of each process
  Vec v1 = d_sphPatchBox.minCorner();
  Vec v2 = d_sphPatchBox.maxCorner();
  Vec vspan = (v2 - v1) / d_mpiProcs;
  Vec lower = v1 + vspan * d_mpiCoords;
  Vec upper = lower + vspan;
  d_patchP = std::make_unique<SPHPatch>(d_cartComm, d_mpiRank, d_mpiCoords,
                                        lower, upper, ghostWidth, dem::EPS);
}

void
SmoothParticleHydro::updatePatch(int iteration, const REAL& ghostWidth)
{
  // determine domain of each process
  Vec v1 = d_sphPatchBox.minCorner();
  Vec v2 = d_sphPatchBox.maxCorner();
  Vec vspan = (v2 - v1) / d_mpiProcs;
  Vec lower = v1 + vspan * d_mpiCoords;
  Vec upper = lower + vspan;
  Box domain(lower, upper);
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
  Vec domainWidth = d_sphPatchBox.maxCorner() - d_sphPatchBox.minCorner();
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

//// the momentum equilibrium equation and state equation are implemented as
// in Monaghan's paper (1994), simulate free surface flow using sph
//// here the neighboring list of SPH particles is searched by the cells,
template <int dim>
void
SmoothParticleHydro::updateParticleInteractions(const Box& spatialDomain,
                                                const REAL& bufferWidth,
                                                const REAL& ghostWidth,
                                                const REAL& kernelSize,
                                                const REAL& smoothLength)
{
  // Initialize the rate quantities that are integrated
  initializeDensityRateAndAcceleration();

  // divide the SPH domain into different cells, each cell may contain SPH
  // particles within it
  assignParticlesToPatchGrid<dim>(spatialDomain, bufferWidth, ghostWidth, kernelSize);

  for (auto cellIndex = 0u; cellIndex < d_sphPatchGrid.size(); ++cellIndex) {

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
                                         kernelSize, smoothLength, spatialDomain);

      } // end for sph_part_b in the same cell

      for (auto& sph_part_b : neighborParticles) {

        // Go to next particle if the distance is too large
        if (sph_part_b->isOutsideInfluenceZone(*sph_part_a, kernelSize)) continue;

        computeParticleInteractions<dim>(sph_part_a, sph_part_b,
                                         kernelSize, smoothLength, spatialDomain);

      } // end for sph_part_b in neighbor cells
    } // end for sph_part_a
  } // end for cellIndex, different cells
} // end updateParticleInteractions()

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
SmoothParticleHydro::assignParticlesToPatchGrid(const Box& domain,
                                                const REAL& bufferWidth,
                                                const REAL& ghostWidth,
                                                const REAL& kernelSize)
{
  REAL small_value = 0.01 * bufferWidth;
  REAL ghostBuffer = ghostWidth + small_value;

  // expand the domain by bufferWidth + small_value
  Box expandedContainer(domain, ghostBuffer);

  // Compute the number of cells in each dimension
  int nx = std::round(expandedContainer.dimX()/kernelSize) + 1;
  int ny = std::round(expandedContainer.dimY()/kernelSize) + 1;
  int nz = std::round(expandedContainer.dimZ()/kernelSize) + 1;
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
    int cellIndex = getCellIndex<dim>(expandedContainer.minCorner(),
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
                                               const IntVec& numGridCells) const
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
                                               const IntVec& numGridCells) const
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
                                               const IntVec& numGridCells) const
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
                                                 const Box& spatialDomain) const
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
      interaction.doDomainBoundaryVelocityCorrection(spatialDomain);
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
    interaction.doDomainBoundaryVelocityCorrection(spatialDomain);
    interaction.updateMomentumExchangeCoeffs(smoothLength, alpha, soundSpeed);
    interaction.updateInteractionBoundaryFree(sph_part_b, sph_part_a, epsilon);

  }
}

template <>
void
SmoothParticleHydro::initializeSPHVelocity<2>(const REAL& delT)
{
  // initial velocity only for free SPH particles based on equation (4.3)
  // and fix the y for all SPH points
  for (auto& particle : d_sphParticleVec) {
    switch (particle->getType()) {
      case SPHParticleType::FREE: // free sph particle
        particle->initialVelocityLeapFrog(delT);
        particle->fixY();
        break;
      case SPHParticleType::BOUNDARY: // boundary sph particle
        particle->fixXYZ();
        break;
      default:
        break;
    } // switch
  }

} // initialSPHVelocity()

template <>
void
SmoothParticleHydro::initializeSPHVelocity<3>(const REAL& delT)
{
  // initial velocity only for free SPH particles based on equation (4.3)
  for (auto& particle : d_sphParticleVec) {
    if (particle->getType() == SPHParticleType::FREE) {
     particle->initialVelocityLeapFrog(delT);
    }
  }

} // initialSPHVelocity()

// update particle position and density based on equation (4.1)
void
SmoothParticleHydro::updateSPHLeapFrogPositionDensity(const REAL& delT)
{ 
  for (auto& particle : d_sphParticleVec) {
    switch (particle->getType()) {
      case SPHParticleType::FREE:
      {
        particle->updatePositionDensityLeapFrog(delT);
        break;
      }
      case SPHParticleType::GHOST:
      {
        // update position and density by dem
        // particle and sph particle density, and velocity
        dem::DEMParticle* dem_particle = particle->getDEMParticle();
        Vec sph_pos_local = particle->getLocalPosition();
        Vec dem_pos_curr = dem_particle->currentPosition();
        Vec sph_pos_curr = dem_particle->localToGlobal(sph_pos_local) +
                           dem_pos_curr;
        Vec sph_dem_rel_pos = sph_pos_curr - dem_pos_curr;
        Vec sph_vel_curr = dem_particle->currentVelocity() +
                           cross(dem_particle->currentAngularVelocity(),
                                 sph_dem_rel_pos);
        particle->setCurrentPositionition(sph_pos_curr);
        particle->setCurrentVelocity(sph_vel_curr);
        particle->updateDensity(delT);
        break;
      }
      case SPHParticleType::BOUNDARY:
      {
        particle->updateDensity(delT);
        break;
      }
      default:
        break;
    } // switch
  }

} // end updateSPHLeapFrogPositionDensity

// update particle velocity only for free SPH particles based on equation (4.2)
void
SmoothParticleHydro::updateSPHLeapFrogVelocity(const REAL& delT)
{ 
  for (auto& particle : d_sphParticleVec) {
    if (particle->getType() == SPHParticleType::FREE) {
      particle->updateVelocity(delT);
    }
  }
} // end updateSPHLeapFrogVelocity


// Create the writer
void 
SmoothParticleHydro::createOutputWriter(const std::string& outputFolder, 
                                 const int& iter) 
{
  bool writeVTK = true;
  if (writeVTK) {
    d_writer = 
      std::make_unique<dem::OutputVTK<SPHParticlePArray> >(outputFolder, iter);
  } else {
    d_writer = 
      std::make_unique<dem::OutputTecplot<SPHParticlePArray> >(outputFolder, iter);
  }
}

void
SmoothParticleHydro::writeParticlesToFile(int frame, REAL time) const
{
  d_writer->writeParticles(&d_allSPHParticleVec, frame, time);
}

void
SmoothParticleHydro::printSPHParticle(const char* str) const
{
  std::ofstream ofs(str);
  if (!ofs) {
    dem::debugInf << "stream error: printSPHParticle" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(dem::OPREC);

  ofs << "Title = \"SPH Particle Information\"" << std::endl;
  ofs << "VARIABLES = \"x\", \"y\",\"z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\""
      << "\"Vy\" \"Vz\" \"Pressure\" \"a_x\" \"a_y\" \"a_z\" \"density_dot\" "
      << "\"density\" "
      << std::endl;

  // Output the coordinates and the array information
  for (const auto& particle: d_allSPHParticleVec) {

    //// not print the most right layer of SPH free particles, 2015/05/19
    // if((*pt)->getInitPosition().getx()==25){
    // continue;
    //}
    if (particle->getType() == SPHParticleType::BOUNDARY)
      continue; // not print out boundary particles

    ofs << std::setw(20) 
        << particle->currentPosition().x() << std::setw(20)
        << particle->currentPosition().y() << std::setw(20)
        << particle->currentPosition().z() << std::setw(20)
        << particle->getDisplacement().x() << std::setw(20)
        << particle->getDisplacement().y() << std::setw(20)
        << particle->getDisplacement().z() << std::setw(20)
        << particle->getVelocity().x() << std::setw(20) 
        << particle->getVelocity().y() << std::setw(20) 
        << particle->getVelocity().z();

    particle->calculatePressure();
    ofs << std::setw(20) << particle->getPressure();

    ofs << std::setw(20) 
        << particle->accelerationeration().x() << std::setw(20)
        << particle->accelerationeration().y() << std::setw(20)
        << particle->accelerationeration().z() << std::setw(20)
        << particle->densityRate() << std::setw(20)
        << particle->density() << std::endl;
  }

  ofs.close();
}


template void
SmoothParticleHydro::updateParticleInteractions<2>(const Box& spatialDomain,
  const REAL& bufferWidth, const REAL& ghostWidth, const REAL& kernelSize,
  const REAL& smoothLength);

template void
SmoothParticleHydro::updateParticleInteractions<3>(const Box& spatialDomain,
  const REAL& bufferWidth, const REAL& ghostWidth, const REAL& kernelSize,
  const REAL& smoothLength);

template void
SmoothParticleHydro::assignParticlesToPatchGrid<2>(const Box& domain,
  const REAL& bufferWidth, const REAL& ghostWidth, const REAL& kernelSize);

template void
SmoothParticleHydro::assignParticlesToPatchGrid<3>(const Box& domain,
  const REAL& bufferWidth, const REAL& ghostWidth, const REAL& kernelSize);

} // end namespace sph
