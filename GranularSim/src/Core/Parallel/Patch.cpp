#include <Core/Parallel/Patch.h>
#include <DiscreteElements/DEMParticle.h>
#include <Peridynamics/PeriContainers.h>
#include <Peridynamics/PeriParticle.h>
#include <SmoothParticleHydro/SPHContainers.h>
#include <SmoothParticleHydro/SPHParticle.h>
#include <Core/Types/IntegerTypes.h>
#include <set>
#include <ostream>


/* Testing error handling*/
#ifdef TEST_MPI_ERR_HANDLER
namespace dem 
{

static int errHandlerCalls = 0;
static int errHandlerErrs = 0;

static void patchCommErrorHandler(MPI_Comm* comm, int* err, ...)
{
  if (*err != MPI_ERR_OTHER) {
    errHandlerErrs++;
    std::cout << "**ERROR** Unexpected MPI error" << std::endl;
  }
  errHandlerCalls++;
}

} // end namespace dem
#endif

using namespace dem;

template<typename TArray>
void 
PatchNeighborComm<TArray>::setNeighbor(MPI_Comm& cartComm, int myRank,
                                      const IntVec& neighborCoords,
                                      PatchBoundary boundaryFlag) 
{
  int neighborRank = -1;

  #ifdef TEST_MPI_ERR_HANDLER
    MPI_Errhandler errHandler;
    MPI_Comm_create_errhandler(dem::patchCommErrorHandler, &errHandler);
    MPI_Comm_set_errhandler(cartComm, errHandler);
  #endif
  MPI_Comm_set_errhandler(cartComm, MPI_ERRORS_RETURN);
  int status = MPI_Cart_rank(cartComm, neighborCoords.data(), &neighborRank);
  if (status != MPI_SUCCESS) {
    //char error_string[1000];
    //int length_of_error_string;
    //MPI_Error_string(status, error_string, &length_of_error_string);
    //std::ostringstream out;
    //out << "rank = " << myRank << ":" << error_string;
    //std::cout << out.str() << "\n";
    d_boundary = boundaryFlag;
  } else {
    d_boundary = PatchBoundary::inside;
  }
  
  d_rank = neighborRank;
  d_coords = neighborCoords;

  //std::ostringstream out;
  //out << "myRank = " << myRank << " neighborCoords: " << neighborCoords 
  //    << " boundaryFlag " << static_cast<int>(d_boundary) << "\n";
  //std::cout << out.str();
}

template<typename TArray>
void 
PatchNeighborComm<TArray>::asyncSendRecv(boost::mpi::communicator& boostWorld,
                                         int myRank, int iteration,
                                         const Box& box, 
                                         const double& tolerance,
                                         const TArray& particles)
{
  //std::ostringstream out;
  //out << "myRank: " << myRank << " neighborRank: " << d_rank 
  //    << " iteration: " << iteration << " box: " << box.minCorner()
  //    << "," << box.maxCorner() 
  //    << " boundaryFlag " << static_cast<int>(d_boundary) << "\n";

  d_sentParticles.clear();
  d_recvParticles.clear();
  findParticlesInBox(box, particles, tolerance, d_sentParticles);
  d_sendRecvReq[0] = 
    boostWorld.isend(d_rank, d_mpiTag, d_sentParticles);
  d_sendRecvReq[1] = 
    boostWorld.irecv(d_rank, d_mpiTag, d_recvParticles);

  //out << " num sent: " << d_sentParticles.size()
  //    << " num recv: " << d_recvParticles.size() << "\n";
  //std::cout << out.str();
}

template<typename TArray>
void 
PatchNeighborComm<TArray>::findParticlesInBox(const Box& box,
                                              const TArray& particles,
                                              const double& tolerance,
                                              TArray& inside)
{
  for (const auto& particle : particles) {
    // it is critical to use EPS
    if (box.inside(particle->currentPosition(), tolerance)) {
      inside.push_back(particle);
    }
  }
}

template<typename TArray>
void 
PatchNeighborComm<TArray>::waitToFinish(int myRank, int iteration) 
{
  boost::mpi::wait_all(d_sendRecvReq, d_sendRecvReq + 2);

  /*
  if (myRank == 0) {
    std::ostringstream out;
    out << "myRank: " << myRank << " neighborRank: " << d_rank 
        << " neighborCoords: " << d_coords
        << " iteration: " << iteration 
        << " num recv: " << d_recvParticles.size() << "\n";
    for (auto& part : d_recvParticles) {
      out << part->getId() << ", ";
    }
    out << "\n";
    std::cout << out.str();
  }
  */
}

template<typename TArray>
void 
PatchNeighborComm<TArray>::combineSentParticles(int myRank, int iteration,
                                                ParticleIDHashMap& sent) 
{
  if (!d_sentParticles.empty()) {
    for (const auto& particle : d_sentParticles) {
      sent.insert(particle->getId());
    }
  }
}

template<typename TArray>
void 
PatchNeighborComm<TArray>::combineReceivedParticles(int myRank, int iteration,
                                                    TArray& received) 
{
  if (!d_recvParticles.empty()) {
    received.insert(received.end(), 
      d_recvParticles.begin(), d_recvParticles.end());
  }
}

template<typename TArray>
Patch<TArray>::Patch(MPI_Comm& cartComm, int rank, const IntVec& mpiCoords, 
                     const Vec& lower, const Vec& upper, double ghostWidth,
                     double tolerance)
{
  d_rank = rank;
  d_patchMPICoords = mpiCoords;
  d_lower = lower;
  d_upper = upper;
  d_ghostWidth = ghostWidth;
  d_tolerance = tolerance;
  setXMinus(cartComm);
  setXPlus(cartComm);
  setYMinus(cartComm);
  setYPlus(cartComm);
  setZMinus(cartComm);
  setZPlus(cartComm);
}

template<typename TArray>
void 
Patch<TArray>::setXMinus(MPI_Comm& cartComm)
{
  IntVec neighborCoords = d_patchMPICoords;
  --neighborCoords.x();
  //std::cout << "x-:" << neighborCoords << "\n";
  d_xMinus.setNeighbor(cartComm, d_rank, neighborCoords, PatchBoundary::xminus);
}

template<typename TArray>
void 
Patch<TArray>::setXPlus(MPI_Comm& cartComm)
{
  IntVec neighborCoords = d_patchMPICoords;
  ++neighborCoords.x();
  //std::cout << "x+:" << neighborCoords << "\n";
  d_xPlus.setNeighbor(cartComm, d_rank, neighborCoords, PatchBoundary::xplus);
}

template<typename TArray>
void 
Patch<TArray>::setYMinus(MPI_Comm& cartComm)
{
  IntVec neighborCoords = d_patchMPICoords;
  --neighborCoords.y();
  //std::cout << "y-:" << neighborCoords << "\n";
  d_yMinus.setNeighbor(cartComm, d_rank, neighborCoords, PatchBoundary::yminus);
}

template<typename TArray>
void 
Patch<TArray>::setYPlus(MPI_Comm& cartComm)
{
  IntVec neighborCoords = d_patchMPICoords;
  ++neighborCoords.y();
  //std::cout << "y+:" << neighborCoords << "\n";
  d_yPlus.setNeighbor(cartComm, d_rank, neighborCoords, PatchBoundary::yplus);
}

template<typename TArray>
void 
Patch<TArray>::setZMinus(MPI_Comm& cartComm)
{
  IntVec neighborCoords = d_patchMPICoords;
  --neighborCoords.z();
  //std::cout << "z-:" << neighborCoords << "\n";
  d_zMinus.setNeighbor(cartComm, d_rank, neighborCoords, PatchBoundary::zminus);
}

template<typename TArray>
void 
Patch<TArray>::setZPlus(MPI_Comm& cartComm)
{
  IntVec neighborCoords = d_patchMPICoords;
  ++neighborCoords.z();
  //std::cout << "z+:" << neighborCoords << "\n";
  d_zPlus.setNeighbor(cartComm, d_rank, neighborCoords, PatchBoundary::zplus);
}

template<typename TArray>
void 
Patch<TArray>::sendRecvGhostXMinus(boost::mpi::communicator& boostWorld, 
                                   int iteration,
                                   const TArray& particles) 
{
  if (d_xMinus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower;
    Vec ghostUpper = d_upper;
    ghostUpper.setX(d_lower.x() + d_ghostWidth);
    Box ghostBox(ghostLower, ghostUpper);
    d_xMinus.asyncSendRecv(boostWorld, d_rank, iteration, ghostBox, d_tolerance,
                           particles);
  }
}

template<typename TArray>
void 
Patch<TArray>::sendRecvGhostXPlus(boost::mpi::communicator& boostWorld, 
                                  int iteration,
                                  const TArray& particles) 
{
  if (d_xPlus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower;
    ghostLower.setX(d_upper.x() - d_ghostWidth);
    Vec ghostUpper = d_upper;
    Box ghostBox(ghostLower, ghostUpper);
    d_xPlus.asyncSendRecv(boostWorld, d_rank, iteration, ghostBox, d_tolerance,
                           particles);
  }
}

template<typename TArray>
void 
Patch<TArray>::sendRecvGhostYMinus(boost::mpi::communicator& boostWorld, 
                                   int iteration,
                                   const TArray& particles) 
{
  if (d_yMinus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower - Vec(d_ghostWidth, 0.0, 0.0);
    Vec ghostUpper = d_upper + Vec(d_ghostWidth, 0.0, 0.0);
    ghostUpper.setY(d_lower.y() + d_ghostWidth);
    Box ghostBox(ghostLower, ghostUpper);
    d_yMinus.asyncSendRecv(boostWorld, d_rank, iteration, ghostBox, d_tolerance,
                           particles);
  }
}

template<typename TArray>
void 
Patch<TArray>::sendRecvGhostYPlus(boost::mpi::communicator& boostWorld, 
                                  int iteration,
                                  const TArray& particles) 
{
  if (d_yPlus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower - Vec(d_ghostWidth, 0.0, 0.0);
    ghostLower.setY(d_upper.y() - d_ghostWidth);
    Vec ghostUpper = d_upper + Vec(d_ghostWidth, 0.0, 0.0);
    Box ghostBox(ghostLower, ghostUpper);
    d_yPlus.asyncSendRecv(boostWorld, d_rank, iteration, ghostBox, d_tolerance,
                           particles);
  }
}

template<typename TArray>
void 
Patch<TArray>::sendRecvGhostZMinus(boost::mpi::communicator& boostWorld, 
                           int iteration,
                           const TArray& particles) 
{
  if (d_zMinus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower - Vec(d_ghostWidth, d_ghostWidth, 0.0);
    Vec ghostUpper = d_upper + Vec(d_ghostWidth, d_ghostWidth, 0.0);
    ghostUpper.setZ(d_lower.z() + d_ghostWidth);
    Box ghostBox(ghostLower, ghostUpper);
    d_zMinus.asyncSendRecv(boostWorld, d_rank, iteration, ghostBox, d_tolerance,
                           particles);
  }
}

template<typename TArray>
void 
Patch<TArray>::sendRecvGhostZPlus(boost::mpi::communicator& boostWorld, 
                          int iteration,
                          const TArray& particles) 
{
  if (d_zPlus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower - Vec(d_ghostWidth, d_ghostWidth, 0.0);
    ghostLower.setZ(d_upper.z() - d_ghostWidth);
    Vec ghostUpper = d_upper + Vec(d_ghostWidth, d_ghostWidth, 0.0);
    Box ghostBox(ghostLower, ghostUpper);
    d_zPlus.asyncSendRecv(boostWorld, d_rank, iteration, ghostBox, d_tolerance,
                           particles);
  }
}

template<typename TArray>
void 
Patch<TArray>::sendRecvMigrateXMinus(boost::mpi::communicator& boostWorld, 
                             int iteration, const Vec& neighborWidth,
                             const TArray& particles) 
{
  if (d_xMinus.d_boundary == PatchBoundary::inside) {
    Vec neighborLower = d_lower - neighborWidth;
    Vec neighborUpper = d_upper + neighborWidth;
    neighborUpper.setX(d_lower.x());
    Box neighborBox(neighborLower, neighborUpper);
    d_xMinus.asyncSendRecv(boostWorld, d_rank, iteration, 
                           neighborBox, d_tolerance,
                           particles);
  }
}

template<typename TArray>
void 
Patch<TArray>::sendRecvMigrateXPlus(boost::mpi::communicator& boostWorld, 
                            int iteration, const Vec& neighborWidth,
                            const TArray& particles) 
{
  if (d_xPlus.d_boundary == PatchBoundary::inside) {
    Vec neighborLower = d_lower - neighborWidth;
    neighborLower.setX(d_upper.x());
    Vec neighborUpper = d_upper + neighborWidth;
    Box neighborBox(neighborLower, neighborUpper);
    d_xPlus.asyncSendRecv(boostWorld, d_rank, iteration, 
                          neighborBox, d_tolerance,
                          particles);
  }
}

template<typename TArray>
void 
Patch<TArray>::sendRecvMigrateYMinus(boost::mpi::communicator& boostWorld, 
                             int iteration, const Vec& neighborWidth,
                             const TArray& particles) 
{
  if (d_yMinus.d_boundary == PatchBoundary::inside) {
    Vec neighborLower = d_lower - neighborWidth;
    Vec neighborUpper = d_upper + neighborWidth;
    neighborUpper.setY(d_lower.y());
    Box neighborBox(neighborLower, neighborUpper);
    d_yMinus.asyncSendRecv(boostWorld, d_rank, iteration, 
                          neighborBox, d_tolerance,
                          particles);
  }
}

template<typename TArray>
void 
Patch<TArray>::sendRecvMigrateYPlus(boost::mpi::communicator& boostWorld, 
                            int iteration, const Vec& neighborWidth,
                            const TArray& particles) 
{
  if (d_yPlus.d_boundary == PatchBoundary::inside) {
    Vec neighborLower = d_lower - neighborWidth;
    neighborLower.setY(d_upper.y());
    Vec neighborUpper = d_upper + neighborWidth;
    Box neighborBox(neighborLower, neighborUpper);
    d_yPlus.asyncSendRecv(boostWorld, d_rank, iteration, 
                          neighborBox, d_tolerance,
                          particles);
  }
}

template<typename TArray>
void 
Patch<TArray>::sendRecvMigrateZMinus(boost::mpi::communicator& boostWorld, 
                             int iteration, const Vec& neighborWidth,
                             const TArray& particles) 
{
  if (d_zMinus.d_boundary == PatchBoundary::inside) {
    Vec neighborLower = d_lower - neighborWidth;
    Vec neighborUpper = d_upper + neighborWidth;
    neighborUpper.setZ(d_lower.z());
    Box neighborBox(neighborLower, neighborUpper);
    d_zMinus.asyncSendRecv(boostWorld, d_rank, iteration, 
                          neighborBox, d_tolerance,
                          particles);
  }
}

template<typename TArray>
void 
Patch<TArray>::sendRecvMigrateZPlus(boost::mpi::communicator& boostWorld, 
                            int iteration, const Vec& neighborWidth,
                            const TArray& particles) 
{
  if (d_zPlus.d_boundary == PatchBoundary::inside) {
    Vec neighborLower = d_lower - neighborWidth;
    neighborLower.setZ(d_upper.z());
    Vec neighborUpper = d_upper + neighborWidth;
    Box neighborBox(neighborLower, neighborUpper);
    d_zPlus.asyncSendRecv(boostWorld, d_rank, iteration, 
                          neighborBox, d_tolerance,
                          particles);
  }
}

template<typename TArray>
void 
Patch<TArray>::waitToFinishX(int iteration) 
{
  if (d_xMinus.d_boundary == PatchBoundary::inside) {
    d_xMinus.waitToFinish(d_rank, iteration);
  }
  if (d_xPlus.d_boundary == PatchBoundary::inside) {
    d_xPlus.waitToFinish(d_rank, iteration);
  }
}

template<typename TArray>
void 
Patch<TArray>::waitToFinishY(int iteration) 
{
  if (d_yMinus.d_boundary == PatchBoundary::inside) {
    d_yMinus.waitToFinish(d_rank, iteration);
  }
  if (d_yPlus.d_boundary == PatchBoundary::inside) {
    d_yPlus.waitToFinish(d_rank, iteration);
  }
}

template<typename TArray>
void 
Patch<TArray>::waitToFinishZ(int iteration) 
{
  if (d_zMinus.d_boundary == PatchBoundary::inside) {
    d_zMinus.waitToFinish(d_rank, iteration);
  }
  if (d_zPlus.d_boundary == PatchBoundary::inside) {
    d_zPlus.waitToFinish(d_rank, iteration);
  }
}

template<typename TArray>
template<typename T>
void 
Patch<TArray>::combineReceivedParticles(int iteration,
                               TArray& received) 
{
  received.clear();
  d_xMinus.combineReceivedParticles(d_rank, iteration, received);
  d_xPlus.combineReceivedParticles(d_rank, iteration, received);
  d_yMinus.combineReceivedParticles(d_rank, iteration, received);
  d_yPlus.combineReceivedParticles(d_rank, iteration, received);
  d_zMinus.combineReceivedParticles(d_rank, iteration, received);
  d_zPlus.combineReceivedParticles(d_rank, iteration, received);
  removeDuplicates<T>(received);
}
 
template<typename TArray>
void 
Patch<TArray>::combineSentParticlesX(int iteration,
                             ParticleIDHashMap& sent) 
{
  d_xMinus.combineSentParticles(d_rank, iteration, sent);
  d_xPlus.combineSentParticles(d_rank, iteration, sent);
}

template<typename TArray>
void 
Patch<TArray>::combineSentParticlesY(int iteration,
                             ParticleIDHashMap& sent) 
{
  d_yMinus.combineSentParticles(d_rank, iteration, sent);
  d_yPlus.combineSentParticles(d_rank, iteration, sent);
}

template<typename TArray>
void 
Patch<TArray>::combineSentParticlesZ(int iteration,
                             ParticleIDHashMap& sent) 
{
  d_zMinus.combineSentParticles(d_rank, iteration, sent);
  d_zPlus.combineSentParticles(d_rank, iteration, sent);
}

template<typename TArray>
void 
Patch<TArray>::combineReceivedParticlesX(int iteration,
                                TArray& received) 
{
  d_xMinus.combineReceivedParticles(d_rank, iteration, received);
  d_xPlus.combineReceivedParticles(d_rank, iteration, received);
}

template<typename TArray>
void 
Patch<TArray>::combineReceivedParticlesY(int iteration,
                                TArray& received) 
{
  d_yMinus.combineReceivedParticles(d_rank, iteration, received);
  d_yPlus.combineReceivedParticles(d_rank, iteration, received);
}

template<typename TArray>
void 
Patch<TArray>::combineReceivedParticlesZ(int iteration,
                                TArray& received) 
{
  d_zMinus.combineReceivedParticles(d_rank, iteration, received);
  d_zPlus.combineReceivedParticles(d_rank, iteration, received);
}

// Taken from: https://stackoverflow.com/questions/12200486/
template<typename TArray>
template<typename T>
void 
Patch<TArray>::removeDuplicates(TArray& input) 
{
  //std::ostringstream out;
  //out << "Removing duplicates: in = " << input.size();

  std::set<ParticleID> seen;
  auto newEnd = 
    std::remove_if(input.begin(), input.end(),
                    [&seen](const T& particle)
                    {
                      if (seen.find(particle->getId()) != std::end(seen)) {
                        return true;
                      }
                      seen.insert(particle->getId());
                      return false;
                    });
  input.erase(newEnd, input.end());

  //out << " out = " << input.size() << "\n";
  //std::cout << out.str();
}

template<typename TArray>
template<typename T>
void
Patch<TArray>::deleteSentParticles(int iteration, 
                           const ParticleIDHashMap& sent,
                           TArray& particles)
{
  //std::ostringstream out;
  //out << "Deleting sent particle ( " << sent.size() 
  //    << ") : in = " << particles.size();

  if (sent.size() > 0) {
    particles.erase(
      std::remove_if(
        particles.begin(), particles.end(),
        [&sent](const T& particle) {
          return (sent.find(particle->getId()) != std::end(sent));
        }),
      std::end(particles));
  }

  //out << " out = " << particles.size() << "\n";
  //std::cout << out.str();
}

template<typename TArray>
void
Patch<TArray>::addReceivedParticles(int iteration, 
                            const TArray& received,
                            TArray& particles)
{
  particles.insert(particles.end(), received.begin(), received.end());
}

template<typename TArray>
template<typename T>
void 
Patch<TArray>::removeParticlesOutsidePatch(TArray& particles)
{
  //std::ostringstream out;
  //out << "Removing outside particles: in = " << particles.size();

  Box box(d_lower, d_upper);
  double epsilon = d_tolerance;
  particles.erase(
    std::remove_if(
      particles.begin(), particles.end(),
      [&box, &epsilon](T& particle) {
        if (box.inside(particle->currentPosition(), epsilon)) {
          return false;
        }
        return true;
      }),
    particles.end());

  //out << " out = " << particles.size() << "\n";
  //std::cout << out.str();

}

template<typename TArray>
void 
Patch<TArray>::update(int iteration,
              const Vec& lower, const Vec& upper, const REAL& ghostWidth)
{
  d_lower = lower;
  d_upper = upper;
  d_ghostWidth = ghostWidth;
}

// Instantiation
namespace dem {
  template struct PatchNeighborComm<DEMParticlePArray>;
  template struct Patch<DEMParticlePArray>;
  template void 
  Patch<DEMParticlePArray>::removeDuplicates<DEMParticleP>(DEMParticlePArray& input);
  template void 
  Patch<DEMParticlePArray>::deleteSentParticles<DEMParticleP>(int iteration, 
    const ParticleIDHashMap& sent, DEMParticlePArray& particles);
  template void 
  Patch<DEMParticlePArray>::removeParticlesOutsidePatch<DEMParticleP>(DEMParticlePArray& particles);

  template struct PatchNeighborComm<pd::PeriParticlePArray>;
  template struct Patch<pd::PeriParticlePArray>;
  template void 
  Patch<pd::PeriParticlePArray>::removeDuplicates<pd::PeriParticleP>(pd::PeriParticlePArray& input);
  template void 
  Patch<pd::PeriParticlePArray>::deleteSentParticles<pd::PeriParticleP>(int iteration, 
    const ParticleIDHashMap& sent, pd::PeriParticlePArray& particles);
  template void 
  Patch<pd::PeriParticlePArray>::removeParticlesOutsidePatch<pd::PeriParticleP>(pd::PeriParticlePArray& particles);

  template struct PatchNeighborComm<sph::SPHParticlePArray>;
  template struct Patch<sph::SPHParticlePArray>;
  template void 
  Patch<sph::SPHParticlePArray>::removeDuplicates<sph::SPHParticleP>(sph::SPHParticlePArray& input);
  template void 
  Patch<sph::SPHParticlePArray>::deleteSentParticles<sph::SPHParticleP>(int iteration, 
    const ParticleIDHashMap& sent, sph::SPHParticlePArray& particles);
  template void 
  Patch<sph::SPHParticlePArray>::removeParticlesOutsidePatch<sph::SPHParticleP>(sph::SPHParticlePArray& particles);
}

