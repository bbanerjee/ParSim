#ifndef DEM_PATCH_H
#define DEM_PATCH_H

#include <DiscreteElements/Containers.h>
#include <Core/Math/IntVec.h>
#include <Core/Math/Vec.h>
#include <Core/Geometry/Box.h>
#include <boost/mpi.hpp>
#include <boost/serialization/shared_ptr.hpp>

namespace dem {

  enum class PatchBoundary : char {
    xminus,
    xplus,
    yminus,
    yplus,
    zminus,
    zplus,
    inside
  };

  struct PatchNeighborComm {
    PatchBoundary d_boundary;   // Whether the patch has a neighbor
    int d_rank;                 // Rank of the neighbor
    int d_mpiTag = 0;
    IntVec d_coords;
    boost::mpi::request d_sendRecvReq[2];
    ParticlePArray d_sentParticles; // For sends to neighbor
    ParticlePArray d_recvParticles; // For receives from neighbor

    void setNeighbor(MPI_Comm& cartComm, int myRank,
                     const IntVec& neighborCoords,
                     PatchBoundary boundaryFlag);

    void asyncSendRecv(boost::mpi::communicator& boostWorld,
                       int myRank, int iteration,
                       const Box& box,
                       const double& tolerance,
                       const ParticlePArray& particles);

    void findParticlesInBox(const Box& box,
                            const ParticlePArray& particles,
                            const double& tolerance,
                            ParticlePArray& inside);

    void waitToFinish(int myRank, int iteration);

    void combineSentParticles(int myRank, int iteration, 
                             ParticleIDHashMap& sent);

    void combineReceivedParticles(int myRank, int iteration, 
                                 ParticlePArray& received);

  };

  struct Patch {
    int d_rank;
    IntVec d_patchMPICoords;
    Vec d_lower;
    Vec d_upper;
    double d_ghostWidth;
    double d_tolerance;
    PatchNeighborComm d_xMinus;
    PatchNeighborComm d_yMinus;
    PatchNeighborComm d_zMinus;
    PatchNeighborComm d_xPlus;
    PatchNeighborComm d_yPlus;
    PatchNeighborComm d_zPlus;

    Patch(MPI_Comm& cartComm, 
          int rank, const IntVec& mpiCoords, const Vec& lower, const Vec& upper,
          double ghostWidth, double tolerance);

    void setXMinus(MPI_Comm& cartComm);

    void setXPlus(MPI_Comm& cartComm);

    void setYMinus(MPI_Comm& cartComm);

    void setYPlus(MPI_Comm& cartComm);

    void setZMinus(MPI_Comm& cartComm);

    void setZPlus(MPI_Comm& cartComm);

    void sendRecvGhostXMinus(boost::mpi::communicator& boostWorld,
                             int iteration,
                             const ParticlePArray& particles);

    void sendRecvGhostXPlus(boost::mpi::communicator& boostWorld,
                            int iteration,
                            const ParticlePArray& particles);

    void sendRecvGhostYMinus(boost::mpi::communicator& boostWorld, 
                             int iteration,
                             const ParticlePArray& particles);

    void sendRecvGhostYPlus(boost::mpi::communicator& boostWorld, 
                            int iteration,
                            const ParticlePArray& particles);

    void sendRecvGhostZMinus(boost::mpi::communicator& boostWorld, 
                             int iteration,
                             const ParticlePArray& particles);

    void sendRecvGhostZPlus(boost::mpi::communicator& boostWorld,
                            int iteration,
                            const ParticlePArray& particles);

    void sendRecvMigrateXMinus(boost::mpi::communicator& boostWorld,
                               int iteration, const Vec& neighborWidth,
                               const ParticlePArray& particles);

    void sendRecvMigrateXPlus(boost::mpi::communicator& boostWorld,
                              int iteration, const Vec& neighborWidth,
                              const ParticlePArray& particles);

    void sendRecvMigrateYMinus(boost::mpi::communicator& boostWorld, 
                               int iteration, const Vec& neighborWidth,
                               const ParticlePArray& particles);

    void sendRecvMigrateYPlus(boost::mpi::communicator& boostWorld, 
                              int iteration, const Vec& neighborWidth,
                              const ParticlePArray& particles);

    void sendRecvMigrateZMinus(boost::mpi::communicator& boostWorld, 
                               int iteration, const Vec& neighborWidth,
                               const ParticlePArray& particles);

    void sendRecvMigrateZPlus(boost::mpi::communicator& boostWorld,
                              int iteration, const Vec& neighborWidth,
                              const ParticlePArray& particles);

    void waitToFinishX(int iteration);

    void waitToFinishY(int iteration);

    void waitToFinishZ(int iteration);

    void combineReceivedParticles(int iteration,
                                 ParticlePArray& received);

    void combineSentParticlesX(int iteration, ParticleIDHashMap& sent);
    void combineSentParticlesY(int iteration, ParticleIDHashMap& sent);
    void combineSentParticlesZ(int iteration, ParticleIDHashMap& sent);
    void combineReceivedParticlesX(int iteration, ParticlePArray& received);
    void combineReceivedParticlesY(int iteration, ParticlePArray& received);
    void combineReceivedParticlesZ(int iteration, ParticlePArray& received);

    // Taken from: https://stackoverflow.com/questions/12200486/
    void removeDuplicates(ParticlePArray& input);

    void deleteSentParticles(int iteration, const ParticleIDHashMap& sent,
                             ParticlePArray& particles);
    void addReceivedParticles(int iteration, const ParticlePArray& received,
                              ParticlePArray& particles);

    void removeParticlesOutsidePatch(ParticlePArray& particles);

    void update(int iteration,
                const Vec& lower, const Vec& upper, const REAL& ghostWidth);

  };
} // end namespace dem

#endif