#ifndef DEM_PATCH_H
#define DEM_PATCH_H

#include <DiscreteElements/DEMContainers.h>
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

  template <typename TArray>
  struct PatchNeighborComm {
    PatchBoundary d_boundary;   // Whether the patch has a neighbor
    int d_rank;                 // Rank of the neighbor
    int d_mpiTag = 0;
    IntVec d_coords;
    boost::mpi::request d_sendRecvReq[2];
    TArray d_sentParticles; // For sends to neighbor
    TArray d_recvParticles; // For receives from neighbor

    void setNeighbor(MPI_Comm& cartComm, int myRank,
                     const IntVec& neighborCoords,
                     PatchBoundary boundaryFlag);

    void asyncSendRecv(boost::mpi::communicator& boostWorld,
                       int myRank, int iteration,
                       const Box& box,
                       const double& tolerance,
                       const TArray& particles);

    void findParticlesInBox(const Box& box,
                            const TArray& particles,
                            const double& tolerance,
                            TArray& inside);

    void waitToFinish(int myRank, int iteration);

    void combineSentParticles(int myRank, int iteration, 
                             ParticleIDHashMap& sent);

    void combineReceivedParticles(int myRank, int iteration, 
                                 TArray& received);

  };

  template <typename TArray>
  struct Patch {
    int d_rank;
    IntVec d_patchMPICoords;
    Vec d_lower;
    Vec d_upper;
    double d_ghostWidth;
    double d_tolerance;
    PatchNeighborComm<TArray> d_xMinus;
    PatchNeighborComm<TArray> d_yMinus;
    PatchNeighborComm<TArray> d_zMinus;
    PatchNeighborComm<TArray> d_xPlus;
    PatchNeighborComm<TArray> d_yPlus;
    PatchNeighborComm<TArray> d_zPlus;

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
                             const TArray& particles);

    void sendRecvGhostXPlus(boost::mpi::communicator& boostWorld,
                            int iteration,
                            const TArray& particles);

    void sendRecvGhostYMinus(boost::mpi::communicator& boostWorld, 
                             int iteration,
                             const TArray& particles);

    void sendRecvGhostYPlus(boost::mpi::communicator& boostWorld, 
                            int iteration,
                            const TArray& particles);

    void sendRecvGhostZMinus(boost::mpi::communicator& boostWorld, 
                             int iteration,
                             const TArray& particles);

    void sendRecvGhostZPlus(boost::mpi::communicator& boostWorld,
                            int iteration,
                            const TArray& particles);

    void sendRecvMigrateXMinus(boost::mpi::communicator& boostWorld,
                               int iteration, const Vec& neighborWidth,
                               const TArray& particles);

    void sendRecvMigrateXPlus(boost::mpi::communicator& boostWorld,
                              int iteration, const Vec& neighborWidth,
                              const TArray& particles);

    void sendRecvMigrateYMinus(boost::mpi::communicator& boostWorld, 
                               int iteration, const Vec& neighborWidth,
                               const TArray& particles);

    void sendRecvMigrateYPlus(boost::mpi::communicator& boostWorld, 
                              int iteration, const Vec& neighborWidth,
                              const TArray& particles);

    void sendRecvMigrateZMinus(boost::mpi::communicator& boostWorld, 
                               int iteration, const Vec& neighborWidth,
                               const TArray& particles);

    void sendRecvMigrateZPlus(boost::mpi::communicator& boostWorld,
                              int iteration, const Vec& neighborWidth,
                              const TArray& particles);

    void waitToFinishX(int iteration);

    void waitToFinishY(int iteration);

    void waitToFinishZ(int iteration);

    template <typename T>
    void combineReceivedParticles(int iteration,
                                  TArray& received);

    void combineSentParticlesX(int iteration, ParticleIDHashMap& sent);
    void combineSentParticlesY(int iteration, ParticleIDHashMap& sent);
    void combineSentParticlesZ(int iteration, ParticleIDHashMap& sent);
    void combineReceivedParticlesX(int iteration, TArray& received);
    void combineReceivedParticlesY(int iteration, TArray& received);
    void combineReceivedParticlesZ(int iteration, TArray& received);

    // Taken from: https://stackoverflow.com/questions/12200486/
    template <typename T>
    void removeDuplicates(TArray& input);

    template <typename T>
    void deleteSentParticles(int iteration, const ParticleIDHashMap& sent,
                             TArray& particles);

    void addReceivedParticles(int iteration, const TArray& received,
                              TArray& particles);

    template <typename T>
    void removeParticlesOutsidePatch(TArray& particles);

    void update(int iteration,
                const Vec& lower, const Vec& upper, const REAL& ghostWidth);

  };
} // end namespace dem

#endif