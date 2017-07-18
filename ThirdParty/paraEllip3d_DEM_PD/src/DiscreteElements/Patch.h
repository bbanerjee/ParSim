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
    int d_mpiTag;
    boost::mpi::request d_sendRecvReq[2];
    ParticlePArray d_sendGhostParticles; // For sends to neighbor
    ParticlePArray d_recvGhostParticles; // For receives from neighbor

    void setNeighbor(MPI_Comm& cartComm, const IntVec& neighborCoords,
                     PatchBoundary boundaryFlag);

    void asynchronousSendRecv(boost::mpi::communicator& boostWorld,
                              int iteration,
                              const ParticlePArray& particles,
                              const Box& ghostBox);

    void findParticlesInGhostBox(const Box& ghostBox,
                                 const ParticlePArray& particles,
                                 ParticlePArray& ghostParticles);

    void waitToFinish(int iteration);

    void insertReceivedParticles(int iteration, ParticlePArray& received);

  };

  struct Patch {
    int d_rank;
    IntVec d_patchMPICoords;
    Vec d_lower;
    Vec d_upper;
    double d_ghostWidth;
    PatchNeighborComm d_xMinus;
    PatchNeighborComm d_yMinus;
    PatchNeighborComm d_zMinus;
    PatchNeighborComm d_xPlus;
    PatchNeighborComm d_yPlus;
    PatchNeighborComm d_zPlus;

    Patch(MPI_Comm& cartComm, 
          int rank, const IntVec& mpiCoords, const Vec& lower, const Vec& upper,
          double ghostWidth);

    void setXMinus(MPI_Comm& cartComm);

    void setXPlus(MPI_Comm& cartComm);

    void setYMinus(MPI_Comm& cartComm);

    void setYPlus(MPI_Comm& cartComm);

    void setZMinus(MPI_Comm& cartComm);

    void setZPlus(MPI_Comm& cartComm);

    void sendRecvXMinus(boost::mpi::communicator& boostWorld,
                        int iteration,
                        const ParticlePArray& particles);

    void sendRecvXPlus(boost::mpi::communicator& boostWorld,
                       int iteration,
                       const ParticlePArray& particles);

    void sendRecvYMinus(boost::mpi::communicator& boostWorld, 
                        int iteration,
                        const ParticlePArray& particles);

    void sendRecvYPlus(boost::mpi::communicator& boostWorld, 
                       int iteration,
                       const ParticlePArray& particles);

    void sendRecvZMinus(boost::mpi::communicator& boostWorld, 
                        int iteration,
                        const ParticlePArray& particles);

    void sendRecvZPlus(boost::mpi::communicator& boostWorld,
                       int iteration,
                       const ParticlePArray& particles);

    void waitToFinishX(int iteration);

    void waitToFinishY(int iteration);

    void waitToFinishZ(int iteration);

    void insertReceivedParticles(int iteration,
                                 ParticlePArray& received);

    // Taken from: https://stackoverflow.com/questions/12200486/
    void removeDuplicates(ParticlePArray& input);

    void update(int iteration,
                const Vec& lower, const Vec& upper, const REAL& ghostWidth);

  };
} // end namespace dem

#endif