#include <DiscreteElements/Patch.h>
#include <DiscreteElements/Particle.h>
#include <set>

using namespace dem;

void 
PatchNeighborComm::setNeighbor(MPI_Comm& cartComm, int myRank,
                               const IntVec& neighborCoords,
                               PatchBoundary boundaryFlag) 
{
  int neighborRank = -1;
  MPI_Cart_rank(cartComm, neighborCoords.data(), &neighborRank);
  d_rank = neighborRank;
  if (neighborRank > -1) {
    d_boundary = PatchBoundary::inside;
  } else {
    d_boundary = boundaryFlag;
  }

  std::ostringstream out;
  out << "myRank = " << myRank << " neighborCoords: " << neighborCoords 
      << " boundaryFlag " << static_cast<int>(boundaryFlag) << "\n";
  std::cout << out.str();
}

void 
PatchNeighborComm::asynchronousSendRecv(boost::mpi::communicator& boostWorld,
                                        int myRank,
                                        int iteration,
                                        const ParticlePArray& particles,
                                        const Box& ghostBox) 
{
  findParticlesInBox(ghostBox, particles, d_sendParticles);
  d_sendRecvReq[0] = 
    boostWorld.isend(d_rank, d_mpiTag, d_sendParticles);
  d_sendRecvReq[1] = 
    boostWorld.irecv(d_rank, d_mpiTag, d_recvParticles);

  std::ostringstream out;
  out << "myRank: " << myRank << " neighborRank: " << d_rank 
      << " iteration: " << iteration << " ghostbox: " << ghostBox
      << " num sent: " << d_sendParticles.size()
      << " num recv: " << d_recvParticles.size() << "\n";
  std::cout << out.str();
}

void 
PatchNeighborComm::findParticlesInBox(const Box& box,
                                      const ParticlePArray& particles,
                                      ParticlePArray& inside)
{
  inside.clear();
  for (const auto& particle : particles) {
    // it is critical to use EPS
    if (box.inside(particle->currentPos(), EPS)) {
      inside.push_back(particle);
    }
  }
}

void 
PatchNeighborComm::waitToFinish(int myRank, int iteration) 
{
  boost::mpi::wait_all(d_sendRecvReq, d_sendRecvReq + 2);

  std::ostringstream out;
  out << "myRank: " << myRank << " neighborRank: " << d_rank 
      << " iteration: " << iteration 
      << " num recv: " << d_recvParticles.size() << "\n";
  std::cout << out.str();
}

void 
PatchNeighborComm::insertReceivedParticles(int myRank, int iteration,
                                           ParticlePArray& received) 
{
  if (!d_recvParticles.empty()) {
    received.insert(received.end(), 
      d_recvParticles.begin(), d_recvParticles.end());
  }
}


Patch::Patch(MPI_Comm& cartComm, int rank, const IntVec& mpiCoords, 
             const Vec& lower, const Vec& upper, double ghostWidth)
{
  d_rank = rank;
  d_patchMPICoords = mpiCoords;
  d_lower = lower;
  d_upper = upper;
  d_ghostWidth = ghostWidth;
  setXMinus(cartComm);
  setXPlus(cartComm);
  setYMinus(cartComm);
  setYPlus(cartComm);
  setZMinus(cartComm);
  setZPlus(cartComm);
}

void 
Patch::setXMinus(MPI_Comm& cartComm)
{
  IntVec neighborCoords = d_patchMPICoords;
  --neighborCoords.x();
  d_xMinus.setNeighbor(cartComm, d_rank, neighborCoords, PatchBoundary::xminus);
}

void 
Patch::setXPlus(MPI_Comm& cartComm)
{
  IntVec neighborCoords = d_patchMPICoords;
  ++neighborCoords.x();
  d_xPlus.setNeighbor(cartComm, d_rank, neighborCoords, PatchBoundary::xplus);
}

void 
Patch::setYMinus(MPI_Comm& cartComm)
{
  IntVec neighborCoords = d_patchMPICoords;
  --neighborCoords.y();
  d_yMinus.setNeighbor(cartComm, d_rank, neighborCoords, PatchBoundary::yminus);
}

void 
Patch::setYPlus(MPI_Comm& cartComm)
{
  IntVec neighborCoords = d_patchMPICoords;
  ++neighborCoords.y();
  d_yPlus.setNeighbor(cartComm, d_rank, neighborCoords, PatchBoundary::yplus);
}

void 
Patch::setZMinus(MPI_Comm& cartComm)
{
  IntVec neighborCoords = d_patchMPICoords;
  --neighborCoords.z();
  d_zMinus.setNeighbor(cartComm, d_rank, neighborCoords, PatchBoundary::zminus);
}

void 
Patch::setZPlus(MPI_Comm& cartComm)
{
  IntVec neighborCoords = d_patchMPICoords;
  ++neighborCoords.z();
  d_zPlus.setNeighbor(cartComm, d_rank, neighborCoords, PatchBoundary::zplus);
}

void 
Patch::sendRecvGhostXMinus(boost::mpi::communicator& boostWorld, 
                           int iteration,
                           const ParticlePArray& particles) 
{
  if (d_xMinus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower - Vec(0.0, d_ghostWidth, d_ghostWidth);
    Vec ghostUpper = d_upper + Vec(0.0, d_ghostWidth, d_ghostWidth);
    ghostUpper.setX(d_lower.x() + d_ghostWidth);
    Box ghostBox(ghostLower, ghostUpper);
    d_xMinus.asynchronousSendRecv(boostWorld, d_rank, iteration, particles, ghostBox);
  }
}

void 
Patch::sendRecvGhostXPlus(boost::mpi::communicator& boostWorld, 
                          int iteration,
                          const ParticlePArray& particles) 
{
  if (d_xPlus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower - Vec(0.0, d_ghostWidth, d_ghostWidth);
    ghostLower.setX(d_upper.x() - d_ghostWidth);
    Vec ghostUpper = d_upper + Vec(0.0, d_ghostWidth, d_ghostWidth);
    Box ghostBox(ghostLower, ghostUpper);
    d_xPlus.asynchronousSendRecv(boostWorld, d_rank, iteration, particles, ghostBox);
  }
}

void 
Patch::sendRecvGhostYMinus(boost::mpi::communicator& boostWorld, 
                           int iteration,
                           const ParticlePArray& particles) 
{
  if (d_yMinus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower - Vec(d_ghostWidth, 0.0, d_ghostWidth);
    Vec ghostUpper = d_upper + Vec(d_ghostWidth, 0.0, d_ghostWidth);
    ghostUpper.setY(d_lower.y() + d_ghostWidth);
    Box ghostBox(ghostLower, ghostUpper);
    d_yMinus.asynchronousSendRecv(boostWorld, d_rank, iteration, particles, ghostBox);
  }
}

void 
Patch::sendRecvGhostYPlus(boost::mpi::communicator& boostWorld, 
                          int iteration,
                          const ParticlePArray& particles) 
{
  if (d_yPlus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower - Vec(d_ghostWidth, 0.0, d_ghostWidth);
    ghostLower.setY(d_upper.y() - d_ghostWidth);
    Vec ghostUpper = d_upper + Vec(d_ghostWidth, 0.0, d_ghostWidth);
    Box ghostBox(ghostLower, ghostUpper);
    d_yPlus.asynchronousSendRecv(boostWorld, d_rank, iteration, particles, ghostBox);
  }
}

void 
Patch::sendRecvGhostZMinus(boost::mpi::communicator& boostWorld, 
                           int iteration,
                           const ParticlePArray& particles) 
{
  if (d_zMinus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower - Vec(d_ghostWidth, d_ghostWidth, 0.0);
    Vec ghostUpper = d_upper + Vec(d_ghostWidth, d_ghostWidth, 0.0);
    ghostUpper.setZ(d_lower.z() + d_ghostWidth);
    Box ghostBox(ghostLower, ghostUpper);
    d_zMinus.asynchronousSendRecv(boostWorld, d_rank, iteration, particles, ghostBox);
  }
}

void 
Patch::sendRecvGhostZPlus(boost::mpi::communicator& boostWorld, 
                          int iteration,
                          const ParticlePArray& particles) 
{
  if (d_zPlus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower - Vec(d_ghostWidth, d_ghostWidth, 0.0);
    ghostLower.setZ(d_upper.z() - d_ghostWidth);
    Vec ghostUpper = d_upper + Vec(d_ghostWidth, d_ghostWidth, 0.0);
    Box ghostBox(ghostLower, ghostUpper);
    d_zPlus.asynchronousSendRecv(boostWorld, d_rank, iteration, particles, ghostBox);
  }
}

void 
Patch::sendRecvMigrateXMinus(boost::mpi::communicator& boostWorld, 
                             int iteration, const REAL& neighborWidth,
                             const ParticlePArray& particles) 
{
  if (d_xMinus.d_boundary == PatchBoundary::inside) {
    Vec neighborLower = d_lower - Vec(neighborWidth, 0.0, 0.0);
    Vec neighborUpper = d_upper;
    neighborUpper.setX(d_lower.x());
    Box neighborBox(neighborLower, neighborUpper);
    d_xMinus.asynchronousSendRecv(boostWorld, d_rank, iteration, 
                                  particles, neighborBox);
  }
}

void 
Patch::sendRecvMigrateXPlus(boost::mpi::communicator& boostWorld, 
                            int iteration, const REAL& neighborWidth,
                            const ParticlePArray& particles) 
{
  if (d_xPlus.d_boundary == PatchBoundary::inside) {
    Vec neighborLower = d_lower;
    neighborLower.setX(d_upper.x());
    Vec neighborUpper = d_upper + Vec(neighborWidth, 0.0, 0.0);
    Box neighborBox(neighborLower, neighborUpper);
    d_yPlus.asynchronousSendRecv(boostWorld, d_rank, iteration, 
                                  particles, neighborBox);
  }
}

void 
Patch::sendRecvMigrateYMinus(boost::mpi::communicator& boostWorld, 
                             int iteration, const REAL& neighborWidth,
                             const ParticlePArray& particles) 
{
  if (d_yMinus.d_boundary == PatchBoundary::inside) {
    Vec neighborLower = d_lower - Vec(0.0, neighborWidth, 0.0);
    Vec neighborUpper = d_upper;
    neighborUpper.setY(d_lower.y());
    Box neighborBox(neighborLower, neighborUpper);
    d_yMinus.asynchronousSendRecv(boostWorld, d_rank, iteration, 
                                  particles, neighborBox);
  }
}

void 
Patch::sendRecvMigrateYPlus(boost::mpi::communicator& boostWorld, 
                            int iteration, const REAL& neighborWidth,
                            const ParticlePArray& particles) 
{
  if (d_yPlus.d_boundary == PatchBoundary::inside) {
    Vec neighborLower = d_lower;
    neighborLower.setY(d_upper.y());
    Vec neighborUpper = d_upper + Vec(0.0, neighborWidth, 0.0);
    Box neighborBox(neighborLower, neighborUpper);
    d_yPlus.asynchronousSendRecv(boostWorld, d_rank, iteration, 
                                  particles, neighborBox);
  }
}

void 
Patch::sendRecvMigrateZMinus(boost::mpi::communicator& boostWorld, 
                             int iteration, const REAL& neighborWidth,
                             const ParticlePArray& particles) 
{
  if (d_zMinus.d_boundary == PatchBoundary::inside) {
    Vec neighborLower = d_lower - Vec(0.0, 0.0, neighborWidth);
    Vec neighborUpper = d_upper;
    neighborUpper.setZ(d_lower.z());
    Box neighborBox(neighborLower, neighborUpper);
    d_zMinus.asynchronousSendRecv(boostWorld, d_rank, iteration, 
                                  particles, neighborBox);
  }
}

void 
Patch::sendRecvMigrateZPlus(boost::mpi::communicator& boostWorld, 
                            int iteration, const REAL& neighborWidth,
                            const ParticlePArray& particles) 
{
  if (d_zPlus.d_boundary == PatchBoundary::inside) {
    Vec neighborLower = d_lower;
    neighborLower.setZ(d_upper.z());
    Vec neighborUpper = d_upper + Vec(0.0, 0.0, neighborWidth);
    Box neighborBox(neighborLower, neighborUpper);
    d_zPlus.asynchronousSendRecv(boostWorld, d_rank, iteration, 
                                  particles, neighborBox);
  }
}

void 
Patch::waitToFinishX(int iteration) 
{
  if (d_xMinus.d_boundary == PatchBoundary::inside) {
    d_xMinus.waitToFinish(d_rank, iteration);
  }
  if (d_xPlus.d_boundary == PatchBoundary::inside) {
    d_xPlus.waitToFinish(d_rank, iteration);
  }
}

void 
Patch::waitToFinishY(int iteration) 
{
  if (d_yMinus.d_boundary == PatchBoundary::inside) {
    d_yMinus.waitToFinish(d_rank, iteration);
  }
  if (d_yPlus.d_boundary == PatchBoundary::inside) {
    d_yPlus.waitToFinish(d_rank, iteration);
  }
}

void 
Patch::waitToFinishZ(int iteration) 
{
  if (d_zMinus.d_boundary == PatchBoundary::inside) {
    d_zMinus.waitToFinish(d_rank, iteration);
  }
  if (d_zPlus.d_boundary == PatchBoundary::inside) {
    d_zPlus.waitToFinish(d_rank, iteration);
  }
}

void 
Patch::insertReceivedParticles(int iteration,
                               ParticlePArray& received) 
{
  received.clear();
  d_xMinus.insertReceivedParticles(d_rank, iteration, received);
  d_xPlus.insertReceivedParticles(d_rank, iteration, received);
  d_yMinus.insertReceivedParticles(d_rank, iteration, received);
  d_yPlus.insertReceivedParticles(d_rank, iteration, received);
  d_zMinus.insertReceivedParticles(d_rank, iteration, received);
  d_zPlus.insertReceivedParticles(d_rank, iteration, received);
  removeDuplicates(received);
}

// Taken from: https://stackoverflow.com/questions/12200486/
void 
Patch::removeDuplicates(ParticlePArray& input) 
{
  std::set<ParticleP> seen;
  auto newEnd = 
    std::remove_if(input.begin(), input.end(),
                    [&seen](const ParticleP& value)
                    {
                      if (seen.find(value) != std::end(seen)) {
                        return true;
                      }
                      seen.insert(value);
                      return false;
                    });
  input.erase(newEnd, input.end());
}

void 
Patch::update(int iteration,
              const Vec& lower, const Vec& upper, const REAL& ghostWidth)
{
  d_lower = lower;
  d_upper = upper;
  d_ghostWidth = ghostWidth;
}
