#include <DiscreteElements/Patch.h>
#include <DiscreteElements/Particle.h>
#include <set>

using namespace dem;

void 
PatchNeighborComm::setNeighbor(MPI_Comm& cartComm, const IntVec& neighborCoords,
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
}

void 
PatchNeighborComm::asynchronousSendRecv(boost::mpi::communicator& boostWorld,
                                        int iteration,
                                        const ParticlePArray& particles,
                                        const Box& ghostBox) 
{
  findParticlesInGhostBox(ghostBox, particles, d_sendGhostParticles);
  d_sendRecvReq[0] = 
    boostWorld.isend(d_rank, d_mpiTag, d_sendGhostParticles);
  d_sendRecvReq[1] = 
    boostWorld.irecv(d_rank, d_mpiTag, d_recvGhostParticles);
}

void 
PatchNeighborComm::findParticlesInGhostBox(const Box& ghostBox,
                                           const ParticlePArray& particles,
                                           ParticlePArray& ghostParticles)
{
  ghostParticles.clear();
  for (const auto& particle : particles) {
    // it is critical to use EPS
    if (ghostBox.inside(particle->currentPos(), EPS)) {
      ghostParticles.push_back(particle);
    }
  }
}

void 
PatchNeighborComm::waitToFinish(int iteration) 
{
  boost::mpi::wait_all(d_sendRecvReq, d_sendRecvReq + 2);
}

void 
PatchNeighborComm::insertReceivedParticles(int iteration,
                                           ParticlePArray& received) 
{
  if (!d_recvGhostParticles.empty()) {
    received.insert(received.end(), 
      d_recvGhostParticles.begin(), d_recvGhostParticles.end());
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
  d_xMinus.setNeighbor(cartComm, neighborCoords, PatchBoundary::xminus);
}

void 
Patch::setXPlus(MPI_Comm& cartComm)
{
  IntVec neighborCoords = d_patchMPICoords;
  ++neighborCoords.x();
  d_xPlus.setNeighbor(cartComm, neighborCoords, PatchBoundary::xplus);
}

void 
Patch::setYMinus(MPI_Comm& cartComm)
{
  IntVec neighborCoords = d_patchMPICoords;
  --neighborCoords.y();
  d_yMinus.setNeighbor(cartComm, neighborCoords, PatchBoundary::yminus);
}

void 
Patch::setYPlus(MPI_Comm& cartComm)
{
  IntVec neighborCoords = d_patchMPICoords;
  ++neighborCoords.y();
  d_yPlus.setNeighbor(cartComm, neighborCoords, PatchBoundary::yplus);
}

void 
Patch::setZMinus(MPI_Comm& cartComm)
{
  IntVec neighborCoords = d_patchMPICoords;
  --neighborCoords.z();
  d_zMinus.setNeighbor(cartComm, neighborCoords, PatchBoundary::zminus);
}

void 
Patch::setZPlus(MPI_Comm& cartComm)
{
  IntVec neighborCoords = d_patchMPICoords;
  ++neighborCoords.z();
  d_zPlus.setNeighbor(cartComm, neighborCoords, PatchBoundary::zplus);
}

void 
Patch::sendRecvXMinus(boost::mpi::communicator& boostWorld, 
                      int iteration,
                      const ParticlePArray& particles) 
{
  if (d_xMinus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower - Vec(0.0, d_ghostWidth, d_ghostWidth);
    Vec ghostUpper = d_upper + Vec(0.0, d_ghostWidth, d_ghostWidth);
    ghostUpper.setX(d_lower.x() + d_ghostWidth);
    Box ghostBox(ghostLower, ghostUpper);
    d_xMinus.asynchronousSendRecv(boostWorld, iteration, particles, ghostBox);
  }
}

void 
Patch::sendRecvXPlus(boost::mpi::communicator& boostWorld, 
                     int iteration,
                     const ParticlePArray& particles) 
{
  if (d_xPlus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower - Vec(0.0, d_ghostWidth, d_ghostWidth);
    ghostLower.setX(d_upper.x() - d_ghostWidth);
    Vec ghostUpper = d_upper + Vec(0.0, d_ghostWidth, d_ghostWidth);
    Box ghostBox(ghostLower, ghostUpper);
    d_xPlus.asynchronousSendRecv(boostWorld, iteration, particles, ghostBox);
  }
}

void 
Patch::sendRecvYMinus(boost::mpi::communicator& boostWorld, 
                      int iteration,
                      const ParticlePArray& particles) 
{
  if (d_yMinus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower - Vec(d_ghostWidth, 0.0, d_ghostWidth);
    Vec ghostUpper = d_upper + Vec(d_ghostWidth, 0.0, d_ghostWidth);
    ghostUpper.setY(d_lower.y() + d_ghostWidth);
    Box ghostBox(ghostLower, ghostUpper);
    d_yMinus.asynchronousSendRecv(boostWorld, iteration, particles, ghostBox);
  }
}

void 
Patch::sendRecvYPlus(boost::mpi::communicator& boostWorld, 
                     int iteration,
                     const ParticlePArray& particles) 
{
  if (d_yPlus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower - Vec(d_ghostWidth, 0.0, d_ghostWidth);
    ghostLower.setY(d_upper.y() - d_ghostWidth);
    Vec ghostUpper = d_upper + Vec(d_ghostWidth, 0.0, d_ghostWidth);
    Box ghostBox(ghostLower, ghostUpper);
    d_yPlus.asynchronousSendRecv(boostWorld, iteration, particles, ghostBox);
  }
}

void 
Patch::sendRecvZMinus(boost::mpi::communicator& boostWorld, 
                      int iteration,
                      const ParticlePArray& particles) 
{
  if (d_zMinus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower - Vec(d_ghostWidth, d_ghostWidth, 0.0);
    Vec ghostUpper = d_upper + Vec(d_ghostWidth, d_ghostWidth, 0.0);
    ghostUpper.setZ(d_lower.z() + d_ghostWidth);
    Box ghostBox(ghostLower, ghostUpper);
    d_zMinus.asynchronousSendRecv(boostWorld, iteration, particles, ghostBox);
  }
}

void 
Patch::sendRecvZPlus(boost::mpi::communicator& boostWorld, 
                     int iteration,
                     const ParticlePArray& particles) 
{
  if (d_zPlus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower - Vec(d_ghostWidth, d_ghostWidth, 0.0);
    ghostLower.setZ(d_upper.z() - d_ghostWidth);
    Vec ghostUpper = d_upper + Vec(d_ghostWidth, d_ghostWidth, 0.0);
    Box ghostBox(ghostLower, ghostUpper);
    d_zPlus.asynchronousSendRecv(boostWorld, iteration, particles, ghostBox);
  }
}

void 
Patch::waitToFinishX(int iteration) 
{
  if (d_xMinus.d_boundary == PatchBoundary::inside) {
    d_xMinus.waitToFinish(iteration);
  }
  if (d_xPlus.d_boundary == PatchBoundary::inside) {
    d_xPlus.waitToFinish(iteration);
  }
}

void 
Patch::waitToFinishY(int iteration) 
{
  if (d_yMinus.d_boundary == PatchBoundary::inside) {
    d_yMinus.waitToFinish(iteration);
  }
  if (d_yPlus.d_boundary == PatchBoundary::inside) {
    d_yPlus.waitToFinish(iteration);
  }
}

void 
Patch::waitToFinishZ(int iteration) 
{
  if (d_zMinus.d_boundary == PatchBoundary::inside) {
    d_zMinus.waitToFinish(iteration);
  }
  if (d_zPlus.d_boundary == PatchBoundary::inside) {
    d_zPlus.waitToFinish(iteration);
  }
}

void 
Patch::insertReceivedParticles(int iteration,
                               ParticlePArray& received) 
{
  received.clear();
  d_xMinus.insertReceivedParticles(iteration, received);
  d_xPlus.insertReceivedParticles(iteration, received);
  d_yMinus.insertReceivedParticles(iteration, received);
  d_yPlus.insertReceivedParticles(iteration, received);
  d_zMinus.insertReceivedParticles(iteration, received);
  d_zPlus.insertReceivedParticles(iteration, received);
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
