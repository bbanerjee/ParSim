#include <DiscreteElements/Patch.h>
#include <DiscreteElements/Particle.h>
#include <Core/Types/integertypes.h>
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
  d_coords = neighborCoords;
  if (neighborRank > -1) {
    d_boundary = PatchBoundary::inside;
  } else {
    d_boundary = boundaryFlag;
  }

  //std::ostringstream out;
  //out << "myRank = " << myRank << " neighborCoords: " << neighborCoords 
  //    << " boundaryFlag " << static_cast<int>(d_boundary) << "\n";
  //std::cout << out.str();
}

void 
PatchNeighborComm::asyncSendRecv(boost::mpi::communicator& boostWorld,
                                 int myRank,
                                 int iteration,
                                 const Box& box, 
                                 const double& tolerance,
                                 const ParticlePArray& particles)
{
  //std::ostringstream out;
  //out << "myRank: " << myRank << " neighborRank: " << d_rank 
  //    << " iteration: " << iteration << " box: " << box.getMinCorner()
  //    << "," << box.getMaxCorner() 
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

void 
PatchNeighborComm::findParticlesInBox(const Box& box,
                                      const ParticlePArray& particles,
                                      const double& tolerance,
                                      ParticlePArray& inside)
{
  for (const auto& particle : particles) {
    // it is critical to use EPS
    if (box.inside(particle->currentPos(), tolerance)) {
      inside.push_back(particle);
    }
  }
}

void 
PatchNeighborComm::waitToFinish(int myRank, int iteration) 
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

void 
PatchNeighborComm::combineSentParticles(int myRank, int iteration,
                                       ParticleIDHashMap& sent) 
{
  if (!d_sentParticles.empty()) {
    for (const auto& particle : d_sentParticles) {
      sent.insert(particle->getId());
    }
  }
}

void 
PatchNeighborComm::combineReceivedParticles(int myRank, int iteration,
                                           ParticlePArray& received) 
{
  if (!d_recvParticles.empty()) {
    received.insert(received.end(), 
      d_recvParticles.begin(), d_recvParticles.end());
  }
}


Patch::Patch(MPI_Comm& cartComm, int rank, const IntVec& mpiCoords, 
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
    Vec ghostLower = d_lower;
    Vec ghostUpper = d_upper;
    ghostUpper.setX(d_lower.x() + d_ghostWidth);
    Box ghostBox(ghostLower, ghostUpper);
    d_xMinus.asyncSendRecv(boostWorld, d_rank, iteration, ghostBox, d_tolerance,
                           particles);
  }
}

void 
Patch::sendRecvGhostXPlus(boost::mpi::communicator& boostWorld, 
                          int iteration,
                          const ParticlePArray& particles) 
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

void 
Patch::sendRecvGhostYMinus(boost::mpi::communicator& boostWorld, 
                           int iteration,
                           const ParticlePArray& particles) 
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

void 
Patch::sendRecvGhostYPlus(boost::mpi::communicator& boostWorld, 
                          int iteration,
                          const ParticlePArray& particles) 
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
    d_zMinus.asyncSendRecv(boostWorld, d_rank, iteration, ghostBox, d_tolerance,
                           particles);
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
    d_zPlus.asyncSendRecv(boostWorld, d_rank, iteration, ghostBox, d_tolerance,
                           particles);
  }
}

void 
Patch::sendRecvMigrateXMinus(boost::mpi::communicator& boostWorld, 
                             int iteration, const Vec& neighborWidth,
                             const ParticlePArray& particles) 
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

void 
Patch::sendRecvMigrateXPlus(boost::mpi::communicator& boostWorld, 
                            int iteration, const Vec& neighborWidth,
                            const ParticlePArray& particles) 
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

void 
Patch::sendRecvMigrateYMinus(boost::mpi::communicator& boostWorld, 
                             int iteration, const Vec& neighborWidth,
                             const ParticlePArray& particles) 
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

void 
Patch::sendRecvMigrateYPlus(boost::mpi::communicator& boostWorld, 
                            int iteration, const Vec& neighborWidth,
                            const ParticlePArray& particles) 
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

void 
Patch::sendRecvMigrateZMinus(boost::mpi::communicator& boostWorld, 
                             int iteration, const Vec& neighborWidth,
                             const ParticlePArray& particles) 
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

void 
Patch::sendRecvMigrateZPlus(boost::mpi::communicator& boostWorld, 
                            int iteration, const Vec& neighborWidth,
                            const ParticlePArray& particles) 
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
Patch::combineReceivedParticles(int iteration,
                               ParticlePArray& received) 
{
  received.clear();
  d_xMinus.combineReceivedParticles(d_rank, iteration, received);
  d_xPlus.combineReceivedParticles(d_rank, iteration, received);
  d_yMinus.combineReceivedParticles(d_rank, iteration, received);
  d_yPlus.combineReceivedParticles(d_rank, iteration, received);
  d_zMinus.combineReceivedParticles(d_rank, iteration, received);
  d_zPlus.combineReceivedParticles(d_rank, iteration, received);
  removeDuplicates(received);
}

void 
Patch::combineSentParticlesX(int iteration,
                            ParticleIDHashMap& sent) 
{
  sent.clear();
  d_xMinus.combineSentParticles(d_rank, iteration, sent);
  d_xPlus.combineSentParticles(d_rank, iteration, sent);
}

void 
Patch::combineSentParticlesY(int iteration,
                            ParticleIDHashMap& sent) 
{
  sent.clear();
  d_yMinus.combineSentParticles(d_rank, iteration, sent);
  d_yPlus.combineSentParticles(d_rank, iteration, sent);
}

void 
Patch::combineSentParticlesZ(int iteration,
                            ParticleIDHashMap& sent) 
{
  sent.clear();
  d_zMinus.combineSentParticles(d_rank, iteration, sent);
  d_zPlus.combineSentParticles(d_rank, iteration, sent);
}

void 
Patch::combineReceivedParticlesX(int iteration,
                                ParticlePArray& received) 
{
  d_xMinus.combineReceivedParticles(d_rank, iteration, received);
  d_xPlus.combineReceivedParticles(d_rank, iteration, received);
}

void 
Patch::combineReceivedParticlesY(int iteration,
                                ParticlePArray& received) 
{
  d_yMinus.combineReceivedParticles(d_rank, iteration, received);
  d_yPlus.combineReceivedParticles(d_rank, iteration, received);
}

void 
Patch::combineReceivedParticlesZ(int iteration,
                                ParticlePArray& received) 
{
  d_zMinus.combineReceivedParticles(d_rank, iteration, received);
  d_zPlus.combineReceivedParticles(d_rank, iteration, received);
}

// Taken from: https://stackoverflow.com/questions/12200486/
void 
Patch::removeDuplicates(ParticlePArray& input) 
{
  //std::ostringstream out;
  //out << "Removing duplicates: in = " << input.size();

  std::set<ParticleID> seen;
  auto newEnd = 
    std::remove_if(input.begin(), input.end(),
                    [&seen](const ParticleP& particle)
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

void
Patch::deleteSentParticles(int iteration, 
                           const ParticleIDHashMap& sent,
                           ParticlePArray& particles)
{
  //std::ostringstream out;
  //out << "Deleting sent particle ( " << sent.size() 
  //    << ") : in = " << particles.size();

  if (sent.size() > 0) {
    particles.erase(
      std::remove_if(
        particles.begin(), particles.end(),
        [&sent](const ParticleP& particle) {
          return (sent.find(particle->getId()) != std::end(sent));
        }),
      std::end(particles));
  }

  //out << " out = " << particles.size() << "\n";
  //std::cout << out.str();
}

void
Patch::addReceivedParticles(int iteration, 
                            const ParticlePArray& received,
                            ParticlePArray& particles)
{
  particles.insert(particles.end(), received.begin(), received.end());
}

void 
Patch::removeParticlesOutsidePatch(ParticlePArray& particles)
{
  //std::ostringstream out;
  //out << "Removing outside particles: in = " << particles.size();

  Box box(d_lower, d_upper);
  double epsilon = d_tolerance;
  particles.erase(
    std::remove_if(
      particles.begin(), particles.end(),
      [&box, &epsilon](ParticleP particle) {
        if (box.inside(particle->currentPos(), epsilon)) {
          return false;
        }
        return true;
      }),
    particles.end());

  //out << " out = " << particles.size() << "\n";
  //std::cout << out.str();

}

void 
Patch::update(int iteration,
              const Vec& lower, const Vec& upper, const REAL& ghostWidth)
{
  d_lower = lower;
  d_upper = upper;
  d_ghostWidth = ghostWidth;
}
