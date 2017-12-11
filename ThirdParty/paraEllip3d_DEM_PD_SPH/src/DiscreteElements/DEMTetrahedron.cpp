#include <DiscreteElements/DEMParticle.h>
#include <DiscreteElements/DEMTetrahedron.h>
#include <iostream>

using namespace dem;

DEMTetrahedron::DEMTetrahedron()
{
  d_nodeIDs.reserve(4);
  d_particles.reserve(4);
  for (auto& nodeID : d_nodeIDs) nodeID = -1;
  d_B.setZero();
  d_velGrad = 0;
  d_defGrad.Identity();
  d_defGradRate = 0;
}

DEMTetrahedron::DEMTetrahedron(NodeID n1, NodeID n2, NodeID n3, NodeID n4,
                                const DEMParticlePArray& particles)
  : DEMTetrahedron()
{
  d_nodeIDs[0] = n1;
  d_nodeIDs[1] = n2;
  d_nodeIDs[2] = n3;
  d_nodeIDs[3] = n4;
  for (const auto& nodeID : d_nodeIDs) {
    DEMParticleP particle = particles[nodeID];
    d_particles.push_back(std::move(particle));
  }
}

REAL 
DEMTetrahedron::volume() const
{
  return std::abs(signedVolume());
}

REAL 
DEMTetrahedron::signedVolume() const
{
  Vec p0 = d_particles[0]->currentPosition();
  Vec p1 = d_particles[1]->currentPosition();
  Vec p2 = d_particles[2]->currentPosition();
  Vec p3 = d_particles[3]->currentPosition();
  return volume(p0, p1, p2, p3);
}

REAL 
DEMTetrahedron::initialVolume() const
{
  return std::abs(signedInitialVolume());
}

REAL 
DEMTetrahedron::signedInitialVolume() const
{
  Vec p0 = d_particles[0]->initialPosition();
  Vec p1 = d_particles[1]->initialPosition();
  Vec p2 = d_particles[2]->initialPosition();
  Vec p3 = d_particles[3]->initialPosition();
  return volume(p0, p1, p2, p3);
}

REAL 
DEMTetrahedron::volume(const Vec& p0, const Vec& p1, 
                       const Vec& p2, const Vec& p3) const
{
  Vec vec01 = p1 - p0;
  Vec vec02 = p2 - p0;
  Vec vec03 = p3 - p0;
  Vec vec01x02 = cross(vec01, vec02);
  return 1.0/6.0*dot(vec01x02, vec03);
}

std::vector<double> 
DEMTetrahedron::basisFunctions(const Vec& point) const
{
  Vec p0 = d_particles[0]->currentPosition();
  Vec p1 = d_particles[1]->currentPosition();
  Vec p2 = d_particles[2]->currentPosition();
  Vec p3 = d_particles[3]->currentPosition();
  Vec rst = transformToRST(point, p0, p1, p2, p3);
  return referenceBasisFunctions(rst[0], rst[1], rst[2]);
}

std::vector<Vec> 
DEMTetrahedron::basisDerivatives(const Vec& point) const
{
  Vec p0 = d_particles[0]->currentPosition();
  Vec p1 = d_particles[1]->currentPosition();
  Vec p2 = d_particles[2]->currentPosition();
  Vec p3 = d_particles[3]->currentPosition();
  Matrix3 A = transformationMatrix(p0, p1, p2, p3);
  Matrix3 Ainv = A.Inverse();
  Vec rst = transformToRST(point, p0, Ainv);
  std::vector<Vec> refDerivs = referenceBasisDerivatives(rst[0], rst[1], rst[2]);
  Vec dN1_dx = Ainv * refDerivs[0];
  Vec dN2_dx = Ainv * refDerivs[1];
  Vec dN3_dx = Ainv * refDerivs[2];
  Vec dN4_dx = Ainv * refDerivs[3];
  return {{dN1_dx, dN2_dx, dN3_dx, dN4_dx}};
}

std::pair<std::vector<double>, std::vector<Vec>> 
DEMTetrahedron::basisFunctionsAndDerivatives(const Vec& point) const
{
  Vec p0 = d_particles[0]->currentPosition();
  Vec p1 = d_particles[1]->currentPosition();
  Vec p2 = d_particles[2]->currentPosition();
  Vec p3 = d_particles[3]->currentPosition();
  Matrix3 A = transformationMatrix(p0, p1, p2, p3);
  Matrix3 Ainv = A.Inverse();
  Vec rst = transformToRST(point, p0, Ainv);
  auto basis = referenceBasisFunctions(rst[0], rst[1], rst[2]);
  std::vector<Vec> refDerivs = referenceBasisDerivatives(rst[0], rst[1], rst[2]);
  Vec dN1_dx = Ainv * refDerivs[0];
  Vec dN2_dx = Ainv * refDerivs[1];
  Vec dN3_dx = Ainv * refDerivs[2];
  Vec dN4_dx = Ainv * refDerivs[3];
  std::vector<Vec> derivs = {{dN1_dx, dN2_dx, dN3_dx, dN4_dx}};
  return std::make_pair(basis, derivs);
}

void DEMTetrahedron::updateBMatrix()
{
  //Vec p0 = d_particles[0]->currentPosition();
  //Vec p1 = d_particles[1]->currentPosition();
  //Vec p2 = d_particles[2]->currentPosition();
  //Vec p3 = d_particles[3]->currentPosition();
}

void DEMTetrahedron::updateVelGrad()
{
}
void DEMTetrahedron::updateDefGrad()
{
}