#include <DiscreteElements/DEMParticle.h>
#include <DiscreteElements/DEMTetrahedron.h>
#include <iostream>

using namespace dem;

DEMTetrahedron::DEMTetrahedron()
{
  d_nodeIDs.resize(4);
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

  if (signedVolume() < 0) {
    updateNodeOrder();
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
  return signedVolume(p0, p1, p2, p3);
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
  return signedVolume(p0, p1, p2, p3);
}

REAL 
DEMTetrahedron::signedVolume(const Vec& p0, const Vec& p1, 
                       const Vec& p2, const Vec& p3) const
{
  Vec vec01 = p1 - p0;
  Vec vec02 = p2 - p0;
  Vec vec03 = p3 - p0;
  Vec vec01x02 = cross(vec01, vec02);
  return 1.0/6.0*dot(vec01x02, vec03);
}

void
DEMTetrahedron::updateNodeOrder()
{
  std::swap(d_nodeIDs[1], d_nodeIDs[2]);
  std::swap(d_particles[1], d_particles[2]);
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
DEMTetrahedron::basisDerivatives() const
{
  Vec p0 = d_particles[0]->currentPosition();
  Vec p1 = d_particles[1]->currentPosition();
  Vec p2 = d_particles[2]->currentPosition();
  Vec p3 = d_particles[3]->currentPosition();
  Matrix3 A = transformationMatrix(p0, p1, p2, p3);
  Matrix3 Ainv = A.Inverse();
  std::vector<Vec> refDerivs = referenceBasisDerivatives();
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
  std::vector<Vec> refDerivs = referenceBasisDerivatives();
  Vec dN1_dx = Ainv * refDerivs[0];
  Vec dN2_dx = Ainv * refDerivs[1];
  Vec dN3_dx = Ainv * refDerivs[2];
  Vec dN4_dx = Ainv * refDerivs[3];
  std::vector<Vec> derivs = {{dN1_dx, dN2_dx, dN3_dx, dN4_dx}};
  return std::make_pair(basis, derivs);
}

void DEMTetrahedron::updateBMatrix()
{
  auto derivs = basisDerivatives();
  d_B(0,0) = derivs[0][0]; d_B(0,1) = derivs[0][1]; d_B(0,2) = derivs[0][2];
  d_B(1,0) = derivs[1][0]; d_B(1,1) = derivs[1][1]; d_B(1,2) = derivs[1][2];
  d_B(2,0) = derivs[2][0]; d_B(2,1) = derivs[2][1]; d_B(2,2) = derivs[2][2];
  d_B(3,0) = derivs[3][0]; d_B(3,1) = derivs[3][1]; d_B(3,2) = derivs[3][2];
}

void DEMTetrahedron::updateVelGrad()
{
  Vec v0 = d_particles[0]->currentVelocity();
  Vec v1 = d_particles[1]->currentVelocity();
  Vec v2 = d_particles[2]->currentVelocity();
  Vec v3 = d_particles[3]->currentVelocity();
  Matrix3x4 vel;
  vel(0,0) = v0[0]; vel(1,0) = v0[1]; vel(2,0) = v0[2];
  vel(0,1) = v1[0]; vel(1,1) = v1[1]; vel(2,1) = v1[2];
  vel(0,2) = v2[0]; vel(1,2) = v2[1]; vel(2,2) = v2[2];
  vel(0,3) = v3[0]; vel(1,3) = v3[1]; vel(2,3) = v3[2];
  auto L = vel * d_B;
  d_velGrad(0,0) = L(0,0); d_velGrad(0,1) = L(0,1); d_velGrad(0,2) = L(0,2);
  d_velGrad(1,0) = L(1,0); d_velGrad(1,1) = L(1,1); d_velGrad(1,2) = L(1,2);
  d_velGrad(2,0) = L(2,0); d_velGrad(2,1) = L(2,1); d_velGrad(2,2) = L(2,2);
}

void DEMTetrahedron::updateDispGrad()
{
  Vec u0 = d_particles[0]->currentPosition() - d_particles[0]->initialPosition();
  Vec u1 = d_particles[1]->currentPosition() - d_particles[1]->initialPosition();
  Vec u2 = d_particles[2]->currentPosition() - d_particles[2]->initialPosition();
  Vec u3 = d_particles[3]->currentPosition() - d_particles[3]->initialPosition();
  Matrix3x4 disp;
  disp(0,0) = u0[0]; disp(1,0) = u0[1]; disp(2,0) = u0[2];
  disp(0,1) = u1[0]; disp(1,1) = u1[1]; disp(2,1) = u1[2];
  disp(0,2) = u2[0]; disp(1,2) = u2[1]; disp(2,2) = u2[2];
  disp(0,3) = u3[0]; disp(1,3) = u3[1]; disp(2,3) = u3[2];
  auto gradu = disp * d_B;
  d_dispGrad(0,0) = gradu(0,0); d_dispGrad(0,1) = gradu(0,1); 
  d_dispGrad(0,2) = gradu(0,2);
  d_dispGrad(1,0) = gradu(1,0); d_dispGrad(1,1) = gradu(1,1); 
  d_dispGrad(1,2) = gradu(1,2);
  d_dispGrad(2,0) = gradu(2,0); d_dispGrad(2,1) = gradu(2,1); 
  d_dispGrad(2,2) = gradu(2,2);
}

void DEMTetrahedron::updateDefGrad()
{
  Matrix3 I; I.Identity();
  d_defGrad = (I - d_dispGrad).Inverse();
  d_defGradRate = d_velGrad * d_defGrad;
}

void DEMTetrahedron::updateGradients()
{
  updateBMatrix();
  updateVelGrad();
  updateDispGrad();
  updateDefGrad();
}