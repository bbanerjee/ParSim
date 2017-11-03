#include <Boundary/CylinderBoundary.h>
#include <DiscreteElements/DEMParticle.h>
// use both pointer to and variable of class DEMParticle

using namespace dem;

CylinderBoundary::CylinderBoundary(Boundary::BoundaryType tp, std::ifstream& ifs)
  : Boundary()
{
  // These are declared in the boundary base class
  b_type = tp;
  ifs >> b_extraNum;
  int id;
  ifs >> id;
  b_id = Boundary::getBoundaryID(id);

  REAL dx, dy, dz, px, py, pz;
  ifs >> dx >> dy >> dz >> px >> py >> pz >> d_radius;
  d_direction = Vec(dx, dy, dz);
  d_point = Vec(px, py, pz);
}

CylinderBoundary::CylinderBoundary(Boundary::BoundaryType tp, BoundaryID id, 
                                   const XMLProblemSpec& ps)
  : Boundary()
{
}

CylinderBoundary::CylinderBoundary(Boundary::BoundaryType tp, BoundaryID id, 
                                   const JsonProblemSpec& ps)
  : Boundary()
{
}

void
CylinderBoundary::findBoundaryContacts(DEMParticlePArray& particles)
{
  b_probableBoundaryParticles.clear();
  b_contacts.clear();

  for (auto& particle : particles) {
    // only process free particles, excluding type 5
    if (particle->getType() == DEMParticle::DEMParticleType::FREE) {
    }
  }
}

void
CylinderBoundary::boundaryForce(BoundaryTangentArrayMap& boundaryTangentMap)
{
  // for each plane boundary, define a temparory variable vtmp to use,
  // better than define a member variable which needs to be cleared.
  // and vtmp is initialized as empty in each iteration.
  BoundaryTangentArray vtmp;

  // for each possible boundary particle
  for (auto it = b_probableBoundaryParticles.begin(); it != b_probableBoundaryParticles.end(); ++it)
    ; // (*it)->cylinderRBForce();

  // checkout tangential forces and displacements after each particle is
  // processed
  boundaryTangentMap[static_cast<size_t>(this->b_id)] = vtmp;

  updateStatForce();
}
