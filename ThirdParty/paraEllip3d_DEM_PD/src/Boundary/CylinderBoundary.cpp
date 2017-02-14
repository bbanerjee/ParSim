#include <Boundary/CylinderBoundary.h>
#include <DiscreteElements/Particle.h>
// use both pointer to and variable of class Particle

using namespace dem;

CylinderBoundary::CylinderBoundary(std::size_t tp, std::ifstream& ifs)
  : Boundary()
{
  // These are declared in the boundary base class
  b_type = tp;
  ifs >> b_extraNum;
  ifs >> b_id;

  REAL dx, dy, dz, px, py, pz;
  ifs >> dx >> dy >> dz >> px >> py >> pz >> radius;
  direc = Vec(dx, dy, dz);
  point = Vec(px, py, pz);
}

CylinderBoundary::CylinderBoundary(BoundaryId id, BoundaryType tp,
                                   const XMLProblemSpec& ps)
  : Boundary()
{
}

CylinderBoundary::CylinderBoundary(BoundaryId id, BoundaryType tp,
                                   const JsonProblemSpec& ps)
  : Boundary()
{
}

void
CylinderBoundary::findBdryContact(ParticlePArray& ptcls)
{
  possParticle.clear();
  contactInfo.clear();

  for (auto& ptcl : ptcls) {
    if (ptcl->getType() == 0) { // only process free particles, excluding type 5
      ;
    }
  }
}

void
CylinderBoundary::boundaryForce(BoundaryTangentArrayMap& boundaryTgtMap)
{
  // for each plane boundary, define a temparory variable vtmp to use,
  // better than define a member variable which needs to be cleared.
  // and vtmp is initialized as empty in each iteration.
  BoundaryTangentArray vtmp;

  // for each possible boundary particle
  for (auto it = possParticle.begin(); it != possParticle.end(); ++it)
    ; // (*it)->cylinderRBForce();

  // checkout tangential forces and displacements after each particle is
  // processed
  boundaryTgtMap[this->b_id] = vtmp;

  updateStatForce();
}
