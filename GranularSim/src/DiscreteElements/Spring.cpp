#include <DiscreteElements/Spring.h>

namespace dem {

Spring::Spring(DEMParticle& ptcl1, DEMParticle& ptcl2, REAL modulus)
  : p1(ptcl1)
  , p2(ptcl2)
  , young(modulus)
{
  init(p1, p2);
}

Spring::Spring(std::vector<DEMParticle*>& ParticleVec, std::size_t id1,
               std::size_t id2, REAL modulus)
  : p1(*ParticleVec[id1])
  , p2(*ParticleVec[id2])
  , young(modulus)
{
  init(p1, p2);
}

Vec
Spring::getDeformation()
{
  Vec dir = p2.currentPosition() - p1.currentPosition();
  REAL length = vnormL2(dir);
  return (length - length0) * normalize(dir);
}

void
Spring::applyForce()
{ // p1 adds this force, p2 subtracts this force
  Vec dir = p2.currentPosition() - p1.currentPosition();
  REAL length = vnormL2(dir);
  if (length - length0 > EPS) {
    Vec force = (length - length0) * normalize(dir);
    force *= ks;
    p1.addForce(force);
    p2.addForce(-force);
  }
}
}
