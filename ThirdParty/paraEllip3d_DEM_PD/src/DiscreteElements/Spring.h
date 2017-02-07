#ifndef SPRING_H
#define SPRING_H

#include <Core/Math/Vec.h>
#include <Core/Types/realtypes.h>
#include <DiscreteElements/Particle.h>
#include <cassert>
#include <cstddef>

namespace dem {

class Spring
{

public:
  Spring(Particle& p1, Particle& p2, REAL young);
  Spring(std::vector<Particle*>& ParticleVec, std::size_t id1, std::size_t id2,
         REAL young);

  REAL getLength0() const { return length0; }
  REAL getLength() const { return vfabs(p2.getCurrPos() - p1.getCurrPos()); }
  Vec getDeformation();
  void applyForce();
  std::size_t getParticleId1() const { return p1.getId(); }
  std::size_t getParticleId2() const { return p2.getId(); }

private:
  Particle& p1;
  Particle& p2;
  REAL young;   // Young's modulus
  REAL ks;      // stiffness
  REAL length0; // equilibrium length

  void init(Particle& p1, Particle& p2)
  {
    length0 = vfabs(p2.getCurrPos() - p1.getCurrPos());
    REAL radius = p1.getA();
    assert(radius == p2.getA());
    ks = young * 4 * radius * radius / length0;
  }
};
}

#endif
