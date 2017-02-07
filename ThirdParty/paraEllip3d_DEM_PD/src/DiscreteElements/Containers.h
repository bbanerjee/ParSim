#ifndef CONTAINERS_H
#define CONTAINERS_H

//#include <DiscreteElements/Particle.h>
//#include <Peridynamics/PeriParticle.h>
#include <memory>
#include <vector>

namespace periDynamics {

class PeriParticle;
class PeriBond;

using PeriParticleUP = std::unique_ptr<PeriParticle>;
using PeriParticleP = std::shared_ptr<PeriParticle>;
using PeriParticlePArray = std::vector<PeriParticleP>;

using PeriBondP = std::shared_ptr<PeriBond>;
using PeriBondPArray = std::vector<PeriBondP>;
}

namespace dem {

class Particle;
class Contact;
class ContactTgt;
class Spring;
class Boundary;
class BoundaryTgt;
class PeriBoundaryBond;
class PeriDEMBond;

using ParticleUP = std::unique_ptr<Particle>;
using ParticleP = std::shared_ptr<Particle>;
using ParticlePArray = std::vector<ParticleP>;

using BoundaryUP = std::unique_ptr<Boundary>;
using BoundaryUPArray = std::vector<BoundaryUP>;

using BoundaryP = std::shared_ptr<Boundary>;
using BoundaryPArray = std::vector<BoundaryP>;

using BoundaryTangentArray = std::vector<BoundaryTgt>;

using MembraneParticlePArray = std::vector<std::vector<ParticlePArray>>;

using ContactArray = std::vector<Contact>;
using ContactTangentArray = std::vector<ContactTgt>;

using SpringUP = std::unique_ptr<Spring>;
using SpringUPArray = std::vector<SpringUP>;

using PeriBoundaryBondP = std::shared_ptr<PeriBoundaryBond>;
using PeriBoundaryBondPArray = std::vector<PeriBoundaryBondP>;

using PeriDEMBondP = std::shared_ptr<PeriDEMBond>;
using PeriDEMBondPArray = std::vector<PeriDEMBondP>;
}

#endif
