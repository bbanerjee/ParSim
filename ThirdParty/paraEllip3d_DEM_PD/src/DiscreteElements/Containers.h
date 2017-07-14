#ifndef CONTAINERS_H
#define CONTAINERS_H

#include <memory>
#include <vector>

namespace dem {

class Particle;
class Contact;
class ContactTgt;
class Spring;

using ParticleUP = std::unique_ptr<Particle>;
using ParticleP = std::shared_ptr<Particle>;
using ParticlePArray = std::vector<ParticleP>;

using MembraneParticlePArray = std::vector<std::vector<ParticlePArray>>;

using ContactArray = std::vector<Contact>;
using ContactTangentArray = std::vector<ContactTgt>;

using SpringUP = std::unique_ptr<Spring>;
using SpringUPArray = std::vector<SpringUP>;
}

#endif
