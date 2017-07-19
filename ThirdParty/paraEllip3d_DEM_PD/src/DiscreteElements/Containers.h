#ifndef CONTAINERS_H
#define CONTAINERS_H

#include <Core/Types/integertypes.h>

#include <memory>
#include <vector>
#include <unordered_set>

namespace dem {

class Particle;
class Contact;
class ContactTgt;
class Spring;

using ParticleUP = std::unique_ptr<Particle>;
using ParticleP = std::shared_ptr<Particle>;
using ParticlePArray = std::vector<ParticleP>;
using ParticleIDHashMap = std::unordered_set<ParticleID>;

using MembraneParticlePArray = std::vector<std::vector<ParticlePArray>>;

using ContactArray = std::vector<Contact>;
using ContactTangentArray = std::vector<ContactTgt>;

using SpringUP = std::unique_ptr<Spring>;
using SpringUPArray = std::vector<SpringUP>;
}

#endif
