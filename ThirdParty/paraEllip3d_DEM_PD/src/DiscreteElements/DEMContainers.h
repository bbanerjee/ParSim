#ifndef CONTAINERS_H
#define CONTAINERS_H

#include <Core/Types/integertypes.h>

#include <memory>
#include <vector>
#include <unordered_set>

namespace dem {

class DEMParticle;
class DEMContact;
class ContactTgt;
class Spring;

using ParticleUP = std::unique_ptr<DEMParticle>;
using ParticleP = std::shared_ptr<DEMParticle>;
using ParticlePArray = std::vector<ParticleP>;
using ParticleIDHashMap = std::unordered_set<ParticleID>;

using MembraneParticlePArray = std::vector<std::vector<ParticlePArray>>;

using ContactArray = std::vector<DEMContact>;
using ContactTangentArray = std::vector<ContactTgt>;

using SpringUP = std::unique_ptr<Spring>;
using SpringUPArray = std::vector<SpringUP>;
}

#endif
