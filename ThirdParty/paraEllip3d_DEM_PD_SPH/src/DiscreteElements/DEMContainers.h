#ifndef CONTAINERS_H
#define CONTAINERS_H

#include <Core/Types/IntegerTypes.h>

#include <memory>
#include <vector>
#include <unordered_set>

namespace dem {

class DEMParticle;
class DEMContact;
class DEMContactTangent;
class Spring;
class DEMTetrahedron;

using DEMParticleUP = std::unique_ptr<DEMParticle>;
using DEMParticleP = std::shared_ptr<DEMParticle>;
using DEMParticlePArray = std::vector<DEMParticleP>;

using ParticleIDHashMap = std::unordered_set<ParticleID>;

using MembraneParticlePArray = std::vector<std::vector<DEMParticlePArray>>;

using ContactArray = std::vector<DEMContact>;
using ContactTangentArray = std::vector<DEMContactTangent>;

using SpringUP = std::unique_ptr<Spring>;
using SpringUPArray = std::vector<SpringUP>;

using DEMTetrahedronUP = std::unique_ptr<DEMTetrahedron>;
using DEMTetrahedronP = std::shared_ptr<DEMTetrahedron>;

using DEMTetrahedronArray = std::vector<DEMTetrahedron>;
using DEMTetrahedronUPArray = std::vector<DEMTetrahedronUP>;
using DEMTetrahedronPArray = std::vector<DEMTetrahedronP>;
}

#endif
