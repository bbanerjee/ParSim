#ifndef SPH_CONTAINERS_H
#define SPH_CONTAINERS_H

#include <memory>
#include <vector>

namespace sph {

class SPHParticle;

using SPHParticleUP = std::unique_ptr<SPHParticle>;
using SPHParticleP = std::shared_ptr<SPHParticle>;
using SPHParticlePArray = std::vector<SPHParticleP>;
}

#endif
