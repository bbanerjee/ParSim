#ifndef PERI_CONTAINERS_H
#define PERI_CONTAINERS_H

#include <Peridynamics/PeriElement.h>

#include <memory>
#include <vector>

namespace pd {

class PeriParticle;
class PeriBond;
class PeriBoundaryBond;
class PeriDEMBond;

using PeriParticleUP = std::unique_ptr<PeriParticle>;
using PeriParticleP = std::shared_ptr<PeriParticle>;
using PeriParticlePArray = std::vector<PeriParticleP>;

using PeriBondP = std::shared_ptr<PeriBond>;
using PeriBondPArray = std::vector<PeriBondP>;

using PeriBoundaryBondP = std::shared_ptr<PeriBoundaryBond>;
using PeriBoundaryBondPArray = std::vector<PeriBoundaryBondP>;

using PeriDEMBondP = std::shared_ptr<PeriDEMBond>;
using PeriDEMBondPArray = std::vector<PeriDEMBondP>;

using PeriElementArray = std::vector<PeriElement>;
}

#endif
