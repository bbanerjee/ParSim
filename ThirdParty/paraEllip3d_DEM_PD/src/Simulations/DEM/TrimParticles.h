#ifndef ELLIP3D_TRIM_PARTICLES_H
#define ELLIP3D_TRIM_PARTICLES_H

#include <Simulations/Command.h>
#include <ostream>

namespace dem {
class TrimParticles : public Command
{
public:
  virtual void execute(DiscreteElements* dem);
  virtual void execute(DiscreteElements* dem, pd::Peridynamics* pd)
  {
    std::cout << "**ERROR** Execute with DEM + PD. "
              << "Should not be called in TrimParticles.\n";
  }
  virtual void execute(DiscreteElements* dem, sph::SmoothParticleHydro* sph)
  {
    std::cout << "**ERROR** Execute with DEM + SPH. "
              << "Should not be called in TrimParticles.\n";
  }
};
}

#endif
