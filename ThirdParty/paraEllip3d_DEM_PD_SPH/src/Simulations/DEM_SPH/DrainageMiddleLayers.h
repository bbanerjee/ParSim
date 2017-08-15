#ifndef ELLIP3D_DRAINAGE__MIDDLE_LAYERS_3D_H
#define ELLIP3D_DRAINAGE__MIDDLE_LAYERS_3D_H

#include <Simulations/Command.h>
#include <ostream>

namespace dem {

class DrainageMiddleLayers : public Command
{
public:
  virtual void execute(DiscreteElements* dem)
  {
    std::cout << "**ERROR** Execute with DEM only. "
              << "Should not be called in DrainageMiddleLayers.\n";
  }
  virtual void execute(DiscreteElements* dem, pd::Peridynamics* pd)
  {
    std::cout << "**ERROR** Execute with DEM-Peridynamics. "
              << "Should not be called in DrainageMiddleLayers.\n";
  }
  virtual void execute(DiscreteElements* dem, sph::SmoothParticleHydro* sph);
};
}

#endif