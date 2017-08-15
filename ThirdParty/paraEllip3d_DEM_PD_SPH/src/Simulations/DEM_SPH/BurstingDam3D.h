#ifndef ELLIP3D_BURSTING_DAM_3D_H
#define ELLIP3D_BURSTING_DAM_3D_H

#include <Simulations/Command.h>
#include <ostream>

namespace dem {

class BurstingDam3D : public Command
{
public:
  virtual void execute(DiscreteElements* dem)
  {
    std::cout << "**ERROR** Execute with DEM only. "
              << "Should not be called in BurstingDam3D.\n";
  }
  virtual void execute(DiscreteElements* dem, pd::Peridynamics* pd)
  {
    std::cout << "**ERROR** Execute with DEM-Peridynamics. "
              << "Should not be called in BurstingDam3D.\n";
  }
  virtual void execute(DiscreteElements* dem, sph::SmoothParticleHydro* sph);
};
}

#endif