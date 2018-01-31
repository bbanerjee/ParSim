#ifndef ELLIP3D_DRAINAGE_3D_H
#define ELLIP3D_DRAINAGE_3D_H

#include <Simulations/Command.h>
#include <ostream>

namespace dem {

class Drainage : public Command
{
public:
  virtual void execute(DiscreteElements* dem)
  {
    std::cout << "**ERROR** Execute with DEM only. "
              << "Should not be called in Drainage.\n";
  }
  virtual void execute(DiscreteElements* dem, pd::Peridynamics* pd)
  {
    std::cout << "**ERROR** Execute with DEM-Peridynamics. "
              << "Should not be called in Drainage.\n";
  }
  virtual void execute(DiscreteElements* dem, sph::SmoothParticleHydro* sph);
};
}

#endif