#ifndef ELLIP3D_BURSTING_DAM_2D_H
#define ELLIP3D_BURSTING_DAM_2D_H

#include <Simulations/Command.h>
#include <ostream>

namespace dem {

class BurstingDam2D: public Command
{
public:
  virtual void execute(DiscreteElements* dem)
  {
    std::cout << "**ERROR** Execute with DEM only. "
              << "Should not be called in BurstingDam2D.\n";
  }
  virtual void execute(DiscreteElements* dem, sph::SmoothParticleHydrodynamics* pd);
};
}

#endif