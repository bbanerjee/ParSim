#ifndef ELLIP3D_PERI_PULLOUT_H
#define ELLIP3D_PERI_PULLOUT_H

#include <Simulations/Command.h>
#include <ostream>

namespace dem {
class PeridynamicsPullOut : public Command
{
public:
  virtual void execute(DiscreteElements* dem)
  {
    std::cout << "**ERROR** Execute with DEM only. "
              << "Should not be called in PeridynamicPullOut.\n";
  }
  virtual void execute(DiscreteElements* dem, pd::Peridynamics* pd);
};
}

#endif