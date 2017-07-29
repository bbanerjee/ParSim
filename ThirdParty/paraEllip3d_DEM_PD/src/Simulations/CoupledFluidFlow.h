#ifndef ELLIP3D_COUPLED_FLUID_FLOW_H
#define ELLIP3D_COUPLED_FLUID_FLOW_H

#include <Simulations/Command.h>
#include <ostream>

namespace dem {
class CoupledFluidFlow : public Command
{
public:
  virtual void execute(DiscreteElements* dem);
  virtual void execute(DiscreteElements* dem, pd::Peridynamics* pd)
  {
    std::cout << "**ERROR** Execute with DEM + PD. "
              << "Should not be called in CoupledFluidFlow.\n";
  }
};
}

#endif