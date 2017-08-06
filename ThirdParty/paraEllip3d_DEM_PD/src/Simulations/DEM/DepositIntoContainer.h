#ifndef ELLIP3D_DEPOSIT_INTO_CONTAINER_H
#define ELLIP3D_DEPOSIT_INTO_CONTAINER_H

#include <Simulations/Command.h>
#include <ostream>

namespace dem {
class DepositIntoContainer : public Command
{
public:
  virtual void execute(DiscreteElements* dem);
  virtual void execute(DiscreteElements* dem, pd::Peridynamics* pd)
  {
    std::cout << "**ERROR** Execute with DEM + PD. "
              << "Should not be called in DepositIntoContainer.\n";
  }
};
}

#endif