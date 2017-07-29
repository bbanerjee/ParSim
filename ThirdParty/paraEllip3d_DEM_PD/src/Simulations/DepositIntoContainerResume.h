#ifndef ELLIP3D_DEPOSIT_INTO_CONTAINER_RESUME_H
#define ELLIP3D_DEPOSIT_INTO_CONTAINER_RESUME_H

#include <Simulations/Command.h>
#include <ostream>

namespace dem {
class DepositIntoContainerResume : public Command
{
public:
  virtual void execute(DiscreteElements* dem);
  virtual void execute(DiscreteElements* dem, pd::Peridynamics* pd)
  {
    std::cout << "**ERROR** Execute with DEM + PD. "
              << "Should not be called in DepositIntoContainerResume.\n";
  }
};
}

#endif