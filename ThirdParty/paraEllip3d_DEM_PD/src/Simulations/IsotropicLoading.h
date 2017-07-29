#ifndef ELLIP3D_ISOTROPIC_LOADING_H
#define ELLIP3D_ISOTROPIC_LOADING_H

#include <Simulations/Command.h>
#include <ostream>

namespace dem {
class IsotropicLoading : public Command
{
public:
  virtual void execute(DiscreteElements* dem);
  virtual void execute(DiscreteElements* dem, pd::Peridynamics* pd)
  {
    std::cout << "**ERROR** Execute with DEM + PD. "
              << "Should not be called in IsotropicLoading.\n";
  }
};
}

#endif