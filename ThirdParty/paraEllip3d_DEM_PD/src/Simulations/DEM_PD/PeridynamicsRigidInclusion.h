#ifndef ELLIP3D_PERI_RIGID_INCLUSION_H
#define ELLIP3D_PERI_RIGID_INCLUSION_H

#include <Simulations/Command.h>
#include <ostream>

namespace dem {
class PeridynamicsRigidInclusion : public Command
{
public:
  virtual void execute(DiscreteElements* dem)
  {
    std::cout << "**ERROR** Execute with DEM only. "
              << "Should not be called in PeridynamicsRigidInclusion.\n";
  }
  virtual void execute(DiscreteElements* dem, pd::Peridynamics* pd);
};
}

#endif