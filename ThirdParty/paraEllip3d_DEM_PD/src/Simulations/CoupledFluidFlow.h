#ifndef ELLIP3D_COUPLED_FLUID_FLOW_H
#define ELLIP3D_COUPLED_FLUID_FLOW_H

#include <Simulations/Command.h>

namespace dem {
class CoupledFluidFlow : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif