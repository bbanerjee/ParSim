#ifndef ELLIP3D_COUPLED_FLUID_FLOW_H
#define ELLIP3D_COUPLED_FLUID_FLOW_H

#include <Commands/Command.h>

namespace dem {
class CoupledFluidFlowCommand : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif