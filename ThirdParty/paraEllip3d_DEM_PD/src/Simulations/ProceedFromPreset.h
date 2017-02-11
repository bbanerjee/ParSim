#ifndef ELLIP3D_PROCEED_FROM_PRESET_COMMAND_H
#define ELLIP3D_PROCEED_FROM_PRESET_COMMAND_H

#include <Simulations/Command.h>
namespace dem {
class ProceedFromPreset : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif