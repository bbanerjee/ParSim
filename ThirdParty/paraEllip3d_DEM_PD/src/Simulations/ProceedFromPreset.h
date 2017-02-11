#ifndef ELLIP3D_PROCEED_FROM_PRESET_COMMAND_H
#define ELLIP3D_PROCEED_FROM_PRESET_COMMAND_H

#include <Commands/Command.h>
namespace dem {
class ProceedFromPresetCommand : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif