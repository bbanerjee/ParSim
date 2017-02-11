#ifndef ELLIP3D_CAVITY_EXPANSION_H
#define ELLIP3D_CAVITY_EXPANSION_H

#include <Commands/Command.h>

namespace dem {

class CavityExpansionCommand : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif