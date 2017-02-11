#ifndef ELLIP3D_CAVITY_EXPANSION_RESUME_H
#define ELLIP3D_CAVITY_EXPANSION_RESUME_H

#include <Commands/Command.h>

namespace dem {

class CavityExpansionResumeCommand : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif