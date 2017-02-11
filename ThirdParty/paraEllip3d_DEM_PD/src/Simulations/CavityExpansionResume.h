#ifndef ELLIP3D_CAVITY_EXPANSION_RESUME_H
#define ELLIP3D_CAVITY_EXPANSION_RESUME_H

#include <Simulations/Command.h>

namespace dem {

class CavityExpansionResume : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif