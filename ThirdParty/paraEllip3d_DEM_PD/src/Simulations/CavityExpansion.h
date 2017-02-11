#ifndef ELLIP3D_CAVITY_EXPANSION_H
#define ELLIP3D_CAVITY_EXPANSION_H

#include <Simulations/Command.h>

namespace dem {

class CavityExpansion : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif