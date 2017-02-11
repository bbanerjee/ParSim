#ifndef ELLIP3D_PERI_PULLOUT_H
#define ELLIP3D_PERI_PULLOUT_H

#include <Simulations/Command.h>
namespace dem {
class PeridynamicsPullOut : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif