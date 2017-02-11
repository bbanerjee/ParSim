#ifndef ELLIP3D_PERI_PULLOUT_H
#define ELLIP3D_PERI_PULLOUT_H

#include <Commands/Command.h>
namespace dem {
class PeridynamicsPullOutCommand : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif