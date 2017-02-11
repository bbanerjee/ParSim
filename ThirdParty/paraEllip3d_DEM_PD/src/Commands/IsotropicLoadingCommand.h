#ifndef ELLIP3D_ISOTROPIC_LOADING_H
#define ELLIP3D_ISOTROPIC_LOADING_H

#include <Commands/Command.h>
namespace dem {
class IsotropicLoadingCommand : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif