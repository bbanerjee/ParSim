#ifndef ELLIP3D_ISOTROPIC_LOADING_H
#define ELLIP3D_ISOTROPIC_LOADING_H

#include <Simulations/Command.h>
namespace dem {
class IsotropicLoading : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif