#ifndef ELLIP3D_TRIAXIAL_LOADING_H
#define ELLIP3D_TRIAXIAL_LOADING_H

#include <Simulations/Command.h>
namespace dem {
class TriaxialLoading : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif