#ifndef ELLIP3D_TRUE_TRIAXIAL_LOADING_H
#define ELLIP3D_TRUE_TRIAXIAL_LOADING_H

#include <Simulations/Command.h>
namespace dem {
class TrueTriaxialLoading : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif