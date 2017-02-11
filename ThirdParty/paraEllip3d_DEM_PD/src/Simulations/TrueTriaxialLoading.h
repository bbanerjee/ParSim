#ifndef ELLIP3D_TRUE_TRIAXIAL_LOADING_H
#define ELLIP3D_TRUE_TRIAXIAL_LOADING_H

#include <Commands/Command.h>
namespace dem {
class TrueTriaxialLoadingCommand : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif