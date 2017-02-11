#ifndef ELLIP3D_TRIAXIAL_LOADING_H
#define ELLIP3D_TRIAXIAL_LOADING_H

#include <Commands/Command.h>
namespace dem {
class TriaxialLoadingCommand : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif