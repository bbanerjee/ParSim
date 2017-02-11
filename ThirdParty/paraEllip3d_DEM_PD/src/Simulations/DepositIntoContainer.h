#ifndef ELLIP3D_DEPOSIT_INTO_CONTAINER_H
#define ELLIP3D_DEPOSIT_INTO_CONTAINER_H

#include <Commands/Command.h>
namespace dem {
class DepositIntoContainerCommand : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif