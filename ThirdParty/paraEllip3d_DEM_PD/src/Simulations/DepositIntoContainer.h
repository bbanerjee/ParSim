#ifndef ELLIP3D_DEPOSIT_INTO_CONTAINER_H
#define ELLIP3D_DEPOSIT_INTO_CONTAINER_H

#include <Simulations/Command.h>
namespace dem {
class DepositIntoContainer : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif