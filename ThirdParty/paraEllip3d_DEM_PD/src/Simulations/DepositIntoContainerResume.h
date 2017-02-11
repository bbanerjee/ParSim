#ifndef ELLIP3D_DEPOSIT_INTO_CONTAINER_RESUME_H
#define ELLIP3D_DEPOSIT_INTO_CONTAINER_RESUME_H

#include <Simulations/Command.h>
namespace dem {
class DepositIntoContainerResume : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif