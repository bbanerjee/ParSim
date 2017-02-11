#ifndef ELLIP3D_DEPOSIT_INTO_CONTAINER_RESUME_H
#define ELLIP3D_DEPOSIT_INTO_CONTAINER_RESUME_H

#include <Commands/Command.h>
namespace dem {
class DepositIntoContainerResumeCommand : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif