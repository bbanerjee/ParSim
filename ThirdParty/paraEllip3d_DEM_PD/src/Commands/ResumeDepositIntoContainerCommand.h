#ifndef ELLIP3D_RESUME_DEPOSIT_INTO_CONTAINER_H
#define ELLIP3D_RESUME_DEPOSIT_INTO_CONTAINER_H

#include <Commands/Command.h>
namespace dem {
class ResumeDepositIntoContainerCommand : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif