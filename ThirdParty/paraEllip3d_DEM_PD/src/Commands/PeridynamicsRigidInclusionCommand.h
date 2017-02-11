#ifndef ELLIP3D_PERI_RIGID_INCLUSION_H
#define ELLIP3D_PERI_RIGID_INCLUSION_H

#include <Commands/Command.h>
namespace dem {
class PeridynamicsRigidInclusionCommand : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif