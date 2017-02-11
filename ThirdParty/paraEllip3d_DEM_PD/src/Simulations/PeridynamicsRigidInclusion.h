#ifndef ELLIP3D_PERI_RIGID_INCLUSION_H
#define ELLIP3D_PERI_RIGID_INCLUSION_H

#include <Simulations/Command.h>
namespace dem {
class PeridynamicsRigidInclusion : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif