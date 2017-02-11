#ifndef ELLIP3D_TRIM_PARTICLES_H
#define ELLIP3D_TRIM_PARTICLES_H

#include <Commands/Command.h>
namespace dem {
class TrimParticlesCommand : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif