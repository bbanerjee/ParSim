#ifndef ELLIP3D_TRIM_PARTICLES_H
#define ELLIP3D_TRIM_PARTICLES_H

#include <Simulations/Command.h>
namespace dem {
class TrimParticles : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif