#ifndef ELLIP3D_COMMAND_H
#define ELLIP3D_COMMAND_H

#include <DiscreteElements/DiscreteElements.h>
#include <Peridynamics/Peridynamics.h>
#include <SmoothParticleHydro/SmoothParticleHydro.h>

namespace dem {

class Command
{
public:
  virtual ~Command(){};
  virtual void execute(DiscreteElements* dem) = 0;
  virtual void execute(DiscreteElements* dem, pd::Peridynamics* pd) = 0;
  virtual void execute(DiscreteElements* dem, sph::SmoothParticleHydro* pd) = 0;
};
}

#endif