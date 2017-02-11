#ifndef ELLIP3D_TUNE_MASS_PERCENTAGE_H
#define ELLIP3D_TUNE_MASS_PERCENTAGE_H

#include <Simulations/Command.h>
namespace dem {
class TuneMassPercentage : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif