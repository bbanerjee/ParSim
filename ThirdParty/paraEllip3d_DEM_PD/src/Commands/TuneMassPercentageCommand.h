#ifndef ELLIP3D_TUNE_MASS_PERCENTAGE_H
#define ELLIP3D_TUNE_MASS_PERCENTAGE_H

#include <Commands/Command.h>
namespace dem {
class TuneMassPercentageCommand : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif