#ifndef ELLIP3D_OEDOMETER_LOADING_H
#define ELLIP3D_OEDOMETER_LOADING_H

#include <Commands/Command.h>
namespace dem {
class OedometerLoadingCommand : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif