#ifndef ELLIP3D_OEDOMETER_LOADING_H
#define ELLIP3D_OEDOMETER_LOADING_H

#include <Simulations/Command.h>
namespace dem {
class OedometerLoading : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif