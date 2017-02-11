#ifndef ELLIP3D_PLANESTRAIN_LOADING_H
#define ELLIP3D_PLANESTRAIN_LOADING_H

#include <Commands/Command.h>
namespace dem {
class PlaneStrainLoadingCommand : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif