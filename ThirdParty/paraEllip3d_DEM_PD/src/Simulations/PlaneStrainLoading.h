#ifndef ELLIP3D_PLANESTRAIN_LOADING_H
#define ELLIP3D_PLANESTRAIN_LOADING_H

#include <Simulations/Command.h>
namespace dem {
class PlaneStrainLoading : public Command
{
public:
  virtual void execute(Assembly* assembly);
};
}

#endif