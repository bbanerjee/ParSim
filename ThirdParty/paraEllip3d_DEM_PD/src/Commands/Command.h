#ifndef ELLIP3D_COMMAND_H
#define ELLIP3D_COMMAND_H

#include <DiscreteElements/Assembly.h>

namespace dem {

class Command
{
public:
  virtual ~Command(){};
  virtual void execute(Assembly* assembly) = 0;
};
}

#endif