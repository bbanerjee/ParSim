#ifndef ELLIP3D_COMMAND_HANDLER_H
#define ELLIP3D_COMMAND_HANDLER_H

#include <memory>

namespace dem {
using CommandP = std::unique_ptr<Command>;

class Command;

class CommandHandler
{
public:
  CommandP handleCommand(int simuType);
};
}

#endif