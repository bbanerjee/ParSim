#include <Commands/Command.h>
#include <Commands/CommandHandler.h>
#include <Commands/ProceedFromPresetCommand.h>

using namespace dem;

CommandP
CommandHandler::handleCommand(int simuType)
{

  switch (simuType) {
    case 001: // proceed from preset state
      return std::make_unique<ProceedFromPresetCommand>();
      break;
    default:
      return nullptr;
  }
  return nullptr;
}