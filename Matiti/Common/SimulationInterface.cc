#include <Common/SimulationInterface.h>
#include <Common/Task.h>

#include <Core/Exceptions/InternalError.h>

using namespace Matiti;

SimulationInterface::SimulationInterface()
{
}

SimulationInterface::~SimulationInterface()
{
}

void
SimulationInterface::scheduleTimeAdvance(SchedulerP&)
{
  throw Uintah::InternalError("no simulation implemented?", __FILE__, __LINE__);
}

double
SimulationInterface::recomputeTimestep(double)
{
  throw Uintah::InternalError("recomputeTimestep not implemented for this component", __FILE__, __LINE__);
}

bool
SimulationInterface::restartableTimesteps()
{
  return false;
}

