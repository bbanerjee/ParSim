#include <Core/SimulationTime.h>

using namespace Matiti;

SimulationTime:: SimulationTime(const Vaango::ProblemSpecP& ps)
{
  problemSetup(ps);
}

SimulationTime::~SimulationTime()
{
}

void 
SimulationTime::problemSetup(const Vaango::ProblemSpec& ps)
{
  Vaango::ProblemSpecP time_ps = ps->findBlock("Time");
  time_ps->require("initial_delta_t", maxInitialDeltaT);
  time_ps->require("min_delta_t", minDeltaT);
  time_ps->require("max_delta_t", maxDeltaT);
  time_ps->require("timestep_multiplier", timestepMultiplier);

  initTime = 0.0;
  time_ps->get("initial_time", initTime);
  
  time_ps->require("max_time", maxTime);
  maxTimeStep = 1000000;
  time_ps->get("max_timesteps", maxTimestep);
}

