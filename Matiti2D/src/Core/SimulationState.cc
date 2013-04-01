#include <Core/SimulationState.h>

SimulationState::SimulationState(Uintah::ProblemSpecP )
{
  current_delt = 0.0;
  prev_delt = 0.0;
  simTime = 0;
  peri_matls = 0;
}

SimulationState::~SimulationState()
{
  peri_matls.clear(); 
}

