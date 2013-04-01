#ifndef MATITI_SIMULATIONINTERFACE_H
#define MATITI_SIMULATIONINTERFACE_H

#include <Vaango/Core/ProblemSpec/ProblemSpecP.h>
#include <Core/SimulationStateP.h>

namespace Matiti {

  class SimulationInterface {

    public:
      SimulationInterface();
      virtual ~SimulationInterface();

      virtual void problemSetup(const ProblemSpecP& ps,
		                SimulationStateP& state) = 0;

    private:

      // Do not allow copying of the interface
      SimulationInterface(const SimulationInterface&);
      SimulationInterface& operator=(const SimulationInterface&);
   };

} // End namespace Matiti

#endif
