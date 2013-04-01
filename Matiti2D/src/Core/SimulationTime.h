
#ifndef MATITI_SIMULATIONTIME_H
#define MATITI_SIMULATIONTIME_H

#include <Vaango/Core/ProblemSpec/ProblemSpecP.h>

namespace Matiti {

  class SimulationTime {

    public:

      double maxTime;
      double initTime;
      double maxDeltaT;
      double minDeltaT;
      double maxInitialDeltaT;
      double timestepMultiplier;

      int maxTimestep;

    public:
      SimulationTime(const Vaango::ProblemSpecP& ps);
      virtual ~SimulationTime();

      void problemSetup(const Vaango::ProblemSpec& ps);

    private:

      // Do not allow copying of the structure
      SimulationTime(const SimulationTime&);
      SimulationTime& operator=(const SimulationTime&);
   };

} // End namespace Matiti

#endif
