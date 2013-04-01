
#ifndef MATITI_SIMULATIONCONTROLLER_H
#define MATITI_SIMULATIONCONTROLLER_H

#include <Vaango/Core/ProblemSpec/ProblemSpecP.h>
#include <Vaango/Core/ProblemSpec/ProblemSpec.h>
#include <Core/SimulationState.h>


namespace Matiti {

  struct SimulationTime;

  //! The main component that controls the execution of the 
  //! entire simulation. 
  class SimulationController {

    public:
      SimulationController(Uintah::ProblemSpecP pspec);
      virtual ~SimulationController();

      //! Execute the simulation
      void run();

    private:

      void problemSetup();

      void calcStartTime   ( void );

      //! Set up, compile, and execute initial timestep
      void doInitialTimestep(DomainP& domain, double& t);

      void scheduleComputeStableTimestep(const DomainP& domain,
                                         SchedulerP&);

      void setStartSimTime ( double t );

      void calcWallTime ( void );

      void executeTimestep(double t, double& delt);

      Uintah::ProblemSpecP d_ups;
      SimulationTime*      d_timeinfo;
      SimulationState*     d_sharedState;

      // Do not allow copying of the controller
      SimulationController(const SimulationController&);
      SimulationController& operator=(const SimulationController&);
   };

} // End namespace Matiti

#endif
