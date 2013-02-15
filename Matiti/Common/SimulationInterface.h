#ifndef MATITI_SIMULATIONINTERFACE_H
#define MATITI_SIMULATIONINTERFACE_H

#include <Common/SerialPort.h>
#include <Common/SchedulerP.h>
#include <Common/Handle.h>
#include <Mesh/DomainP.h>
#include <Mesh/SimulationStateP.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/OS/Dir.h>

namespace Uintah {

  class DataWarehouse;

  class SimulationInterface : public SerialPort {

    public:
      SimulationInterface();
      virtual ~SimulationInterface();
      
      virtual void problemSetup(const Uintah::ProblemSpecP& params, 
                                const Uintah::ProblemSpecP& restart_prob_spec,
                                DomainP& grid, SimulationStateP& state) = 0;

      virtual void outputProblemSpec(Uintah::ProblemSpecP& ps) {}
      virtual void outputPS(SCIRun::Dir& dir) {}
      
      virtual void scheduleInitialize(SchedulerP&) = 0;

      //////////
      // restartInitialize() is called once and only once if and when a simulation is restarted.
      // This allows the simulation component to handle initializations that are necessary when
      // a simulation is restarted.
      // 
      virtual void restartInitialize() {}

      virtual void switchInitialize(SchedulerP&) {}
      
      virtual void scheduleComputeStableTimestep(SchedulerP&) = 0;
      
      virtual void scheduleTimeAdvance(SchedulerP&);

      // this is for wrapping up a timestep when it can't be done in scheduleTimeAdvance.
      virtual void scheduleFinalizeTimestep(SchedulerP&) {}

      // Redo a timestep if current time advance is not converging.
      // Returned time is the new dt to use.
      virtual double recomputeTimestep(double delt);
      virtual bool restartableTimesteps();

       //////////
       // ask the component if it needs to be recompiled
       virtual bool needRecompile(double /*time*/, double /*dt*/,
                                  const DomainP& /*grid*/) {return false;}

     private:
       SimulationInterface(const SimulationInterface&);
       SimulationInterface& operator=(const SimulationInterface&);
   };
} // End namespace Uintah
   


#endif
