#ifndef MATITI_OUTPUT_H
#define MATITI_OUTPUT_H

#include <Common/SerialPort.h>
#include <Common/SchedulerP.h>
#include <Mesh/MeshP.h>

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/OS/Dir.h>

#include <string>


namespace Matiti {

  class SimulationState;

  class Output : public SerialPort {

    public:
      Output();
      virtual ~Output();
      
      virtual void problemSetup(const Uintah::ProblemSpecP& params,
                                SimulationState* state) = 0;

      virtual void initializeOutput(const Uintah::ProblemSpecP& params) = 0;

      //////////
      // Call this when restarting from a checkpoint after calling
      // problemSetup.
      virtual void restartSetup(SCIRun::Dir& restartFromDir, int startTimestep,
                                int timestep, double time, bool fromScratch,
                                bool removeOldDir) = 0;

      virtual bool needRecompile(double time, double delt,
                                 const MeshP& grid) = 0;

      //////////
      // Call this after all other tasks have been added to the scheduler
      virtual void finalizeTimestep(double t, double delt, const MeshP&,
                                    SchedulerP&, bool recompile = false,
                                    int addMaterial = 0) = 0;

      //////////
      // Call this after a timestep restart to make sure we still
      // have an output timestep
      virtual void reEvaluateOutputTimestep(double old_delt, double new_delt)=0;

      //////////
      // Call this after the timestep has been executed.
      virtual void executedTimestep(double delt, const MeshP&) = 0;
     
      virtual const std::string getOutputLocation() const = 0;

      virtual int getCurrentTimestep() = 0;

      virtual double getCurrentTime() = 0;

      // Get the time the next output will occur
      virtual double getNextOutputTime() = 0;

      // Get the timestep the next output will occur
      virtual int getNextOutputTimestep() = 0;

      // Get the time the next checkpoint will occur
      virtual double getNextCheckpointTime() = 0;

      // Get the timestep the next checkpoint will occur
      virtual int getNextCheckpointTimestep() = 0;
      
      // Returns true if data will be output this timestep
      virtual bool isOutputTimestep() = 0;

      // Returns true if data will be checkpointed this timestep
      virtual bool isCheckpointTimestep() = 0;

      // Returns true if the label is being saved
      virtual bool isLabelSaved(std::string label) = 0;

      //////////
      // Get the directory of the current time step for outputting info.
      virtual const std::string& getLastTimestepOutputLocation() const = 0;

    private:
      Output(const Output&);
      Output& operator=(const Output&);
  };

} // End namespace Matiti

#endif
