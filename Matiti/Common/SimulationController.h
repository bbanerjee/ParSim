
#ifndef MATITI_SIMULATIONCONTROLLER_H
#define MATITI_SIMULATIONCONTROLLER_H

#include <Core/Parallel/SerialComponent.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <Core/Mesh/MeshP.h>
#include <Core/Mesh/LevelP.h>
#include <Core/Mesh/SimulationStateP.h>
#include <Core/Mesh/SimulationState.h>
#include <CCA/Ports/SchedulerP.h>
#include <CCA/Ports/Scheduler.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/ProblemSpec/ProblemSpec.h>


namespace Matiti {

  class SimulationInterface;
  class Output;
  struct SimulationTime;
  class DataArchive;


  //! The main component that controls the execution of the 
  //! entire simulation. 
  class SimulationController : public SerialComponent {

    public:
      SimulationController(Uintah::ProblemSpecP pspec);
      virtual ~SimulationController();

      //! Notifies (before calling run) the SimulationController
      //! that this is simulation is a restart.
      void doRestart(std::string restartFromDir, int timestep);

      //! Execute the simulation
      void run();

    private:

      void preMeshSetup();

      MeshP meshSetup();

      void postMeshSetup(MeshP& mesh, double& t);

      void calcStartTime   ( void );

      //! Set up, compile, and execute initial timestep
      void doInitialTimestep(MeshP& mesh, double& t);

      void scheduleComputeStableTimestep(const MeshP& mesh,
                                         SchedulerP&);

      void setStartSimTime ( double t );

      void initSimulationStatsVars ( void );

      //! Asks a variety of components if one of them needs the taskgraph
      //! to recompile.
      bool needRecompile(double t, double delt, const MeshP& level);

      void recompile(double t, double delt, MeshP& currentMesh);

      void calcWallTime ( void );

      void printSimulationStats ( int timestep, double delt, double time );

      void executeTimestep(double t, double& delt, MeshP& currentMesh);

      //! adjust delt based on timeinfo and other parameters
      //    'first' is whether this is the first time adjustDelT is called.
      void adjustDelT(double& delt, double prev_delt, bool first, double t);

      double getWallTime     ( void );

      double getStartTime    ( void );

      ProblemSpecP         d_ups;
      ProblemSpecP         d_mesh_ps;         // Problem Spec for the Mesh
      SimulationStateP     d_sharedState;
      SchedulerP           d_scheduler;
      Output*              d_output;
      SimulationTime*      d_timeinfo;
      SimulationInterface* d_sim;
      DataArchive*         d_archive;

      /* for restarting */
      bool d_restarting;
      std::string d_fromDir;
      int d_restartTimestep;
      int d_restartIndex;

      int d_lastRecompileTimestep;

      // If d_restartFromScratch is true then don't copy or move any of
      // the old timesteps or dat files from the old directory.  Run as
      // as if it were running from scratch but with initial conditions
      // given by the restart checkpoint.
      bool d_restartFromScratch;

      int    d_n;
      double d_wallTime;              // current wall time
      double d_startTime;             // starting wall time
      double d_startSimTime;          // starting sim time
      double d_prevWallTime;
     
      // this is for calculating an exponential moving average
      double d_movingAverage;

      // Do not allow copying of the controller
      SimulationController(const SimulationController&);
      SimulationController& operator=(const SimulationController&);
   };

} // End namespace Matiti

#endif
