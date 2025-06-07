/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef VAANGO_CCA_COMPONENTS_SIMULATIONCONTROLLER_SIMULATIONCONTROLLER_H
#define VAANGO_CCA_COMPONENTS_SIMULATIONCONTROLLER_SIMULATIONCONTROLLER_H

#include <Core/Parallel/UintahParallelComponent.h>

#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Util/InfoMapper.h>
#include <Core/Util/Timers/Timers.hpp>

#include <CCA/Ports/DataWarehouseP.h>
#include <CCA/Ports/Scheduler.h>
#include <CCA/Ports/SchedulerP.h>

#include <CCA/Components/Schedulers/RuntimeStatsEnum.h>

// Window size for the overhead calculation
#define OVERHEAD_WINDOW 40

// Window size for the exponential moving average
#define AVERAGE_WINDOW 10

namespace Uintah {

class SimulationInterface;
class Output;
class LoadBalancer;
class Regridder;
class DataArchive;

/*! WallTimers: Utility class to manage the Wall Time */
class WallTimers
{
public:
  WallTimers()
  {
    d_num_samples = 0;
    d_wall_timer.start();
  };

public:
  Timers::Simple Timestep;         // Total time for all time steps
  Timers::Simple ExpMovingAverage; // Execution exponential moving average
                                   // for N time steps.
  Timers::Simple InSitu;           // In-situ time for previous time step

  auto
  getWindow() -> int
  {
    return AVERAGE_WINDOW;
  };
  void
  resetWindow()
  {
    d_num_samples = 0;
  };

  auto
  updateExpMovingAverage() -> Timers::nanoseconds
  {

    Timers::nanoseconds laptime = Timestep.lap();

    // Ignore the first sample as that is for initialization.
    if (d_num_samples) {
      // Calculate the exponential moving average for this time step.
      // Multiplier: (2 / (Time periods + 1) )
      // EMA: {current - EMA(previous)} x multiplier + EMA(previous).

      double mult =
        2.0 / ((double)std::min(d_num_samples, AVERAGE_WINDOW) + 1.0);

      ExpMovingAverage = mult * laptime + (1.0 - mult) * ExpMovingAverage();
    } else {
      ExpMovingAverage = laptime;
    }

    ++d_num_samples;

    return laptime;

  } // end Timers::nanoseconds

  auto
  GetWallTime() -> double
  {
    return d_wall_timer().seconds();
  };

private:
  int d_num_samples{ 0 }; // Number of samples for the moving average
  Timers::Simple d_wall_timer{};
};

/*! SimulationController: Abstract baseclass for the SimulationControllers. */
//! This is the main component that controls the execution of the
//! entire simulation.
class SimulationController : public UintahParallelComponent
{

public:
  SimulationController(const ProcessorGroup* myworld,
                       ProblemSpecP pspec,
                       const std::string& input_ups_dir);

  ~SimulationController() override = default;

  // eliminate copy, assignment and move
  SimulationController(const SimulationController&) = delete;
  SimulationController(SimulationController&&)      = delete;

  auto
  operator=(const SimulationController&) -> SimulationController& = delete;
  auto
  operator=(SimulationController&&) -> SimulationController& = delete;

  // Methods for managing the components attached via the ports.
  void
  setComponents([[maybe_unused]] UintahParallelComponent* comp) override{};

  void
  getComponents() override;

  void
  releaseComponents() override;

  //! Notifies (before calling run) the SimulationController
  //! that this is simulation is a restart.
  void
  doRestart(const std::string& restartFromDir,
            int index,
            bool fromScratch,
            bool removeOldDir);

  //! Execute the simulation
  virtual void
  run() = 0;

  //! Set simulation controller flags
  void
  setPostProcessFlags();

  auto
  getProblemSpecP() -> ProblemSpecP
  {
    return d_ups;
  }

  auto
  getGridProblemSpecP() -> ProblemSpecP
  {
    return d_grid_ps;
  }

  auto
  getSchedulerP() -> SchedulerP
  {
    return d_scheduler;
  }

  auto
  getLoadBalancer() -> LoadBalancer*
  {
    return d_loadBalancer;
  }

  auto
  getOutput() -> Output*
  {
    return d_output;
  }

  auto
  getSimulationInterface() -> SimulationInterface*
  {
    return d_simulator;
  }

  auto
  getRegridder() -> Regridder*
  {
    return d_regridder;
  }

  auto
  getWallTimers() -> WallTimers*
  {
    return &d_wall_timers;
  }

  [[nodiscard]] auto
  getRecompileTaskGraph() const -> bool
  {
    return d_recompile_taskgraph;
  }

  void
  setRecompileTaskGraph(bool val)
  {
    d_recompile_taskgraph = val;
  }

  void
  ScheduleReportStats(bool header);

  void
  ReportStats(const ProcessorGroup*,
              const PatchSubset*,
              const MaterialSubset*,
              DataWarehouse*,
              DataWarehouse*,
              bool header);

  auto
  getRuntimeStats() -> ReductionInfoMapper<RuntimeStatsEnum, double>&
  {
    return d_runtime_stats;
  };

protected:
  void
  restartArchiveSetup();

  void
  outputSetup();

  void
  gridSetup();

  void
  regridderSetup();

  void
  loadBalancerSetup();

  void
  simulatorSetup();

  void
  schedulerSetup();

  void
  timeStateSetup();

  void
  finalSetup();

  void
  resetStats();

  void
  getMemoryStats(bool create = false);

  ProblemSpecP d_ups{ nullptr };
  ProblemSpecP d_grid_ps{ nullptr };    // Problem Spec for the Grid
  ProblemSpecP d_restart_ps{ nullptr }; // Problem Spec for the Grid

  GridP d_current_gridP{ nullptr };

  const std::string d_input_ups_dir{ "" };

  SimulationInterface* d_simulator{ nullptr };
  SchedulerP d_scheduler{ nullptr };
  Output* d_output{ nullptr };
  LoadBalancer* d_loadBalancer{ nullptr };
  Regridder* d_regridder{ nullptr };

  // Only used when restarting: Data from checkpoint UDA.
  std::shared_ptr<DataArchive> d_restart_archive{ nullptr };

  bool d_do_multi_taskgraphing{ false };

  WallTimers d_wall_timers;

  // Used when restarting.
  bool d_restarting{ false };
  std::string d_from_dir;
  int d_restart_timestep{ 0 };
  int d_restart_index{ 0 };
  int d_last_recompile_timeStep{ 0 };
  bool d_post_process_uda{ false };

  // d_restartFromScratch is true then don't copy or move any of
  // the old timesteps or dat files from the old directory.  Run as
  // as if it were running from scratch but with initial conditions
  // given by the restart checkpoint.
  bool d_restart_from_scratch{ false };

  // If !restartFromScratch, then this indicates whether to move
  // or copy the old timesteps.
  bool d_restart_remove_old_dir{ false };

  bool d_recompile_taskgraph{ false };

  // Runtime stat mappers.
  ReductionInfoMapper<RuntimeStatsEnum, double> d_runtime_stats;

private:
  // For reporting stats Frequency > OnTimestep
  unsigned int d_reportStatsFrequency{ 1 };
  unsigned int d_reportStatsOnTimestep{ 0 };

  // Percent time in overhead samples
  double d_overhead_values[OVERHEAD_WINDOW];
  double d_overhead_weights[OVERHEAD_WINDOW];
  int d_overhead_index{ 0 }; // Next sample for writing
  int d_num_samples{ 0 };
};

} // End namespace Uintah

#endif
