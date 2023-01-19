/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#include <CCA/Components/SimulationController/SimulationController.h>

#include <CCA/Components/Schedulers/MPIScheduler.h>

#include <Core/DataArchive/DataArchive.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/OS/Dir.h>
#include <Core/OS/ProcessInfo.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Util/DOUT.hpp>
#include <Core/Util/Timers/Timers.hpp>

#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Output.h>
#include <CCA/Ports/ProblemSpecInterface.h>
#include <CCA/Ports/Regridder.h>
#include <CCA/Ports/Scheduler.h>
#include <CCA/Ports/SimulationInterface.h>

#include <sci_defs/malloc_defs.h>

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sys/param.h>
#include <vector>

#include <pwd.h>

namespace {

Uintah::Dout g_sim_stats("SimulationStats",
                         "SimulationController",
                         "Simulation general stats",
                         true);
Uintah::Dout g_sim_stats_mem("SimulationStatsMem",
                             "SimulationController",
                             "Simulation memory stats",
                             true);

Uintah::Dout g_comp_stats("ComponentStats",
                          "SimulationController",
                          "Aggregated component stats",
                          false);
Uintah::Dout g_comp_node_stats("ComponentNodeStats",
                               "SimulationController",
                               "Aggregated node component stats",
                               false);
Uintah::Dout g_comp_indv_stats("ComponentIndividualStats",
                               "SimulationController",
                               "Individual component stats",
                               false);

Uintah::Dout g_app_stats("SimulationStats",
                         "SimulationController",
                         "Aggregated simulation stats",
                         false);
Uintah::Dout g_app_node_stats("SimulationNodeStats",
                              "SimulationController",
                              "Aggregated node simulation stats",
                              false);
Uintah::Dout g_app_indv_stats("SimulationIndividualStats",
                              "SimulationController",
                              "Individual simulation stats",
                              false);

}
namespace Uintah {

SimulationController::SimulationController(const ProcessorGroup* myworld,
                                           ProblemSpecP ps)
  : UintahParallelComponent(myworld)
  , d_ups(ps)
{
  // initialize the overhead percentage
  for (int i = 0; i < OVERHEAD_WINDOW; ++i) {
    double x              = (double)i / (double)(OVERHEAD_WINDOW / 2);
    d_overhead_values[i]  = 0;
    d_overhead_weights[i] = 8.0 - x * x * x;
  }

  d_grid_ps = d_ups->findBlock("Grid");

  ProblemSpecP simController_ps = d_ups->findBlock("SimulationController");

  if (simController_ps) {
    ProblemSpecP runtimeStats_ps = simController_ps->findBlock("RuntimeStats");

    if (runtimeStats_ps) {
      runtimeStats_ps->get("frequency", d_reportStatsFrequency);
      runtimeStats_ps->get("onTimeStep", d_reportStatsOnTimeStep);

      if (d_reportStatsOnTimeStep >= d_reportStatsFrequency) {
        proc0cout << "Error: the frequency of reporting the run time stats "
                  << d_reportStatsFrequency << " "
                  << "is less than or equal to the time step ordinality "
                  << d_reportStatsOnTimeStep << " "
                  << ". Resetting the ordinality to ";

        if (d_reportStatsFrequency > 1) {
          d_reportStatsOnTimeStep = 1;
        } else {
          d_reportStatsOnTimeStep = 0;
        }
        proc0cout << d_reportStatsOnTimeStep << std::endl;
      }
    }
  }

  std::string timeStr("seconds");
  std::string bytesStr("MBytes");

  d_runtime_stats.insert(CompilationTime, std::string("Compilation"), timeStr);
  d_runtime_stats.insert(RegriddingTime, std::string("Regridding"), timeStr);
  d_runtime_stats.insert(RegriddingCompilationTime,
                         std::string("RegriddingCompilation"),
                         timeStr);
  d_runtime_stats.insert(RegriddingCopyDataTime,
                         std::string("RegriddingCopyData"),
                         timeStr);
  d_runtime_stats.insert(LoadBalancerTime,
                         std::string("LoadBalancer"),
                         timeStr);

  d_runtime_stats.insert(TaskExecTime, std::string("TaskExec"), timeStr);
  d_runtime_stats.insert(TaskLocalCommTime,
                         std::string("TaskLocalComm"),
                         timeStr);
  d_runtime_stats.insert(TaskWaitCommTime,
                         std::string("TaskWaitCommTime"),
                         timeStr);
  d_runtime_stats.insert(TaskReduceCommTime,
                         std::string("TaskReduceCommTime"),
                         timeStr);
  d_runtime_stats.insert(TaskWaitThreadTime,
                         std::string("TaskWaitThread"),
                         timeStr);

  d_runtime_stats.insert(XMLIOTime, std::string("XMLIO"), timeStr);
  d_runtime_stats.insert(OutputIOTime, std::string("OutputIO"), timeStr);
  d_runtime_stats.insert(OutputGlobalIOTime,
                         std::string("OutputGlobalIO"),
                         timeStr);
  d_runtime_stats.insert(CheckpointIOTime,
                         std::string("CheckpointIO"),
                         timeStr);
  d_runtime_stats.insert(CheckpointGlobalIOTime,
                         std::string("CheckpointGlobalIO"),
                         timeStr);
  d_runtime_stats.insert(TotalIOTime, std::string("TotalIO"), timeStr);

  d_runtime_stats.insert(OutputIORate,
                         std::string("OutputIORate"),
                         "MBytes/sec");
  d_runtime_stats.insert(OutputGlobalIORate,
                         std::string("OutputGlobalIORate"),
                         "MBytes/sec");
  d_runtime_stats.insert(CheckpointIORate,
                         std::string("CheckpointIORate"),
                         "MBytes/sec");
  d_runtime_stats.insert(CheckpointGlobalIORate,
                         std::string("CheckpointGlobalIORate"),
                         "MBytes/sec");

  d_runtime_stats.insert(NumTasks, std::string("NumberOfTasks"), "tasks");
  d_runtime_stats.insert(NumPatches, std::string("NumberOfPatches"), "patches");
  d_runtime_stats.insert(NumCells, std::string("NumberOfCells"), "cells");
  d_runtime_stats.insert(NumParticles,
                         std::string("NumberOfParticles"),
                         "paticles");

  d_runtime_stats.insert(SCIMemoryUsed, std::string("SCIMemoryUsed"), bytesStr);
  d_runtime_stats.insert(SCIMemoryMaxUsed,
                         std::string("SCIMemoryMaxUsed"),
                         bytesStr);
  d_runtime_stats.insert(SCIMemoryHighwater,
                         std::string("SCIMemoryHighwater"),
                         bytesStr);
  d_runtime_stats.insert(MemoryUsed, std::string("MemoryUsed"), bytesStr);
  d_runtime_stats.insert(MemoryResident,
                         std::string("MemoryResident"),
                         bytesStr);

  d_runtime_stats.calculateRankMinimum(true);
  d_runtime_stats.calculateRankStdDev(true);

  if (g_comp_node_stats) {
    d_runtime_stats.calculateNodeSum(true);
    d_runtime_stats.calculateNodeMinimum(true);
    d_runtime_stats.calculateNodeAverage(true);
    d_runtime_stats.calculateNodeMaximum(true);
    d_runtime_stats.calculateNodeStdDev(true);
  }

} // end SimulationController constructor

void
SimulationController::setPostProcessFlags()
{
  d_post_process_uda = true;
}

void
SimulationController::getComponents(void)
{
  d_simulator = dynamic_cast<SimulationInterface*>(getPort("simulator"));
  if (!d_simulator) {
    throw InternalError("dynamic_cast of 'd_simulator' failed!",
                        __FILE__,
                        __LINE__);
  }

  d_loadBalancer = dynamic_cast<LoadBalancer*>(getPort("load balancer"));
  if (!d_loadBalancer) {
    throw InternalError("dynamic_cast of 'd_loadBalancer' failed!",
                        __FILE__,
                        __LINE__);
  }

  d_output = dynamic_cast<Output*>(getPort("output"));
  if (!d_output) {
    throw InternalError("dynamic_cast of 'd_output' failed!",
                        __FILE__,
                        __LINE__);
  }

  d_regridder = dynamic_cast<Regridder*>(getPort("regridder"));
  if (d_simulator->isDynamicRegridding() && !d_regridder) {
    throw InternalError("dynamic_cast of 'd_regridder' failed!",
                        __FILE__,
                        __LINE__);
  }

  d_scheduler = dynamic_cast<Scheduler*>(getPort("scheduler"));
  if (!d_scheduler) {
    throw InternalError("dynamic_cast of 'd_scheduler' failed!",
                        __FILE__,
                        __LINE__);
  }
}

void
SimulationController::releaseComponents(void)
{
  releasePort("simulator");
  releasePort("load balancer");
  releasePort("output");
  releasePort("regridder");
  releasePort("scheduler");

  d_simulator    = nullptr;
  d_loadBalancer = nullptr;
  d_output       = nullptr;
  d_regridder    = nullptr;
  d_scheduler    = nullptr;
}

void
SimulationController::doRestart(const std::string& restartFromDir,
                                int index,
                                bool fromScratch,
                                bool removeOldDir)
{
  d_restarting             = true;
  d_from_dir               = restartFromDir;
  d_restart_index          = index;
  d_restart_from_scratch   = fromScratch;
  d_restart_remove_old_dir = removeOldDir;
}

void
SimulationController::restartArchiveSetup()
{
  // Set up the restart archive now as it is need by the output.
  if (d_restarting) {
    // Create the DataArchive here, and store it, as it is used a few
    // times. The grid needs to be read before the ProblemSetup, and
    // not all of the data can be read until after ProblemSetup, so
    // DataArchive operations are needed.

    Dir restartFromDir(d_from_dir);
    Dir checkpointRestartDir = restartFromDir.getSubdir("checkpoints");

    d_restart_archive =
      std::make_shared<DataArchive>(checkpointRestartDir.getName(),
                                    d_myworld->myRank(),
                                    d_myworld->nRanks());

    std::vector<int> indices;
    std::vector<double> times;

    try {
      d_restart_archive->queryTimesteps(indices, times);
    } catch (InternalError& ie) {
      std::cerr << "\n";
      std::cerr << "An internal error was caught while trying to restart:\n";
      std::cerr << "\n";
      std::cerr << ie.message() << "\n";
      std::cerr << "This most likely means that the simulation UDA that you "
                   "have specified\n";
      std::cerr << "to use for the restart does not have any checkpoint data "
                   "in it.  Look\n";
      std::cerr << "in <uda>/checkpoints/ for timestep directories (t#####/) "
                   "to verify.\n";
      std::cerr << "\n";
      Parallel::exitAll(1);
    }

    // Find the right checkpoint timestep to query the grid
    if (indices.size() == 0) {
      std::ostringstream message;
      message << "No restart checkpoints found.";
      throw InternalError(message.str(), __FILE__, __LINE__);
    } else if (d_restart_index < 0) {
      d_restart_index = (unsigned int)(indices.size() - 1);
    } else if (d_restart_index >= static_cast<int>(indices.size())) {
      std::ostringstream message;
      message << "Invalid restart checkpoint index " << d_restart_index << ". "
              << "Found " << indices.size() << " checkpoints";
      throw InternalError(message.str(), __FILE__, __LINE__);
    }

    d_restart_timestep = indices[d_restart_index];

    // Do this call before calling DataArchive::restartInitialize,
    // because problemSetup() creates VarLabels the DataArchive needs.
    d_restart_ps =
      d_restart_archive->getTimestepDocForComponent(d_restart_index);

    proc0cout << "Restart directory: \t'" << restartFromDir.getName() << "'\n"
              << "Restart time step: \t" << d_restart_timestep << "\n";
  }
}

void
SimulationController::outputSetup()
{
  // Set up the output - needs to be done before the application is setup.
  d_output->setRuntimeStats(&d_runtime_stats);

  d_output->problemSetup(d_ups,
                         d_restart_ps,
                         d_simulator->getMaterialManagerP());
}

void
SimulationController::gridSetup()
{
  GridP grid;

  // Set up the grid.
  if (d_restarting) {
    // tsaad & bisaac: At this point, and during a restart, there not
    // legitimate load balancer. This means that the grid obtained
    // from the data archiver will global domain BCs on every MPI Rank
    // - i.e. every rank will have knowledge of ALL OTHER patches and
    // their boundary conditions.  This leads to a noticeable and
    // unacceptable increase in memory usage especially when hundreds
    // of boundaries (and boundary conditions) are present. That being
    // said, we query the grid WITHOUT requiring boundary
    // conditions. Once that is done, a legitimate load balancer will
    // be created later on - after which we use said balancer and
    // assign BCs to the grid.  NOTE the "false" argument below.
    d_current_gridP =
      d_restart_archive->queryGrid(d_restart_index, d_ups, false);
  } else /* if( !m_restarting ) */ {
    d_current_gridP = scinew Grid();

    // The call to preGridProblemSetup() by the simulation interface allows
    // for a call to grid->setExtraCells() to be made before the grid
    // problemSetup() so that if there are no extra cells defined in the
    // ups file, all levels will use the grid extra cell value.

    // For instance, Wasatch does not allow users to specify extra
    // cells through the input file. Instead, Wasatch wants to specify
    // it internally. This call gives the option to do just that though
    // it does not follow normal paradigm of calling problemSetup
    // immediately after a component or other object is created.
    d_simulator->preGridProblemSetup(d_ups, d_current_gridP);

    // Now that the simulation interface has made its changes do the normal grid
    // problemSetup()
    d_current_gridP->problemSetup(d_ups, d_myworld, d_simulator->isAMR());
  }

  if (d_current_gridP->numLevels() == 0) {
    throw InternalError("No problem (no levels in grid) specified.",
                        __FILE__,
                        __LINE__);
  }

  // Print out metadata
  if (d_myworld->myRank() == 0) {
    d_current_gridP->printStatistics();
  }
}

void
SimulationController::regridderSetup(void)
{
  // Set up the regridder.
  // Do this step before fully setting up the application interface so that the
  // Switcher (being an application) can reset the state of the regridder.
  if (d_regridder) {
    d_regridder->problemSetup(d_ups,
                              d_current_gridP,
                              d_simulator->getMaterialManagerP());
  }
}

void
SimulationController::schedulerSetup(void)
{
  // Now that the grid is completely set up, set up the scheduler.
  d_scheduler->setRuntimeStats(&d_runtime_stats);

  d_scheduler->problemSetup(d_ups, d_simulator->getMaterialManagerP());

  // Additional set up calls.
  d_scheduler->setInitTimestep(true);
  d_scheduler->setRestartInitTimestep(d_restarting);
  d_scheduler->initialize(1, 1);
  d_scheduler->clearTaskMonitoring();

  d_scheduler->advanceDataWarehouse(d_current_gridP, true);
}

void
SimulationController::loadBalancerSetup(void)
{
  // Set up the load balancer.
  d_loadBalancer->setRuntimeStats(&d_runtime_stats);

  //  Set the dimensionality of the problem.
  IntVector low, high, size;
  d_current_gridP->getLevel(0)->findCellIndexRange(low, high);

  size = high - low -
         d_current_gridP->getLevel(0)->getExtraCells() * IntVector(2, 2, 2);

  d_loadBalancer->setDimensionality(size[0] > 1, size[1] > 1, size[2] > 1);

  // In addition, do this step after regridding setup as the minimum
  // patch size that the regridder will create will be known.
  d_loadBalancer->problemSetup(d_ups,
                               d_current_gridP,
                               d_simulator->getMaterialManagerP());
}

void
SimulationController::simulatorSetup(void)
{
  // Pass the m_restart_ps to the component's problemSetup.  For
  // restarting, pull the <MaterialProperties> from the m_restart_ps.
  // If the properties are not available, then pull the properties
  // from the m_ups instead.  This step needs to be done before
  // DataArchive::restartInitialize.
  d_simulator->problemSetup(d_ups, d_restart_ps, d_current_gridP);

  // Finalize the materials
  d_simulator->getMaterialManagerP()->finalizeMaterials();

  d_simulator->setRestartTimeStep(d_restarting);
}

void
SimulationController::timeStateSetup()
{
  // Restarting so initialize time state using the archive data.
  if (d_restarting) {

    double simTime;

    d_restart_archive->restartInitialize(d_restart_index,
                                         d_current_gridP,
                                         d_scheduler->get_dw(1),
                                         d_loadBalancer,
                                         &simTime);

    // Set the time step to the restart time step which is immediately
    // written to the DW.
    d_simulator->setTimeStep(d_restart_timestep);

    // Set the simulation time to the restart simulation time which is
    // immediately written to the DW.
    d_simulator->setSimTime(simTime);

    // Set the next delta T which is immediately written to the DW.

    // Note the old delta t is the delta t used for that time step.
    d_simulator->setNextDelT(d_restart_archive->getOldDelt(d_restart_index),
                             d_restarting);

    // Tell the scheduler the generation of the re-started simulation.
    // (Add +1 because the scheduler will be starting on the next
    // time step.)
    d_scheduler->setGeneration(d_restart_timestep + 1);

    // This delete is an enigma. If it is called then memory is not
    // leaked, but sometimes if is called, then everything segfaults.
    // delete d_restart_archive;
  } else {
    // Set the default time step to which is immediately written to the DW.
    d_simulator->setTimeStep(d_simulator->getTimeStep());

    // Set the default simulation time which is immediately written to the DW.
    d_simulator->setSimTime(d_simulator->getSimTime());

    // Note the above seems back asswards but the initial sim time
    // must be set in the UPS file, the time step will always default
    // to 0. However, it is read in the problem setup stage and the
    // data warehouse is not yet available. So the above gets the
    // values in the data warehouse.
  }
}

void
SimulationController::finalSetup()
{
  // This step is done after the call to d_simulator->problemSetup to get
  // the defaults set by the simulation interface into the input.xml,
  // which the output writes along with index.xml
  d_output->initializeOutput(d_ups, d_current_gridP);

  // This step is done after the output is initialized so that global
  // reduction output vars are copied to the new uda. Further, it must
  // be called after timeStateSetup() is call so that checkpoints are
  // copied to the new uda as well.
  if (d_restarting) {
    Dir dir(d_from_dir);
    d_output->restartSetup(dir,
                           0,
                           d_restart_timestep,
                           d_simulator->getSimTime(),
                           d_restart_from_scratch,
                           d_restart_remove_old_dir);
  }

  // Miscellaneous initializations.
  ProblemSpecP amr_ps = d_ups->findBlock("AMR");
  if (amr_ps) {
    amr_ps->get("doMultiTaskgraphing", d_do_multi_taskgraphing);
  }
}

void
SimulationController::resetStats()
{
  d_runtime_stats.reset(0);
  d_simulator->resetSimulationStats(0);
}

void
SimulationController::ScheduleReportStats(bool header)
{
  Task* task = scinew Task("SimulationController::ReportStats",
                           this,
                           &SimulationController::ReportStats,
                           header);

  task->setType(Task::OncePerProc);

  // Require delta T so that the task gets scheduled
  // correctly. Otherwise the scheduler/taskgraph will toss an error :
  // Caught std exception: map::at: key not found
  task->requires(Task::NewDW, d_simulator->getDelTLabel());

  d_scheduler->addTask(task,
                       d_loadBalancer->getPerProcessorPatchSet(d_current_gridP),
                       d_simulator->getMaterialManagerP()->allMaterials());
}

void
SimulationController::ReportStats(const ProcessorGroup*,
                                  const PatchSubset*,
                                  const MaterialSubset*,
                                  DataWarehouse*,
                                  DataWarehouse*,
                                  bool header)
{
  bool reportStats = false;

  // If the reporting frequency is greater than 1 check to see if output is
  // needed.
  if (d_reportStatsFrequency == 1) {
    reportStats = true;
  }

  // Note: this check is split up so to be assured that the call to
  // isLastTimeStep is last and called only when needed. Unfortunately,
  // the greater the frequency the more often it will be called.
  else {
    if (header) {
      reportStats = true;
    } else if (d_simulator->getTimeStep() % d_reportStatsFrequency ==
               d_reportStatsOnTimeStep) {
      reportStats = true;
    } else {
      // Get the wall time if is needed, otherwise ignore it.
      double walltime;

      if (d_simulator->getWallTimeMax() > 0) {
        walltime = d_wall_timers.GetWallTime();
      } else {
        walltime = 0;
      }
      reportStats = d_simulator->isLastTimeStep(walltime);
    }
  }

  // Get and reduce the performance runtime stats
  getMemoryStats();

#ifdef HAVE_VISIT
  bool reduce = getVisIt();
#else
  bool reduce = false;
#endif

  // Reductions are only need if these are true.
  if ((d_regridder && d_regridder->useDynamicDilation()) || g_sim_stats_mem ||
      g_comp_stats || g_comp_node_stats || reduce) {

    d_runtime_stats.reduce(d_regridder && d_regridder->useDynamicDilation(),
                           d_myworld);

    // Reduce the MPI runtime stats.
    MPIScheduler* mpiScheduler =
      dynamic_cast<MPIScheduler*>(d_scheduler.get_rep());

    if (mpiScheduler) {
      mpiScheduler->d_mpi_info.reduce(d_regridder &&
                                        d_regridder->useDynamicDilation(),
                                      d_myworld);
    }
  }

  if (g_app_stats || g_app_node_stats || reduce) {
    d_simulator->getSimulationStats().calculateNodeSum(true);
    d_simulator->getSimulationStats().calculateNodeMinimum(true);
    d_simulator->getSimulationStats().calculateNodeAverage(true);
    d_simulator->getSimulationStats().calculateNodeMaximum(true);
    d_simulator->getSimulationStats().calculateNodeStdDev(true);

    d_simulator->reduceSimulationStats(d_regridder &&
                                         d_regridder->useDynamicDilation(),
                                       d_myworld);
  }

  // Update the moving average and get the wall time for this time step.
  // Timers::nanoseconds timeStepTime =
  d_wall_timers.updateExpMovingAverage();

  // Print the stats for this time step
  if (d_myworld->myRank() == 0 && g_sim_stats) {
    std::ostringstream message;

    if (header) {
      message << std::endl
              << "Simulation and run time stats are reported "
              << "at the end of each time step" << std::endl
              << "EMA == Wall time as an exponential moving average "
              << "using a window of the last " << d_wall_timers.getWindow()
              << " time steps\n"
              << "ETC == Estimated time to completion [HH:MM:SS]\n"
              << "_____________________________________________________________"
                 "_________"
              << std::endl;
    }

    //__________________________________
    // compute ETC  Estimated time to completion
    double nTimesteps =
      trunc((d_simulator->getSimTimeMax() - d_simulator->getSimTime()) /
            d_simulator->getNextDelT());
    int ETC_sec =
      trunc(d_wall_timers.ExpMovingAverage().seconds() * nTimesteps);

    using namespace std::chrono;
    seconds secs(ETC_sec);
    auto mins = duration_cast<minutes>(secs); // sec -> minutes
    secs -= duration_cast<seconds>(mins);     //
    auto hrs = duration_cast<hours>(mins);    // min -> hrs
    mins -= duration_cast<minutes>(hrs);

    int hrs_setw =
      std::fmax(2, floor(log10(hrs.count()))); // hrs can be a large number

    std::ios orgFormat(NULL); // keep track of
    orgFormat.copyfmt(message);

    message
      << std::left << "Timestep " << std::setw(8) << d_simulator->getTimeStep()
      << "Time=" << std::setw(12) << d_simulator->getSimTime()
      << "Next delT=" << std::setw(12) << d_simulator->getNextDelT()
      << "Wall Time=" << std::setw(10)
      << d_wall_timers.GetWallTime()
      //          << "Net Wall Time=" << std::setw(10) << timeStepTime.seconds()
      << "EMA=" << std::setw(12) << d_wall_timers.ExpMovingAverage().seconds()
      << std::right << "ETC=" << std::setfill('0') << std::setw(hrs_setw)
      << hrs.count() << ":" << std::setfill('0') << std::setw(2) << mins.count()
      << ":" << std::setfill('0') << std::setw(2) << secs.count() << "  ";
    message.copyfmt(orgFormat);
    message << std::left;

    // Report on the memory used.
    if (g_sim_stats_mem) {
      // With the sum reduces, use double, since with memory it is possible that
      // it will overflow
      double avg_memused        = d_runtime_stats.getRankAverage(SCIMemoryUsed);
      unsigned long max_memused = d_runtime_stats.getRankMaximum(SCIMemoryUsed);
      int max_memused_rank = d_runtime_stats.getRankForMaximum(SCIMemoryUsed);

      double avg_highwater = d_runtime_stats.getRankAverage(SCIMemoryHighwater);
      unsigned long max_highwater =
        d_runtime_stats.getRankMaximum(SCIMemoryHighwater);
      int max_highwater_rank =
        d_runtime_stats.getRankForMaximum(SCIMemoryHighwater);

      if (avg_memused == max_memused && avg_highwater == max_highwater) {
        message << "Memory Use=" << std::setw(8)
                << ProcessInfo::toHumanUnits((unsigned long)avg_memused);

        if (avg_highwater) {
          message << "    Highwater Memory Use=" << std::setw(8)
                  << ProcessInfo::toHumanUnits((unsigned long)avg_highwater);
        }
      } else {
        message << "Memory Used=" << std::setw(10)
                << ProcessInfo::toHumanUnits((unsigned long)avg_memused)
                << " (avg) " << std::setw(10)
                << ProcessInfo::toHumanUnits(max_memused)
                << " (max on rank: " << std::setw(6) << max_memused_rank << ")";

        if (avg_highwater) {
          message << "    Highwater Memory Used=" << std::setw(10)
                  << ProcessInfo::toHumanUnits((unsigned long)avg_highwater)
                  << " (avg) " << std::setw(10)
                  << ProcessInfo::toHumanUnits(max_highwater)
                  << " (max on rank: " << std::setw(6) << max_highwater_rank
                  << ")";
        }
      }
    } else {
      double memused   = d_runtime_stats[SCIMemoryUsed];
      double highwater = d_runtime_stats[SCIMemoryHighwater];

      message << "Memory Use=" << std::setw(8)
              << ProcessInfo::toHumanUnits((unsigned long)memused);

      if (highwater) {
        message << "    Highwater Memory Use=" << std::setw(8)
                << ProcessInfo::toHumanUnits((unsigned long)highwater);
      }

      message << " (on rank 0 only)";
    }

    DOUT(reportStats, message.str());
    std::cout.flush();
  }

  // Variable for calculating the percentage of time spent in overhead.
  double percent_overhead = 0;

  if ((d_regridder && d_regridder->useDynamicDilation()) || g_comp_stats) {

    // Sum up the average time for overhead related components.
    double overhead_time =
      (d_runtime_stats.getRankAverage(CompilationTime) +
       d_runtime_stats.getRankAverage(RegriddingTime) +
       d_runtime_stats.getRankAverage(RegriddingCompilationTime) +
       d_runtime_stats.getRankAverage(RegriddingCopyDataTime) +
       d_runtime_stats.getRankAverage(LoadBalancerTime));

    // Sum up the average times for simulation components.
    double total_time =
      (overhead_time + d_runtime_stats.getRankAverage(TaskExecTime) +
       d_runtime_stats.getRankAverage(TaskLocalCommTime) +
       d_runtime_stats.getRankAverage(TaskWaitCommTime) +
       d_runtime_stats.getRankAverage(TaskReduceCommTime) +
       d_runtime_stats.getRankAverage(TaskWaitThreadTime));

    // Calculate percentage of time spent in overhead.
    percent_overhead = overhead_time / total_time;
  }

  double overheadAverage = 0;

  // Set the overhead percentage. Ignore the first sample as that is
  // for initialization.
  if (d_num_samples) {
    d_overhead_values[d_overhead_index] = percent_overhead;

    double overhead = 0;
    double weight   = 0;

    int sample_size = std::min(d_num_samples, OVERHEAD_WINDOW);

    // Calculate total weight by incrementing through the overhead
    // sample array backwards and multiplying samples by the weights
    for (int i = 0; i < sample_size; ++i) {
      unsigned int index =
        (d_overhead_index - i + OVERHEAD_WINDOW) % OVERHEAD_WINDOW;
      overhead += d_overhead_values[index] * d_overhead_weights[i];
      weight += d_overhead_weights[i];
    }

    // Increment the overhead index
    d_overhead_index = (d_overhead_index + 1) % OVERHEAD_WINDOW;

    overheadAverage = overhead / weight;

    if (d_regridder) {
      d_regridder->setOverheadAverage(overheadAverage);
    }
  }

  // Ignore the first sample as that is for initialization.
  if (reportStats && d_num_samples) {

    // Infrastructure proc runtime performance stats.
    if (g_comp_stats && d_myworld->myRank() == 0) {
      d_runtime_stats.reportRankSummaryStats("Runtime Summary ",
                                             "",
                                             d_myworld->myRank(),
                                             d_myworld->nRanks(),
                                             d_simulator->getTimeStep(),
                                             d_simulator->getSimTime(),
                                             BaseInfoMapper::Dout,
                                             true);

      // Report the overhead percentage.
      if (!std::isnan(overheadAverage)) {
        std::ostringstream message;
        message << "  Percentage of time spent in overhead : "
                << overheadAverage * 100.0;

        // This code is here in case one wants to write to disk the
        // stats. Currently theses are written via Dout.
        if (1) {
          DOUT(true, message.str());
        }
        // else if( 1 ) {
        //   std::ofstream fout;
        //   std::string filename = "Runtime Summary " +
        //     (nRanks != -1 ? "." + std::to_string(nRanks)   : "") +
        //     (rank   != -1 ? "." + std::to_string(rank)     : "") +
        //     (oType == Write_Separate ? "." + std::to_string(timeStep) : "");

        //   if( oType == Write_Append )
        //     fout.open(filename, std::ofstream::out | std::ofstream::app);
        //   else
        //     fout.open(filename, std::ofstream::out);

        //   fout << message.str() << std::endl;
        //   fout.close();
        // }
      }
    }

    // Infrastructure per node runtime performance stats.
    if (g_comp_node_stats && d_myworld->myNode_myRank() == 0) {
      d_runtime_stats.reportNodeSummaryStats(
        ("Runtime Node " + d_myworld->myNodeName()).c_str(),
        "",
        d_myworld->myNode_myRank(),
        d_myworld->myNode_nRanks(),
        d_myworld->myNode(),
        d_myworld->nNodes(),
        d_simulator->getTimeStep(),
        d_simulator->getSimTime(),
        BaseInfoMapper::Dout,
        true);
    }

    // Infrastructure per proc runtime performance stats
    if (g_comp_indv_stats) {
      d_runtime_stats.reportIndividualStats("Runtime",
                                            "",
                                            d_myworld->myRank(),
                                            d_myworld->nRanks(),
                                            d_simulator->getTimeStep(),
                                            d_simulator->getSimTime(),
                                            BaseInfoMapper::Dout);
    }

    // Simulation proc runtime performance stats.
    if (g_app_stats && d_myworld->myRank() == 0) {
      d_simulator->getSimulationStats().reportRankSummaryStats(
        "Simulation Summary",
        "",
        d_myworld->myRank(),
        d_myworld->nRanks(),
        d_simulator->getTimeStep(),
        d_simulator->getSimTime(),
        BaseInfoMapper::Dout,
        false);
    }

    // Simulation per node runtime performance stats.
    if (g_app_node_stats && d_myworld->myNode_myRank() == 0) {
      d_simulator->getSimulationStats().reportNodeSummaryStats(
        ("Simulation Node " + d_myworld->myNodeName()).c_str(),
        "",
        d_myworld->myNode_myRank(),
        d_myworld->myNode_nRanks(),
        d_myworld->myNode(),
        d_myworld->nNodes(),
        d_simulator->getTimeStep(),
        d_simulator->getSimTime(),
        BaseInfoMapper::Dout,
        false);
    }

    // Simulation per proc runtime performance stats
    if (g_app_indv_stats) {
      d_simulator->getSimulationStats().reportIndividualStats(
        "Simulation",
        "",
        d_myworld->myRank(),
        d_myworld->nRanks(),
        d_simulator->getTimeStep(),
        d_simulator->getSimTime(),
        BaseInfoMapper::Dout);
    }
  }

  ++d_num_samples;

} // end reportStats()

void
SimulationController::getMemoryStats(bool create /* = false */)
{
  unsigned long memUsed;
  unsigned long highwater;
  unsigned long maxMemUsed;

  d_scheduler->checkMemoryUse(memUsed, highwater, maxMemUsed);

  d_runtime_stats[SCIMemoryUsed]      = memUsed;
  d_runtime_stats[SCIMemoryMaxUsed]   = maxMemUsed;
  d_runtime_stats[SCIMemoryHighwater] = highwater;

  if (ProcessInfo::isSupported(ProcessInfo::MEM_SIZE)) {
    d_runtime_stats[MemoryUsed] = ProcessInfo::getMemoryUsed();
  }

  if (ProcessInfo::isSupported(ProcessInfo::MEM_RSS)) {
    d_runtime_stats[MemoryResident] = ProcessInfo::getMemoryResident();
  }

  // Get memory stats for each proc if MALLOC_PERPROC is in the environment.
  if (getenv("MALLOC_PERPROC")) {
    std::ostream* mallocPerProcStream = nullptr;
    char* filenamePrefix              = getenv("MALLOC_PERPROC");

    // provide a default filename if none provided
    if (!filenamePrefix || strlen(filenamePrefix) == 0) {
      filenamePrefix = (char*)"malloc.log";
    }

    char filename[256];
    sprintf(filename, "%s.%d", filenamePrefix, d_myworld->myRank());

    if (create) {
      mallocPerProcStream =
        scinew std::ofstream(filename, std::ios::out | std::ios::trunc);
    } else {
      mallocPerProcStream =
        scinew std::ofstream(filename, std::ios::out | std::ios::app);
    }

    *mallocPerProcStream << "Proc " << d_myworld->myRank() << "   ";
    *mallocPerProcStream << "TimeStep " << d_simulator->getTimeStep()
                         << "   ";

    if (ProcessInfo::isSupported(ProcessInfo::MEM_SIZE)) {
      *mallocPerProcStream << "Size " << ProcessInfo::getMemoryUsed() << "   ";
    }

    if (ProcessInfo::isSupported(ProcessInfo::MEM_RSS)) {
      *mallocPerProcStream << "RSS " << ProcessInfo::getMemoryResident()
                           << "   ";
    }

    *mallocPerProcStream << "Sbrk "
                         << (char*)sbrk(0) - d_scheduler->getStartAddr()
                         << "   ";

#ifndef DISABLE_SCI_MALLOC
    *mallocPerProcStream << "Sci_Malloc_MemUsed " << memUsed << "   ";
    *mallocPerProcStream << "Sci_Malloc_MaxMemUsed " << maxMemUsed << "   ";
    *mallocPerProcStream << "Sci_Malloc_Highwater " << highwater;
#endif

    *mallocPerProcStream << std::endl;

    if (mallocPerProcStream) {
      delete mallocPerProcStream;
    }
  }
}

} // namespace Uintah