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

#include <CCA/Components/SimulationController/AMRSimulationController.h>

#include <CCA/Components/PostProcessUda/PostProcessUda.h>

#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Output.h>
#include <CCA/Ports/ProblemSpecInterface.h>
#include <CCA/Ports/Regridder.h>
#include <CCA/Ports/Scheduler.h>
#include <CCA/Ports/SimulationInterface.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/PerPatch.h>
#include <Core/Grid/Variables/ReductionVariable.h>
#include <Core/Grid/Variables/SoleVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarLabelMatl.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Math/MiscMath.h>
#include <Core/OS/ProcessInfo.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Util/DOUT.hpp>
#include <Core/Util/DebugStream.h>
#include <Core/Util/Timers/Timers.hpp>

#include <sci_defs/gperftools_defs.h>
#include <sci_defs/malloc_defs.h>

#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace {

Uintah::Dout amrout("AMR",
                    "AMRSimulationController",
                    "AMR - report patch layout",
                    false);
Uintah::Dout dbg("AMRSimulationController",
                 "AMRSimulationController",
                 "Task/Cycle debug stream",
                 false);
Uintah::Dout dbg_barrier("MPIBarriers",
                         "AMRSimulationController",
                         "MPIBarriers debug stream",
                         false);
Uintah::Dout dbg_dwmem("LogDWMemory",
                       "AMRSimulationController",
                       "LogDWMemory debug stream",
                       false);
Uintah::Dout gprofiler("CPUProfiler",
                       "AMRSimulationController",
                       "Google Prof CPUProfiler",
                       false);
Uintah::Dout gheapprofiler("HeapProfiler",
                           "AMRSimulationController",
                           "Google Prof HeapProfiler",
                           false);
Uintah::Dout gheapchecker("HeapChecker",
                          "AMRSimulationController",
                          "Google Prof HeapChecker",
                          false);

} // namespace

namespace Uintah {

AMRSimulationController::AMRSimulationController(
  const ProcessorGroup* myworld,
  ProblemSpecP pspec,
  const std::string& input_ups_dir)
  : SimulationController(myworld, pspec, input_ups_dir)
{
}

void
AMRSimulationController::run()
{
  bool first = true;

#ifdef USE_GPERFTOOLS
  if (gprofile.active()) {
    char gprofname[512];
    sprintf(gprofname, "cpuprof-rank%d", d_myworld->myRank());
    ProfilerStart(gprofname);
  }
  if (gheapprofile.active()) {
    char gheapprofname[512];
    sprintf(gheapprofname, "heapprof-rank%d", d_myworld->myRank());
    HeapProfilerStart(gheapprofname);
  }

  HeapLeakChecker* heap_checker = nullptr;
  if (gheapchecker.active()) {
    if (!gheapprofile.active()) {
      char gheapchkname[512];
      sprintf(gheapchkname, "heapchk-rank%d", d_myworld->myRank());
      heap_checker = new HeapLeakChecker(gheapchkname);
    } else {
      std::cout << "HEAPCHECKER: Cannot start with heapprofiler"
                << "\n";
    }
  }
#endif

  //********************************************************************
  // Begin the zero time step. Which is either initialization or restart.
  //********************************************************************

  // Start the wall timer for the initialization time step
  d_wall_timers.Timestep.reset(true);

  //---------------------------------------------------------------------
  // The order of the setup is important. See the individual component
  // setup calls for more details.
  // --------------------------------------------------------------------

  // Setup the restart archive first as the output needs it.
  restartArchiveSetup();

  // Setup the output as the simulator interface needs it.
  outputSetup();

  // Setup the grid using the restart archive and simulator interface.
  gridSetup();

  // Setup the regridder using the grid.
  regridderSetup();

  // Setup the scheduler.
  schedulerSetup();

  // Setup the load balancer using the grid.
  loadBalancerSetup();

  // Setup the simulator using the restart archive and under the
  // hood the output and scheduler.
  simulatorSetup();

  // Setup the time state using the restart archive, grid, scheduler,
  // and load balancer.
  timeStateSetup();

  // Setup the final bits including the output.
  finalSetup();

  // Once the grid is set up pass it on to the GPU.
#ifdef HAVE_CUDA
  GpuUtilities::assignPatchesToGpus(d_current_gridP);
#endif

  // Setup, compile, and run the taskgraph for the initialization time step
  doInitialTimestep();

  // Update the profiler weights
  d_loadBalancer->finalizeContributions(d_current_gridP);
  d_loadBalancer->resetCostForecaster();

  // Done with all the initialization.
  d_scheduler->setInitTimestep(false);

  //*******************************************************************
  // End the zero time step. Which is either initialization or restart.
  //*******************************************************************

#ifndef DISABLE_SCI_MALLOC
  AllocatorSetDefaultTagLineNumber(d_simulator->getTimestep());
#endif

  // Set the timer for the main loop. This timer is sync'ed with the
  // simulation time to get a measurement of the simulation to wall time.
  d_wall_timers.Timestep.reset(true);

  double walltime = d_wall_timers.GetWallTime();

  //**************************************************************************
  // The main time loop; here the specified problem is actually getting solved
  //**************************************************************************
  while (!d_simulator->isLastTimestep(walltime)) {

    // Perform a bunch of housekeeping operations at the top of the
    // loop. Performing them here assures that everything is ready
    // after the initial time step. It also reduces duplicate code.

    // Before any work is done including incrementing the time step
    // check to see if this iteration may be the last one. The
    // DataArchiver uses it for determining whether to output or
    // checkpoint the last time step.
    // "maybeLastTimestep" uses the wall time and is sync'd across all ranks.

    // Get the wall time if is needed, otherwise ignore it.
    double predictedWalltime;

    // The predicted time is a best guess at what the wall time will be
    // when the time step is finished. It is currently used only for
    // outputting and checkpointing. Both of which typically take much
    // longer than the simulation calculation.
    if (d_simulator->getWallTimeMax() > 0) {
      predictedWalltime =
        walltime + 1.5 * d_wall_timers.ExpMovingAverage().seconds();
    } else {
      predictedWalltime = 0;
    }

    if (d_simulator->maybeLastTimestep(predictedWalltime)) {
      d_output->maybeLastTimestep(true);
    } else {
      d_output->maybeLastTimestep(false);
    }

    // Set the current wall time for this rank (i.e. this value is
    // sync'd across all ranks). The Data Archive uses it for
    // determining when to output or checkpoint.
    d_output->setElapsedWallTime(walltime);

    // Get the next output checkpoint time step. This step is not done
    // in d_output->beginOutputTimestep because the original values
    // are needed to compare with if there is a timestep restart so it
    // is performed here.

    // NOTE: It is called BEFORE d_simulator->prepareForNextTimestep
    // because at this point the delT, nextDelT, time step, sim time,
    // and all wall times are all in sync.
    d_output->findNext_OutputCheckPointTimestep(first && d_restarting,
                                                d_current_gridP);

    // Reset the runtime performance stats.
    resetStats();

    // Reset memory use tracking variable.
    d_scheduler->resetMaxMemValue();

    // Clear the task monitoring.
    d_scheduler->clearTaskMonitoring();

    // Increment (by one) the current time step number so components
    // know what time step they are on and get the delta T that will
    // be used.
    d_simulator->prepareForNextTimestep();

    /* Ready for next timestep */

#ifdef USE_GPERFTOOLS
    if (gheapprofile.active()) {
      char heapename[512];
      sprintf(heapename, "Timestep %d", d_simulator->getTimestep());
      HeapProfilerDump(heapename);
    }
#endif

    if (dbg_barrier.active()) {
      for (double& barrier_time : m_barrier_times) {
        barrier_time = 0;
      }
    }

    // Regridding
    if (d_regridder) {

      d_simulator->setRegridTimestep(false);

      // If not the first time step or restarting check for regridding
      if ((!first || d_restarting) &&
          d_regridder->needsToReGrid(d_current_gridP)) {

        proc0cout << " Need to regrid for next time step "
                  << d_simulator->getTimestep() << " "
                  << "at current sim time " << d_simulator->getSimTime()
                  << std::endl;

        doRegridding(false);
      }

      // Covers single-level regridder case (w/ restarts)
      else if (d_regridder->doRegridOnce() && d_regridder->isAdaptive()) {
        proc0cout << " Regridding once for next time step "
                  << d_simulator->getTimestep() << " "
                  << "at current sim time " << d_simulator->getSimTime()
                  << std::endl;

        d_scheduler->setRestartInitTimestep(false);
        doRegridding(false);
        d_regridder->setAdaptivity(false);
      }
    }

    // Compute number of dataWarehouses - multiplies by the time refinement
    // ratio for each level you increase
    int totalFine = 1;
    if (!d_simulator->isLockstepAMR()) {
      for (int i = 1; i < d_current_gridP->numLevels(); i++) {
        totalFine *= d_current_gridP->getLevel(i)->getRefinementRatioMaxDim();
      }
    }

    if (dbg_dwmem.active()) {
      // Remember, this isn't logged if DISABLE_SCI_MALLOC is set (so
      // usually in optimized mode this will not be run.)
      d_scheduler->logMemoryUse();
      std::ostringstream fn;
      fn << "alloc." << std::setw(5) << std::setfill('0') << d_myworld->myRank()
         << ".out";
      std::string filename(fn.str());

#ifndef DISABLE_SCI_MALLOC
      DumpAllocator(DefaultAllocator(), filename.c_str());
#endif
    }

    if (dbg_barrier.active()) {
      m_barrier_timer.reset(true);
      Uintah::MPI::Barrier(d_myworld->getComm());
      m_barrier_times[2] += m_barrier_timer().seconds();
    }

    // This step is a hack but it is the only way to get a new grid
    // from postProcessUda and needs to be done before
    // advanceDataWarehouse is called.
    if (d_post_process_uda) {
      d_current_gridP =
        static_cast<PostProcessUda*>(d_simulator)->getGrid(d_current_gridP);
    }

    // After one step (either time step or initialization) and the
    // updating of delta T finalize the old time step, e.g. finalize
    // and advance the DataWarehouse
    d_scheduler->advanceDataWarehouse(d_current_gridP);

#ifndef DISABLE_SCI_MALLOC
    AllocatorSetDefaultTagLineNumber(d_simulator->getTimestep());
#endif

    // Various components can request a recompile including the in
    // situ which will set the d_recompile_taskgraph flag directly.
    d_recompile_taskgraph |=
      (d_simulator->needRecompile(d_current_gridP) ||
       d_output->needRecompile(d_current_gridP) ||
       d_loadBalancer->needRecompile(d_current_gridP) ||
       (d_regridder && d_regridder->needRecompile(d_current_gridP)));

    if (d_recompile_taskgraph || first) {

      // Recompile taskgraph, re-assign BCs, reset recompile flag.
      if (d_recompile_taskgraph) {
        d_current_gridP->assignBCS(d_grid_ps, d_loadBalancer);
        d_current_gridP->performConsistencyCheck();
        d_recompile_taskgraph = false;
      }

      d_scheduler->setRestartInitTimestep(false);

      compileTaskGraph(totalFine);
    } else {
      // This is not correct if we have switched to a different
      // component, since the delT will be wrong
      d_output->finalizeTimestep(d_current_gridP, d_scheduler, false);
    }

    if (dbg_barrier.active()) {
      m_barrier_timer.reset(true);
      Uintah::MPI::Barrier(d_myworld->getComm());
      m_barrier_times[3] += m_barrier_timer().seconds();
    }

    // Execute the current time step, restarting if necessary.
    executeTimestep(totalFine);

    // If debugging, output the barrier times.
    if (dbg_barrier.active()) {
      m_barrier_timer.reset(true);
      Uintah::MPI::Barrier(d_myworld->getComm());
      m_barrier_times[4] += m_barrier_timer().seconds();

      double avg[5];
      Uintah::MPI::Reduce(
        m_barrier_times, avg, 5, MPI_DOUBLE, MPI_SUM, 0, d_myworld->getComm());

      std::ostringstream mesg;
      if (d_myworld->myRank() == 0) {
        mesg << "Barrier Times: ";
        for (double& val : avg) {
          val /= d_myworld->nRanks();
          mesg << "[" << val << "]"
               << "  ";
        }
        DOUT(dbg_barrier, mesg.str())
      }
    }

    // ARS - CAN THIS BE SCHEDULED??
    d_output->writeto_xml_files(d_current_gridP);

    // ARS - FIX ME - SCHEDULE INSTEAD
    ReportStats(nullptr, nullptr, nullptr, nullptr, nullptr, false);

    // Update the profiler weights
    d_loadBalancer->finalizeContributions(d_current_gridP);

    // Done with the first time step.
    if (first) {
      d_scheduler->setRestartInitTimestep(false);
      d_simulator->setRestartTimestep(false);

      first = false;
    }

    // The wall time is needed at the top of the loop in the while
    // conditional. So get it here.
    walltime = d_wall_timers.GetWallTime();

  } // end while ( time )

#ifdef USE_GPERFTOOLS
  if (gprofile.active()) {
    ProfilerStop();
  }
  if (gheapprofile.active()) {
    HeapProfilerStop();
  }
  if (gheapchecker.active() && !gheapprofile.active()) {
    if (heap_checker && !heap_checker->NoLeaks()) {
      std::cout << "HEAPCHECKER: MEMORY LEACK DETECTED!"
                << "\n";
    }
    delete heap_checker;
  }
#endif

} // end run()

void
AMRSimulationController::doInitialTimestep()
{
  d_scheduler->mapDataWarehouse(Task::OldDW, 0);
  d_scheduler->mapDataWarehouse(Task::NewDW, 1);
  d_scheduler->mapDataWarehouse(Task::CoarseOldDW, 0);
  d_scheduler->mapDataWarehouse(Task::CoarseNewDW, 1);

  Timers::Simple taskGraphTimer; // Task graph time

  if (d_restarting) {

    // for dynamic lb's, set up restart patch config
    d_loadBalancer->possiblyDynamicallyReallocate(d_current_gridP,
                                                  LoadBalancer::RESTART_LB);

    // tsaad & bisaac: At this point, during a restart, a grid does
    // NOT have knowledge of the boundary conditions.  (See other
    // comments in SimulationController.cc for why that is the
    // case). Here, and given a legitimate load balancer, we can
    // assign the BCs to the grid in an efficient manner.
    d_current_gridP->assignBCS(d_grid_ps, d_loadBalancer);
    d_current_gridP->performConsistencyCheck();

    // Initialize the system var (time step and simulation time). Must
    // be done before all other simulator tasks as they may need
    // these values.
    d_simulator->scheduleInitializeSystemVars(
      d_current_gridP,
      d_loadBalancer->getPerProcessorPatchSet(d_current_gridP),
      d_scheduler);

    for (int i = d_current_gridP->numLevels() - 1; i >= 0; i--) {
      d_simulator->scheduleRestartInitialize(d_current_gridP->getLevel(i),
                                             d_scheduler);
    }

    // Report all of the stats before doing any possible in-situ work
    // as that effects the lap timer for the time steps.
    // ScheduleReportStats( true );

    // If compiled with VisIt check the in-situ status for work.
    // ScheduleCheckInSitu( true );

    taskGraphTimer.reset(true);
    d_scheduler->compile();
    taskGraphTimer.stop();

    d_runtime_stats[CompilationTime] += taskGraphTimer().seconds();

    proc0cout << "Done with taskgraph compile (" << taskGraphTimer().seconds()
              << " seconds)" << std::endl;

    // No scrubbing for initial step
    d_scheduler->get_dw(1)->setScrubbing(DataWarehouse::ScrubNone);
    d_scheduler->execute();

    // Now we know we're done with any additions to the new DW - finalize it
    d_scheduler->get_dw(1)->finalize();

    if (d_regridder && d_regridder->isAdaptive()) {
      // On restart:
      //   we must set up the tasks (but not compile) so we can have the
      //   initial OldDW Requirements in order to regrid straight away
      for (int i = d_current_gridP->numLevels() - 1; i >= 0; i--) {
        d_simulator->scheduleTimeAdvance(d_current_gridP->getLevel(i),
                                         d_scheduler);
      }
    }

    // Monitoring tasks must be scheduled last!!
    for (int i = d_current_gridP->numLevels() - 1; i >= 0; i--) {
      d_scheduler->scheduleTaskMonitoring(d_current_gridP->getLevel(i));
    }
  } else /* if (!d_restarting) */ {
    // for dynamic lb's, set up initial patch config
    d_loadBalancer->possiblyDynamicallyReallocate(d_current_gridP,
                                                  LoadBalancer::INIT_LB);

    d_current_gridP->assignBCS(d_grid_ps, d_loadBalancer);
    d_current_gridP->performConsistencyCheck();

    bool needNewLevel = false;

    do {
      proc0cout << "\nCompiling initialization taskgraph." << std::endl;

      // Initialize the system var (time step and simulation
      // time). Must be done before all other simulator tasks as
      // they may need these values.
      d_simulator->scheduleInitializeSystemVars(
        d_current_gridP,
        d_loadBalancer->getPerProcessorPatchSet(d_current_gridP),
        d_scheduler);

      // Initialize the CFD and/or MPM data
      for (int i = d_current_gridP->numLevels() - 1; i >= 0; i--) {
        d_simulator->scheduleInitialize(d_current_gridP->getLevel(i),
                                        d_scheduler);

        if (d_regridder) {
          // So we can initially regrid
          d_regridder->scheduleInitializeErrorEstimate(
            d_current_gridP->getLevel(i));
          d_simulator->scheduleInitialErrorEstimate(
            d_current_gridP->getLevel(i), d_scheduler);

          // We don't use error estimates if we don't make another
          // level, so don't dilate.
          if (i < d_regridder->maxLevels() - 1) {
            d_regridder->scheduleDilation(d_current_gridP->getLevel(i),
                                          d_simulator->isLockstepAMR());
          }
        }
      }

      // Compute the next time step.
      scheduleComputeStableTimestep();

      // ARS COMMENT THESE TASKS WILL BE SCHEDULED FOR EACH
      // LEVEL. THAT MAY OR MAY NOT BE REASONABLE.

      // NOTE ARS - FIXME before the output so the values can be saved.
      // Monitoring tasks must be scheduled last!!
      for (int i = 0; i < d_current_gridP->numLevels(); i++) {
        d_scheduler->scheduleTaskMonitoring(d_current_gridP->getLevel(i));
      }

      // Output tasks
      const bool recompile = true;

      d_output->finalizeTimestep(d_current_gridP, d_scheduler, recompile);

      d_output->sched_allOutputTasks(d_current_gridP, d_scheduler, recompile);

      // Report all of the stats before doing any possible in-situ work
      // as that effects the lap timer for the time steps.
      // ScheduleReportStats( true );

      // If compiled with VisIt check the in-situ status for work.
      // ScheduleCheckInSitu( true );

      taskGraphTimer.reset(true);
      d_scheduler->compile();
      taskGraphTimer.stop();

      d_runtime_stats[CompilationTime] += taskGraphTimer().seconds();

      proc0cout << "Done with taskgraph compile (" << taskGraphTimer().seconds()
                << " seconds)" << std::endl;

      // No scrubbing for initial step
      d_scheduler->get_dw(1)->setScrubbing(DataWarehouse::ScrubNone);
      d_scheduler->execute();

      needNewLevel = (d_regridder && d_regridder->isAdaptive() &&
                      d_current_gridP->numLevels() < d_regridder->maxLevels() &&
                      doRegridding(true));

      if (needNewLevel) {
        d_scheduler->initialize(1, 1);
        d_scheduler->advanceDataWarehouse(d_current_gridP, true);
      }

    } while (needNewLevel);

    d_output->writeto_xml_files(d_current_gridP);
  }

  // ARS - FIX ME - SCHEDULE INSTEAD
  ReportStats(nullptr, nullptr, nullptr, nullptr, nullptr, true);

} // end doInitialTimestep()

void
AMRSimulationController::executeTimestep(int totalFine)
{
  // If the time step needs to be recomputed, this loop will execute
  // multiple times.
  bool success = false;

  int tg_index = d_simulator->getTaskGraphIndex();

  // Execute at least once.
  while (!success) {
    d_simulator->setDelTForAllLevels(d_scheduler, d_current_gridP, totalFine);

    // Standard data warehouse scrubbing.
    if (m_scrub_datawarehouse && d_loadBalancer->getNthRank() == 1) {
      if (d_simulator->activeReductionVariable(recomputeTimestep_name)) {
        d_scheduler->get_dw(0)->setScrubbing(DataWarehouse::ScrubNonPermanent);

        static bool reported = false;

        if (!reported) {

          proc0cout << "The ability to recompute a time step is turned on. "
                       "This will affect the data warehouse scrubbing for all "
                       "time steps "
                       "regardless if a time step is recomputed or not. "
                       "There may be a memory increase."
                    << std::endl;
        }

        reported = true;
      } else {
        d_scheduler->get_dw(0)->setScrubbing(DataWarehouse::ScrubComplete);
      }

      // The other data warehouse as well as those for other levels.
      for (int i = 1; i <= totalFine; ++i) {
        d_scheduler->get_dw(i)->setScrubbing(DataWarehouse::ScrubNonPermanent);
      }
    }
    // If not scrubbing or getNthRank requires the variables after
    // they would have been scrubbed so turn off all scrubbing.
    else { // if (!d_scrub_datawarehouse || d_loadBalancer->getNthRank() > 1)
      for (int i = 0; i <= totalFine; ++i) {
        d_scheduler->get_dw(i)->setScrubbing(DataWarehouse::ScrubNone);
      }
    }

    if (d_do_multi_taskgraphing) {
      subCycleExecute(0, totalFine, 0, true);
    }
    // TG index set by component that requested temporal scheduling
    //   (multiple primary task graphs) this is passed to
    //   scheduler->execute(), default index is 0
    else {
      int iteration =
        (d_last_recompile_timeStep == d_simulator->getTimestep()) ? 0 : 1;

      d_scheduler->execute(tg_index, iteration);
    }

    //  If time step is to be recomputed adjust the delta T and recompute.
    if (d_simulator->getReductionVariable(recomputeTimestep_name)) {

      for (int i = 1; i <= totalFine; ++i) {
        d_scheduler->replaceDataWarehouse(i, d_current_gridP);
      }

      // Recompute the delta T.
      d_simulator->recomputeDelT();

      // Use the recomputed DelT and check the need for performing an
      // output and checkpoint time step.
      d_output->recompute_OutputCheckPointTimestep();

      success = false;
    } else if (d_simulator->getReductionVariable(abortTimestep_name)) {
      proc0cout << "Time step aborted and cannot recompute it, "
                // << "outputing and checkpointing the time step. "
                << "Ending the simulation." << std::endl;

      d_simulator->setReductionVariable(
        d_scheduler->getLastDW(), abortTimestep_name, true);

      // This should be a for the previous time step.
      // d_output->setOutputTimestep( true, d_current_gridP );
      // d_output->setCheckpointTimestep( true, d_current_gridP );

      success = true;
    } else {
      success = true;
    }
  }
}

auto
AMRSimulationController::doRegridding(bool initialTimestep) -> bool
{
  Timers::Simple regriddingTimer; // Regridding time

  regriddingTimer.start();

  bool retVal = false;

  if (!initialTimestep) {
    proc0cout << "_____________________________________________________________"
                 "_________\n";
  }

  GridP oldGrid = d_current_gridP;
  d_current_gridP =
    d_regridder->regrid(oldGrid.get_rep(), d_simulator->getTimestep());

  if (dbg_barrier.active()) {
    m_barrier_timer.reset(true);
    Uintah::MPI::Barrier(d_myworld->getComm());
    m_barrier_times[0] += m_barrier_timer().seconds();
  }

  regriddingTimer.stop();

  d_runtime_stats[RegriddingTime] += regriddingTimer().seconds();

  d_simulator->setRegridTimestep(false);

  int lbstate =
    initialTimestep ? LoadBalancer::INIT_LB : LoadBalancer::REGRID_LB;

  if (d_current_gridP != oldGrid) {

    d_simulator->setRegridTimestep(true);

    d_loadBalancer->possiblyDynamicallyReallocate(d_current_gridP, lbstate);

    if (dbg_barrier.active()) {
      m_barrier_timer.reset(true);
      Uintah::MPI::Barrier(d_myworld->getComm());
      m_barrier_times[1] += m_barrier_timer().seconds();
    }

    d_current_gridP->assignBCS(d_grid_ps, d_loadBalancer);
    d_current_gridP->performConsistencyCheck();

    //__________________________________
    //  output regridding stats
    if (d_myworld->myRank() == 0) {

      proc0cout << "  REGRIDDING:";

      // amrout << "---------- OLD GRID ----------" << std::endl <<
      // *(oldGrid.get_rep());
      for (int i = 0; i < d_current_gridP->numLevels(); i++) {
        proc0cout << " Level " << i << " has "
                  << d_current_gridP->getLevel(i)->numPatches()
                  << " patch(es).";
      }
      proc0cout << std::endl;

      if (amrout.active()) {
        DOUT(true,
             "---------- NEW GRID ----------\n"
               << "Grid has " << d_current_gridP->numLevels() << " level(s)");

        for (int levelIndex = 0; levelIndex < d_current_gridP->numLevels();
             levelIndex++) {
          LevelP level = d_current_gridP->getLevel(levelIndex);

          DOUT(true,
               "  Level " << level->getID() << ", indx: " << level->getIndex()
                          << " has " << level->numPatches() << " patch(es).");

          for (auto patchIter = level->patchesBegin();
               patchIter < level->patchesEnd();
               patchIter++) {
            const Patch* patch = *patchIter;
            DOUT(true,
                 "(Patch " << patch->getID() << " proc "
                           << d_loadBalancer->getPatchwiseProcessorAssignment(
                                patch)
                           << ": box=" << patch->getExtraBox()
                           << ", lowIndex=" << patch->getExtraCellLowIndex()
                           << ", highIndex=" << patch->getExtraCellHighIndex()
                           << ")");
          }
        }
      }
    } // rank 0

    Timers::Simple schedulerTimer;

    if (!initialTimestep) {
      schedulerTimer.start();
      d_scheduler->scheduleAndDoDataCopy(d_current_gridP);
      schedulerTimer.stop();
    }

    proc0cout << "Done regridding for next time step "
              << d_simulator->getTimestep() << " "
              << "at current sim time " << d_simulator->getSimTime() << ", "
              << "total time took "
              << regriddingTimer().seconds() + schedulerTimer().seconds()
              << " seconds, "
              << "regridding took " << regriddingTimer().seconds()
              << " seconds";

    if (!initialTimestep) {
      proc0cout << ", scheduling and copying took "
                << schedulerTimer().seconds() << " seconds";
    }

    proc0cout << "." << std::endl;

    retVal = true;
  } // grid != oldGrid

  if (!initialTimestep) {
    proc0cout << "_____________________________________________________________"
                 "_________\n";
  }

  return retVal;
}

void
AMRSimulationController::compileTaskGraph(int totalFine)
{
  Timers::Simple taskGraphTimer;

  taskGraphTimer.start();

  proc0cout << "Compiling taskgraph(s)"
            << "..." << std::endl;

  d_output->recompile(d_current_gridP);

  d_last_recompile_timeStep = d_simulator->getTimestep();

  d_scheduler->initialize(1, totalFine);
  d_scheduler->fillDataWarehouses(d_current_gridP);

  // Set up new DWs, DW mappings.
  d_scheduler->clearMappings();
  d_scheduler->mapDataWarehouse(Task::OldDW, 0);
  d_scheduler->mapDataWarehouse(Task::NewDW, totalFine);
  d_scheduler->mapDataWarehouse(Task::CoarseOldDW, 0);
  d_scheduler->mapDataWarehouse(Task::CoarseNewDW, totalFine);

  int my_rank = d_myworld->myRank();
  if (d_do_multi_taskgraphing) {
    for (int i = 0; i < d_current_gridP->numLevels(); i++) {
      // taskgraphs 0-numlevels-1
      if (i > 0) {
        // we have the first one already
        d_scheduler->addTaskGraph(Scheduler::NormalTaskGraph);
      }
      DOUT(dbg, my_rank << "   Creating level " << i << " task graph");

      d_simulator->scheduleTimeAdvance(d_current_gridP->getLevel(i),
                                       d_scheduler);
    }

    for (int i = 0; i < d_current_gridP->numLevels(); i++) {
      if (d_simulator->isAMR() && d_current_gridP->numLevels() > 1) {
        DOUT(dbg,
             my_rank << "   Doing Intermediate TG level " << i
                     << " task graph");
        // taskgraphs numlevels-2*numlevels-1
        d_scheduler->addTaskGraph(Scheduler::IntermediateTaskGraph);
      }

      // schedule a coarsen from the finest level to this level
      for (int j = d_current_gridP->numLevels() - 2; j >= i; j--) {
        DOUT(dbg, my_rank << "   schedule coarsen on level " << j);
        d_simulator->scheduleCoarsen(d_current_gridP->getLevel(j), d_scheduler);
      }

      d_simulator->scheduleFinalizeTimestep(d_current_gridP->getLevel(i),
                                            d_scheduler);

      // schedule a refineInterface from this level to the finest level
      for (int j = i; j < d_current_gridP->numLevels(); j++) {
        if (j != 0) {
          DOUT(dbg,
               my_rank << "   schedule RI on level " << j << " for tg " << i
                       << " coarseold " << (j == i) << " coarsenew " << true);
          d_simulator->scheduleRefineInterface(
            d_current_gridP->getLevel(j), d_scheduler, j == i, true);
        }
      }
    }
    // for the final error estimate and stable timestep tasks
    d_scheduler->addTaskGraph(Scheduler::IntermediateTaskGraph);
  } else /* if ( !d_do_multi_taskgraphing ) */ {
    subCycleCompile(0, totalFine, 0, 0);

    d_scheduler->clearMappings();
    d_scheduler->mapDataWarehouse(Task::OldDW, 0);
    d_scheduler->mapDataWarehouse(Task::NewDW, totalFine);
  }

  // If regridding schedule error estimates
  for (int i = d_current_gridP->numLevels() - 1; i >= 0; i--) {
    DOUT(dbg, my_rank << "   final TG " << i);

    if (d_regridder) {
      d_regridder->scheduleInitializeErrorEstimate(
        d_current_gridP->getLevel(i));
      d_simulator->scheduleErrorEstimate(d_current_gridP->getLevel(i),
                                         d_scheduler);

      if (i < d_regridder->maxLevels() -
                1) { // we don't use error estimates if we don't make another
                     // level, so don't dilate
        d_regridder->scheduleDilation(d_current_gridP->getLevel(i),
                                      d_simulator->isLockstepAMR());
      }
    }
  }

  // After all tasks are done schedule the on-the-fly and other analysis.
  for (int i = 0; i < d_current_gridP->numLevels(); i++) {
    d_simulator->scheduleAnalysis(d_current_gridP->getLevel(i), d_scheduler);
  }

  // Compute the next time step.
  scheduleComputeStableTimestep();

  // NOTE ARS - FIXME before the output so the values can be saved.
  // Monitoring tasks must be scheduled last!!
  for (int i = 0; i < d_current_gridP->numLevels(); i++) {
    d_scheduler->scheduleTaskMonitoring(d_current_gridP->getLevel(i));
  }

  // Output tasks
  d_output->finalizeTimestep(d_current_gridP, d_scheduler, true);

  d_output->sched_allOutputTasks(d_current_gridP, d_scheduler, true);

  // Update the system var (time step and simulation time). Must be
  // done after the output and after scheduleComputeStableTimestep.
  d_simulator->scheduleUpdateSystemVars(
    d_current_gridP,
    d_loadBalancer->getPerProcessorPatchSet(d_current_gridP),
    d_scheduler);

  // Report all of the stats before doing any possible in-situ work
  // as that effects the lap timer for the time steps.
  // ScheduleReportStats( false );

  // If compiled with VisIt check the in-situ status for work.
  // ScheduleCheckInSitu( false );

  d_scheduler->compile();

  taskGraphTimer.stop();

  d_runtime_stats[CompilationTime] += taskGraphTimer().seconds();

  proc0cout << "Done with taskgraph compile (" << taskGraphTimer().seconds()
            << " seconds)" << std::endl;

} // end compileTaskGraph()

void
AMRSimulationController::subCycleCompile(int startDW,
                                         int dwStride,
                                         int numLevel,
                                         int step)
{
  // amrout << "Start AMRSimulationController::subCycleCompile, level=" <<
  // numLevel << '\n';

  // We are on (the fine) level numLevel
  LevelP fineLevel = d_current_gridP->getLevel(numLevel);
  LevelP coarseLevel;
  int coarseStartDW;
  int coarseDWStride;
  int numCoarseSteps; // how many steps between this level and the coarser
  int numFineSteps;   // how many steps between this level and the finer
  if (numLevel > 0) {
    numCoarseSteps =
      d_simulator->isLockstepAMR() ? 1 : fineLevel->getRefinementRatioMaxDim();
    coarseLevel    = d_current_gridP->getLevel(numLevel - 1);
    coarseDWStride = dwStride * numCoarseSteps;
    coarseStartDW  = (startDW / coarseDWStride) * coarseDWStride;
  } else {
    coarseDWStride = dwStride;
    coarseStartDW  = startDW;
    numCoarseSteps = 0;
  }

  ASSERT(dwStride > 0 && numLevel < d_current_gridP->numLevels())
  d_scheduler->clearMappings();
  d_scheduler->mapDataWarehouse(Task::OldDW, startDW);
  d_scheduler->mapDataWarehouse(Task::NewDW, startDW + dwStride);
  d_scheduler->mapDataWarehouse(Task::CoarseOldDW, coarseStartDW);
  d_scheduler->mapDataWarehouse(Task::CoarseNewDW,
                                coarseStartDW + coarseDWStride);

  d_simulator->scheduleTimeAdvance(fineLevel, d_scheduler);

  if (d_simulator->isAMR()) {
    if (numLevel + 1 < d_current_gridP->numLevels()) {
      numFineSteps  = d_simulator->isLockstepAMR()
                        ? 1
                        : fineLevel->getFinerLevel()->getRefinementRatioMaxDim();
      int newStride = dwStride / numFineSteps;

      for (int substep = 0; substep < numFineSteps; substep++) {
        subCycleCompile(
          startDW + substep * newStride, newStride, numLevel + 1, substep);
      }

      // Coarsen and then refine_CFI at the end of the W-cycle
      d_scheduler->clearMappings();
      d_scheduler->mapDataWarehouse(Task::OldDW, 0);
      d_scheduler->mapDataWarehouse(Task::NewDW, startDW + dwStride);
      d_scheduler->mapDataWarehouse(Task::CoarseOldDW, startDW);
      d_scheduler->mapDataWarehouse(Task::CoarseNewDW, startDW + dwStride);
      d_simulator->scheduleCoarsen(fineLevel, d_scheduler);
    }
  }

  d_scheduler->clearMappings();
  d_scheduler->mapDataWarehouse(Task::OldDW, startDW);
  d_scheduler->mapDataWarehouse(Task::NewDW, startDW + dwStride);
  d_scheduler->mapDataWarehouse(Task::CoarseOldDW, coarseStartDW);
  d_scheduler->mapDataWarehouse(Task::CoarseNewDW,
                                coarseStartDW + coarseDWStride);
  d_simulator->scheduleFinalizeTimestep(fineLevel, d_scheduler);

  // do refineInterface after the freshest data we can get; after the
  // finer level's coarsen completes do all the levels at this point
  // in time as well, so all the coarsens go in order, and then the
  // refineInterfaces
  if (d_simulator->isAMR() && (step < numCoarseSteps - 1 || numLevel == 0)) {

    for (int i = fineLevel->getIndex(); i < fineLevel->getGrid()->numLevels();
         i++) {
      if (i == 0) {
        continue;
      }
      if (i == fineLevel->getIndex() && numLevel != 0) {
        d_scheduler->mapDataWarehouse(Task::CoarseOldDW, coarseStartDW);
        d_scheduler->mapDataWarehouse(Task::CoarseNewDW,
                                      coarseStartDW + coarseDWStride);
        d_simulator->scheduleRefineInterface(
          fineLevel, d_scheduler, true, true);
      } else {
        // look in the NewDW all the way down
        d_scheduler->mapDataWarehouse(Task::CoarseOldDW, 0);
        d_scheduler->mapDataWarehouse(Task::CoarseNewDW, startDW + dwStride);
        d_simulator->scheduleRefineInterface(
          fineLevel->getGrid()->getLevel(i), d_scheduler, false, true);
      }
    }
  }
}

void
AMRSimulationController::subCycleExecute(int startDW,
                                         int dwStride,
                                         int levelNum,
                                         [[maybe_unused]] bool rootCycle)
{
  // there are 2n+1 taskgraphs, n for the basic timestep, n for intermediate
  // timestep work, and 1 for the errorEstimate and stableTimestep, where n
  // is the number of levels.

  // amrout << "Start AMRSimulationController::subCycleExecute, level=" <<
  // numLevel << '\n';

  // We are on (the fine) level numLevel
  int numSteps;
  if (levelNum == 0 || d_simulator->isLockstepAMR()) {
    numSteps = 1;
  } else {
    numSteps = d_current_gridP->getLevel(levelNum)->getRefinementRatioMaxDim();
  }

  int newDWStride = dwStride / numSteps;

  DataWarehouse::ScrubMode oldScrubbing =
    (d_simulator->activeReductionVariable(
      recomputeTimestep_name) /*|| d_loadBalancer->isDynamic()*/)
      ? DataWarehouse::ScrubNonPermanent
      : DataWarehouse::ScrubComplete;

  int curDW = startDW;
  for (int step = 0; step < numSteps; step++) {

    if (step > 0) {
      curDW += newDWStride; // can't increment at the end, or the FINAL tg for
                            // L0 will use the wrong DWs
    }

    d_scheduler->clearMappings();
    d_scheduler->mapDataWarehouse(Task::OldDW, curDW);
    d_scheduler->mapDataWarehouse(Task::NewDW, curDW + newDWStride);
    d_scheduler->mapDataWarehouse(Task::CoarseOldDW, startDW);
    d_scheduler->mapDataWarehouse(Task::CoarseNewDW, startDW + dwStride);

    // we really only need to pass in whether the current DW is mapped to 0 or
    // not
    // TODO - fix inter-Taskgraph scrubbing
    d_scheduler->get_dw(curDW)->setScrubbing(oldScrubbing); // OldDW
    d_scheduler->get_dw(curDW + newDWStride)
      ->setScrubbing(DataWarehouse::ScrubNonPermanent);       // NewDW
    d_scheduler->get_dw(startDW)->setScrubbing(oldScrubbing); // CoarseOldDW
    d_scheduler->get_dw(startDW + dwStride)
      ->setScrubbing(DataWarehouse::ScrubNonPermanent); // CoarseNewDW

    // we need to unfinalize because execute finalizes all new DWs,
    // and we need to write into them still (even if we finalized only
    // the NewDW in execute, we will still need to write into that DW)
    d_scheduler->get_dw(curDW + newDWStride)->unfinalize();

    // iteration only matters if it's zero or greater than 0
    int iteration =
      curDW + (d_last_recompile_timeStep == d_simulator->getTimestep() ? 0 : 1);

    DOUT(dbg,
         d_myworld->myRank()
           << "   Executing TG on level " << levelNum << " with old DW "
           << curDW << "=" << d_scheduler->get_dw(curDW)->getID() << " and new "
           << curDW + newDWStride << "="
           << d_scheduler->get_dw(curDW + newDWStride)->getID()
           << " CO-DW: " << startDW << " CNDW " << startDW + dwStride
           << " on iteration " << iteration);

    d_scheduler->execute(levelNum, iteration);

    if (levelNum + 1 < d_current_gridP->numLevels()) {
      ASSERT(newDWStride > 0);
      subCycleExecute(curDW, newDWStride, levelNum + 1, false);
    }

    if (d_simulator->isAMR() && d_current_gridP->numLevels() > 1 &&
        (step < numSteps - 1 || levelNum == 0)) {
      // Since the execute of the intermediate is time-based,
      // execute the intermediate TG relevant to this level, if we are in the
      // middle of the subcycle or at the end of level 0.
      // the end of the cycle will be taken care of by the parent level sybcycle
      d_scheduler->clearMappings();
      d_scheduler->mapDataWarehouse(Task::OldDW, curDW);
      d_scheduler->mapDataWarehouse(Task::NewDW, curDW + newDWStride);
      d_scheduler->mapDataWarehouse(Task::CoarseOldDW, startDW);
      d_scheduler->mapDataWarehouse(Task::CoarseNewDW, startDW + dwStride);

      d_scheduler->get_dw(curDW)->setScrubbing(oldScrubbing); // OldDW
      d_scheduler->get_dw(curDW + newDWStride)
        ->setScrubbing(DataWarehouse::ScrubNonPermanent);       // NewDW
      d_scheduler->get_dw(startDW)->setScrubbing(oldScrubbing); // CoarseOldDW
      d_scheduler->get_dw(startDW + dwStride)
        ->setScrubbing(DataWarehouse::ScrubNonPermanent); // CoarseNewDW

      DOUT(dbg,
           d_myworld->myRank()
             << "   Executing INT TG on level " << levelNum << " with old DW "
             << curDW << "=" << d_scheduler->get_dw(curDW)->getID()
             << " and new " << curDW + newDWStride << "="
             << d_scheduler->get_dw(curDW + newDWStride)->getID()
             << " CO-DW: " << startDW << " CNDW " << startDW + dwStride
             << " on iteration " << iteration);

      d_scheduler->get_dw(curDW + newDWStride)->unfinalize();
      d_scheduler->execute(levelNum + d_current_gridP->numLevels(), iteration);
    }

    if (curDW % dwStride != 0) {
      // the currentDW(old datawarehouse) should no longer be needed - in the
      // case of NonPermanent OldDW scrubbing
      d_scheduler->get_dw(curDW)->clear();
    }
  }

  if (levelNum == 0) {
    // execute the final TG
    DOUT(dbg,
         d_myworld->myRank()
           << "   Executing Final TG on level " << levelNum << " with old DW "
           << curDW << " = " << d_scheduler->get_dw(curDW)->getID()
           << " and new " << curDW + newDWStride << " = "
           << d_scheduler->get_dw(curDW + newDWStride)->getID());

    d_scheduler->get_dw(curDW + newDWStride)->unfinalize();
    d_scheduler->execute(d_scheduler->getNumTaskGraphs() - 1, 1);
  }
} // end subCycleExecute()

void
AMRSimulationController::scheduleComputeStableTimestep()
{
  // Schedule the simulator to compute the next time step on a per
  // patch basis.
  for (int i = 0; i < d_current_gridP->numLevels(); i++) {
    d_simulator->scheduleComputeStableTimestep(d_current_gridP->getLevel(i),
                                               d_scheduler);
  }

  // Schedule the reduction of the time step and other variables on a
  // per patch basis to a per rank basis.
  d_simulator->scheduleReduceSystemVars(
    d_current_gridP,
    d_loadBalancer->getPerProcessorPatchSet(d_current_gridP),
    d_scheduler);
}

} // namespace Uintah