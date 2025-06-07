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

/*
 *  vaango.cc: Vaango - an extension of the uintah simulation system
 *
 */

#include <StandAlone/Utils/vaango_options.h>
#include <StandAlone/Utils/vaango_utils.h>

#include <CCA/Components/DataArchiver/DataArchiver.h>
#include <CCA/Components/LoadBalancers/LoadBalancerFactory.h>
#include <CCA/Components/Models/ModelFactory.h>
#include <CCA/Components/Parent/ComponentFactory.h>
#include <CCA/Components/ProblemSpecification/ProblemSpecReader.h>
#include <CCA/Components/Regridder/RegridderFactory.h>
#include <CCA/Components/Schedulers/SchedulerFactory.h>
#include <CCA/Components/SimulationController/AMRSimulationController.h>
#include <CCA/Components/Solvers/SolverFactory.h>

#ifdef HAVE_CUDA
#include <CCA/Components/Schedulers/UnifiedScheduler.h>
#endif

#include <CCA/Ports/DataWarehouse.h>

#include <Core/Exceptions/Exception.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidGrid.h>
#include <Core/Exceptions/ProblemSetupException.h>

#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/ProcessorGroup.h>

#include <Core/Util/DOUT.hpp>
#include <Core/Util/DebugStream.h>
#include <Core/Util/Environment.h>
#include <Core/Util/FileUtils.h>
#include <Core/Util/StringUtil.h>
#include <Core/Util/Timers/Timers.hpp>

#include <sci_defs/cuda_defs.h>
#include <sci_defs/hypre_defs.h>
#include <sci_defs/malloc_defs.h>
#include <sci_defs/uintah_defs.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#if HAVE_IEEEFP_H
#include <ieeefp.h>
#endif
#if 0
#include <fenv.h>
#endif

#include <cstdio>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

Uintah::MasterLock cerr_mutex{};

Uintah::Dout g_stack_debug("ExceptionStack",
                           "vaango",
                           "vaango exception stack debug stream",
                           true);
Uintah::Dout g_wait_for_debugger(
  "WaitForDebugger",
  "vaango",
  "halt program, print out pid and attach a debugger",
  false);
Uintah::Dout g_show_env(
  "ShowEnv",
  "vaango",
  "vaango show environment (the SCI env that was built up)",
  false);

} // namespace

int
main(int argc, char* argv[], char* env[])
{
  // sanity check
  Vaango::Utils::check_malloc();

  std::string oldTag;

#if HAVE_IEEEFP_H
  fpsetmask(FP_X_OFL | FP_X_DZ | FP_X_INV);
#endif
#if 0
  feenableexcept(FE_INVALID|FE_OVERFLOW|FE_DIVBYZERO);
#endif

  // Parse arguments
  Vaango::Utils::Options::parse(argc, argv);

  // Set threads
  Uintah::Parallel::setNumThreads(Vaango::Utils::Options::num_threads());
  Uintah::Parallel::setNumPartitions(Vaango::Utils::Options::num_partitions());
  Uintah::Parallel::setThreadsPerPartition(
    Vaango::Utils::Options::threads_per_partition());

#ifdef HAVE_CUDA
  // Set gpus
  if (Vaango::Utils::Options::use_gpu()) {
    Uintah::Parallel::setUsingDevice(true);
  }
#endif

  // Pass the env into the sci env so it can be used there...
  Uintah::create_sci_environment(env, nullptr, true);

  if (g_wait_for_debugger) {
    Uintah::TURN_ON_WAIT_FOR_DEBUGGER();
  }

  char* start_addr     = (char*)sbrk(0);
  bool thrownException = false;

  try {

    // Initialize after parsing the args
    Uintah::Parallel::initializeManager(argc, argv);

#if defined(MALLOC_TRACE)
    std::ostringstream traceFilename;
    traceFilename << "mallocTrace-" << Uintah::Parallel::getMPIRank();
    MALLOC_TRACE_LOG_FILE(traceFilename.str().c_str());
    // mallocTraceInfo.setTracingState( false );
#endif

    char* st = getenv("INITIAL_SLEEP_TIME");
    if (st != nullptr) {
      char name[256];
      gethostname(name, 256);
      int sleepTime = atoi(st);
      if (Uintah::Parallel::getMPIRank() == 0) {
        std::cout << "SLEEPING FOR " << sleepTime
                  << " SECONDS TO ALLOW DEBUGGER ATTACHMENT\n";
      }
      std::cout << "PID for rank " << Uintah::Parallel::getMPIRank() << " ("
                << name << ") is " << getpid() << "\n";
      std::cout.flush();

      struct timespec ts;
      ts.tv_sec  = (int)sleepTime;
      ts.tv_nsec = (int)(1.e9 * (sleepTime - ts.tv_sec));

      nanosleep(&ts, &ts);
    }

    // Read input file
    Uintah::ProblemSpecP ups = nullptr;

    auto filename = Vaango::Utils::Options::uda_filename();
    Vaango::Utils::set_input_ups_path(filename);

    try {
      ups = Uintah::ProblemSpecReader().readInputFile(
        filename, Vaango::Utils::Options::validate_ups());
    } catch (const Uintah::ProblemSetupException& err) {
      proc0cout << "\nERROR caught while parsing UPS file: "
                << Vaango::Utils::Options::uda_filename()
                << "\nDetails follow.\n"
                << err.message() << "\n";
      Vaango::Utils::stop_mpi_and_exit(0);
    } catch (...) {
      // Bulletproofing.  Catches the case where a user accidentally specifies a
      // UDA directory instead of a UPS file.
      proc0cout << "\n";
      proc0cout << "ERROR - Failed to parse UPS file: " << filename << ".\n";

      if (Uintah::validDir(filename)) {
        proc0cout << "ERROR - Note: '" << filename
                  << "' is a directory! Did you mistakenly specify a UDA "
                     "instead of an UPS file?\n";
      }
      proc0cout << "\n";
      Vaango::Utils::stop_mpi_and_exit(0);
    }

    if (Vaango::Utils::Options::only_validate_ups()) {
      std::cout << "\nValidation of .ups File finished... good bye.\n\n";
      ups = nullptr; // This cleans up memory held by the 'ups'.
      Vaango::Utils::stop_mpi_and_exit(0);
    }

    const Uintah::ProcessorGroup* world =
      Uintah::Parallel::getRootProcessorGroup();

    std::unique_ptr<Uintah::SimulationController> simController =
      std::make_unique<Uintah::AMRSimulationController>(
        world, ups, Vaango::Utils::get_input_ups_path());

    if (Vaango::Utils::Options::postprocess_uda()) {
      simController->setPostProcessFlags();
    }

    std::unique_ptr<Uintah::UintahParallelComponent> simComponent =
      Uintah::ComponentFactory::create(
        ups, world, nullptr, Vaango::Utils::Options::uda_dir());

    Uintah::SimulationInterface* simulator =
      dynamic_cast<Uintah::SimulationInterface*>(simComponent.get());

    // Read the UPS file to get the general application details.
    simulator->problemSetup(ups);

    simController->attachPort("simulator", simulator);

    // Can not do a postProcess uda with AMR
    if (Vaango::Utils::Options::postprocess_uda() && simulator->isAMR()) {
      Vaango::Utils::Options::usage(
        "You cannot use '-postprocess_uda' for an AMR simulation",
        "-postprocess_uda",
        argv[0]);
    }

    // Solver
    std::shared_ptr<Uintah::SolverInterface> solver =
      Uintah::SolverFactory::create(
        ups, world, Vaango::Utils::Options::solver_name());

    Uintah::UintahParallelComponent* solverComponent =
      dynamic_cast<Uintah::UintahParallelComponent*>(solver.get());

    simComponent->attachPort("solver", solver.get());
    solverComponent->attachPort("simulator", simulator);

    // Load balancer
    std::unique_ptr<Uintah::LoadBalancerCommon> loadBalancer =
      Uintah::LoadBalancerFactory::create(ups, world);

    loadBalancer->attachPort("simulator", simulator);
    simController->attachPort("load balancer", loadBalancer.get());
    simComponent->attachPort("load balancer", loadBalancer.get());

    // Scheduler
    Uintah::SchedulerCommon* scheduler =
      Uintah::SchedulerFactory::create(ups, world);

    scheduler->attachPort("load balancer", loadBalancer.get());
    scheduler->attachPort("simulator", simulator);

    simComponent->attachPort("scheduler", scheduler);
    simController->attachPort("scheduler", scheduler);
    loadBalancer->attachPort("scheduler", scheduler);

    scheduler->setStartAddr(start_addr);
    scheduler->addReference();

    if (Vaango::Utils::Options::emit_graphs()) {
      scheduler->doEmitTaskGraphDocs();
    }

    // Output
    std::unique_ptr<Uintah::DataArchiver> dataArchiver =
      std::make_unique<Uintah::DataArchiver>(
        world, Vaango::Utils::Options::uda_suffix());

    dataArchiver->attachPort("simulator", simulator);
    dataArchiver->attachPort("load balancer", loadBalancer.get());

    dataArchiver->setUseLocalFileSystems(
      Vaango::Utils::Options::local_filesystem());

    simController->attachPort("output", dataArchiver.get());
    simComponent->attachPort("output", dataArchiver.get());
    scheduler->attachPort("output", dataArchiver.get());

    // Regridder - optional
    std::unique_ptr<Uintah::RegridderCommon> regridder = nullptr;

    if (simulator->isAMR()) {
      regridder = Uintah::RegridderFactory::create(ups, world);

      if (regridder) {
        regridder->attachPort("scheduler", scheduler);
        regridder->attachPort("load balancer", loadBalancer.get());
        regridder->attachPort("simulator", simulator);

        simController->attachPort("regridder", regridder.get());
        simComponent->attachPort("regridder", regridder.get());

        loadBalancer->attachPort("regridder", regridder.get());
      }
    }

    // Get all the components.
    if (regridder) {
      regridder->getComponents();
    }

    scheduler->getComponents();
    loadBalancer->getComponents();
    solverComponent->getComponents();
    dataArchiver->getComponents();

    simComponent->getComponents();
    simController->getComponents();

    /*
     * Start the simulation controller
     */
    if (Vaango::Utils::Options::restart()) {
      simController->doRestart(
        Vaango::Utils::Options::uda_dir(),
        Vaango::Utils::Options::restart_checkpoint_index(),
        Vaango::Utils::Options::restart_from_scratch(),
        Vaango::Utils::Options::restart_remove_old_dir());
    }

    // This gives memory held by the 'ups' back before the simulation
    // starts... Assuming no one else is holding on to it...
    ups = nullptr;

    simController->run();

    // Clean up release all the components.
    if (regridder) {
      regridder->releaseComponents();
    }

    dataArchiver->releaseComponents();
    scheduler->releaseComponents();
    loadBalancer->releaseComponents();
    solverComponent->releaseComponents();
    simComponent->releaseComponents();
    simController->releaseComponents();

    scheduler->removeReference();
    delete scheduler;
  } catch (const Uintah::ProblemSetupException& e) {
    // Don't show a stack trace in the case of ProblemSetupException.
    std::lock_guard<Uintah::MasterLock> cerr_guard(cerr_mutex);
    std::cerr << "\n\n(Proc: " << Uintah::Parallel::getMPIRank()
              << ") Caught: " << e.message() << "\n\n";
    thrownException = true;
  } catch (const Uintah::Exception& e) {
    std::lock_guard<Uintah::MasterLock> cerr_guard(cerr_mutex);
    std::cerr << "\n\n(Proc " << Uintah::Parallel::getMPIRank()
              << ") Caught exception: " << e.message() << "\n\n";
    if (e.stackTrace()) {
      DOUT(g_stack_debug, "Stack trace: " << e.stackTrace());
    }
    thrownException = true;
  } catch (const std::bad_alloc& e) {
    std::lock_guard<Uintah::MasterLock> cerr_guard(cerr_mutex);
    std::cerr << Uintah::Parallel::getMPIRank()
              << " Caught std exception 'bad_alloc': " << e.what() << '\n';
    thrownException = true;
  } catch (const std::bad_exception& e) {
    std::lock_guard<Uintah::MasterLock> cerr_guard(cerr_mutex);
    std::cerr << Uintah::Parallel::getMPIRank()
              << " Caught std exception: 'bad_exception'" << e.what() << '\n';
    thrownException = true;
  } catch (const std::ios_base::failure& e) {
    std::lock_guard<Uintah::MasterLock> cerr_guard(cerr_mutex);
    std::cerr << Uintah::Parallel::getMPIRank()
              << " Caught std exception 'ios_base::failure': " << e.what()
              << '\n';
    thrownException = true;
  } catch (const std::runtime_error& e) {
    std::lock_guard<Uintah::MasterLock> cerr_guard(cerr_mutex);
    std::cerr << Uintah::Parallel::getMPIRank()
              << " Caught std exception 'runtime_error': " << e.what() << '\n';
    thrownException = true;
  } catch (const std::exception& e) {
    std::lock_guard<Uintah::MasterLock> cerr_guard(cerr_mutex);
    std::cerr << Uintah::Parallel::getMPIRank()
              << " Caught std exception: " << e.what() << '\n';
    thrownException = true;
  } catch (...) {
    std::lock_guard<Uintah::MasterLock> cerr_guard(cerr_mutex);
    std::cerr << Uintah::Parallel::getMPIRank()
              << " Caught unknown exception\n";
    thrownException = true;
  }

  Uintah::TypeDescription::deleteAll();

  /*
   * Finalize MPI
   */
  Uintah::Parallel::finalizeManager(thrownException
                                      ? Uintah::Parallel::Abort
                                      : Uintah::Parallel::NormalShutdown);

  if (thrownException) {
    if (Uintah::Parallel::getMPIRank() == 0) {
      std::cout << "\n\nAN EXCEPTION WAS THROWN... Goodbye.\n\n";
    }
    Uintah::Parallel::exitAll(1);
  }

  if (Uintah::Parallel::getMPIRank() == 0) {
    std::cout << "Vaango: going down successfully\n";
  }

  // use exitAll(0) since return does not work
  Uintah::Parallel::exitAll(0);
  return 0;

} // end main()
