#include <TauProfilerForSCIRun.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Exceptions/InvalidGrid.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/Tracker/TrackerClient.h>

#include <CCA/Components/ProblemSpecification/ProblemSpecReader.h>
#include <CCA/Components/SimulationController/AMRSimulationController.h>
#include <CCA/Components/Models/ModelFactory.h>
#include <CCA/Components/Solvers/CGSolver.h>
#include <CCA/Components/Solvers/DirectSolve.h>

#ifdef HAVE_HYPRE
#include <CCA/Components/Solvers/HypreSolver.h>
#endif

#include <CCA/Components/PatchCombiner/PatchCombiner.h>
#include <CCA/Components/PatchCombiner/UdaReducer.h>
#include <CCA/Components/DataArchiver/DataArchiver.h>
#include <CCA/Components/Solvers/SolverFactory.h>
#include <CCA/Components/Regridder/RegridderFactory.h>
#include <CCA/Components/LoadBalancers/LoadBalancerFactory.h>
#include <CCA/Components/Schedulers/SchedulerFactory.h>
#include <CCA/Components/Parent/ComponentFactory.h>
#include <CCA/Ports/DataWarehouse.h>

#include <Core/Exceptions/Exception.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Thread/Mutex.h>
#include <Core/Thread/Time.h>
#include <Core/Thread/Thread.h>
#include <Core/Util/DebugStream.h>
#include <Core/Util/Environment.h>
#include <Core/Util/FileUtils.h>

#include <sci_defs/hypre_defs.h>
#include <sci_defs/malloc_defs.h>
#include <sci_defs/mpi_defs.h>
#include <sci_defs/uintah_defs.h>
#include <sci_defs/cuda_defs.h>

#include <svn_info.h>

#include <Core/Malloc/Allocator.h>

#ifdef USE_VAMPIR
#  include <Core/Parallel/Vampir.h>
#endif

#if HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif
#if 0
#  include <fenv.h>
#endif

#ifdef _WIN32
#  include <process.h>
#  include <winsock2.h>
#endif

#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <stdexcept>
#include <sys/stat.h>

#include <time.h>

using namespace Matiti;
using namespace std;


// Debug: Used to sync cerr so it is readable (when output by
// multiple threads at the same time)
// Mutex cerrLock( "cerr lock" );
// DebugStream mixedDebug( "MixedScheduler Debug Output Stream", false );
// DebugStream fullDebug( "MixedScheduler Full Debug", false );

extern Mutex cerrLock;
extern DebugStream mixedDebug;
extern DebugStream fullDebug;
static DebugStream stackDebug("ExceptionStack", true);
static DebugStream dbgwait("WaitForDebugger", false);

static
void
quit( const std::string & msg = "" )
{
  if (msg != "") {
    cerr << msg << "\n";
  }
  Matiti::Parallel::finalizeManager();
  Thread::exitAll( 2 );
}

static
void
usage( const std::string & message,
       const std::string& badarg,
       const std::string& progname)
{
#ifndef HAVE_MPICH_OLD
  int argc = 0;
  char **argv;
  argv = 0;

  // Initialize MPI so that "usage" is only printed by proc 0.
  // (If we are using MPICH, then MPI_Init() has already been called.)
  Matiti::Parallel::initializeManager( argc, argv );
#endif

  if( Matiti::Parallel::getMPIRank() == 0 ) {
      cerr << "\n";
      if(badarg != "") {
        cerr << "Error parsing argument: " << badarg << '\n';
      }
      cerr << "\n";
      cerr << message << "\n";
      cerr << "\n";
      cerr << "Usage: " << progname << " [options] <input_file_name>\n\n";
      cerr << "Valid options are:\n";
      cerr << "-h[elp]              : This usage information.\n";
      cerr << "-AMR                 : use AMR simulation controller\n";
#ifdef HAVE_CUDA
      cerr << "-gpu                 : use available GPU devices, requires a multi-threaded GPU scheduler \n";
#endif
      cerr << "-nthreads <#>        : number of threads per MPI process, requires a multi-threaded scheduler\n";
      cerr << "-layout NxMxO        : Eg: 2x1x1.  MxNxO must equal number\n";
      cerr << "                           of boxes you are using.\n";
      cerr << "-emit_taskgraphs     : Output taskgraph information\n";
      cerr << "-restart             : Give the checkpointed uda directory as the input file\n";
      cerr << "-combine_patches     : Give a uda directory as the input file\n";  
      cerr << "-reduce_uda          : Reads <uda-dir>/input.xml file and removes unwanted labels (see FAQ).\n";
      cerr << "-uda_suffix <number> : Make a new uda dir with <number> as the default suffix\n";      
      cerr << "-t <timestep>        : Restart timestep (last checkpoint is default,\n\t\t\tyou can use -t 0 for the first checkpoint)\n";
      cerr << "-svnDiff             : runs svn diff <src/...../Packages/Matiti \n";
      cerr << "-svnStat             : runs svn stat -u & svn info <src/...../Packages/Matiti \n";
      cerr << "-copy                : Copy from old uda when restarting\n";
      cerr << "-move                : Move from old uda when restarting\n";
      cerr << "-nocopy              : Default: Don't copy or move old uda timestep when\n\t\t\trestarting\n";
      cerr << "-validate            : Verifies the .ups file is valid and quits!\n";
      cerr << "-do_not_validate     : Skips .ups file validation! Please avoid this flag if at all possible.\n";
      cerr << "-track               : Turns on (external) simulation tracking... continues w/o tracking if connection fails.\n";
      cerr << "-TRACK               : Turns on (external) simulation tracking... dies if connection fails.\n";
      cerr << "\n\n";
    }
  quit();
}


#include <iomanip>
int
main( int argc, char *argv[], char *env[] )
{
  string oldTag;
  bool restart = false;
  bool restartFromScratch = true;
  bool restartRemoveOldDir = false;
  bool validateUps = true;
  string filename;

#if HAVE_IEEEFP_H
  fpsetmask(FP_X_OFL|FP_X_DZ|FP_X_INV);
#endif
#if 0
  feenableexcept(FE_INVALID|FE_OVERFLOW|FE_DIVBYZERO);
#endif

  //  Parse arguments
  for(int ii=1; ii<argc; ii++) {
    string arg = argv[ii];
    if( (arg == "-help") || (arg == "-h") ) {
      usage( "", "", argv[0]);
    } elseif (arg == "-restart") {
      restart=true;
    } else if(arg == "-copy") {
      restartFromScratch = false;
      restartRemoveOldDir = false;
    } else if(arg == "-move") {
      restartFromScratch = false;
      restartRemoveOldDir = true;
    } else if(arg == "-t") {
      if (ii < argc-1) {
        restartTimestep = atoi(argv[++ii]);
      }
    } else {
      if( filename != "" ) {
        usage("", arg, argv[0]);
      }
      else if( argv[i][0] == '-' ) { // Don't allow 'filename' to begin with '-'.
        usage("Error!  It appears that the filename you specified begins with a '-'.\n"
              "        This is not allowed.  Most likely there is problem with your\n"
              "        command line.",
              argv[i], argv[0]);        
      } 
      else {
        filename = argv[i];
      }
    }
  }
 
  if( filename == "" ) {
    usage("No input file specified", "", argv[0]);
  }

  if (restart) {
    // check if state.xml is present
    // if not do normal
    udaDir = filename;
    filename = filename + "/input.xml";

    // If restarting (etc), make sure that the uda specified is not a symbolic link to an Uda.
    // This is because the sym link can (will) be updated to point to a new uda, thus creating
    // an inconsistency.  Therefore it is just better not to use the sym link in the first place.
    if( isSymLink( udaDir.c_str() ) ) {
      cout << "\n";
      cout << "Error: " + udaDir + " is a symbolic link.  Please use the full name of the UDA.\n";
      cout << "\n";
    }
  }

  bool thrownException = false;

  try {

    time_t t = time(NULL) ;
    string time_string(ctime(&t));
    char name[256];
    gethostname(name, 256);
    cout << "Date:    " << time_string; // has its own newline
    cout << "Machine: " << name << "\n";
    cout << "CFLAGS: " << CFLAGS << "\n";

    //__________________________________
    // Read input file
    Uintah::ProblemSpecP ups =  Uintah::ProblemSpecReader().readInputFile( filename, validateUps );

    SimulationController* ctl = 
      new AMRSimulationController(world, do_AMR, ups);

    RegridderCommon* reg = 0;
    if(do_AMR) {
      reg = RegridderFactory::create(ups, world);
      if (reg) {
        ctl->attachPort("regridder", reg);
      }
    }

    //______________________________________________________________________
    // Create the components

    //__________________________________
    // Component
    // try to make it from the command line first, then look in ups file
    MatitiParallelComponent* comp = ComponentFactory::create(ups, world, do_AMR, udaDir);
    SimulationInterface* sim = dynamic_cast<SimulationInterface*>(comp);

    ctl->attachPort("sim", sim);
    comp->attachPort("solver", solve);
    comp->attachPort("regridder", reg);
    
    //__________________________________
    // Load balancer
    LoadBalancerCommon* lbc = LoadBalancerFactory::create(ups, world);
    lbc->attachPort("sim", sim);
    if(reg) {
      reg->attachPort("load balancer", lbc);
      lbc->attachPort("regridder",reg);
    }
    
    //__________________________________
    // Output
    DataArchiver* dataarchiver = new DataArchiver(world, udaSuffix);
    Output* output = dataarchiver;
    ctl->attachPort("output", dataarchiver);
    dataarchiver->attachPort("load balancer", lbc);
    comp->attachPort("output", dataarchiver);
    dataarchiver->attachPort("sim", sim);
    
    //__________________________________
    // Scheduler
    SchedulerCommon* sched = SchedulerFactory::create(ups, world, output);
    sched->attachPort("load balancer", lbc);
    ctl->attachPort("scheduler", sched);
    lbc->attachPort("scheduler", sched);
    comp->attachPort("scheduler", sched);

    sched->setStartAddr( start_addr );
    
    if (reg) {
      reg->attachPort("scheduler", sched);
    }
    sched->addReference();
    
    if (emit_graphs) {
      sched->doEmitTaskGraphDocs();
    }
    
    MALLOC_TRACE_TAG(oldTag);
    /*
     * Start the simulation controller
     */
    if (restart) {
      ctl->doRestart(udaDir, restartTimestep, restartFromScratch, restartRemoveOldDir);
    }
    
    // This gives memory held by the 'ups' back before the simulation starts... Assuming
    // no one else is holding on to it...
    ups = 0;

    ctl->run();
    delete ctl;

    sched->removeReference();
    delete sched;
    if (reg) {
      delete reg;
    }
    delete lbc;
    delete sim;
    delete solve;
    delete output;

  } catch (ProblemSetupException& e) {
    // Don't show a stack trace in the case of ProblemSetupException.
    cerrLock.lock();
    cout << "\n\n" << Matiti::Parallel::getMPIRank() << " Caught exception: " << e.message() << "\n\n";
    cerrLock.unlock();
    thrownException = true;
  } catch (Exception& e) {
    cerrLock.lock();
    cout << "\n\n" << Matiti::Parallel::getMPIRank() << " Caught exception: " << e.message() << "\n\n";
    if(e.stackTrace())
      stackDebug << "Stack trace: " << e.stackTrace() << '\n';
    cerrLock.unlock();
    thrownException = true;
  } catch (std::bad_alloc& e) {
    cerrLock.lock();
    cerr << Matiti::Parallel::getMPIRank() << " Caught std exception 'bad_alloc': " << e.what() << '\n';
    cerrLock.unlock();
    thrownException = true;
  } catch (std::bad_exception& e) {
    cerrLock.lock();
    cerr << Matiti::Parallel::getMPIRank() << " Caught std exception: 'bad_exception'" << e.what() << '\n';
    cerrLock.unlock();
    thrownException = true;
  } catch (std::ios_base::failure& e) {
    cerrLock.lock();
    cerr << Matiti::Parallel::getMPIRank() << " Caught std exception 'ios_base::failure': " << e.what() << '\n';
    cerrLock.unlock();
    thrownException = true;
  } catch (std::runtime_error& e) {
    cerrLock.lock();
    cerr << Matiti::Parallel::getMPIRank() << " Caught std exception 'runtime_error': " << e.what() << '\n';
    cerrLock.unlock();
    thrownException = true;
  } catch (std::exception& e) {
    cerrLock.lock();
    cerr << Matiti::Parallel::getMPIRank() << " Caught std exception: " << e.what() << '\n';
    cerrLock.unlock();
    thrownException = true;
  } catch(...) {
    cerrLock.lock();
    cerr << Matiti::Parallel::getMPIRank() << " Caught unknown exception\n";
    cerrLock.unlock();
    thrownException = true;
  }
  
  Matiti::TypeDescription::deleteAll();
  
  /*
   * Finalize MPI
   */
  Matiti::Parallel::finalizeManager( thrownException ?
                                        Matiti::Parallel::Abort : Matiti::Parallel::NormalShutdown);

  if (thrownException) {
    if( Matiti::Parallel::getMPIRank() == 0 ) {
      cout << "\n\nAN EXCEPTION WAS THROWN... Goodbye.\n\n";
    }
    Thread::exitAll(1);
  }
  
  if( Matiti::Parallel::getMPIRank() == 0 ) {
    cout << "Sus: going down successfully\n";
  }

  // use exitAll(0) since return does not work
  Thread::exitAll(0);
  return 0;

} // end main()

