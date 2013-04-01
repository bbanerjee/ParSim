#include <Vaango/src/Core/Util/FileUtils.h>
#include <string>
#include <iostream>
#include <iomanip>

using namespace Matiti;

static void quit( const std::string & msg = "" )
{
  if (msg != "") {
    std::cerr << msg << "\n";
  }
}

static void usage(const std::string & message,
                  const std::string& badarg,
                  const std::string& progname)
{
  std::cerr << "\n";
  if(badarg != "") {
    std::cerr << "Error parsing argument: " << badarg << '\n';
  }
  std::cerr << "\n";
  std::cerr << message << "\n";
  std::cerr << "\n";
  std::cerr << "Usage: " << progname << " [options] <input_file_name>\n\n";
  std::cerr << "Valid options are:\n";
  std::cerr << "-h[elp]        : This usage information.\n";
  std::cerr << "\n\n";
  quit();
}

int
main( int argc, char *argv[], char *env[] )
{
  string filename;       // Input file name

  //  Parse arguments
  for(int ii=1; ii<argc; ii++) {
    string arg = argv[ii];
    if( (arg == "-help") || (arg == "-h") ) {
      usage( "", "", argv[0]);
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

  bool thrownException = false;

  try {

    time_t t = time(NULL) ;
    std::string time_string(ctime(&t));

    // Read input file
    Uintah::ProblemSpecP ups =  Uintah::ProblemSpecReader().readInputFile( filename, false );

    // Set up simulation controller
    SimulationController* ctl = new SimulationController(ups);

    //______________________________________________________________________
    // Create the components

    //__________________________________
    // Component
    // try to make it from the command line first, then look in ups file
    SerialComponent* comp = ComponentFactory::create(ups, world, do_AMR, udaDir);
    SimulationInterface* sim = dynamic_cast<SimulationInterface*>(comp);

    ctl->attachPort("sim", sim);
    comp->attachPort("solver", solve);
    
    //__________________________________
    // Output
    DataArchiver* dataarchiver = new DataArchiver(world, udaSuffix);
    Output* output = dataarchiver;
    ctl->attachPort("output", dataarchiver);
    comp->attachPort("output", dataarchiver);
    dataarchiver->attachPort("sim", sim);
    
    //__________________________________
    // Scheduler
    SchedulerCommon* sched = SchedulerFactory::create(ups, world, output);
    ctl->attachPort("scheduler", sched);
    comp->attachPort("scheduler", sched);

    sched->setStartAddr( start_addr );
    
    if (reg) {
      reg->attachPort("scheduler", sched);
    }
    sched->addReference();
    
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
    delete sim;
    delete solve;
    delete output;

  } catch (ProblemSetupException& e) {
    // Don't show a stack trace in the case of ProblemSetupException.
    std::cout << "\n\n Caught exception: " << e.message() << "\n\n";
    thrownException = true;
  } catch (Exception& e) {
    std::cout << "\n\n Caught exception: " << e.message() << "\n\n";
    if(e.stackTrace())
      stackDebug << "Stack trace: " << e.stackTrace() << '\n';
    thrownException = true;
  } catch (std::bad_alloc& e) {
    std::cerr << " Caught std exception 'bad_alloc': " << e.what() << '\n';
    thrownException = true;
  } catch (std::bad_exception& e) {
    std::cerr << " Caught std exception: 'bad_exception'" << e.what() << '\n';
    thrownException = true;
  } catch (std::ios_base::failure& e) {
    std::cerr <<  " Caught std exception 'ios_base::failure': " << e.what() << '\n';
    thrownException = true;
  } catch (std::runtime_error& e) {
    std::cerr << " Caught std exception 'runtime_error': " << e.what() << '\n';
    thrownException = true;
  } catch (std::exception& e) {
    std::cerr <<  " Caught std exception: " << e.what() << '\n';
    thrownException = true;
  } catch(...) {
    std::cerr << " Caught unknown exception\n";
    thrownException = true;
  }
  
  Matiti::TypeDescription::deleteAll();
  
  if (thrownException) {
    std::cout << "\n\nAN EXCEPTION WAS THROWN... Goodbye.\n\n";
  }
  
  std::cout << "Sus: going down successfully\n";

  return 0;

} // end main()

