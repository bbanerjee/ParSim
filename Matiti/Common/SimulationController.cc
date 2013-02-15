domaininclude <sci_defs/malloc_defs.h>

#include <Common/SimulationController.h>
#include <Core/Mesh/SimulationState.h>
#include <Core/Mesh/SimulationTime.h>
#include <Core/Mesh/Domain.h>
#include <Core/Mesh/Variables/VarTypes.h>
#include <Core/DataArchive/DataArchive.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <CCA/Ports/ProblemSpecInterface.h>
#include <CCA/Ports/SimulationInterface.h>
#include <CCA/Ports/Output.h>
#include <CCA/Ports/Scheduler.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/OS/ProcessInfo.h>
#include <Core/OS/Dir.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstring>
#include <map>

#define SECONDS_PER_MINUTE 60.0
#define SECONDS_PER_HOUR   3600.0
#define SECONDS_PER_DAY    86400.0
#define SECONDS_PER_WEEK   604800.0
#define SECONDS_PER_YEAR   31536000.0

#define AVERAGE_WINDOW 10
using namespace std;

namespace Matiti {

  struct double_int
  {
    double val;
    int loc;
    double_int(double val, int loc): val(val), loc(loc) {}
    double_int(): val(0), loc(-1) {}
  };

  double stdDeviation(double sum_of_x, double sum_of_x_squares, int n)
  {
    return sqrt((n*sum_of_x_squares - sum_of_x*sum_of_x)/(n*n));
  }

  string toHumanUnits( unsigned long value )
  {
    char tmp[64];
  
    sprintf( tmp, "%.2lf", value / 1000000.0 );
    return tmp;
  }

}

using namespace Matiti;

SimulationController::SimulationController(Uintah::ProblemSpecP pspec) : d_ups(pspec)
{
  d_n = 0;
  d_wallTime = 0;
  d_startTime = 0;
  d_prevWallTime = 0;
  d_movingAverage=0;

  d_restarting = false;
  d_archive = NULL;
  d_sim = 0;
  d_timeinfo = NULL;

  d_domain_ps=d_ups->findBlock("Domain");

} // end SimulationController constructor

SimulationController::~SimulationController()
{
  delete d_timeinfo;
}

void 
SimulationController::doRestart(std::string restartFromDir, int timestep)
{
  d_restarting = true;
  d_fromDir = restartFromDir;
  d_restartTimestep = timestep;
  d_restartFromScratch = true;
  d_restartRemoveOldDir = false;
}

double barrier_times[5]={0};
void
SimulationController::run()
{
  // sets up sharedState, timeinfo, output, scheduler, lb
  preMeshSetup();

  // create domain
  DomainP currentDomain = domainSetup();

  d_scheduler->initialize(1, 1);
  d_scheduler->advanceDataWarehouse(currentDomain, true);
    
  double time;

  // set up sim and finalize sharedState
  // also reload from the DataArchive on restart
  postDomainSetup( currentDomain, time );

  calcStartTime();

  // setup, compile, and run the taskgraph for the initialization timestep
  doInitialTimestep( currentDomain, time );

  setStartSimTime( time );
  initSimulationStatsVars();

  ////////////////////////////////////////////////////////////////////////////
  // The main time loop; here the specified problem is actually getting solved
   
  bool   first = true;
  int    iterations = d_sharedState->getCurrentTopLevelTimeStep();
  double delt = 0;

  double start;
  
  while( ( time < d_timeinfo->maxTime ) &&
          ( iterations < d_timeinfo->maxTimestep ) && 
          ( d_timeinfo->max_wall_time == 0 || getWallTime() < d_timeinfo->max_wall_time )  ) {
  
     TrackerClient::trackEvent( Tracker::TIMESTEP_STARTED, time );

     // Compute number of dataWarehouses - multiplies by the time refinement
     // ratio for each level you increase
     int totalFine=1;
     
     d_sharedState->d_prev_delt = delt;
     iterations++;
 
     // get delt and adjust it
     delt_vartype delt_var;
     DataWarehouse* newDW = d_scheduler->getLastDW();
     newDW->get(delt_var, d_sharedState->get_delt_label());

     delt = delt_var;

     // delt adjusted based on timeinfo parameters
     adjustDelT( delt, d_sharedState->d_prev_delt, first, time );
     newDW->override(delt_vartype(delt), d_sharedState->get_delt_label());

     // After one step (either timestep or initialization) and correction
     // the delta we can finally, finalize our old timestep, eg. 
     // finalize and advance the Datawarehouse
     d_scheduler->advanceDataWarehouse(currentDomain);

     // Put the current time into the shared state so other components
     // can access it.  Also increment (by one) the current time step
     // number so components can tell what timestep they are on. 
     d_sharedState->setElapsedTime( time );
     d_sharedState->incrementCurrentTopLevelTimeStep();

     // Each component has their own init_delt specified.  On a switch
     // from one component to the next, we need to adjust the delt to
     // that specified in the input file.  To detect the switch of components,
     // we compare the old_init_delt before the needRecompile() to the 
     // new_init_delt after the needRecompile().  

     double old_init_delt = d_timeinfo->max_initial_delt;
     double new_init_delt = 0.;

     bool nr;
     if( (nr=needRecompile( time, delt, currentDomain )) || first ){
       if(nr)
       {
         //if needRecompile returns true it has reload balanced and thus we need 
         //to assign the boundary conditions.
          currentDomain->assignBCS(d_domain_ps,d_lb);
          currentDomain->performConsistencyCheck();
       }
       new_init_delt = d_timeinfo->max_initial_delt;
       if (new_init_delt != old_init_delt) {
         // writes to the DW in the next section below
         delt = new_init_delt;
       }
       first = false;
       recompile( time, delt, currentDomain, totalFine );
     }
     else {
       if (d_output){
         // This is not correct if we have switched to a different
         // component, since the delt will be wrong 
         d_output->finalizeTimestep( time, delt, currentDomain, d_scheduler, 0 );
       }
     }

     // adjust the delt for each level and store it in all applicable dws.
     double delt_fine = delt;
     int skip=totalFine;
     for(int i=0;i<currentDomain->numLevels();i++){
       const Level* level = currentDomain->getLevel(i).get_rep();
	DataWarehouse* dw = d_scheduler->get_dw(idw);
	dw->override(delt_vartype(delt_fine), d_sharedState->get_delt_label(),
		      level);
     }
     
     // override for the global level as well (which only matters on dw 0)
     d_scheduler->get_dw(0)->override(delt_vartype(delt),
                                      d_sharedState->get_delt_label());

     calcWallTime();

     printSimulationStats( d_sharedState->getCurrentTopLevelTimeStep()-1, delt, time );

     // Execute the current timestep, restarting if necessary
     d_sharedState->d_current_delt = delt;
     executeTimestep( time, delt, currentDomain, totalFine );
     
     if(d_output){
       d_output->executedTimestep(delt, currentDomain);
     }
     time += delt;
   } // end while ( time )

   // print for the final timestep, as the one above is in the middle of a while loop - get new delt, and set walltime first
   delt_vartype delt_var;
   d_scheduler->getLastDW()->get(delt_var, d_sharedState->get_delt_label());
   delt = delt_var;
   adjustDelT( delt, d_sharedState->d_prev_delt, d_sharedState->getCurrentTopLevelTimeStep(), time );
   calcWallTime();
   printSimulationStats( d_sharedState->getCurrentTopLevelTimeStep(), delt, time );

} // end run()

void 
SimulationController::preDomainSetup( void )
{
  d_sharedState = new SimulationState(d_ups);
    
  d_output = dynamic_cast<Output*>(getPort("output"));
    
  Scheduler* sched = dynamic_cast<Scheduler*>(getPort("scheduler"));
  sched->problemSetup(d_ups, d_sharedState);
  d_scheduler = sched;
    
  if( !d_output ){
    cout << "dynamic_cast of 'd_output' failed!\n";
    throw InternalError("dynamic_cast of 'd_output' failed!", __FILE__, __LINE__);
  }
  d_output->problemSetup(d_ups, d_sharedState.get_rep());

  // Parse time struct
  d_timeinfo = new SimulationTime(d_ups);
  d_sharedState->d_simTime = d_timeinfo;
}

DomainP 
SimulationController::domainSetup( void ) 
{
    DomainP domain;

    if (d_restarting) {
      // create the DataArchive here, and store it, as we use it a few times...
      // We need to read the domain before ProblemSetup, and we can't load all
      // the data until after problemSetup, so we have to do a few 
      // different DataArchive operations

      Dir restartFromDir(d_fromDir);
      Dir checkpointRestartDir = restartFromDir.getSubdir("checkpoints");
      d_archive = new DataArchive(checkpointRestartDir.getName(),
                          d_myworld->myrank(), d_myworld->size());

      vector<int> indices;
      vector<double> times;

      try {
        d_archive->queryTimesteps(indices, times);
      } catch( InternalError & ie ) {
        cerr << "\n";
        cerr << "An internal error was caught while trying to restart:\n";
        cerr << "\n";
        cerr << ie.message() << "\n";
        cerr << "This most likely means that the simulation UDA that you have specified\n";
        cerr << "to use for the restart does not have any checkpoint data in it.  Look\n";
        cerr << "in <uda>/checkpoints/ for timestep directories (t#####/) to verify.\n";
        cerr << "\n";
        Thread::exitAll(1);
      }

      // find the right time to query the domain
      if (d_restartTimestep == 0) {
        d_restartIndex = 0; // timestep == 0 means use the first timestep
        // reset d_restartTimestep to what it really is
        d_restartTimestep = indices[0];
      }
      else if (d_restartTimestep == -1 && indices.size() > 0) {
        d_restartIndex = (unsigned int)(indices.size() - 1); 
        // reset d_restartTimestep to what it really is
        d_restartTimestep = indices[indices.size() - 1];
      }
      else {
        for (int index = 0; index < (int)indices.size(); index++)
          if (indices[index] == d_restartTimestep) {
            d_restartIndex = index;
            break;
          }
      }
      
      if (d_restartIndex == (int) indices.size()) {
        // timestep not found
        ostringstream message;
        message << "Timestep " << d_restartTimestep << " not found";
        throw InternalError(message.str(), __FILE__, __LINE__);
      }
    }

    if (!d_restarting) {
      domain = new Domain;
      domain->problemSetup(d_ups, d_myworld, d_doAMR);
    }
    else {
      domain = d_archive->queryDomain(d_restartIndex, d_ups.get_rep());
    }
    if(domain->numLevels() == 0){
      throw InternalError("No problem (no levels in domain) specified.", __FILE__, __LINE__);
    }
   
    // Print out meta data
    if (d_myworld->myrank() == 0){
      domain->printStatistics();
      amrout << "Restart domain\n" << *domain.get_rep() << endl;
    }

    // set the dimensionality of the problem.
    IntVector low, high, size;
    domain->getLevel(0)->findCellIndexRange(low, high);
    size = high-low - domain->getLevel(0)->getExtraCells()*IntVector(2,2,2);
    d_sharedState->setDimensionality(size[0] > 1, size[1] > 1, size[2] > 1);

    return domain;
}

void 
SimulationController::postDomainSetup( DomainP& domain, double& t)
{
    
    // Initialize the peridynamics components
    d_sim = dynamic_cast<SimulationInterface*>(getPort("sim"));
    if(!d_sim)
      throw InternalError("No simulation component", __FILE__, __LINE__);

    ProblemSpecP restart_prob_spec = 0;

    if (d_restarting) {
      // do these before calling archive->restartInitialize, since problemSetup creates VarLabes the DA needs
      restart_prob_spec = d_archive->getTimestepDoc(d_restartIndex);
    }

    // Pass the restart_prob_spec to the problemSetup.  For restarting, 
    // pull the <MaterialProperties> from the restart_prob_spec.  If it is not
    // available, then we will pull the properties from the d_ups instead.
    // Needs to be done before DataArchive::restartInitialize
    d_sim->problemSetup(d_ups, restart_prob_spec, domain, d_sharedState);

    if (d_restarting) {
      simdbg << "Restarting... loading data\n";    
      d_archive->restartInitialize(d_restartIndex, domain, d_scheduler->get_dw(1), d_lb, &t);
      

      // set prevDelt to what it was in the last simulation.  If in the last 
      // sim we were clamping delt based on the values of prevDelt, then
      // delt will be off if it doesn't match.
      ProblemSpecP timeSpec = restart_prob_spec->findBlock("Time");
      if (timeSpec) {
        d_sharedState->d_prev_delt = 0.0;
        if (!timeSpec->get("oldDelt", d_sharedState->d_prev_delt))
          // the delt is deprecated since it is misleading, but older udas may have it...
          timeSpec->get("delt", d_sharedState->d_prev_delt);
      }

      d_sharedState->setCurrentTopLevelTimeStep( d_restartTimestep );
      // Tell the scheduler the generation of the re-started simulation.
      // (Add +1 because the scheduler will be starting on the next
      // timestep.)
      d_scheduler->setGeneration( d_restartTimestep+1 );
      
      // just in case you want to change the delt on a restart....
      if (d_timeinfo->override_restart_delt != 0) {
        double newdelt = d_timeinfo->override_restart_delt;
        if (d_myworld->myrank() == 0)
          cout << "Overriding restart delt with " << newdelt << endl;
        d_scheduler->get_dw(1)->override(delt_vartype(newdelt), 
                                        d_sharedState->get_delt_label());
        double delt_fine = newdelt;
        for(int i=0;i<domain->numLevels();i++){
          const Level* level = domain->getLevel(i).get_rep();
          if(i != 0 && !d_sharedState->isLockstepAMR()) {
            delt_fine /= level->getRefinementRatioMaxDim();
          }
          d_scheduler->get_dw(1)->override(delt_vartype(delt_fine), d_sharedState->get_delt_label(),
                                          level);
        }
      }
      d_scheduler->get_dw(1)->finalize();
      
      // don't need it anymore...
      delete d_archive;
    }

    // Finalize the shared state/materials
    d_sharedState->finalizeMaterials();
    
    // done after the sim->problemSetup to get defaults into the
    // input.xml, which it writes along with index.xml
    d_output->initializeOutput(d_ups);

    if (d_restarting) {
      Dir dir(d_fromDir);
      d_output->restartSetup(dir, 0, d_restartTimestep, t,
                             d_restartFromScratch, d_restartRemoveOldDir);
    }

}
  
//______________________________________________________________________
void
SimulationController::doInitialTimestep(DomainP& domain, double& t)
{
  double start = Time::currentSeconds();
  d_scheduler->mapDataWarehouse(Task::OldDW, 0);
  d_scheduler->mapDataWarehouse(Task::NewDW, 1);
  d_scheduler->mapDataWarehouse(Task::CoarseOldDW, 0);
  d_scheduler->mapDataWarehouse(Task::CoarseNewDW, 1);
  
  if(d_restarting){
    domain->performConsistencyCheck();
    d_sim->restartInitialize();
  } else {
    d_sharedState->setCurrentTopLevelTimeStep( 0 );
    domain->assignBCS(d_domain_ps,d_lb);
    domain->performConsistencyCheck();
    t = d_timeinfo->initTime;

    cout << "Compiling initialization taskgraph...\n";

    // Initialize the peridynamics data
    d_sim->scheduleInitialize(domain->getLevel(0), d_scheduler);

    scheduleComputeStableTimestep(domain,d_scheduler);

    if(d_output) d_output->finalizeTimestep(t, 0, domain, d_scheduler, 1);
      
    d_scheduler->compile();

    double end = Time::currentSeconds() - start;
    cout << "done taskgraph compile (" << end << " seconds)\n";

    // No scrubbing for initial step
    d_scheduler->get_dw(1)->setScrubbing(DataWarehouse::ScrubNone);
    d_scheduler->execute();

    if(d_output) d_output->executedTimestep(0, domain);
  }
}

//______________________________________________________________________
bool
SimulationController::needRecompile(double time, double delt,
				    const DomainP& domain)
{
  // Currently, d_output, d_sim, can request a recompile.  --bryan
  bool recompile = false;
  
  // do it this way so everybody can have a chance to maintain their state
  recompile |= (d_output && d_output->needRecompile(time, delt, domain));
  recompile |= (d_sim && d_sim->needRecompile(time, delt, domain));
  return recompile;
}

//______________________________________________________________________
void
SimulationController::recompile(double t, double delt, DomainP& currentDomain, int totalFine)
{
  cout << "Compiling taskgraph...\n";
  d_lastRecompileTimestep = d_sharedState->getCurrentTopLevelTimeStep();
  double start = Time::currentSeconds();
  
  d_scheduler->initialize(1, totalFine);
  d_scheduler->fillDataWarehouses(currentDomain);
  
  // Set up new DWs, DW mappings.
  d_scheduler->clearMappings();
  d_scheduler->mapDataWarehouse(Task::OldDW, 0);
  d_scheduler->mapDataWarehouse(Task::NewDW, totalFine);
  d_scheduler->mapDataWarehouse(Task::CoarseOldDW, 0);
  d_scheduler->mapDataWarehouse(Task::CoarseNewDW, totalFine);  
  
  d_scheduler->clearMappings();
  d_scheduler->mapDataWarehouse(Task::OldDW, 0);
  d_scheduler->mapDataWarehouse(Task::NewDW, totalFine);
    
  scheduleComputeStableTimestep(currentDomain, d_scheduler);

  if(d_output){
    d_output->finalizeTimestep(t, delt, currentDomain, d_scheduler, true, d_sharedState->needAddMaterial());
  }
  
  d_scheduler->compile();
 
  double dt=Time::currentSeconds() - start;
  cout << "DONE TASKGRAPH RE-COMPILE (" << dt << " seconds)\n";
  d_sharedState->compilationTime += dt;
}

//______________________________________________________________________
void
SimulationController::executeTimestep(double t, double& delt, DomainP& currentDomain, int totalFine)
{
  // If the timestep needs to be
  // restarted, this loop will execute multiple times.
  bool success = true;
  double orig_delt = delt;
  do {
    bool restartable = d_sim->restartableTimesteps();
    d_scheduler->setRestartable(restartable);
    if (restartable)
      d_scheduler->get_dw(0)->setScrubbing(DataWarehouse::ScrubNonPermanent);
    else
      d_scheduler->get_dw(0)->setScrubbing(DataWarehouse::ScrubComplete);
    
    d_scheduler->get_dw(1)->setScrubbing(DataWarehouse::ScrubNonPermanent);
    
    d_scheduler->execute(0, d_lastRecompileTimestep == d_sharedState->getCurrentTopLevelTimeStep() ? 0 : 1);
    
    //__________________________________
    //  If timestep has been restarted
    if(d_scheduler->get_dw(totalFine)->timestepRestarted()){
      
      // Figure out new delt
      double new_delt = d_sim->recomputeTimestep(delt);

      cout << "Restarting timestep at " << t << ", changing delt from "
           << delt << " to " << new_delt << '\n';
     
      // bulletproofing
      if(new_delt < d_timeinfo->delt_min || new_delt <= 0 ){
        ostringstream warn;
        warn << "The new delT (" << new_delt << ") is either less than delT_min (" << d_timeinfo->delt_min
             << ") or equal to 0";
        throw InternalError(warn.str(), __FILE__, __LINE__);
      }
      
      
      d_output->reEvaluateOutputTimestep(orig_delt, new_delt);
      delt = new_delt;
      
      d_scheduler->get_dw(0)->override(delt_vartype(new_delt),
                                       d_sharedState->get_delt_label());

      double delt_fine = delt;
      int skip=totalFine;
      for(int i=0;i<currentDomain->numLevels();i++){
        const Level* level = currentDomain->getLevel(i).get_rep();
        
        for(int idw=0;idw<totalFine;idw+=skip){
          DataWarehouse* dw = d_scheduler->get_dw(idw);
          dw->override(delt_vartype(delt_fine), d_sharedState->get_delt_label(),
                       level);
        }
      }
      success = false;
      
    } else {
      success = true;
      if(d_scheduler->get_dw(1)->timestepAborted()){
        throw InternalError("Execution aborted, cannot restart timestep\n", __FILE__, __LINE__);
      }
    }
  } while(!success);
} // end executeTimestep()

void
SimulationController::scheduleComputeStableTimestep(const DomainP& domain,
                                                    SchedulerP& sched )
{
  d_sim->scheduleComputeStableTimestep(domain->getLevel(0), sched);
}

void 
SimulationController::adjustDelT(double& delt, double prev_delt, bool first, double t) 
{
    delt *= d_timeinfo->delt_factor;
      
    if(delt < d_timeinfo->delt_min){
      cout << "WARNING: raising delt from " << delt
             << " to minimum: " << d_timeinfo->delt_min << '\n';
      delt = d_timeinfo->delt_min;
    }
    if( !first && 
        d_timeinfo->max_delt_increase < 1.e90 &&
        delt > (1+d_timeinfo->max_delt_increase)*prev_delt) {
      cout << "WARNING (a): lowering delt from " << delt 
             << " to maxmimum: " << (1+d_timeinfo->max_delt_increase)*prev_delt
             << " (maximum increase of " << d_timeinfo->max_delt_increase
             << ")\n";
      delt = (1+d_timeinfo->max_delt_increase)*prev_delt;
    }
    if( t <= d_timeinfo->initial_delt_range && delt > d_timeinfo->max_initial_delt ) {
      cout << "WARNING (b): lowering delt from " << delt 
             << " to maximum: " << d_timeinfo->max_initial_delt
             << " (for initial timesteps)\n";
      delt = d_timeinfo->max_initial_delt;
    }
    if( delt > d_timeinfo->delt_max ) {
      cout << "WARNING (c): lowering delt from " << delt 
             << " to maximum: " << d_timeinfo->delt_max << '\n';
      delt = d_timeinfo->delt_max;
    }
    // clamp timestep to output/checkpoint
    if( d_timeinfo->timestep_clamping && d_output ) {
      double orig_delt = delt;
      double nextOutput = d_output->getNextOutputTime();
      double nextCheckpoint = d_output->getNextCheckpointTime();
      if (nextOutput != 0 && t + delt > nextOutput) {
        delt = nextOutput - t;       
      }
      if (nextCheckpoint != 0 && t + delt > nextCheckpoint) {
        delt = nextCheckpoint - t;
      }
      if (delt != orig_delt) {
        if(d_myworld->myrank() == 0)
          cout << "WARNING (d): lowering delt from " << orig_delt 
               << " to " << delt
               << " to line up with output/checkpoint time\n";
      }
    }
    if (d_timeinfo->end_on_max_time && t + delt > d_timeinfo->maxTime){
       delt = d_timeinfo->maxTime - t;
    }
}

double 
SimulationController::getWallTime  ( void )
{
    return d_wallTime;
}

void 
SimulationController::calcWallTime ( void )
{
    d_wallTime = Time::currentSeconds() - d_startTime;
}

double 
SimulationController::getStartTime ( void )
{
    return d_startTime;
}

void 
SimulationController::calcStartTime ( void )
{
    d_startTime = Time::currentSeconds();
}

void 
SimulationController::setStartSimTime ( double t )
{
    d_startSimTime = t;
}

void 
SimulationController::initSimulationStatsVars ( void )
{
    d_n = 0;
    d_wallTime = 0;
    d_prevWallTime = Time::currentSeconds();
}

void
SimulationController::printSimulationStats ( int timestep, double delt, double time )
{
  unsigned long memuse, highwater, maxMemUse;
  d_scheduler->checkMemoryUse( memuse, highwater, maxMemUse );

  // with the sum reduces, use double, since with memory it is possible that
  // it will overflow
  double avg_memuse = memuse;
  unsigned long max_memuse = memuse;
  int max_memuse_loc = -1;
  double avg_highwater = highwater;
  unsigned long max_highwater = highwater;

  // a little ugly, but do it anyway so we only have to do one reduce for sum and
  // one reduce for max
  std::vector<double> toReduce, avgReduce;
  std::vector<double_int> toReduceMax;
  std::vector<double_int> maxReduce;
  std::vector<const char*> statLabels;
  int rank=d_myworld->myrank();
  double total_time=0, overhead_time=0, percent_overhead=0;
  toReduce.push_back(memuse);
  toReduceMax.push_back(double_int(memuse,rank));
  toReduce.push_back(d_sharedState->compilationTime);
  toReduceMax.push_back(double_int(d_sharedState->compilationTime,rank));
  toReduce.push_back(d_sharedState->taskExecTime);
  toReduceMax.push_back(double_int(d_sharedState->taskExecTime,rank));
  toReduce.push_back(d_sharedState->taskGlobalCommTime);
  toReduceMax.push_back(double_int(d_sharedState->taskGlobalCommTime,rank));
  toReduce.push_back(d_sharedState->taskLocalCommTime);
  toReduceMax.push_back(double_int(d_sharedState->taskLocalCommTime,rank));
  toReduce.push_back(d_sharedState->taskWaitCommTime);
  toReduceMax.push_back(double_int(d_sharedState->taskWaitCommTime,rank));
  toReduce.push_back(d_sharedState->outputTime);
  toReduceMax.push_back(double_int(d_sharedState->outputTime,rank));
  toReduce.push_back(d_sharedState->taskWaitThreadTime);
  toReduceMax.push_back(double_int(d_sharedState->taskWaitThreadTime,rank));
  statLabels.push_back("Mem usage");
  statLabels.push_back("Recompile");
  statLabels.push_back("TaskExec");
  statLabels.push_back("TaskGlobalComm");
  statLabels.push_back("TaskLocalComm");
  statLabels.push_back("TaskWaitCommTime");
  statLabels.push_back("Output");
  statLabels.push_back("TaskWaitThreadTime");

  if (highwater) { // add highwater to the end so we know where everything else is (as highwater is conditional)
    toReduce.push_back(highwater);
  }
  avgReduce.resize(toReduce.size());
  maxReduce.resize(toReduce.size());


  //sum up the times for simulation components
  total_time=d_sharedState->compilationTime
    +d_sharedState->taskExecTime;

  //sum up the average time for overhead related components
  overhead_time=d_sharedState->compilationTime;

  //calculate percentage of time spent in overhead
  percent_overhead=overhead_time/total_time;

  //set the overhead sample
  if(d_n>2)  //ignore the first 3 samples, they are not good samples
  {
    d_sharedState->overhead[d_sharedState->overheadIndex]=percent_overhead;
    //increment the overhead index

    double overhead=0;
    double weight=0;

    int t=min(d_n-2,OVERHEAD_WINDOW);
    //calcualte total weight by incrementing through the overhead sample array backwards and multiplying samples by the weights
    for(int i=0;i<t;i++)
    {
      overhead+=d_sharedState->overhead[(d_sharedState->overheadIndex+OVERHEAD_WINDOW-i)%OVERHEAD_WINDOW]*d_sharedState->overheadWeights[i];
      weight+=d_sharedState->overheadWeights[i];
    }
    d_sharedState->overheadAvg=overhead/weight; 

    d_sharedState->overheadIndex=(d_sharedState->overheadIndex+1)%OVERHEAD_WINDOW;
    //increase overhead size if needed
  } 
  d_sharedState->clearStats();

  // calculate mean/std dev
  //double stdDev = 0;
  double mean = 0;
  double walltime = d_wallTime-d_prevWallTime;

  if (d_n > 2) { // ignore times 0,1,2
    //walltimes.push_back();
    //d_sumOfWallTimes += (walltime);
    //d_sumOfWallTimeSquares += pow(walltime,2);

    //alpha=2/(N+1)
    float alpha=2.0/(min(d_n-2,AVERAGE_WINDOW)+1);  
    d_movingAverage=alpha*walltime+(1-alpha)*d_movingAverage;
    mean=d_movingAverage;

  }
  /*
     if (d_n > 3) {
// divide by n-2 and not n, because we wait till n>2 to keep track
// of our stats
stdDev = stdDeviation(d_sumOfWallTimes, d_sumOfWallTimeSquares, d_n-2);
//mean = d_sumOfWallTimes / (d_n-2);
//         ofstream timefile("avg_elapsed_wallTime.txt");
//         timefile << mean << " +- " << stdDev << "\n";
}
*/

// output timestep statistics

if (istats.active()) {
  for (unsigned i = 1; i < statLabels.size(); i++) { // index 0 is memuse
    if (toReduce[i] > 0)
      istats << "rank: " << d_myworld->myrank() << " " << statLabels[i] << " avg: " << toReduce[i] << "\n";
  }
} 

  char walltime[96];
  if (d_n > 3) {
    //sprintf(walltime, ", elap T = %.2lf, mean: %.2lf +- %.3lf", d_wallTime, mean, stdDev);
    sprintf(walltime, ", elap T = %.2lf, mean: %.2lf", d_wallTime, mean);
  }
  else {
    sprintf(walltime, ", elap T = %.2lf", d_wallTime);
  }
  ostringstream message;

  message << "Time="         << time
    << " (timestep "  << timestep 
    << "), delT="     << delt
    << walltime;
  message << ", Mem Use (MB)= ";
  if (avg_memuse == max_memuse && avg_highwater == max_highwater) {
    message << toHumanUnits((unsigned long) avg_memuse);
    if(avg_highwater) {
      message << "/" << toHumanUnits((unsigned long) avg_highwater);
    }

  dbg << message.str() << "\n";
  dbg.flush();
  cout.flush();

  if (stats.active()) {
    for (unsigned i = 1; i < statLabels.size(); i++) { // index 0 is memuse
      if (toReduce[i] > 0)
        stats << statLabels[i] << " avg: " << toReduce[i] << " max: " << toReduce[i] << " maxloc:" << 0
            << " LIB%: " << 0 << "\n";
    }

    if(d_n>2 && !isnan(d_sharedState->overheadAvg))
      stats << "Percent Time in overhead:" << d_sharedState->overheadAvg*100 <<  "\n";
  } 


  if ( d_n > 0 ) {
    double realSecondsNow = (d_wallTime - d_prevWallTime)/delt;
    double realSecondsAvg = (d_wallTime - d_startTime)/(time-d_startSimTime);

    dbgTime << "1 sim second takes ";

    dbgTime << left << showpoint << setprecision(3) << setw(4);

    if (realSecondsNow < SECONDS_PER_MINUTE) {
      dbgTime << realSecondsNow << " seconds (now), ";
    } else if ( realSecondsNow < SECONDS_PER_HOUR ) {
      dbgTime << realSecondsNow/SECONDS_PER_MINUTE << " minutes (now), ";
    } else if ( realSecondsNow < SECONDS_PER_DAY  ) {
      dbgTime << realSecondsNow/SECONDS_PER_HOUR << " hours (now), ";
    } else if ( realSecondsNow < SECONDS_PER_WEEK ) {
      dbgTime << realSecondsNow/SECONDS_PER_DAY << " days (now), ";
    } else if ( realSecondsNow < SECONDS_PER_YEAR ) {
      dbgTime << realSecondsNow/SECONDS_PER_WEEK << " weeks (now), ";
    } else {
      dbgTime << realSecondsNow/SECONDS_PER_YEAR << " years (now), ";
    }

    dbgTime << setw(4);

    if (realSecondsAvg < SECONDS_PER_MINUTE) {
      dbgTime << realSecondsAvg << " seconds (avg) ";
    } else if ( realSecondsAvg < SECONDS_PER_HOUR ) {
      dbgTime << realSecondsAvg/SECONDS_PER_MINUTE << " minutes (avg) ";
    } else if ( realSecondsAvg < SECONDS_PER_DAY  ) {
      dbgTime << realSecondsAvg/SECONDS_PER_HOUR << " hours (avg) ";
    } else if ( realSecondsAvg < SECONDS_PER_WEEK ) {
      dbgTime << realSecondsAvg/SECONDS_PER_DAY << " days (avg) ";
    } else if ( realSecondsAvg < SECONDS_PER_YEAR ) {
      dbgTime << realSecondsAvg/SECONDS_PER_WEEK << " weeks (avg) ";
    } else {
      dbgTime << realSecondsAvg/SECONDS_PER_YEAR << " years (avg) ";
    }

    dbgTime << "to calculate." << "\n";
  }

  d_prevWallTime = d_wallTime;
}
d_n++;

// Reset mem use tracking variable for next iteration
d_scheduler->resetMaxMemValue();

} // end printSimulationStats()


  
