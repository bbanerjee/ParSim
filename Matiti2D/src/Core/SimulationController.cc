#include <iostream>

using namespace std;
using namespace Matiti;

SimulationController::SimulationController(Uintah::ProblemSpecP pspec) : d_ups(pspec)
{
  d_timeinfo = NULL;

} // end SimulationController constructor

SimulationController::~SimulationController()
{
  delete d_timeinfo;
}

void
SimulationController::run()
{
  // set up global data (time, dataarchive, domain)
  problemSetup();

  // do initialization timestep (set up grid, materials etc.)
  doInitialTimestep();

  bool   first = true;
  double time = 0.0;
  double delt = 0.0;
  double start = 0.0;
  int iterations = 0;
  
  while( ( time < d_timeinfo->maxTime ) && ( iterations < d_timeinfo->maxTimestep )) {
  
     d_timeinfo->d_prev_delt = delt;
     iterations++;
 
     // get delt and adjust it
     delt = d_timeinfo->get_delt();

     // delt adjusted based on timeinfo parameters
     adjustDelT(delt, d_timeinfo->d_prev_delt, first, time );

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
SimulationController::problemSetup()
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

      domain = new Domain;
      domain->problemSetup(d_ups, d_myworld, d_doAMR);
   
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

    // Pass the restart_prob_spec to the problemSetup.  For restarting, 
    // pull the <MaterialProperties> from the restart_prob_spec.  If it is not
    // available, then we will pull the properties from the d_ups instead.
    // Needs to be done before DataArchive::restartInitialize
    d_sim->problemSetup(d_ups, restart_prob_spec, domain, d_sharedState);

    // Finalize the shared state/materials
    d_sharedState->finalizeMaterials();
    
    // done after the sim->problemSetup to get defaults into the
    // input.xml, which it writes along with index.xml
    d_output->initializeOutput(d_ups);
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

