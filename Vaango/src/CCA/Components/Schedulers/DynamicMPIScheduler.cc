/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2022    Parresia Research Limited, New Zealand
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

#include <CCA/Components/Schedulers/DynamicMPIScheduler.h>
#include <CCA/Components/Schedulers/OnDemandDataWarehouse.h>
#include <CCA/Components/Schedulers/TaskGraph.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Thread/Time.h>
#include <Core/Thread/Mutex.h>

#include   <cstring>

using namespace Uintah;
using namespace Uintah;

// Debug: Used to sync cerr so it is readable (when output by
// multiple threads at the same time)  From sus.cc:
extern Uintah::Mutex coutLock;
extern Uintah::Mutex cerrLock;

extern DebugStream taskdbg;
extern DebugStream taskorder;

static DebugStream dbg("DynamicMPI_DBG", false);
static DebugStream timeout("DynamicMPI_TimingsOut", false);
static DebugStream queuelength("DynamicMPI_QueueLength", false);

DynamicMPIScheduler::DynamicMPIScheduler(const ProcessorGroup* myworld,
                                         const Output* oport,
                                         DynamicMPIScheduler* parentScheduler)
  : MPIScheduler( myworld, oport, parentScheduler)
{
  taskQueueAlg_ =  MostMessages;

  if (timeout.active()) {
    char filename[64];
    sprintf(filename, "timingStats.%d", d_myworld->myRank());
    timingStats.open(filename);
    if (d_myworld->myRank() == 0) {
      sprintf(filename, "timingStats.avg");
      avgStats.open(filename);
      sprintf(filename, "timingStats.max");
      maxStats.open(filename);
    }
  }
}

DynamicMPIScheduler::~DynamicMPIScheduler()
{
  if (timeout.active()) {
    timingStats.close();
    if (d_myworld->myRank() == 0) {
      avgStats.close();
      maxStats.close();
    }
  }
}

void
DynamicMPIScheduler::problemSetup(const ProblemSpecP& prob_spec,
                                  SimulationStateP& state)
{
  std::string taskQueueAlg = "";

  ProblemSpecP params = prob_spec->findBlock("Scheduler");
  if(params){
    params->get("taskReadyQueueAlg", taskQueueAlg);
  }
  if (taskQueueAlg == "") 
    taskQueueAlg = "MostMessages"; //default taskReadyQueueAlg

  if (taskQueueAlg == "FCFS") 
    taskQueueAlg_ =  FCFS;
  else if (taskQueueAlg == "Random")
    taskQueueAlg_ =  Random;
  else if (taskQueueAlg == "Stack")
    taskQueueAlg_ =  Stack;
  else if (taskQueueAlg == "MostChildren")
    taskQueueAlg_ =  MostChildren;
  else if (taskQueueAlg == "LeastChildren")
    taskQueueAlg_ =  LeastChildren;
  else if (taskQueueAlg == "MostAllChildren")
    taskQueueAlg_ =  MostChildren;
  else if (taskQueueAlg == "LeastAllChildren")
    taskQueueAlg_ =  LeastChildren;
  else if (taskQueueAlg == "MostL2Children")
    taskQueueAlg_ =  MostL2Children;
  else if (taskQueueAlg == "LeastL2Children")
    taskQueueAlg_ =  LeastL2Children;
  else if (taskQueueAlg == "MostMessages")
    taskQueueAlg_ =  MostMessages;
  else if (taskQueueAlg == "LeastMessages")
    taskQueueAlg_ =  LeastMessages;
  else if (taskQueueAlg == "PatchOrder")
    taskQueueAlg_ =  PatchOrder;
  else if (taskQueueAlg == "PatchOrderRandom")
    taskQueueAlg_ =  PatchOrderRandom;
  else {
    throw ProblemSetupException("Unknown task ready queue algorithm", __FILE__, __LINE__);
  }
  log.problemSetup(prob_spec);
  SchedulerCommon::problemSetup(prob_spec, state);
}

SchedulerP
DynamicMPIScheduler::createSubScheduler()
{
  DynamicMPIScheduler* newsched = scinew DynamicMPIScheduler(d_myworld, m_outPort, this);
  newsched->d_sharedState = d_sharedState;
  UintahParallelPort* lbp = getPort("load balancer");
  newsched->attachPort("load balancer", lbp);
  newsched->d_sharedState=d_sharedState;
  return newsched;
}


void
DynamicMPIScheduler::execute(int tgnum /*=0*/, int iteration /*=0*/)
{
  if (d_sharedState->isCopyDataTimestep()) {
    MPIScheduler::execute(tgnum, iteration);
    return;
  }

  MALLOC_TRACE_TAG_SCOPE("DynamicMPIScheduler::execute");

  ASSERTRANGE(tgnum, 0, (int)d_graphs.size());
  TaskGraph* tg = d_graphs[tgnum];
  tg->setIteration(iteration);
  d_currentTG = tgnum;

  if (d_graphs.size() > 1) {
    // tg model is the multi TG model, where each graph is going to need to
    // have its dwmap reset here (even with the same tgnum)
    tg->remapTaskDWs(d_dwmap);
  }

  DetailedTasks* dts = tg->getDetailedTasks();

  if(dts == 0) {
    if (d_myworld->myRank() == 0) {
      std::cerr << "DynamicMPIScheduler skipping execute, no tasks\n";
    }
    return;
  }

  int ntasks = dts->numLocalTasks();
  dts->initializeScrubs(d_dws, d_dwmap);
  dts->initTimestep();

  for (int i = 0; i < ntasks; i++) {
    dts->localTask(i)->resetDependencyCounts();
  }

  if(timeout.active()) {
    d_labels.clear();
    d_times.clear();
    //emitTime("time since last execute");
  }

  int me = d_myworld->myRank();
  makeTaskGraphDoc(dts, me);

  //if(timeout.active())
  //emitTime("taskGraph output");

  mpi_info_.totalreduce = 0;
  mpi_info_.totalsend = 0;
  mpi_info_.totalrecv = 0;
  mpi_info_.totaltask = 0;
  mpi_info_.totalreducempi = 0;
  mpi_info_.totalsendmpi = 0;
  mpi_info_.totalrecvmpi = 0;
  mpi_info_.totaltestmpi = 0;
  mpi_info_.totalwaitmpi = 0;

  int numTasksDone = 0;
  bool abort=false;
  int abort_point = 987654;

  int i = 0;

  if (d_reloc_new_posLabel && d_dws[d_dwmap[Task::OldDW]] != 0) {
    d_dws[d_dwmap[Task::OldDW]]->exchangeParticleQuantities(dts, getLoadBalancer(), d_reloc_new_posLabel, iteration);
  }

#if 0
  // hook to post all the messages up front
  if (useDynamicScheduling_ && !d_sharedState->isCopyDataTimestep()) {
    // post the receives in advance
    for (int i = 0; i < ntasks; i++)
      initiateTask( dts->localTask(i), abort, abort_point, iteration );
  }
#endif

  int currphase=0;
  std::map<int, int> phaseTasks;
  std::map<int, int> phaseTasksDone;
  std::map<int,  DetailedTask *> phaseSyncTask;
  dts->setTaskPriorityAlg(taskQueueAlg_ );

  for (int i = 0; i < ntasks; i++) {
    phaseTasks[dts->localTask(i)->getTask()->d_phase]++;
  }
  
  if( dbg.active()) {
    
    cerrLock.lock();
    dbg << me << " Executing " << dts->numTasks() << " tasks (" 
        << ntasks << " local)";
    for (auto it=phaseTasks.begin() ; it != phaseTasks.end(); it++ ) {
      dbg << ", phase["<< (*it).first << "] = " << (*it).second;
    }
    dbg << endl;
    cerrLock.unlock();
  }

  static std::vector<int> histogram;
  static int totaltasks;
  std::set<DetailedTask*> pending_tasks;

  while( numTasksDone < ntasks) {
    i++;

    // 
    // The following checkMemoryUse() is commented out to allow for
    // maintaining the same functionality as before this commit...
    // In other words, so that memory highwater checking is only done
    // at the end of a timestep, and not between tasks... Once the
    // RT settles down we will uncomment this section and then
    // memory use checks will occur before every task.
    //
    // Note, the results (memuse, highwater, maxMemUse) from the following
    // checkMemoryUse call are not used... the call, however, records
    // the maxMemUse for future reference, and that is why we are calling
    // it.
    //
    //unsigned long memuse, highwater, maxMemUse;
    //checkMemoryUse( memuse, highwater, maxMemUse );

    DetailedTask * task = 0;

    // if we have an internally-ready task, initiate its recvs
    while(dts->numInternalReadyTasks() > 0) { 
      DetailedTask * task = dts->getNextInternalReadyTask();

      if ((task->getTask()->getType() == Task::Reduction) || (task->getTask()->usesMPI())) {  //save the reduction task for later
        phaseSyncTask[task->getTask()->d_phase] = task;
        taskdbg << d_myworld->myRank() << " Task Reduction ready " << *task << " deps needed: " << task->getExternalDepCount() << endl;
      } else {
        initiateTask(task, abort, abort_point, iteration);
        task->markInitiated();
        task->checkExternalDepCount();
        taskdbg << d_myworld->myRank() << " Task internal ready " << *task << " deps needed: " << task->getExternalDepCount() << endl;
        // if MPI has completed, it will run on the next iteration
        pending_tasks.insert(task);
      }
    }

    if (dts->numExternalReadyTasks() > 0) {
      // run a task that has its communication complete
      // tasks get in this queue automatically when their receive count hits 0
      //   in DependencyBatch::received, which is called when a message is delivered.
      if(queuelength.active())
      {
        if((int)histogram.size()<dts->numExternalReadyTasks()+1)
          histogram.resize(dts->numExternalReadyTasks()+1);
        histogram[dts->numExternalReadyTasks()]++;
      }
     
      DetailedTask * task = dts->getNextExternalReadyTask();
      if (taskdbg.active()) {
        cerrLock.lock();
        taskdbg << d_myworld->myRank() << " Running task " << *task << "(" << dts->numExternalReadyTasks() <<"/"<< pending_tasks.size() <<" tasks in queue)"<<endl;
        cerrLock.unlock();
      }

      pending_tasks.erase(pending_tasks.find(task));
      ASSERTEQ(task->getExternalDepCount(), 0);
      runTask(task, iteration);
      numTasksDone++;
      if (taskorder.active()) {
        if (d_myworld->myRank() == d_myworld->nRanks() / 2) {
          cerrLock.lock();
          taskorder << d_myworld->myRank() << " Running task static order: " << task->getStaticOrder() << " , scheduled order: "
                    << numTasksDone << std::endl;
          cerrLock.unlock();
        }
      }
      phaseTasksDone[task->getTask()->d_phase]++;
      //cout << d_myworld->myRank() << " finished task(0) " << *task << " scheduled in phase: " << task->getTask()->d_phase << ", tasks finished in that phase: " <<  phaseTasksDone[task->getTask()->d_phase] << " current phase:" << currphase << endl; 
    } 

    if ((phaseSyncTask.find(currphase)!= phaseSyncTask.end()) && (phaseTasksDone[currphase] == phaseTasks[currphase]-1)){ //if it is time to run the reduction task
      if(queuelength.active())
      {
        if((int)histogram.size()<dts->numExternalReadyTasks()+1)
          histogram.resize(dts->numExternalReadyTasks()+1);
        histogram[dts->numExternalReadyTasks()]++;
      }
      DetailedTask *reducetask = phaseSyncTask[currphase];
      if (reducetask->getTask()->getType() == Task::Reduction){
        if(!abort) {
          cerrLock.lock();
          taskdbg << d_myworld->myRank() << " Running Reduce task " << reducetask->getTask()->getName() << endl;
          cerrLock.unlock();
        }
        initiateReduction(reducetask);
      }
      else { // Task::OncePerProc task
        ASSERT(reducetask->getTask()->usesMPI());
        initiateTask( reducetask, abort, abort_point, iteration );
        reducetask->markInitiated();
        ASSERT(reducetask->getExternalDepCount() == 0);
        runTask(reducetask, iteration);
        if (taskdbg.active()) {
          cerrLock.lock();
          taskdbg << d_myworld->myRank() << " Runnding OPP task:  \t";
          printTask(taskdbg, reducetask); 
          taskdbg << '\n';
          cerrLock.unlock();
        }
      }
      ASSERT(reducetask->getTask()->d_phase==currphase);

      numTasksDone++;
      if (taskorder.active()) {
        if (d_myworld->myRank() == d_myworld->nRanks() / 2) {
          taskorder << d_myworld->myRank() << " Running task static order: " << reducetask->getStaticOrder()
                    << " , scheduled order: " << numTasksDone << std::endl;
        }
      }
      phaseTasksDone[reducetask->getTask()->d_phase]++;
    }

    if (numTasksDone < ntasks){
      if (phaseTasks[currphase] == phaseTasksDone[currphase]) {
        currphase++;
      } else if(dts->numExternalReadyTasks()>0 || dts->numInternalReadyTasks()>0 ||
                (phaseSyncTask.find(currphase)!= phaseSyncTask.end() && 
                 phaseTasksDone[currphase] == phaseTasks[currphase]-1) ) //if there is work to do
      {
        processMPIRecvs(TEST);  // receive what is ready and do not block
      }
      else
      {
        // we have nothing to do, so wait until we get something
        processMPIRecvs(WAIT_ONCE); //There is no other work to do so block until some receives are completed
      }
    }

    if(!abort && d_dws[d_dws.size()-1] && d_dws[d_dws.size()-1]->timestepAborted()){
      // TODO - abort might not work with external queue...
      abort = true;
      abort_point = task->getTask()->getSortedOrder();
      dbg << "Aborting timestep after task: " << *task->getTask() << '\n';
    }
  } // end while( numTasksDone < ntasks )

  if(queuelength.active())
  {
    float lengthsum=0;
    totaltasks += ntasks;
    for (unsigned int i=1; i<histogram.size(); i++)
    {
      lengthsum = lengthsum + i*histogram[i];
    }

    float queuelength = lengthsum/totaltasks;
    float allqueuelength = 0;

    MPI_Reduce(&queuelength, &allqueuelength, 1 , MPI_FLOAT, MPI_SUM, 0, d_myworld->getComm());
    proc0cout  << "average queue length:" << allqueuelength/d_myworld->nRanks() << std::endl;
  }
  
  if(timeout.active()){
    emitTime("MPI send time", mpi_info_.totalsendmpi);
    emitTime("MPI Testsome time", mpi_info_.totaltestmpi);
    emitTime("Total send time", 
             mpi_info_.totalsend - mpi_info_.totalsendmpi - mpi_info_.totaltestmpi);
    emitTime("MPI recv time", mpi_info_.totalrecvmpi);
    emitTime("MPI wait time", mpi_info_.totalwaitmpi);
    emitTime("Total recv time", 
             mpi_info_.totalrecv - mpi_info_.totalrecvmpi - mpi_info_.totalwaitmpi);
    emitTime("Total task time", mpi_info_.totaltask);
    emitTime("Total MPI reduce time", mpi_info_.totalreducempi);
    emitTime("Total reduction time", 
             mpi_info_.totalreduce - mpi_info_.totalreducempi);
    emitTime("Total comm time", 
             mpi_info_.totalrecv + mpi_info_.totalsend + mpi_info_.totalreduce);

    double time      = Time::currentSeconds();
    double totalexec = time - d_lasttime;
    
    d_lasttime = time;

    emitTime("Other excution time", totalexec - mpi_info_.totalsend -
             mpi_info_.totalrecv - mpi_info_.totaltask - mpi_info_.totalreduce);
  }

  if (d_sharedState != 0) { // subschedulers don't have a sharedState
    d_sharedState->taskExecTime += mpi_info_.totaltask - d_sharedState->outputTime; // don't count output time...
    d_sharedState->taskLocalCommTime += mpi_info_.totalrecv + mpi_info_.totalsend;
    d_sharedState->taskWaitCommTime += mpi_info_.totalwaitmpi;
    d_sharedState->taskGlobalCommTime += mpi_info_.totalreduce;
  }

  // Don't need to lock sends 'cause all threads are done at this point.
  sends_[0].waitall(d_myworld);
  ASSERT(sends_[0].numRequests() == 0);

  if(d_restartable && tgnum == (int) d_graphs.size() -1) {
    // Copy the restart flag to all processors
    int myrestart = d_dws[d_dws.size()-1]->timestepRestarted();
    int netrestart;
    MPI_Allreduce(&myrestart, &netrestart, 1, MPI_INT, MPI_LOR,
                  d_myworld->getComm());
    if(netrestart) {
      d_dws[d_dws.size()-1]->restartTimestep();
      if (d_dws[0])
        d_dws[0]->setRestarted();
    }
  }

  finalizeTimestep();
  log.finishTimestep();

  if( timeout.active() && !parentScheduler_ ) {  // only do on toplevel scheduler
    outputTimingStats("DynamicMPIScheduler");
  }

  if( dbg.active()) {
    coutLock.lock();
    dbg << me << " DynamicMPIScheduler finished\n";
    coutLock.unlock();
  }

}

