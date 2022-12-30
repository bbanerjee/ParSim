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

#include <CCA/Components/Schedulers/MPIScheduler.h>
#include <CCA/Components/Schedulers/OnDemandDataWarehouse.h>
#include <CCA/Components/Schedulers/SendState.h>
#include <CCA/Components/Schedulers/CommRecMPI.h>
#include <CCA/Components/Schedulers/DetailedTasks.h>
#include <CCA/Components/Schedulers/TaskGraph.h>
#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Output.h>

#include <Core/Parallel/ProcessorGroup.h>
#include <Core/Grid/Variables/ParticleSubset.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Thread/Time.h>
#include <Core/Thread/Mutex.h>
#include <Core/Util/DebugStream.h>
#include <Core/Util/FancyAssert.h>

#include <sci_defs/mpi_defs.h> // For MPIPP_H on SGI

#include   <sstream>
#include   <iomanip>
#include   <map>
#include   <cstring>

// Pack data into a buffer before sending -- testing to see if this
// works better and avoids certain problems possible when you allow
// tasks to modify data that may have a pending send.
#define USE_PACKING

using namespace Uintah;

// Debug: Used to sync cerr so it is readable (when output by
// multiple threads at the same time)  From sus.cc:
extern Uintah::Mutex coutLock;
extern Uintah::Mutex cerrLock;

static DebugStream dbg("MPIScheduler", false);
static DebugStream dbgst("SendTiming", false);
static DebugStream timeout("MPIScheduler.timings", false);
static DebugStream reductionout("ReductionTasks",  false);

DebugStream taskorder("TaskOrder", false);
DebugStream waitout("WaitTimes", false);
DebugStream execout("ExecTimes", false);
DebugStream taskdbg("TaskDBG", false);
DebugStream taskLevel_dbg("TaskLevel", false);
DebugStream mpidbg("MPIDBG", false);

static double CurrentWaitTime=0;

std::map<std::string,double> waittimes;
std::map<std::string,double> exectimes;


MPIScheduler::MPIScheduler( const ProcessorGroup * myworld,
                            const Output         * oport,
                            MPIScheduler   * inParentScheduler) :
  SchedulerCommon( myworld, oport ),
  parentScheduler_( inParentScheduler ), 
  log( myworld, oport ), 
  oport_(oport),
  numMessages_(0),
  messageVolume_(0),
  recvLock("MPI receive lock"),
  sendLock("MPI send lock"),
  dlbLock("loadbalancer lock"),
  waittimesLock("waittimes lock")
{
  d_lasttime=Time::currentSeconds();
  d_reloc_new_posLabel=0;

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


void
MPIScheduler::problemSetup(const ProblemSpecP& prob_spec,
                           SimulationStateP& state)
{
  log.problemSetup(prob_spec);
  SchedulerCommon::problemSetup(prob_spec, state);
}

MPIScheduler::~MPIScheduler()
{
  if (timeout.active()) {
    timingStats.close();
    if (d_myworld->myRank() == 0) {
      avgStats.close();
      maxStats.close();
    }
  }
}

SchedulerP
MPIScheduler::createSubScheduler()
{
  MPIScheduler* newsched = scinew MPIScheduler(d_myworld, m_outPort, this);
  UintahParallelPort* lbp = getPort("load balancer");
  newsched->attachPort("load balancer", lbp);
  newsched->d_sharedState=d_sharedState;
  return newsched;
}

void
MPIScheduler::verifyChecksum()
{
#if SCI_ASSERTION_LEVEL >= 3

    // Compute a simple checksum to make sure that all processes
    // are trying to execute the same graph.  We should do two
    // things in the future:
    //  - make a flag to turn this off
    //  - make the checksum more sophisticated
    int checksum = 0;
    int numSpatialTasks = 0;
    for (unsigned i = 0; i < d_graphs.size(); i++) {
      checksum += d_graphs[i]->getTasks().size();

      // This begins addressing the issue of making the global checksum more sophisticated:
      //   check if any tasks were spatially scheduled - TaskType::Spatial, meaning no computes, requires or modifies
      //     e.g. RMCRT radiometer task, which is not scheduled on all patches
      //          these Spatial tasks won't count toward the global checksum
      std::vector<Task*> tasks = d_graphs[i]->getTasks();
      std::vector<Task*>::const_iterator tasks_iter = tasks.begin();
      for (; tasks_iter != tasks.end(); ++tasks_iter) {
        Task* task = *tasks_iter;
        if (task->getType() == Task::Spatial) {
          numSpatialTasks++;
        }
      }
    }

    // Spatial tasks don't count against the global checksum
    checksum -= numSpatialTasks;

    if (mpidbg.active()) {
      coutLock.lock();
      mpidbg << d_myworld->myRank() << " (MPI_Allreduce) Checking checksum of " << checksum << '\n';
      coutLock.unlock();
    }

    int result_checksum;
    MPI_Allreduce(&checksum, &result_checksum, 1, MPI_INT, MPI_MIN,
                  d_myworld->getComm());

    if(checksum != result_checksum){
      std::cerr << "Failed task checksum comparison! Not all processes are executing the same taskgraph \n";
      std::cerr << " Rank- " << d_myworld->myRank() << " of "
                << d_myworld->nRanks() - 1 << ": has sum " << checksum
                << " and global is " << result_checksum << '\n';
      MPI_Abort(d_myworld->getComm(), 1);
    }

    if (mpidbg.active()) {
      coutLock.lock();
      mpidbg << d_myworld->myRank() << " (MPI_Allreduce) Check succeeded\n";
      coutLock.unlock();
    }
#endif
}

void
MPIScheduler::initiateTask( DetailedTask * task,
                            bool           only_old_recvs, 
                            int            abort_point, 
                            int            iteration )
{
  postMPIRecvs(task, only_old_recvs, abort_point, iteration);

  if(only_old_recvs) {
    return;
  }
} // end initiateTask()

void
MPIScheduler::initiateReduction( DetailedTask* task )
{
  if (reductionout.active() && d_myworld->myRank() == 0) {
    coutLock.lock();
    reductionout << "Running Reduction Task: " << task->getName() << std::endl;
    coutLock.unlock();
  }

  double reducestart = Time::currentSeconds();
    
  runReductionTask(task);

  double reduceend = Time::currentSeconds();

  emitNode(task, reducestart, reduceend - reducestart, 0);

  mpi_info_.totalreduce += reduceend - reducestart;
  mpi_info_.totalreducempi += reduceend - reducestart;
}

void
MPIScheduler::runTask( DetailedTask* task, 
                       int           iteration,
                       int           thread_id /*=0*/ )
{
  if(waitout.active()) {
    waittimesLock.lock();
    waittimes[task->getTask()->getName()] += CurrentWaitTime;
    CurrentWaitTime = 0;
    waittimesLock.unlock();
  }

  double taskstart = Time::currentSeconds();

  if (d_trackingVarsPrintLocation & SchedulerCommon::PRINT_BEFORE_EXEC) {
    printTrackedVars(task, SchedulerCommon::PRINT_BEFORE_EXEC);
  }
  std::vector<DataWarehouseP> plain_old_dws(d_dws.size());
  for (int i = 0; i < (int)d_dws.size(); i++) {
    plain_old_dws[i] = d_dws[i].get_rep();
  }

  task->doit(d_myworld, d_dws, plain_old_dws);

  if (d_trackingVarsPrintLocation & SchedulerCommon::PRINT_AFTER_EXEC) {
    printTrackedVars(task, SchedulerCommon::PRINT_AFTER_EXEC);
  }

  double total_task_time = Time::currentSeconds() - taskstart;

  dlbLock.lock();
  {
    if (execout.active()) {
      exectimes[task->getTask()->getName()] += total_task_time;
    }

    // if I do not have a sub scheduler
    if (!task->getTask()->getHasSubScheduler()) {
      //add my task time to the total time
      mpi_info_.totaltask += total_task_time;
      if (!d_sharedState->isCopyDataTimestep() && task->getTask()->getType() != Task::Output) {
        //add contribution for patchlist
        getLoadBalancer()->addContribution(task, total_task_time);
      }
    }
  }
  dlbLock.unlock();

  postMPISends( task, iteration, thread_id );

  task->done(d_dws); // should this be timed with taskstart? - BJW

  double teststart = Time::currentSeconds();

  sends_[thread_id].testsome( d_myworld );

  mpi_info_.totaltestmpi += Time::currentSeconds() - teststart;
 
  // Add subscheduler timings to the parent scheduler and reset subscheduler timings
  if(parentScheduler_) {
    parentScheduler_->mpi_info_.totaltask+=mpi_info_.totaltask;
    parentScheduler_->mpi_info_.totaltestmpi+=mpi_info_.totaltestmpi;
    parentScheduler_->mpi_info_.totalrecv+=mpi_info_.totalrecv;
    parentScheduler_->mpi_info_.totalsend+=mpi_info_.totalsend;
    parentScheduler_->mpi_info_.totalwaitmpi+=mpi_info_.totalwaitmpi;
    parentScheduler_->mpi_info_.totalreduce+=mpi_info_.totalreduce;
    mpi_info_.totalreduce = 0;
    mpi_info_.totalsend = 0;
    mpi_info_.totalrecv = 0;
    mpi_info_.totaltask = 0;
    mpi_info_.totalreducempi = 0;
    mpi_info_.totaltestmpi = 0;
    mpi_info_.totalwaitmpi = 0;
  }

  emitNode(task, taskstart, total_task_time, 0);
  
} // end runTask()

void
MPIScheduler::runReductionTask( DetailedTask* task )
{
  const Task::Dependency* mod = task->getTask()->getModifies();
  ASSERT(!mod->next);
  
  OnDemandDataWarehouse* dw = d_dws[mod->mapDataWarehouse()].get_rep();
  ASSERT(task->getTask()->d_comm>=0);
  dw->reduceMPI(mod->var, mod->reductionLevel, mod->matls, task->getTask()->d_comm);
  task->done(d_dws);
}

void
MPIScheduler::postMPISends( DetailedTask* task, 
                            int            iteration,
                            int            thread_id /*=0*/ )
{
  double sendstart = Time::currentSeconds();
  bool dbg_active = dbg.active();

  int me = d_myworld->myRank();
  if (dbg_active) {
    cerrLock.lock();
    dbg << "Rank-" << me << " postMPISends - task " << *task << '\n';
    cerrLock.unlock();
  }

  int numSend = 0;
  int volSend = 0;

  // Send data to dependendents
  for(DependencyBatch* batch = task->getComputes();
      batch != 0; batch = batch->comp_next){

    // Prepare to send a message
#ifdef USE_PACKING
    PackBufferInfo mpibuff;
#else
    BufferInfo mpibuff;
#endif
    // Create the MPI type
    int to = batch->toTasks.front()->getAssignedResourceIndex();
    ASSERTRANGE(to, 0, d_myworld->nRanks());

    std::ostringstream ostr;
    ostr.clear();
    for(DetailedDep* req = batch->head; req != 0; req = req->next){

      ostr << *req << ' '; // for CommRecMPI::add()

      if ((req->condition == DetailedDep::FirstIteration && iteration > 0) || 
          (req->condition == DetailedDep::SubsequentIterations && iteration == 0)||
          (d_notCopyDataVars.count(req->req->var->getName()) > 0  )){
        // See comment in DetailedDep about CommCondition
        if (dbg_active) {
          cerrLock.lock();
          dbg << "Rank-" << me << "   Ignoring conditional send for " << *req << "\n";
          cerrLock.unlock();
        }
        continue;
      }

      // if we send/recv to an output task, don't send/recv if not an output timestep
      if (req->toTasks.front()->getTask()->getType() == Task::Output && 
          !oport_->isOutputTimestep() && !oport_->isCheckpointTimestep()) {
        if (dbg_active) {
          cerrLock.lock();
          dbg << "Rank-" << me << "   Ignoring non-output-timestep send for " << *req << "\n";
          cerrLock.unlock();
        }
        continue;
      }

      OnDemandDataWarehouse* dw = d_dws[req->req->mapDataWarehouse()].get_rep();

      if (dbg_active) {
        cerrLock.lock();
        {
          dbg << "Rank-" << me << " --> sending " << *req << ", ghost type: " << "\""
              << Ghost::getGhostTypeName(req->req->gtype) << "\", " << "num req ghost "
              << Ghost::getGhostTypeName(req->req->gtype) << ": " << req->req->numGhostCells
              << ", Ghost::direction: " << Ghost::getGhostTypeDir(req->req->gtype)
              << ", from dw " << dw->getID() << '\n';
        }
        cerrLock.unlock();
      }

      // the load balancer is used to determine where data was in the old dw on the prev timestep -
      // pass it in if the particle data is on the old dw
      const VarLabel* posLabel;
      OnDemandDataWarehouse* posDW;
      LoadBalancer* lb = 0;

      if(!d_reloc_new_posLabel && parentScheduler_) {
        posDW = d_dws[req->req->task->mapDataWarehouse(Task::ParentOldDW)].get_rep();
        posLabel = parentScheduler_->d_reloc_new_posLabel;
      } else {
        // on an output task (and only on one) we require particle variables from the NewDW
        if (req->toTasks.front()->getTask()->getType() == Task::Output)
          posDW = d_dws[req->req->task->mapDataWarehouse(Task::NewDW)].get_rep();
        else {
          posDW = d_dws[req->req->task->mapDataWarehouse(Task::OldDW)].get_rep();
          lb = getLoadBalancer();
        }
        posLabel = d_reloc_new_posLabel;
      }
      MPIScheduler* top = this;
      while(top->parentScheduler_) top = top->parentScheduler_;

      dw->sendMPI(batch, posLabel, mpibuff, posDW, req, lb);
    }

    // Post the send
    if(mpibuff.count()>0){
      ASSERT(batch->messageTag > 0);
      double start = Time::currentSeconds();
      void* buf;
      int count;
      MPI_Datatype datatype;
     
#ifdef USE_PACKING
      mpibuff.get_type(buf, count, datatype, d_myworld->getComm());
      mpibuff.pack(d_myworld->getComm(), count);
#else
      mpibuff.get_type(buf, count, datatype);
#endif

      //only send message if size is greather than zero
      //we need this empty message to enforce modify after read dependencies 
      //if(count>0)
      //{
      if (mpidbg.active()) {
        cerrLock.lock();
        mpidbg << "Rank-" << me << " Posting send for message number " << batch->messageTag << " to   rank-" << to << ", length: " << count
               << " (bytes)\n";
        cerrLock.unlock();
      }

      numMessages_++;
      numSend++;
      int typeSize;

      MPI_Type_size(datatype,&typeSize);
      messageVolume_+=count*typeSize;
      volSend += count * typeSize;

      MPI_Request requestid;
      MPI_Isend(buf, count, datatype, to, batch->messageTag,
                d_myworld->getComm(), &requestid);
      int bytes = count;

      sendLock.writeLock();
      sends_[thread_id].add( requestid, bytes, mpibuff.takeSendlist(), ostr.str(), batch->messageTag );
      sendLock.writeUnlock();

      mpi_info_.totalsendmpi += Time::currentSeconds() - start;
      //}
    }
  } // end for (DependencyBatch * batch = task->getComputes() )

  double dsend = Time::currentSeconds()-sendstart;
  mpi_info_.totalsend += dsend;

  if (dbgst.active() && numSend > 0) {
    if (d_myworld->myRank() == d_myworld->nRanks() / 2) {
      if (dbgst.active()) {
        cerrLock.lock();
        dbgst << d_myworld->myRank() << " Time: " << Time::currentSeconds() << " , NumSend= " << numSend << " , VolSend: "
              << volSend << std::endl;
        cerrLock.unlock();
      }
    }
  } 

} // end postMPISends();


struct CompareDep {
  bool operator()(DependencyBatch* a, DependencyBatch* b)
    {
      return a->messageTag < b->messageTag;
    }
};

void
MPIScheduler::postMPIRecvs( DetailedTask* task, 
                            bool          only_old_recvs, 
                            int           abort_point, 
                            int           iteration)
{
  double recvstart = Time::currentSeconds();
  bool dbg_active = dbg.active();


  if( dbg_active) {
    cerrLock.lock();
    dbg << "Rank-" << d_myworld->myRank() << " postMPIRecvs - task " << *task << '\n';
    cerrLock.unlock();
  }

  if (d_trackingVarsPrintLocation & SchedulerCommon::PRINT_BEFORE_COMM)
    printTrackedVars(task, SchedulerCommon::PRINT_BEFORE_COMM);

  // sort the requires, so in case there is a particle send we receive it with
  // the right message tag

  std::vector<DependencyBatch*> sorted_reqs;

  for(auto iter= task->getRequires().begin() ; iter != task->getRequires().end(); iter++) {
    sorted_reqs.push_back(iter->first);
  }

  CompareDep comparator;
  std::sort(sorted_reqs.begin(), sorted_reqs.end(), comparator);
  //std::vector<DependencyBatch*>::iterator sorted_iter = sorted_reqs.begin();

  // Receive any of the foreign requires
  recvLock.writeLock();
  {
    for(auto sorted_iter = sorted_reqs.begin(); sorted_iter != sorted_reqs.end(); sorted_iter++) {
      DependencyBatch* batch = *sorted_iter;

      // The first thread that calls this on the batch will return true
      // while subsequent threads calling this will block and wait for
      // that first thread to receive the data.

      task->incrementExternalDepCount();
      if( !batch->makeMPIRequest() ) {
        if (dbg_active) {
          cerrLock.lock();
          dbg << "Someone else already receiving it\n";
          cerrLock.unlock();
        }
        continue;
      }

      if (only_old_recvs){
        if (dbg_active){
          dbg << "abort analysis: " << batch->fromTask->getTask()->getName()
              << ", so=" << batch->fromTask->getTask()->getSortedOrder()
              << ", abort_point=" << abort_point << '\n';
          if(batch->fromTask->getTask()->getSortedOrder() <= abort_point)
            dbg << "posting MPI recv for pre-abort message " 
                << batch->messageTag << '\n';
        }
        if(!(batch->fromTask->getTask()->getSortedOrder() <= abort_point)) {
          continue;
        }
      }

      // Prepare to receive a message
      BatchReceiveHandler* pBatchRecvHandler = scinew BatchReceiveHandler(batch);
      PackBufferInfo* p_mpibuff = 0;
#ifdef USE_PACKING
      p_mpibuff = scinew PackBufferInfo();
      PackBufferInfo& mpibuff = *p_mpibuff;
#else
      BufferInfo mpibuff;
#endif

      std::ostringstream ostr;
      ostr.clear();
      // Create the MPI type
      for(DetailedDep* req = batch->head; req != 0; req = req->next){

        ostr << *req << ' ';  // for CommRecMPI::add()

        OnDemandDataWarehouse* dw = d_dws[req->req->mapDataWarehouse()].get_rep();
        if ((req->condition == DetailedDep::FirstIteration && iteration > 0) || 
            (req->condition == DetailedDep::SubsequentIterations && iteration == 0) ||
            (d_notCopyDataVars.count(req->req->var->getName()) > 0  )){

          // See comment in DetailedDep about CommCondition
          if (dbg_active) {
            cerrLock.lock();
            dbg << "Rank-" << d_myworld->myRank() << "   Ignoring conditional receive for " << *req << std::endl;
            cerrLock.unlock();
          }
          continue;
        }

        // if we send/recv to an output task, don't send/recv if not an output timestep
        if (req->toTasks.front()->getTask()->getType() == Task::Output && !oport_->isOutputTimestep() 
            && !oport_->isCheckpointTimestep()) {
          cerrLock.lock();
          dbg << "Rank-" << d_myworld->myRank() << "   Ignoring non-output-timestep receive for " << *req << std::endl;
          cerrLock.unlock();
          continue;
        }
        if (dbg_active) {
          cerrLock.lock();
          {
            dbg << "Rank-" << d_myworld->myRank() << " <-- receiving " << *req << ", ghost type: " << "\""
                << Ghost::getGhostTypeName(req->req->gtype) << "\", " << "num req ghost "
                << Ghost::getGhostTypeName(req->req->gtype) << ": " << req->req->numGhostCells
                << ", Ghost::direction: " << Ghost::getGhostTypeDir(req->req->gtype)
                << ", into dw " << dw->getID() << '\n';
          }
          cerrLock.unlock();
        }

        OnDemandDataWarehouse* posDW;
        //const VarLabel* posLabel;

        // the load balancer is used to determine where data was in the old dw on the prev timestep
        // pass it in if the particle data is on the old dw
        LoadBalancer* lb = 0;
        if(!d_reloc_new_posLabel && parentScheduler_){
          posDW = d_dws[req->req->task->mapDataWarehouse(Task::ParentOldDW)].get_rep();
          //posLabel = parentScheduler_->d_reloc_new_posLabel;
        } else {
          // on an output task (and only on one) we require particle variables from the NewDW
          if (req->toTasks.front()->getTask()->getType() == Task::Output) {
            posDW = d_dws[req->req->task->mapDataWarehouse(Task::NewDW)].get_rep();
          } else {
            posDW = d_dws[req->req->task->mapDataWarehouse(Task::OldDW)].get_rep();
            lb = getLoadBalancer();
          }
          //posLabel = d_reloc_new_posLabel;
        }

        MPIScheduler* top = this;
        while(top->parentScheduler_) top = top->parentScheduler_;

        dw->recvMPI(batch, mpibuff, posDW, req, lb);

        if (!req->isNonDataDependency()) {
          d_graphs[d_currentTG]->getDetailedTasks()->setScrubCount(req->req, req->matl, req->fromPatch, d_dws);
        }

      }

      // Post the receive
      if(mpibuff.count()>0){

        ASSERT(batch->messageTag > 0);
        double start = Time::currentSeconds();
        void* buf;
        int count;
        MPI_Datatype datatype;

#ifdef USE_PACKING
        mpibuff.get_type(buf, count, datatype, d_myworld->getComm());
#else
        mpibuff.get_type(buf, count, datatype);
#endif

        //only receive message if size is greater than zero
        //we need this empty message to enforce modify after read dependencies 
        //if(count>0)
        //{
        int from = batch->fromTask->getAssignedResourceIndex();
        ASSERTRANGE(from, 0, d_myworld->nRanks());
        MPI_Request requestid;

        if (mpidbg.active()) {
          cerrLock.lock();
          mpidbg << "Rank-" << d_myworld->myRank() << " Posting recv for message number " << batch->messageTag << " from rank-" << from
                 << ", length: " << count << " (bytes)\n";
          cerrLock.unlock();
        }

        MPI_Irecv(buf, count, datatype, from, batch->messageTag,
                  d_myworld->getComm(), &requestid);
        int bytes = count;
        recvs_.add(requestid, bytes,
                   scinew ReceiveHandler(p_mpibuff, pBatchRecvHandler),
                   ostr.str(), batch->messageTag);
        mpi_info_.totalrecvmpi += Time::currentSeconds() - start;
        /*}
          else
          {
          //no message was sent so clean up buffer and handler
          delete p_mpibuff;
          delete pBatchRecvHandler;
          }*/
      }
      else {
        // Nothing really need to be received, but let everyone else know
        // that it has what is needed (nothing).
        batch->received(d_myworld);
#ifdef USE_PACKING
        // otherwise, these will be deleted after it receives and unpacks
        // the data.
        delete p_mpibuff;
        delete pBatchRecvHandler;
#endif        
      }
    } // end for
  }
  recvLock.writeUnlock();

  double drecv = Time::currentSeconds() - recvstart;
  mpi_info_.totalrecv += drecv;

} // end postMPIRecvs()

void
MPIScheduler::processMPIRecvs(int how_much)
{
  // Should only have external receives in the MixedScheduler version which
  // shouldn't use this function.
  // ASSERT(outstandingExtRecvs.empty());
  if (recvs_.numRequests() == 0) return;
  double start = Time::currentSeconds();

  recvLock.writeLock();
  {
    switch (how_much) {
      case TEST :
        recvs_.testsome(d_myworld);
        break;
      case WAIT_ONCE :
        coutLock.lock();
        mpidbg << "Rank-" << d_myworld->myRank() << " Start waiting once (WAIT_ONCE)...\n";
        coutLock.unlock();

      recvs_.waitsome(d_myworld);

        coutLock.lock();
        mpidbg << "Rank-" << d_myworld->myRank() << " Done waiting once (WAIT_ONCE)...\n";
        coutLock.unlock();
        break;
      case WAIT_ALL :
        // This will allow some receives to be "handled" by their
        // AfterCommincationHandler while waiting for others.
        coutLock.lock();
        mpidbg << "Rank-" << d_myworld->myRank() << "  Start waiting (WAIT_ALL)...\n";
        coutLock.unlock();

      while( (recvs_.numRequests() > 0)) {
        bool keep_waiting = recvs_.waitsome(d_myworld);
        if (!keep_waiting) break;
      }

        coutLock.lock();
        mpidbg << "Rank-" << d_myworld->myRank() << "  Done waiting (WAIT_ALL)...\n";
        coutLock.unlock();
        break;
    } // end switch
  }
  recvLock.writeUnlock();

  mpi_info_.totalwaitmpi+=Time::currentSeconds()-start;
  CurrentWaitTime+=Time::currentSeconds()-start;

} // end processMPIRecvs()

void
MPIScheduler::execute(int tgnum     /*=0*/, 
                      int iteration /*=0*/)
{

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

  if(dts == 0){
    if (d_myworld->myRank() == 0) {
      cerrLock.lock();
      std::cerr << "MPIScheduler skipping execute, no tasks\n";
      cerrLock.unlock();
    }
    return;
  }

  //ASSERT(pg_ == 0);
  //pg_ = pg;
  
  int ntasks = dts->numLocalTasks();
  dts->initializeScrubs(d_dws, d_dwmap);
  dts->initTimestep();

  for (int i = 0; i < ntasks; i++)
    dts->localTask(i)->resetDependencyCounts();

  if(timeout.active()){
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

  if( dbg.active()) {
    cerrLock.lock();
    dbg << me << " Executing " << dts->numTasks() << " tasks (" 
        << ntasks << " local)\n";
    cerrLock.unlock();
  }

  bool abort=false;
  int abort_point = 987654;


  if (d_reloc_new_posLabel && d_dws[d_dwmap[Task::OldDW]] != 0) {
    d_dws[d_dwmap[Task::OldDW]]->exchangeParticleQuantities(dts, getLoadBalancer(), d_reloc_new_posLabel, iteration);
  }

  int i = 0;
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

    DetailedTask * task = dts->getNextInternalReadyTask();

    numTasksDone++;

    if (taskorder.active()) {
      taskorder << d_myworld->myRank() << " Running task static order: " << task->getStaticOrder() << " , scheduled order: "
                << numTasksDone << std::endl;
    }

    if (taskdbg.active()) {
      taskdbg << me << " Initiating task:  \t";
      printTask(taskdbg, task);
      taskdbg << '\n';
    }

    if (task->getTask()->getType() == Task::Reduction) {
      if (!abort) initiateReduction(task);
    } else {
      initiateTask( task, abort, abort_point, iteration );
      processMPIRecvs(WAIT_ALL);
      ASSERT(recvs_.numRequests() == 0);
      runTask(task, iteration);

      //  output to help application developers
      if (taskdbg.active()) {
        taskdbg << d_myworld->myRank() << " Completed task:  \t";
        printTask(taskdbg, task);
        taskdbg << '\n';
        printTaskLevels(d_myworld, taskLevel_dbg, task);
      }
    }
  
    if(!abort && d_dws[d_dws.size()-1] && d_dws[d_dws.size()-1]->timestepAborted()){
      abort = true;
      abort_point = task->getTask()->getSortedOrder();
      dbg << "Aborting timestep after task: " << *task->getTask() << '\n';
    }
  } // end while( numTasksDone < ntasks )

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

    emitTime("Other execution time", totalexec - mpi_info_.totalsend -
             mpi_info_.totalrecv - mpi_info_.totaltask - mpi_info_.totalreduce);
  }

  if (!parentScheduler_) {//if this schedule is root scheduler
    d_sharedState->taskExecTime += mpi_info_.totaltask - d_sharedState->outputTime; // don't count output time...
    d_sharedState->taskLocalCommTime += mpi_info_.totalrecv + mpi_info_.totalsend;
    d_sharedState->taskWaitCommTime += mpi_info_.totalwaitmpi;
    d_sharedState->taskGlobalCommTime += mpi_info_.totalreduce;
  }

  // Don't need to lock sends 'cause all threads are done at this point.
  sends_[0].waitall(d_myworld);

  ASSERT(sends_[0].numRequests() == 0);
  //if(timeout.active())
  //emitTime("final wait");
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
  if ( !parentScheduler_ && (execout.active() || timeout.active() || waitout.active()) ) {  // only do on toplevel scheduler
    outputTimingStats("MPIScheduler");
  }

  if (dbg.active()) {
    coutLock.lock();
    dbg << me << " MPIScheduler finished\n";
    coutLock.unlock();
  }
}

void
MPIScheduler::emitTime(const char* label)
{
  double time = Time::currentSeconds();
  emitTime(label, time-d_lasttime);
  d_lasttime=time;
}

void
MPIScheduler::emitTime(const char* label, double dt)
{
  d_labels.push_back(label);
  d_times.push_back(dt);
}
//______________________________________________________________________
//
void
MPIScheduler::outputTimingStats(const char* label)
{
  // add number of cells, patches, and particles
  int numCells = 0, numParticles = 0;
  OnDemandDataWarehouseP dw = d_dws[d_dws.size() - 1];
  const GridP grid(const_cast<Grid*>(dw->getGrid()));
  const PatchSubset* myPatches = getLoadBalancer()->getPerProcessorPatchSet(grid)->getSubset(d_myworld->myRank());
  for (int p = 0; p < myPatches->size(); p++) {
    const Patch* patch = myPatches->get(p);
    IntVector range = patch->getExtraCellHighIndex() - patch->getExtraCellLowIndex();
    numCells += range.x() * range.y() * range.z();

    // go through all materials since getting an MPMMaterial correctly would depend on MPM
    for (int m = 0; m < d_sharedState->getNumMatls(); m++) {
      if (dw->haveParticleSubset(m, patch))
        numParticles += dw->getParticleSubset(m, patch)->numParticles();
    }
  }

  emitTime("NumPatches", myPatches->size());
  emitTime("NumCells", numCells);
  emitTime("NumParticles", numParticles);
  std::vector<double> d_totaltimes(d_times.size());
  std::vector<double> d_maxtimes(d_times.size());
  std::vector<double> d_avgtimes(d_times.size());
  double avgTask = -1, maxTask = -1;
  double avgComm = -1, maxComm = -1;
  double avgCell = -1, maxCell = -1;

  MPI_Comm comm = d_myworld->getComm();
  MPI_Reduce(&d_times[0], &d_totaltimes[0], static_cast<int>(d_times.size()), MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(&d_times[0], &d_maxtimes[0],   static_cast<int>(d_times.size()), MPI_DOUBLE, MPI_MAX, 0, comm);

  double total = 0, avgTotal = 0, maxTotal = 0;
  for (int i = 0; i < (int)d_totaltimes.size(); i++) {
    d_avgtimes[i] = d_totaltimes[i] / d_myworld->nRanks();
    if (strcmp(d_labels[i], "Total task time") == 0) {
      avgTask = d_avgtimes[i];
      maxTask = d_maxtimes[i];
    }
    else if (strcmp(d_labels[i], "Total comm time") == 0) {
      avgComm = d_avgtimes[i];
      maxComm = d_maxtimes[i];
    }
    else if (strncmp(d_labels[i], "Num", 3) == 0) {
      if (strcmp(d_labels[i], "NumCells") == 0) {
        avgCell = d_avgtimes[i];
        maxCell = d_maxtimes[i];
      }
      // these are independent stats - not to be summed
      continue;
    }

    total    += d_times[i];
    avgTotal += d_avgtimes[i];
    maxTotal += d_maxtimes[i];
  }

  // to not duplicate the code
  std::vector<std::ofstream*> files;
  std::vector<std::vector<double>*> data;
  files.push_back(&timingStats);
  data.push_back(&d_times);

  int me = d_myworld->myRank();

  if (me == 0) {
    files.push_back(&avgStats);
    files.push_back(&maxStats);
    data.push_back(&d_avgtimes);
    data.push_back(&d_maxtimes);
  }

  for (unsigned file = 0; file < files.size(); file++) {
    std::ofstream& out = *files[file];
    out << "Timestep " << d_sharedState->getCurrentTopLevelTimeStep() << std::endl;
    for (int i = 0; i < (int)(*data[file]).size(); i++) {
      out << label << ": " << d_labels[i] << ": ";
      int len = static_cast<int>(strlen(d_labels[i]) + strlen("MPIScheduler: ") + strlen(": "));
      for (int j = len; j < 55; j++)
        out << ' ';
      double percent;
      if (strncmp(d_labels[i], "Num", 3) == 0) {
        percent = d_totaltimes[i] == 0 ? 100 : (*data[file])[i] / d_totaltimes[i] * 100;
      }
      else {
        percent = (*data[file])[i] / total * 100;
      }
      out << (*data[file])[i] << " (" << percent << "%)\n";
    }
    out << std::endl << std::endl;
  }

  if (me == 0) {
    timeout << "  Avg. exec: " << avgTask << ", max exec: " << maxTask << " = " << (1 - avgTask / maxTask) * 100 << " load imbalance (exec)%\n";
    timeout << "  Avg. comm: " << avgComm << ", max comm: " << maxComm << " = " << (1 - avgComm / maxComm) * 100 << " load imbalance (comm)%\n";
    timeout << "  Avg.  vol: " << avgCell << ", max  vol: " << maxCell << " = " << (1 - avgCell / maxCell) * 100 << " load imbalance (theoretical)%\n";
  }

  double time = Time::currentSeconds();
  d_lasttime = time;

  if (execout.active()) {
    static int count = 0;

    // only output the exec times every 10 timesteps
    if (++count % 10 == 0) {
      std::ofstream fout;
      char filename[100];
      sprintf(filename, "exectimes.%d.%d", d_myworld->nRanks(), d_myworld->myRank());
      fout.open(filename);

      // Report which timesteps TaskExecTime values have been accumulated over
      fout << "Reported values are cumulative over 10 timesteps ("
           << d_sharedState->getCurrentTopLevelTimeStep()-9
           << " through "
           << d_sharedState->getCurrentTopLevelTimeStep()
           << ")" << std::endl;

      for (auto iter = exectimes.begin(); iter != exectimes.end(); iter++) {
        fout << std::fixed << d_myworld->myRank() << ": TaskExecTime(s): " << iter->second << " Task:" << iter->first << std::endl;
      }
      fout.close();
      exectimes.clear();
    }
  }

  if (waitout.active()) {
    static int count = 0;

    // only output the exec times every 10 timesteps
    if (++count % 10 == 0) {

      if (d_myworld->myRank() == 0 || d_myworld->myRank() == d_myworld->nRanks() / 2 || d_myworld->myRank() == d_myworld->nRanks() - 1) {

        std::ofstream wout;
        char fname[100];
        sprintf(fname, "waittimes.%d.%d", d_myworld->nRanks(), d_myworld->myRank());
        wout.open(fname);

        for (auto iter = waittimes.begin(); iter != waittimes.end(); iter++) {
          wout << std::fixed << d_myworld->myRank() << ":   TaskWaitTime(TO): " << iter->second << " Task:" << iter->first << std::endl;
        }

        for (auto iter = DependencyBatch::waittimes.begin(); iter != DependencyBatch::waittimes.end();
             iter++) {
          wout << std::fixed << d_myworld->myRank() << ": TaskWaitTime(FROM): " << iter->second << " Task:" << iter->first << std::endl;
        }
        wout.close();
      }

      waittimes.clear();
      DependencyBatch::waittimes.clear();
    }
  }
}

