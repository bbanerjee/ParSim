/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#include <CCA/Components/Schedulers/CommunicationList.h>
#include <CCA/Components/Schedulers/RuntimeStats.h>
#include <CCA/Components/Schedulers/SendState.h>
#include <CCA/Components/Schedulers/TaskGraph.h>

#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Output.h>
#include <CCA/Ports/SimulationInterface.h>

#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/ParticleSubset.h>

#include <Core/Malloc/Allocator.h>

#include <Core/Parallel/ProcessorGroup.h>

#include <Core/Util/DOUT.hpp>
#include <Core/Util/FancyAssert.h>
#include <Core/Util/Timers/Timers.hpp>

#include <sci_defs/kokkos_defs.h>

#ifdef UINTAH_ENABLE_KOKKOS
#include <Kokkos_Core.hpp>
#endif // UINTAH_ENABLE_KOKKOS

#include <cstring>
#include <iomanip>
#include <map>
#include <sstream>

// Pack data into a buffer before sending -- testing to see if this
// works better and avoids certain problems possible when you allow
// tasks to modify data that may have a pending send.
#define USE_PACKING

namespace {

Uintah::MasterLock g_lb_mutex{};   // load balancer lock
Uintah::MasterLock g_recv_mutex{}; // for postMPIRecvs

Uintah::MasterLock g_msg_vol_mutex{}; // to report thread-safe msg volume info
Uintah::MasterLock
  g_send_time_mutex{}; // for reporting thread-safe MPI send times
Uintah::MasterLock
  g_recv_time_mutex{}; // for reporting thread-safe MPI recv times
Uintah::MasterLock
  g_wait_time_mutex{}; // for reporting thread-safe MPI wait times

Uintah::Dout g_dbg("MPIScheduler_DBG",
                   "MPIScheduler",
                   "general dbg info for MPIScheduler",
                   false);
Uintah::Dout g_send_stats("MPISendStats",
                          "MPIScheduler",
                          "MPI send statistics, num_sends, send volume",
                          false);
Uintah::Dout g_reductions("ReductionTasks",
                          "MPIScheduler",
                          "rank-0 reports each reduction task",
                          false);
Uintah::Dout g_time_out(
  "MPIScheduler_TimingsOut",
  "MPIScheduler",
  "write MPI timing files: timingstats.avg, timingstats.max",
  false);
Uintah::Dout g_task_level(
  "TaskLevel",
  "MPIScheduler",
  "output task name and each level's beginning patch when done",
  false);
} // namespace

namespace Uintah {

bool do_task_exec_stats = false;

// These are used externally, keep them visible outside this unit
Dout g_task_order("TaskOrder",
                  "MPIScheduler",
                  "task order debug stream",
                  false);
Dout g_task_dbg(
  "TaskDBG",
  "MPIScheduler",
  "output each task name as it begins/ends or when threaded, ready",
  false);
Dout g_task_run("TaskRun",
                "MPIScheduler",
                "output each task name as it runs",
                false);
Dout g_mpi_dbg("MPIDBG", "MPIScheduler", "MPI debug stream", false);
Dout g_exec_out("ExecOut", "MPIScheduler", "exec debug stream", false);

MPIScheduler::MPIScheduler(const ProcessorGroup* myworld,
                           MPIScheduler* inParentScheduler)
  : SchedulerCommon(myworld)
  , d_parent_scheduler(inParentScheduler)
{
#ifdef UINTAH_ENABLE_KOKKOS
  Kokkos::initialize();
#endif // UINTAH_ENABLE_KOKKOS

  if (g_time_out) {
    char filename[64];
    if (d_myworld->myRank() == 0) {
      sprintf(filename, "timingstats.avg");
      d_avg_stats.open(filename);
      sprintf(filename, "timingstats.max");
      d_max_stats.open(filename);
    }
  }

  std::string timeStr("seconds");

  d_mpi_info.insert(TotalSend, std::string("TotalSend"), timeStr);
  d_mpi_info.insert(TotalRecv, std::string("TotalRecv"), timeStr);
  d_mpi_info.insert(TotalTest, std::string("TotalTest"), timeStr);
  d_mpi_info.insert(TotalWait, std::string("TotalWait"), timeStr);
  d_mpi_info.insert(TotalReduce, std::string("TotalReduce"), timeStr);
  d_mpi_info.insert(TotalTask, std::string("TotalTask"), timeStr);

  d_task_info.setKeyName("Task");
  d_task_info.insert(ExecTime, std::string("ExecTime"), timeStr);
  d_task_info.insert(WaitTime, std::string("WaitTime"), timeStr);

  d_task_info.calculateSum(true);
  d_task_info.calculateAverage(true);
  d_task_info.calculateMinimum(true);
  d_task_info.calculateMaximum(true);
  d_task_info.calculateStdDev(true);
  m_num_schedulers += 1;
}

MPIScheduler::~MPIScheduler()
{
  if ((g_time_out) && (d_myworld->myRank() == 0)) {
    d_avg_stats.close();
    d_max_stats.close();
  }

#ifdef UINTAH_ENABLE_KOKKOS
  Kokkos::finalize();
#endif // UINTAH_ENABLE_KOKKOS
}

void
MPIScheduler::problemSetup(const ProblemSpecP& prob_spec,
                           const MaterialManagerP& mat_manager)
{
  SchedulerCommon::problemSetup(prob_spec, mat_manager);
}

SchedulerP
MPIScheduler::createSubScheduler()
{
  MPIScheduler* newsched = scinew MPIScheduler(d_myworld, this);

  newsched->setComponents(this);
  newsched->d_materialManager = d_materialManager;

  newsched->m_num_schedulers += 1;
  m_num_schedulers += 1;

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
  int checksum        = 0;
  int numSpatialTasks = 0;
  for (unsigned i = 0; i < d_task_graphs.size(); i++) {
    checksum += d_task_graphs[i]->getTasks().size();

    // This begins addressing the issue of making the global checksum more
    // sophisticated:
    //   check if any tasks were spatially scheduled - TaskType::Spatial,
    //   meaning no computes, requires or modifies
    //     e.g. RMCRT radiometer task, which is not scheduled on all patches
    //          these Spatial tasks won't count toward the global checksum
    std::vector<std::shared_ptr<Task>> tasks = d_task_graphs[i]->getTasks();
    for (auto& task : tasks) {
      if (task->getType() == Task::Spatial) {
        numSpatialTasks++;
      }
    }
  }

  // Spatial tasks don't count against the global checksum
  checksum -= numSpatialTasks;

  int my_rank = d_myworld->myRank();
  DOUTR(g_mpi_dbg,
        " (Uintah::MPI::Allreduce) Checking checksum of " << checksum);

  int result_checksum;
  Uintah::MPI::Allreduce(
    &checksum, &result_checksum, 1, MPI_INT, MPI_MIN, d_myworld->getComm());

  if (checksum != result_checksum) {
    std::cerr << "Failed task checksum comparison! Not all processes are "
                 "executing the same taskgraph \n";
    std::cerr << " Rank- " << my_rank << " of " << d_myworld->nRanks() - 1
              << ": has sum " << checksum << " and global is "
              << result_checksum << '\n';
    Uintah::MPI::Abort(d_myworld->getComm(), 1);
  }

  DOUTR(g_mpi_dbg, " (Uintah::MPI::Allreduce) Check succeeded");
#endif
}

void
MPIScheduler::initiateTask(DetailedTask* task,
                           bool only_old_recvs,
                           int abort_point,
                           int iteration)
{
  if (only_old_recvs) {
    return;
  }

  postMPIRecvs(task, only_old_recvs, abort_point, iteration);
}

void
MPIScheduler::initiateReduction(DetailedTask* dtask)
{
  bool ans = (g_reductions || g_task_run);
  DOUTR(ans, " Running Reduction Task: " << dtask->getName());

  Timers::Simple timer;

  timer.start();

  runReductionTask(dtask);

  timer.stop();

  d_mpi_info[TotalReduce] += timer().seconds();

  if (g_exec_out || do_task_exec_stats) {
    d_task_info[dtask->getTask()->getName()][TaskStatsEnum::ExecTime] +=
      timer().seconds();
    d_task_info[dtask->getTask()->getName()][TaskStatsEnum::WaitTime] = 0.0;
  }
}

void
MPIScheduler::runTask(DetailedTask* dtask, int iteration)
{
  if (d_tracking_vars_print_location & SchedulerCommon::PRINT_BEFORE_EXEC) {
    printTrackedVars(dtask, SchedulerCommon::PRINT_BEFORE_EXEC);
  }

  std::vector<DataWarehouseP> plain_old_dws(d_dws.size());
  for (size_t i = 0; i < d_dws.size(); i++) {
    plain_old_dws[i] = d_dws[i].get_rep();
  }

  DOUTR(g_task_run, " Running task:   " << *dtask);

  dtask->doit(d_myworld, d_dws, plain_old_dws);

  if (d_tracking_vars_print_location & SchedulerCommon::PRINT_AFTER_EXEC) {
    printTrackedVars(dtask, SchedulerCommon::PRINT_AFTER_EXEC);
  }

  postMPISends(dtask, iteration);

  dtask->done(d_dws);

  g_lb_mutex.lock();
  {
    // Do the global and local per task monitoring
    sumTaskMonitoringValues(dtask);

    double total_task_time = dtask->task_exec_time();
    if (g_exec_out || do_task_exec_stats) {
      d_task_info[dtask->getTask()->getName()][TaskStatsEnum::ExecTime] +=
        total_task_time;
      d_task_info[dtask->getTask()->getName()][TaskStatsEnum::WaitTime] +=
        dtask->task_wait_time();
    }

    // if I do not have a sub scheduler
    if (!dtask->getTask()->getHasSubScheduler()) {
      // add my task time to the total time
      d_mpi_info[TotalTask] += total_task_time;
      if (!d_is_copy_data_timestep &&
          dtask->getTask()->getType() != Task::Output) {
        // add contribution for patchlist
        d_load_balancer->addContribution(dtask, total_task_time);
      }
    }
  }
  g_lb_mutex.unlock();

  //---------------------------------------------------------------------------
  // New way of managing single MPI requests - avoids Uintah::MPI::Waitsome &
  // MPI_Donesome - APH 07/20/16
  // ---------------------------------------------------------------------------
  // test a pending request
  auto ready_request = [](CommRequest const& r) -> bool { return r.test(); };
  CommRequestPool::iterator comm_sends_iter = d_sends.find_any(ready_request);
  if (comm_sends_iter) {
    MPI_Status status;
    comm_sends_iter->finishedCommunication(d_myworld, status);
    d_sends.erase(comm_sends_iter);
  }
  //-----------------------------------

  // Add subscheduler timings to the parent scheduler and reset subscheduler
  // timings
  if (d_parent_scheduler) {
    size_t num_elems = d_mpi_info.size();
    for (size_t i = 0; i < num_elems; ++i) {
      d_parent_scheduler->d_mpi_info[i] += d_mpi_info[i];
    }
    d_mpi_info.reset(0);
  }
} // end runTask()

void
MPIScheduler::runReductionTask(DetailedTask* task)
{
  const Task::Dependency* mod = task->getTask()->getModifies();
  ASSERT(!mod->next);

  OnDemandDataWarehouse* dw = d_dws[mod->mapDataWarehouse()].get_rep();
  ASSERT(task->getTask()->m_comm >= 0);
  dw->reduceMPI(
    mod->var, mod->reduction_level, mod->matls, task->getTask()->m_comm);
  task->done(d_dws);
}

void
MPIScheduler::postMPISends(DetailedTask* task, int iteration)
{
  Timers::Simple send_timer;
  send_timer.start();

  MPI_Comm my_comm = d_myworld->getComm();

  DOUTR(g_dbg, " postMPISends - task " << *task);

  // Send data to dependendents
  for (DependencyBatch* batch = task->getComputes(); batch != nullptr;
       batch                  = batch->comp_next) {

    // Prepare to send a message
#ifdef USE_PACKING
    PackBufferInfo mpibuff;
#else
    BufferInfo mpibuff;
#endif

    // Create the MPI type
    int to = batch->to_tasks.front()->getAssignedResourceIndex();
    ASSERTRANGE(to, 0, d_myworld->nRanks());

    for (DetailedDep* req = batch->head; req != nullptr; req = req->m_next) {

      if ((req->m_comm_condition == DetailedDep::FirstIteration &&
           iteration > 0) ||
          (req->m_comm_condition == DetailedDep::SubsequentIterations &&
           iteration == 0) ||
          (d_not_copy_data_vars.count(req->m_req->var->getName()) > 0)) {

        // See comment in DetailedDep about CommCondition
        DOUTR(g_dbg, "   Ignoring conditional send for " << *req);
        continue;
      }

      // if we send/recv to an output task, don't send/recv if not an output
      // timestep

      // ARS NOTE: Outputing and Checkpointing may be done out of snyc
      // now, i.e. turned on just before it happens rather than turned
      // on before the task graph execution.  As such, one should also
      // be checking:

      // d_simulator->activeReductionVariable( "outputInterval" );
      // d_simulator->activeReductionVariable( "checkpointInterval" );

      // However, if active the code below would be called regardless
      // if an output or checkpoint time step or not. Not sure that is
      // desired but not sure of the effect of not calling it and doing
      // an out of sync output or checkpoint.
      if (req->m_to_tasks.front()->getTask()->getType() == Task::Output &&
          !d_output->isOutputTimestep() && !d_output->isCheckpointTimestep()) {
        DOUTR(g_dbg, "   Ignoring non-output-timestep send for " << *req);
        continue;
      }

      OnDemandDataWarehouse* dw = d_dws[req->m_req->mapDataWarehouse()].get_rep();

      DOUTR(g_dbg, " --> sending " << *req);
      DOUTR(g_dbg, "     From task: " << batch->from_task->getName());
      DOUTR(g_dbg,
            "     To task:   " << req->m_to_tasks.front()->getTask()->getName()
                               << " and  rank " << to);
      DOUTR(g_dbg,
            "     ghost type: "
              << Ghost::getGhostTypeName(req->m_req->gtype)
              << ", num req ghost: " << req->m_req->num_ghost_cells
              << ", Ghost::direction: "
              << Ghost::getGhostTypeDir(req->m_req->gtype) << ", from dw "
              << dw->getID());

      // the load balancer is used to determine where data was in the old dw on
      // the prev timestep - pass it in if the particle data is on the old dw
      const VarLabel* posLabel;
      OnDemandDataWarehouse* posDW;

      if (!d_reloc_new_pos_label && d_parent_scheduler) {
        posDW =
          d_dws[req->m_req->task->mapDataWarehouse(Task::ParentOldDW)].get_rep();
        posLabel = d_parent_scheduler->d_reloc_new_pos_label;
      } else {
        // on an output task (and only on one) we require particle variables
        // from the NewDW
        if (req->m_to_tasks.front()->getTask()->getType() == Task::Output) {
          posDW = d_dws[req->m_req->task->mapDataWarehouse(Task::NewDW)].get_rep();
        } else {
          posDW = d_dws[req->m_req->task->mapDataWarehouse(Task::OldDW)].get_rep();
        }
        posLabel = d_reloc_new_pos_label;
      }

      MPIScheduler* top = this;
      while (top->d_parent_scheduler) {
        top = top->d_parent_scheduler;
      }

      dw->sendMPI(batch, posLabel, mpibuff, posDW, req, d_load_balancer);
    }

    // Post the send
    if (mpibuff.count() > 0) {
      ASSERT(batch->message_tag > 0);
      void* buf{ nullptr };
      int count;
      MPI_Datatype datatype;

#ifdef USE_PACKING
      mpibuff.get_type(buf, count, datatype, d_myworld->getComm());
      mpibuff.pack(d_myworld->getComm(), count);
#else
      mpibuff.get_type(buf, count, datatype);
#endif

      if (!buf) {
        printf("postMPISends() - ERROR, the send MPI buffer is nullptr\n");
        SCI_THROW(
          InternalError("The send MPI buffer is null", __FILE__, __LINE__));
      }
      DOUTR(g_mpi_dbg,
            " Posting send for message number "
              << batch->message_tag << " to   rank-" << to
              << ", length: " << count << " (bytes)");

      d_num_messages++;
      int typeSize;

      Uintah::MPI::Type_size(datatype, &typeSize);

      {
        std::lock_guard<Uintah::MasterLock> msg_vol_lock(g_msg_vol_mutex);
        d_message_volume += count * typeSize;
      }

      //---------------------------------------------------------------------------
      // New way of managing single MPI requests - avoids Uintah::MPI::Waitsome
      // & MPI_Donesome - APH 07/20/16
      //---------------------------------------------------------------------------
      CommRequestPool::iterator comm_sends_iter =
        d_sends.emplace(new SendHandle(mpibuff.takeSendlist()));
      Uintah::MPI::Isend(buf,
                         count,
                         datatype,
                         to,
                         batch->message_tag,
                         my_comm,
                         comm_sends_iter->request());
      comm_sends_iter.clear();
      //---------------------------------------------------------------------------
    }
  } // end for (DependencyBatch * batch = task->getComputes() )

  send_timer.stop();
  {
    std::lock_guard<Uintah::MasterLock> send_time_lock(g_send_time_mutex);
    d_mpi_info[TotalSend] += send_timer().seconds();
  }
} // end postMPISends();

struct CompareDep
{
  bool
  operator()(DependencyBatch* a, DependencyBatch* b)
  {
    return a->message_tag < b->message_tag;
  }
};

void
MPIScheduler::postMPIRecvs(DetailedTask* task,
                           bool only_old_recvs,
                           int abort_point,
                           int iteration)
{
  Timers::Simple recv_timer;
  recv_timer.start();

  MPI_Comm my_comm = d_myworld->getComm();

  DOUTR(g_dbg, " postMPIRecvs - task " << *task);

  if (d_tracking_vars_print_location & SchedulerCommon::PRINT_BEFORE_COMM) {
    printTrackedVars(task, SchedulerCommon::PRINT_BEFORE_COMM);
  }

  // sort the requires, so in case there is a particle send we receive it with
  // the right message tag
  std::vector<DependencyBatch*> sorted_reqs;
  for (auto iter = task->getRequires().cbegin();
       iter != task->getRequires().cend();
       iter++) {
    sorted_reqs.push_back(iter->first);
  }

  CompareDep comparator;
  std::sort(sorted_reqs.begin(), sorted_reqs.end(), comparator);

  // Need this until race condition on foreign variables is resolved - APH,
  // 09/19/17
  std::lock_guard<Uintah::MasterLock> recv_lock(g_recv_mutex);
  {
    for (auto sorted_iter = sorted_reqs.cbegin();
         sorted_iter != sorted_reqs.cend();
         sorted_iter++) {
      DependencyBatch* batch = *sorted_iter;

      task->incrementExternalDepCount();

      // The first thread that calls this on the batch will return true
      // while subsequent threads calling this will block and wait for
      // that first thread to receive the data.
      if (!batch->makeMPIRequest()) {
        DOUTR(g_dbg, "Someone else already receiving it");
        continue;
      }

      if (only_old_recvs) {
        DOUTR(g_dbg,
              " abort analysis: "
                << batch->from_task->getTask()->getName()
                << ", so=" << batch->from_task->getTask()->getSortedOrder()
                << ", abort_point=" << abort_point);

        if (batch->from_task->getTask()->getSortedOrder() <= abort_point) {
          DOUTR(g_dbg,
                "posting MPI recv for pre-abort message "
                  << batch->message_tag);
        }
        if (!(batch->from_task->getTask()->getSortedOrder() <= abort_point)) {
          continue;
        }
      }

      // Prepare to receive a message
      BatchReceiveHandler* pBatchRecvHandler =
        scinew BatchReceiveHandler(batch);
      PackBufferInfo* p_mpibuff = nullptr;

#ifdef USE_PACKING
      p_mpibuff               = scinew PackBufferInfo();
      PackBufferInfo& mpibuff = *p_mpibuff;
#else
      BufferInfo mpibuff;
#endif

      // Create the MPI type
      for (DetailedDep* req = batch->head; req != 0; req = req->m_next) {

        OnDemandDataWarehouse* dw = d_dws[req->m_req->mapDataWarehouse()].get_rep();
        if ((req->m_comm_condition == DetailedDep::FirstIteration &&
             iteration > 0) ||
            (req->m_comm_condition == DetailedDep::SubsequentIterations &&
             iteration == 0) ||
            (d_not_copy_data_vars.count(req->m_req->var->getName()) > 0)) {

          // See comment in DetailedDep about CommCondition
          DOUTR(g_dbg, "   Ignoring conditional receive for " << *req);
          continue;
        }

        // if we send/recv to an output task, don't send/recv if not an output
        // timestep

        // ARS NOTE: Outputing and Checkpointing may be done out of
        // snyc now. I.e. turned on just before it happens rather than
        // turned on before the task graph execution.  As such, one
        // should also be checking:

        // d_simulator->activeReductionVariable( "outputInterval" );
        // d_simulator->activeReductionVariable( "checkpointInterval" );

        // However, if active the code below would be called regardless
        // if an output or checkpoint time step or not. Not sure that is
        // desired but not sure of the effect of not calling it and doing
        // an out of sync output or checkpoint.
        if (req->m_to_tasks.front()->getTask()->getType() == Task::Output &&
            !d_output->isOutputTimestep() &&
            !d_output->isCheckpointTimestep()) {
          DOUTR(g_dbg, "   Ignoring non-output-timestep receive for " << *req);
          continue;
        }

        DOUTR(g_dbg, " <-- receiving " << *req);
        DOUTR(g_dbg, "     From task: " << batch->from_task->getName());
        DOUTR(g_dbg,
              "     ghost type: " << Ghost::getGhostTypeName(req->m_req->gtype)
                                  << ", num req ghost "
                                  << req->m_req->num_ghost_cells
                                  << ", Ghost::direction: "
                                  << Ghost::getGhostTypeDir(req->m_req->gtype)
                                  << ", into dw " << dw->getID());

        OnDemandDataWarehouse* posDW;

        // the load balancer is used to determine where data was in the old dw
        // on the prev timestep pass it in if the particle data is on the old dw
        if (!d_reloc_new_pos_label && d_parent_scheduler) {
          posDW =
            d_dws[req->m_req->task->mapDataWarehouse(Task::ParentOldDW)].get_rep();
        } else {
          // on an output task (and only on one) we require particle variables
          // from the NewDW
          if (req->m_to_tasks.front()->getTask()->getType() == Task::Output) {
            posDW =
              d_dws[req->m_req->task->mapDataWarehouse(Task::NewDW)].get_rep();
          } else {
            posDW =
              d_dws[req->m_req->task->mapDataWarehouse(Task::OldDW)].get_rep();
          }
        }

        MPIScheduler* top = this;
        while (top->d_parent_scheduler) {
          top = top->d_parent_scheduler;
        }

        dw->recvMPI(batch, mpibuff, posDW, req, d_load_balancer);

        if (!req->isNonDataDependency()) {
          d_task_graphs[d_current_task_graph]
            ->getDetailedTasks()
            ->setScrubCount(req->m_req, req->m_matl, req->m_from_patch, d_dws);
        }
      }

      // Post the receive
      if (mpibuff.count() > 0) {

        ASSERT(batch->message_tag > 0);
        void* buf{ nullptr };
        int count;
        MPI_Datatype datatype;

#ifdef USE_PACKING
        mpibuff.get_type(buf, count, datatype, d_myworld->getComm());
#else
        mpibuff.get_type(buf, count, datatype);
#endif

        if (!buf) {
          printf("postMPIRecvs() - ERROR, the receive MPI buffer is nullptr\n");
          SCI_THROW(InternalError(
            "The receive MPI buffer is nullptr", __FILE__, __LINE__));
        }

        int from = batch->from_task->getAssignedResourceIndex();
        ASSERTRANGE(from, 0, d_myworld->nRanks());

        DOUTR(g_mpi_dbg,
              " Posting recv for message number "
                << batch->message_tag << " from rank-" << from
                << ", length: " << count << " (bytes)");

        //---------------------------------------------------------------------------
        // New way of managing single MPI requests - avoids
        // Uintah::MPI::Waitsome & MPI_Donesome - APH 07/20/16
        //---------------------------------------------------------------------------
        CommRequestPool::iterator comd_d_recvsiter =
          d_recvs.emplace(new RecvHandle(p_mpibuff, pBatchRecvHandler));
        Uintah::MPI::Irecv(buf,
                           count,
                           datatype,
                           from,
                           batch->message_tag,
                           my_comm,
                           comd_d_recvsiter->request());
        comd_d_recvsiter.clear();
        //---------------------------------------------------------------------------

      } else {
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

    recv_timer.stop();
  }

  {
    std::lock_guard<Uintah::MasterLock> recv_time_lock(g_recv_time_mutex);
    d_mpi_info[TotalRecv] += recv_timer().seconds();
  }

} // end postMPIRecvs()

void
MPIScheduler::processMPIRecvs(int how_much)
{
  if (d_recvs.size() == 0u) {
    return;
  }

  Timers::Simple process_recv_timer;
  process_recv_timer.start();

  //---------------------------------------------------------------------------
  // New way of managing single MPI requests - avoids Uintah::MPI::Waitsome &
  // MPI_Donesome - APH 07/20/16
  //---------------------------------------------------------------------------
  auto test_request = [](CommRequest const& n) -> bool { return n.test(); };
  auto wait_request = [](CommRequest const& n) -> bool { return n.wait(); };

  CommRequestPool::iterator comm_iter;

  switch (how_much) {
    case TEST: {
      RuntimeStats::TestTimer mpi_test_timer;
      comm_iter = d_recvs.find_any(test_request);
      if (comm_iter) {
        MPI_Status status;
        comm_iter->finishedCommunication(d_myworld, status);
        d_recvs.erase(comm_iter);
      }
      break;
    }

    case WAIT_ONCE: {
      RuntimeStats::WaitTimer mpi_wait_timer;
      comm_iter = d_recvs.find_any(wait_request);
      if (comm_iter) {
        MPI_Status status;
        comm_iter->finishedCommunication(d_myworld, status);
        d_recvs.erase(comm_iter);
      }
      break;
    }

    case WAIT_ALL: {
      RuntimeStats::WaitTimer mpi_wait_timer;
      while (d_recvs.size() != 0u) {
        comm_iter = d_recvs.find_any(wait_request);
        if (comm_iter) {
          MPI_Status status;
          comm_iter->finishedCommunication(d_myworld, status);
          d_recvs.erase(comm_iter);
        }
      }
      break;
    }
  } // end switch

  process_recv_timer.stop();

  {
    std::lock_guard<Uintah::MasterLock> wait_time_lock(g_wait_time_mutex);
    d_mpi_info[TotalWait] += process_recv_timer().seconds();
  }

} // end processMPIRecvs()

void
MPIScheduler::execute(int tgnum /*=0*/, int iteration /*=0*/)
{
  // Track total scheduler execution time across timesteps.
  d_exec_timer.reset(true);

  // If doing in situ monitoring clear the times before each time step
  // otherwise the times are accumulated over N time steps.
  if (do_task_exec_stats) {
    d_task_info.reset(0);
  }

  // create the various timers
  RuntimeStats::initialize_timestep(m_num_schedulers, d_task_graphs);

  ASSERTRANGE(tgnum, 0, (int)d_task_graphs.size());
  TaskGraph* tg = d_task_graphs[tgnum].get();
  tg->setIteration(iteration);
  d_current_task_graph = tgnum;

  if (d_task_graphs.size() > 1) {
    // tg model is the multi TG model, where each graph is going to need to
    // have its dwmap reset here (even with the same tgnum)
    tg->remapTaskDWs(d_dw_map);
  }

  DetailedTasks* dts = tg->getDetailedTasks();

  if (dts == nullptr) {
    proc0cout << "MPIScheduler skipping execute, no tasks" << std::endl;
    return;
  }

  int ntasks = dts->numLocalTasks();

  if (d_runtime_stats) {
    (*d_runtime_stats)[NumTasks] += ntasks;
  }

  dts->initializeScrubs(d_dws, d_dw_map);
  dts->initTimestep();

  for (int i = 0; i < ntasks; i++) {
    dts->localTask(i)->resetDependencyCounts();
  }

  int my_rank = d_myworld->myRank();

  makeTaskGraphDoc(dts, my_rank);

  d_mpi_info.reset(0);

  DOUTR(g_dbg,
        ", MPI Scheduler executing taskgraph: "
          << tgnum << ", timestep: " << d_simulator->getTimestep() << " with "
          << dts->numTasks() << " tasks (" << ntasks << " local)");

  if (d_reloc_new_pos_label && d_dws[d_dw_map[Task::OldDW]] != nullptr) {
    d_dws[d_dw_map[Task::OldDW]]->exchangeParticleQuantities(
      dts, d_load_balancer, d_reloc_new_pos_label, iteration);
  }

  bool abort       = false;
  int abort_point  = 987654;
  int numTasksDone = 0;
  int i            = 0;

  while (numTasksDone < ntasks) {
    i++;

    DetailedTask* task = dts->getNextInternalReadyTask();

    numTasksDone++;

    if (g_task_order && d_myworld->myRank() == d_myworld->nRanks() / 2) {
      std::ostringstream task_name;
      task_name << "  Running task: \"" << task->getTask()->getName() << "\" ";

      std::ostringstream task_type;
      task_type << "(" << task->getTask()->getType() << ") ";

      // task ordering debug info - please keep this here, APH 05/30/18
      DOUTR(true,
            std::setw(60) << std::left << task_name.str() << std::setw(14)
                          << std::left << task_type.str() << std::setw(15)
                          << " static order: " << std::setw(3) << std::left
                          << task->getStaticOrder() << std::setw(18)
                          << " scheduled order: " << std::setw(3) << std::left
                          << numTasksDone);
    }

    if (task->getTask()->getType() == Task::Reduction) {
      if (!abort) {
        initiateReduction(task);
        DOUTR(g_task_dbg, " Completed task:   " << *task)
      }
    } else {
      initiateTask(task, abort, abort_point, iteration);
      task->markInitiated();
      processMPIRecvs(WAIT_ALL);
      ASSERT(d_recvs.size() == 0u);
      runTask(task, iteration);

      DOUTR(g_task_dbg, " Completed task:   " << *task);
      printTaskLevels(d_myworld, g_task_level, task);
    }

    // ARS - FIXME CHECK THE WAREHOUSE
    OnDemandDataWarehouse* dw = d_dws[d_dws.size() - 1].get_rep();
    if (!abort && dw && dw->abortTimestep()) {
      // TODO - abort might not work with external queue...
      abort       = true;
      abort_point = task->getTask()->getSortedOrder();

      DOUTR(true,
            "  WARNING: Aborting time step after task: "
              << task->getTask()->getName());
    }
  } // end while( numTasksDone < ntasks )

  //---------------------------------------------------------------------------
  // New way of managing single MPI requests - avoids Uintah::MPI::Waitsome &
  // MPI_Donesome - APH 07/20/16
  // ---------------------------------------------------------------------------
  // wait on all pending requests
  auto ready_request = [](CommRequest const& r) -> bool { return r.wait(); };
  CommRequestPool::handle find_handle;
  while (d_sends.size() != 0u) {
    CommRequestPool::iterator comm_sends_iter;
    if ((comm_sends_iter = d_sends.find_any(find_handle, ready_request))) {
      find_handle = comm_sends_iter;
      d_sends.erase(comm_sends_iter);
    } else {
      // TODO - make this a sleep? APH 07/20/16
    }
  }
  //---------------------------------------------------------------------------

  ASSERT(d_sends.size() == 0u);
  ASSERT(d_recvs.size() == 0u);

  finalizeTimestep();

  d_exec_timer.stop();

  // compute the net timings
  computeNetRuntimeStats();

  // only do on top-level scheduler
  if (d_parent_scheduler == nullptr) {

    // This seems like the best place to collect and save these runtime stats.
    // They are reported in outputTimingStats.
    if (d_runtime_stats) {
      int numCells = 0, numParticles = 0;
      OnDemandDataWarehouse* dw = d_dws[d_dws.size() - 1].get_rep();
      const GridP grid(const_cast<Grid*>(dw->getGrid()));
      const PatchSubset* myPatches =
        d_load_balancer->getPerProcessorPatchSet(grid)->getSubset(my_rank);

      for (auto p = 0; p < myPatches->size(); p++) {
        const Patch* patch = myPatches->get(p);
        IntVector range =
          patch->getExtraCellHighIndex() - patch->getExtraCellLowIndex();
        numCells += range.x() * range.y() * range.z();

        // Go through all materials since getting an MPMMaterial
        // correctly would depend on MPM
        for (unsigned int m = 0; m < d_materialManager->getNumMaterials();
             m++) {
          if (dw->haveParticleSubset(m, patch)) {
            numParticles += dw->getParticleSubset(m, patch)->numParticles();
          }
        }
      }

      (*d_runtime_stats)[NumPatches]   = myPatches->size();
      (*d_runtime_stats)[NumCells]     = numCells;
      (*d_runtime_stats)[NumParticles] = numParticles;
    }

    outputTimingStats("MPIScheduler");
  }

  RuntimeStats::report(d_myworld->getComm());
}

//  Take the various timers and compute the net results
void
MPIScheduler::computeNetRuntimeStats()
{
  if (d_runtime_stats) {
    // don't count output time
    (*d_runtime_stats)[TaskExecTime] +=
      d_mpi_info[TotalTask] - (*d_runtime_stats)[TotalIOTime];
    (*d_runtime_stats)[TaskLocalCommTime] +=
      d_mpi_info[TotalRecv] + d_mpi_info[TotalSend];
    (*d_runtime_stats)[TaskWaitCommTime] += d_mpi_info[TotalWait];
    (*d_runtime_stats)[TaskReduceCommTime] += d_mpi_info[TotalReduce];
  }
}

void
MPIScheduler::emitTime(const char* label, double dt)
{
  d_labels.push_back(label);
  d_times.push_back(dt);
}

void
MPIScheduler::outputTimingStats(const char* label)
{
  int my_rank      = d_myworld->myRank();
  int my_comm_size = d_myworld->nRanks();
  MPI_Comm my_comm = d_myworld->getComm();

  // for ExecTimes
  if (g_exec_out) {
    static int accumulate = 10;
    static int count      = 0;

    ++count;

    // Only output the exec times every N timesteps.
    if (d_simulator->getTimestep() % accumulate == 0) {

      // Report which timesteps the values have been accumulated
      // over. If doing in situ monitoring the values will be reset to
      // zero in the "execute" method so there is no accumulation. As
      // such the times will be for current time step ONLY otherwise
      // the times are accumulated over N time steps.
      std::ostringstream preamble;

      if (do_task_exec_stats) {
        preamble << "Reported values are for timestep : "
                 << d_simulator->getTimestep() << " ONLY";
      } else {
        preamble << "# Reported values are cumulative over " << count
                 << " timesteps (" << d_simulator->getTimestep() - (count - 1)
                 << " through " << d_simulator->getTimestep() << ")\n"
                 << "# Tasks run inside a subscheduler are not included";
      }

      // Report the stats for each task. Writing over any previous
      // files.
      d_task_info.reportIndividualStats("TaskStats",
                                        preamble.str(),
                                        my_rank,
                                        my_comm_size,
                                        d_simulator->getTimestep(),
                                        d_simulator->getSimTime(),
                                        BaseInfoMapper::Write_Last);

      // Not clear if writing reductions is useful so commented out for now.

      // m_task_info.reduce( false );

      // m_task_info.reportSummaryStats( "TaskStatsSummary", preamble.str(),
      //                                 my_rank, my_comm_size,
      //                                 d_simulator->getTimestep(),
      //                                 d_simulator->getSimTime(),
      //                                 BaseInfoMapper::Write_Last, false );

      count = 0;

      // If doing in situ monitoring do not reset the values to zero
      // as they will be reset in the "execute" method.
      if (!do_task_exec_stats) {
        d_task_info.reset(0);
      }
    }
  }

  // for file-based MPI timings
  if (g_time_out) {

    d_labels.clear();
    d_times.clear();

    double totalexec = d_exec_timer().seconds();

    if (d_runtime_stats) {
      emitTime("NumPatches", (*d_runtime_stats)[NumPatches]);
      emitTime("NumCells", (*d_runtime_stats)[NumCells]);
      emitTime("NumParticles", (*d_runtime_stats)[NumParticles]);
    }

    emitTime("Total send time", d_mpi_info[TotalSend]);
    emitTime("Total recv time", d_mpi_info[TotalRecv]);
    emitTime("Total test time", d_mpi_info[TotalTest]);
    emitTime("Total wait time", d_mpi_info[TotalWait]);
    emitTime("Total reduce time", d_mpi_info[TotalReduce]);
    emitTime("Total task time", d_mpi_info[TotalTask]);
    emitTime("Total comm time",
             d_mpi_info[TotalSend] + d_mpi_info[TotalRecv] +
               d_mpi_info[TotalTest] + d_mpi_info[TotalWait] +
               d_mpi_info[TotalReduce]);

    emitTime("Total execution time", totalexec);
    emitTime("Non-comm execution time",
             totalexec - d_mpi_info[TotalSend] - d_mpi_info[TotalRecv] -
               d_mpi_info[TotalTest] - d_mpi_info[TotalWait] -
               d_mpi_info[TotalReduce]);

    std::vector<double> totaltimes(d_times.size());
    std::vector<double> maxtimes(d_times.size());
    std::vector<double> avgtimes(d_times.size());
    double avgTask = -1, maxTask = -1;
    double avgComm = -1, maxComm = -1;
    double avgCell = -1, maxCell = -1;

    MPI_Comm comm = d_myworld->getComm();
    Uintah::MPI::Reduce(&d_times[0],
                        &totaltimes[0],
                        static_cast<int>(d_times.size()),
                        MPI_DOUBLE,
                        MPI_SUM,
                        0,
                        comm);
    Uintah::MPI::Reduce(&d_times[0],
                        &maxtimes[0],
                        static_cast<int>(d_times.size()),
                        MPI_DOUBLE,
                        MPI_MAX,
                        0,
                        comm);

    double total    = 0;
    double avgTotal = 0;
    double maxTotal = 0;
    for (size_t i = 0; i < totaltimes.size(); i++) {
      avgtimes[i] = totaltimes[i] / my_comm_size;
      if (strcmp(d_labels[i], "Total task time") == 0) {
        avgTask = avgtimes[i];
        maxTask = maxtimes[i];
      } else if (strcmp(d_labels[i], "Total comm time") == 0) {
        avgComm = avgtimes[i];
        maxComm = maxtimes[i];
      } else if (strncmp(d_labels[i], "Num", 3) == 0) {
        if (strcmp(d_labels[i], "NumCells") == 0) {
          avgCell = avgtimes[i];
          maxCell = maxtimes[i];
        }
        // these are independent stats - not to be summed
        continue;
      }

      total += d_times[i];
      avgTotal += avgtimes[i];
      maxTotal += maxtimes[i];
    }

    // to not duplicate the code
    std::vector<std::ofstream*> files;
    std::vector<std::vector<double>*> data;
    data.push_back(&d_times);

    if (my_rank == 0) {
      files.push_back(&d_avg_stats);
      files.push_back(&d_max_stats);
      data.push_back(&avgtimes);
      data.push_back(&maxtimes);
    }

    for (size_t file = 0; file < files.size(); ++file) {
      std::ofstream& out = *files[file];
      out << "Timestep " << d_simulator->getTimestep() << std::endl;
      for (size_t i = 0; i < (*data[file]).size(); i++) {
        out << label << ": " << d_labels[i] << ": ";
        int len = static_cast<int>(strlen(d_labels[i]) +
                                   strlen("MPIScheduler: ") + strlen(": "));
        for (int j = len; j < 55; j++) {
          out << ' ';
        }
        double percent;
        if (strncmp(d_labels[i], "Num", 3) == 0) {
          percent =
            totaltimes[i] == 0 ? 100 : (*data[file])[i] / totaltimes[i] * 100;
        } else {
          percent = (*data[file])[i] / total * 100;
        }
        out << (*data[file])[i] << " (" << percent << "%)\n";
      }
      out << std::endl << std::endl;
    }

    if (my_rank == 0) {
      std::ostringstream message;
      message << "\n";
      message << "  avg exec: " << std::setw(12) << avgTask
              << ",   max exec: " << std::setw(12) << maxTask
              << "    load imbalance (exec)%:        " << std::setw(6)
              << (1 - avgTask / maxTask) * 100 << "\n";
      message << "  avg comm: " << std::setw(12) << avgComm
              << ",   max comm: " << std::setw(12) << maxComm
              << "    load imbalance (comm)%:        " << std::setw(6)
              << (1 - avgComm / maxComm) * 100 << "\n";
      message << "  avg  vol: " << std::setw(12) << avgCell
              << ",   max  vol: " << std::setw(12) << maxCell
              << "    load imbalance (theoretical)%: " << std::setw(6)
              << (1 - avgCell / maxCell) * 100 << "\n";
      DOUT(g_time_out, message.str());
    }
  } // end g_time_out

  // for MPISendStats
  if (g_send_stats) {
    unsigned int total_messages;
    unsigned int max_messages;
    double total_volume;
    double max_volume;

    // do SUM and MAX reduction for m_num_messages and m_message_volume
    Uintah::MPI::Reduce(
      &d_num_messages, &total_messages, 1, MPI_UNSIGNED, MPI_SUM, 0, my_comm);
    Uintah::MPI::Reduce(
      &d_message_volume, &total_volume, 1, MPI_DOUBLE, MPI_SUM, 0, my_comm);
    Uintah::MPI::Reduce(
      &d_num_messages, &max_messages, 1, MPI_UNSIGNED, MPI_MAX, 0, my_comm);
    Uintah::MPI::Reduce(
      &d_message_volume, &max_volume, 1, MPI_DOUBLE, MPI_MAX, 0, my_comm);

    if (my_rank == 0) {
      std::ostringstream message;
      message << "MPISendStats: Num Send Messages   (avg): " << std::setw(12)
              << total_messages / (static_cast<double>(my_comm_size))
              << "    (max):" << std::setw(12) << max_messages << "\n";
      message << "MPISendStats: Send Message Volume (avg): " << std::setw(12)
              << total_volume / (static_cast<double>(my_comm_size))
              << "    (max):" << std::setw(12) << max_volume << "\n";
      DOUT(g_send_stats, message.str());
    }
  }
}

} // namespace Uintah