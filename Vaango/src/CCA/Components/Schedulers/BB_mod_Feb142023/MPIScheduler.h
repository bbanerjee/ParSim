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

#ifndef VAANGO_CCA_COMPONENTS_SCHEDULERS_MPISCHEDULER_H
#define VAANGO_CCA_COMPONENTS_SCHEDULERS_MPISCHEDULER_H

#include <CCA/Components/Schedulers/CommunicationList.h>
#include <CCA/Components/Schedulers/DetailedTasks.h>
#include <CCA/Components/Schedulers/OnDemandDataWarehouseP.h>
#include <CCA/Components/Schedulers/SchedulerCommon.h>

#include <Core/Grid/MaterialManagerP.h>

#include <Core/Util/InfoMapper.h>
#include <Core/Util/Timers/Timers.hpp>

#include <fstream>
#include <map>
#include <vector>

namespace Uintah {

class MPIScheduler : public SchedulerCommon
{
public:
  MPIScheduler(const ProcessorGroup* myworld,
               MPIScheduler* inParentScheduler = 0);

  virtual ~MPIScheduler();

  // eliminate copy, assignment and move
  MPIScheduler(const MPIScheduler&) = delete;
  MPIScheduler&
  operator=(const MPIScheduler&) = delete;
  MPIScheduler(MPIScheduler&&)   = delete;
  MPIScheduler&
  operator=(MPIScheduler&&) = delete;

  virtual void
  problemSetup(const ProblemSpecP& prob_spec, const MaterialManagerP& mat_manager);

  virtual void
  execute(int tgnum = 0, int iteration = 0);

  virtual SchedulerP
  createSubScheduler();

  virtual void
  processMPIRecvs(int test_type);

  void
  postMPISends(DetailedTask* task, int iteration);

  void
  postMPIRecvs(DetailedTask* task,
               bool only_old_recvs,
               int abort_point,
               int iteration);

  void
  runTask(DetailedTask* task, int iteration);

  void
  runReductionTask(DetailedTask* task);

  void
  compile()
  {
    d_num_messages   = 0;
    d_message_volume = 0;
    SchedulerCommon::compile();
  }

  // Performs the reduction task. (In threaded, Unified scheduler, a single
  // worker thread will execute this.)
  virtual void
  initiateReduction(DetailedTask* dtask);

  void
  computeNetRuntimeStats();

public:
  // timing statistics for Uintah infrastructure overhead
  enum TimingStatsEnum
  {
    TotalSend = 0,
    TotalRecv,
    TotalTest,
    TotalWait,
    TotalReduce,
    TotalTask
  };

  enum
  {
    TEST,
    WAIT_ONCE,
    WAIT_ALL
  };

  ReductionInfoMapper<TimingStatsEnum, double> d_mpi_info;
  MapInfoMapper<std::string, TaskStatsEnum, double> d_task_info;

  MPIScheduler* d_parent_scheduler{ nullptr };

protected:
  virtual void
  initiateTask(DetailedTask* task,
               bool only_old_recvs,
               int abort_point,
               int iteration);

  virtual void
  verifyChecksum();

  void
  emitTime(const char* label, double time);

  void
  outputTimingStats(const char* label);

protected:
  CommRequestPool d_sends{};
  CommRequestPool d_recvs{};

  std::vector<const char*> d_labels;
  std::vector<double> d_times;

  std::ofstream d_max_stats;
  std::ofstream d_avg_stats;

  std::atomic<unsigned int> d_num_messages{ 0 };
  double d_message_volume{ 0.0 };

  Timers::Simple d_exec_timer;
};

} // End namespace Uintah

#endif
