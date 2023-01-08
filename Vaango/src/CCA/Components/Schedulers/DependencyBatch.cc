/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2023 Biswajit Banerjee
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

#include <CCA/Components/Schedulers/DependencyBatch.h>

#include <Core/Parallel/MasterLock.h>
#include <Core/Util/DOUT.hpp>

#include <sstream>

namespace Uintah {

namespace {

Dout g_received_dbg("DependencyBatch",
                    "DependencyBatch",
                    "report when a DependencyBatch is received",
                    false);

Uintah::MasterLock g_dep_batch_mutex{};

}

DependencyBatch::DependencyBatch(int to,
                                 DetailedTask* fromTask,
                                 DetailedTask* toTask)
  : from_task(fromTask)
  , to_rank(to)
{
  to_tasks.push_back(toTask);
}

DependencyBatch::~DependencyBatch() {}

void
DependencyBatch::reset()
{
  d_received = false;
  d_made_mpi_request.store(false, std::memory_order_relaxed);
}

bool
DependencyBatch::makeMPIRequest()
{
  bool expected_val = false;
  return d_made_mpi_request.compare_exchange_strong(
    expected_val, true, std::memory_order_seq_cst);
}

void
DependencyBatch::received(const ProcessorGroup* pg)
{
  std::lock_guard<Uintah::MasterLock> dep_batch_lock(g_dep_batch_mutex);

  d_received = true;

  // set all the toVars to valid, meaning the MPI has been completed
  for (auto& to_var : d_to_vars) {
    to_var->setValid();
  }

  // prepare for placement into the external ready queue
  for (auto& to_task : to_tasks) {
    // if the count is 0, the task will add itself to the external ready queue
    to_task->decrementExternalDepCount();
    to_task->checkExternalDepCount();
  }

  // clear the variables that have outstanding MPI as they are completed now.
  d_to_vars.clear();
}

void
DependencyBatch::addVar(Variable* var)
{
  std::lock_guard<Uintah::MasterLock> dep_batch_lock(g_dep_batch_mutex);

  d_to_vars.push_back(var);
}

} // namespace Uintah
