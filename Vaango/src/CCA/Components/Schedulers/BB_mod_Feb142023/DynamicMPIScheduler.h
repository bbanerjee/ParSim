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

#ifndef VAANGO_CCA_COMPONENTS_SCHEDULERS_DMPISCHEDULER_H
#define VAANGO_CCA_COMPONENTS_SCHEDULERS_DMPISCHEDULER_H

#include <CCA/Components/Schedulers/MPIScheduler.h>

namespace Uintah {

class Task;

class DynamicMPIScheduler : public MPIScheduler
{

public:
  DynamicMPIScheduler(const ProcessorGroup* myworld,
                      DynamicMPIScheduler* parentScheduler = 0);

  virtual ~DynamicMPIScheduler();

  // eliminate copy, assignment and move
  DynamicMPIScheduler(const DynamicMPIScheduler&) = delete;
  DynamicMPIScheduler&
  operator=(const DynamicMPIScheduler&)      = delete;
  DynamicMPIScheduler(DynamicMPIScheduler&&) = delete;
  DynamicMPIScheduler&
  operator=(DynamicMPIScheduler&&) = delete;

  virtual void
  problemSetup(const ProblemSpecP& prob_spec, const MaterialManagerP& state);

  virtual SchedulerP
  createSubScheduler();

  virtual void
  execute(int tgnum = 0, int iteration = 0);

  virtual bool
  useInternalDeps()
  {
    return !d_is_copy_data_timestep;
  }

private:

  QueueAlg d_task_queue_algo{ MostMessages };
};

} // End namespace Uintah

#endif
