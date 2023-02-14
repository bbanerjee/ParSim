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

#ifndef CCA_COMPONENTS_SCHEDULERS_DETAILED_DEP_H
#define CCA_COMPONENTS_SCHEDULERS_DETAILED_DEP_H

#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>

#include <list>

namespace Uintah {

class DetailedTask;

struct DetailedDep
{
public:
  enum CommCondition
  {
    Always,
    FirstIteration,
    SubsequentIterations
  };

  DetailedDep(DetailedDep* next_in,
              Task::Dependency* comp_in,
              Task::Dependency* req_in,
              DetailedTask* toTask_in,
              const Patch* fromPatch_in,
              int matl_in,
              const IntVector& low_in,
              const IntVector& high_in,
              CommCondition cond_in)
    : m_next(next_in)
    , m_comp(comp_in)
    , m_req(req_in)
    , m_from_patch(fromPatch_in)
    , m_low(low_in)
    , m_high(high_in)
    , m_matl(matl_in)
    , m_comm_condition(cond_in)
    , m_patch_low(low_in)
    , m_patch_high(high_in)
  {
    ASSERT(Min(high_in - low_in, IntVector(1, 1, 1)) == IntVector(1, 1, 1));

    USE_IF_ASSERTS_ON(Patch::VariableBasis basis = Patch::translateTypeToBasis(
                        req_in->var->typeDescription()->getType(), true);)

    ASSERT(
      fromPatch_in == 0 ||
      (Min(
         m_low,
         fromPatch_in->getExtraLowIndex(basis, req_in->var->getBoundaryLayer())) ==
       fromPatch_in->getExtraLowIndex(basis, req_in->var->getBoundaryLayer())));

    ASSERT(
      fromPatch_in == 0 ||
      (Max(
         m_high,
         fromPatch_in->getExtraHighIndex(basis, req_in->var->getBoundaryLayer())) ==
       fromPatch_in->getExtraHighIndex(basis, req_in->var->getBoundaryLayer())));

    m_to_tasks.push_back(toTask_in);
  }

  // eliminate copy, assignment and move
  DetailedDep(const DetailedDep&) = delete;
  DetailedDep& operator=(const DetailedDep&) = delete;
  DetailedDep(DetailedDep&&) = delete;
  DetailedDep& operator=(DetailedDep&&) = delete;

  // As an arbitrary convention, non-data dependency have a nullptr fromPatch.
  // These types of dependencies exist between a modifying task and any task
  // that requires the data (from ghost cells in particular) before it is
  // modified, preventing the possibility of modifying data while it is being
  // used.
  bool isNonDataDependency() const { return (m_from_patch == nullptr); }

  DetailedDep* m_next;
  Task::Dependency* m_comp;
  Task::Dependency* m_req;
  std::list<DetailedTask*> m_to_tasks;
  const Patch* m_from_patch;
  IntVector m_low;
  IntVector m_high;
  int m_matl;

  // this is to satisfy a need created by the DynamicLoadBalancer.  To keep it
  // unrestricted on when it can perform, and to avoid a costly second recompile
  // on the next timestep, we add a comm condition which will send/recv data
  // based on whether some condition is met at run time - in this case whether
  // it is the first execution or not.
  CommCondition m_comm_condition;

  // for SmallMessages - if we don't copy the complete patch, we need to know
  // the range so we can store all segments properly
  IntVector m_patch_low;
  IntVector m_patch_high;
};

std::ostream&
operator<<(std::ostream& out, const Uintah::DetailedDep& task);

} // namespace Uintah

#endif // CCA_COMPONENTS_SCHEDULERS_DETAILED_DEP_H
