/*
 * The MIT License
 *
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

#include <CCA/Components/Schedulers/TaskGraphCompTable.h>

#include <CCA/Components/Schedulers/DetailedTask.h>

namespace {
Uintah::Dout g_find_computes_dbg(
  "ComputeTableFindComputes",
  "ComputeTable",
  "info on computing task for particular requires",
  false);
Uintah::Dout g_detailed_deps_dbg("ComputeTableDetailedDeps",
                                 "ComputeTable",
                                 "detailed dep info for each DetailedTask",
                                 false);
}

namespace Uintah {

void
CompTable::remembercomp(DetailedTask* task,
                        Task::Dependency* comp,
                        const PatchSubset* patches,
                        const MaterialSubset* matls,
                        const ProcessorGroup* pg)
{
  if (patches && matls) {
    for (int p = 0; p < patches->size(); p++) {
      const Patch* patch = patches->get(p);
      for (int m = 0; m < matls->size(); m++) {
        int matl      = matls->get(m);
        Data* newData = scinew Data(task, comp, patch, matl);
        remembercomp(newData, pg);
      }
    }
  } else if (matls) {
    for (int m = 0; m < matls->size(); m++) {
      int matl      = matls->get(m);
      Data* newData = scinew Data(task, comp, 0, matl);
      remembercomp(newData, pg);
    }
  } else if (patches) {
    for (int p = 0; p < patches->size(); p++) {
      const Patch* patch = patches->get(p);
      Data* newData      = scinew Data(task, comp, patch, 0);
      remembercomp(newData, pg);
    }
  } else {
    Data* newData = scinew Data(task, comp, nullptr, 0);
    remembercomp(newData, pg);
  }
}

void
CompTable::remembercomp(Data* newData, const ProcessorGroup* pg)
{
  if (g_detailed_deps_dbg) {
    std::ostringstream message;
    message << "Rank-" << pg->myRank() << " remembercomp: " << *newData->m_comp
            << ", matl=" << newData->m_matl;
    if (newData->m_patch) {
      message << ", patch=" << *newData->m_patch;
    }
    DOUT(true, message.str());
  }

  TypeDescription::Type vartype =
    newData->m_comp->var->typeDescription()->getType();

  // can't have two computes for the same variable (need modifies)

  // ARS A VarLabel can have the following condition added:
  // allowMultipleComputes that allows multiple computes. As such why
  // do we check specifically for a ReductionVariable var and skip it?
  // Seems like we should use the conditional that is part of the
  // label to skip.

  // ARS - Treat sole vars the same as reduction vars??
  if (newData->m_comp->dep_type != Task::Modifies &&
      vartype != TypeDescription::Type::ReductionVariable &&
      vartype != TypeDescription::Type::SoleVariable) {
    if (m_data.lookup(newData)) {
      std::cout << "Multiple compute found:\n";
      std::cout << "  matl: " << newData->m_matl << "\n";
      if (newData->m_patch) {
        std::cout << "  patch: " << *newData->m_patch << "\n";
      }
      std::cout << "  " << *newData->m_comp << "\n";
      std::cout << "  " << *newData->m_dtask << "\n\n";
      std::cout << "  It was originally computed by the following task(s):\n";
      for (Data* old = m_data.lookup(newData); old != nullptr;
           old       = m_data.nextMatch(newData, old)) {
        std::cout << "  " << *old->m_dtask << std::endl;
        old->m_comp->task->displayAll(std::cout);
      }
      SCI_THROW(InternalError("Multiple computes for variable: " +
                                newData->m_comp->var->getName(),
                              __FILE__,
                              __LINE__));
    }
  }
  m_data.insert(newData);
}

bool
CompTable::findcomp(Task::Dependency* req,
                    const Patch* patch,
                    int matlIndex,
                    DetailedTask*& dt,
                    Task::Dependency*& comp,
                    const ProcessorGroup* pg)
{
  DOUTR(g_find_computes_dbg,
        "Rank-" << pg->myRank() << ": Finding comp of req: " << *req
                << " for task: " << *req->task << "/");

  Data key(nullptr, req, patch, matlIndex);
  Data* result = nullptr;
  for (Data* p = m_data.lookup(&key); p != nullptr;
       p       = m_data.nextMatch(&key, p)) {

    DOUT(g_find_computes_dbg,
         "Rank-" << pg->myRank()
                 << ": Examining comp from: " << p->m_comp->task->getName()
                 << ", order=" << p->m_comp->task->getSortedOrder());

    // TODO - fix why this assert is tripped when the gold standard,
    // MPM/ARL/NanoPillar2D_FBC_sym.ups is run using a non-optimized build.
    // On a debug, inputs/MPMdisks_complex.ups also hits this.
    // Clue: This assertion is tripped if there are two modifies() in a single
    // task.
    // ASSERT(!result || p->m_comp->task->getSortedOrder() !=
    // result->m_comp->task->getSortedOrder());

    if (p->m_comp->task->getSortedOrder() < req->task->getSortedOrder()) {
      if (!result || p->m_comp->task->getSortedOrder() >
                       result->m_comp->task->getSortedOrder()) {

        DOUT(g_find_computes_dbg,
             "Rank-" << pg->myRank() << ": New best is comp from: "
                     << p->m_comp->task->getName()
                     << ", order=" << p->m_comp->task->getSortedOrder());

        result = p;
      }
    }
  }

  if (result) {

    DOUT(g_find_computes_dbg,
         "Rank-" << pg->myRank()
                 << ": Found comp at: " << result->m_comp->task->getName()
                 << ", order=" << result->m_comp->task->getSortedOrder());

    dt   = result->m_dtask;
    comp = result->m_comp;
    return true;
  } else {
    return false;
  }
}

bool
CompTable::findReductionComps(Task::Dependency* req,
                              const Patch* patch,
                              int matlIndex,
                              std::vector<DetailedTask*>& creators,
                              const ProcessorGroup* pg)
{
  // reduction variables for each level can be computed by several tasks (once
  // per patch) return the list of all tasks nearest the req

  Data key(nullptr, req, patch, matlIndex);
  int bestSortedOrder = -1;
  for (Data* p = m_data.lookup(&key); p != nullptr;
       p       = m_data.nextMatch(&key, p)) {
    DOUT(g_detailed_deps_dbg,
         "Rank-" << pg->myRank()
                 << " Examining comp from: " << p->m_comp->task->getName()
                 << ", order=" << p->m_comp->task->getSortedOrder() << " ("
                 << req->task->getName()
                 << " order: " << req->task->getSortedOrder() << ")");

    if (p->m_comp->task->getSortedOrder() < req->task->getSortedOrder() &&
        p->m_comp->task->getSortedOrder() >= bestSortedOrder) {
      if (p->m_comp->task->getSortedOrder() > bestSortedOrder) {
        creators.clear();
        bestSortedOrder = p->m_comp->task->getSortedOrder();
        DOUT(g_detailed_deps_dbg,
             "Rank-" << pg->myRank()
                     << "    New best sorted order: " << bestSortedOrder);
      }
      DOUT(g_detailed_deps_dbg,
           "Rank-" << pg->myRank()
                   << " Adding comp from: " << p->m_comp->task->getName()
                   << ", order=" << p->m_comp->task->getSortedOrder());
      creators.push_back(p->m_dtask);
    }
  }
  return (creators.size() > 0);
}

} // namespace Uintah
