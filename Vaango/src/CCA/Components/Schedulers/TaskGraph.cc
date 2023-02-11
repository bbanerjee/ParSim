/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#include <CCA/Components/Schedulers/TaskGraph.h>

#include <CCA/Components/Schedulers/DetailedTasks.h>
#include <CCA/Components/Schedulers/OnDemandDataWarehouse.h>
#include <CCA/Components/Schedulers/SchedulerCommon.h>
#include <CCA/Components/Schedulers/TaskGraphCompTable.h>

#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/LoadBalancer.h>

#include <Core/Containers/FastHashTable.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/TypeMismatchException.h>

#include <Core/Grid/Grid.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>

#include <Core/Malloc/Allocator.h>

#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/ProcessorGroup.h>

#include <Core/Util/DOUT.hpp>
#include <Core/Util/FancyAssert.h>
#include <Core/Util/ProgressiveWarning.h>

#include <iostream>
#include <map>
#include <memory>
#include <sstream>

namespace {

Uintah::Dout g_tg_phase_dbg(
  "TaskGraphPhases",
  "TaskGraph",
  "task phase assigned to each task by the task graph",
  false);
Uintah::Dout g_proc_neighborhood_dbg(
  "ProcNeighborhood",
  "TaskGraph",
  "info on local or distal patch neighborhoods",
  false);
Uintah::Dout g_find_computes_dbg(
  "FindComputes",
  "TaskGraph",
  "info on computing task for particular requires",
  false);
Uintah::Dout g_detailed_deps_dbg("TaskGraphDetailedDeps",
                                 "TaskGraph",
                                 "detailed dep info for each DetailedTask",
                                 false);
Uintah::Dout g_add_task_dbg(
  "TaskGraphAddTask",
  "TaskGraph",
  "report task name, computes and requires: every task",
  false);
Uintah::Dout g_detailed_task_dbg("TaskGraphDetailedTasks",
                                 "TaskGraph",
                                 "high-level info on creation of DetailedTasks",
                                 false);
Uintah::Dout g_topological_deps_dbg(
  "TopologicalDetailedDeps",
  "TaskGraph",
  "topologiocal sort detailed dependnecy info",
  false);

} // namespace

namespace Uintah {

TaskGraph::TaskGraph(SchedulerCommon* sc,
                     const ProcessorGroup* pg,
                     Scheduler::tgType type,
                     [[maybe_unused]] int index)
  : d_scheduler(sc)
  , d_proc_group(pg)
  , d_type(type)

{
  d_load_balancer =
    dynamic_cast<LoadBalancer*>(d_scheduler->getPort("load balancer"));
}

TaskGraph::~TaskGraph()
{
  initialize(); // Frees all of the memory...
}

void
TaskGraph::initialize()
{
  if (d_detailed_tasks) {
    delete d_detailed_tasks;
  }

  d_tasks.clear();

  d_num_task_phases   = 0;
  d_current_iteration = 0;
}

void
TaskGraph::nullSort(std::vector<Task*>& tasks)
{
  DOUTR(g_detailed_task_dbg, " TaskGraph::nullSort ");

  // No longer going to sort them... let the UnifiedScheduler (threaded) take
  // care of calling the tasks when all dependencies are satisfied. Sorting the
  // tasks causes problem because now tasks (actually task groups) run in
  // different orders on different MPI processes.
  int n = 0;
  for (auto task : d_tasks) {

    std::ostringstream mesg;
    mesg << "      " << std::left << std::setw(50) << task->getName();

    // For all reduction tasks filtering out the one that is not in
    // LevelReductionTasksMap
    if (task->getType() == Task::Reduction) {

      for (auto& [varlabel_matl, reduction_task] :
           d_scheduler->d_reduction_tasks) {
        if (task == reduction_task) {
          task->setSortedOrder(n++);
          tasks.push_back(task.get());
          mesg << "     Added to sorted tasks (reduction task)";
          break;
        }
      }
    } else {
      mesg << "     Added to sorted tasks (normal task)";
      task->setSortedOrder(n++);
      tasks.push_back(task.get());
    }
    DOUTR(g_detailed_task_dbg, mesg.str());
  }
}

void
TaskGraph::addTask(std::shared_ptr<Task> task,
                   const PatchSet* patchset,
                   const MaterialSet* matlset)
{
  task->setSets(patchset, matlset);

  // empty tasks?
  if ((patchset && patchset->totalsize() == 0) ||
      (matlset && matlset->totalsize() == 0)) {
    task.reset(); // release ownership of managed std::shared_ptr<Task>
    DOUTR(g_topological_deps_dbg, " Killing empty task: " << *task);
  } else {
    d_tasks.push_back(task);

    DOUTR(g_add_task_dbg,
          " TG[" << d_index << "] adding task: " << task->getName());
    task->displayAll_DOUT(g_add_task_dbg);

#if 0
    // This snippet will find all the tasks that require a label
    for (Task::Dependency* m_req = task->getRequires(); m_req != nullptr; m_req = m_req->m_next) {
      const VarLabel* label = m_req->m_var;
      std::string name = label->getName();
      std::string search_name("abskg");
      if (name == search_name) {
        std::cout << "\n" << "Rank-" << std::setw(4) << std::left << Parallel::getMPIRank() << " "
                  << task->getName() << " requires label " << "\"" << search_name << "\"";
        task->displayAll_DOUT(g_add_task_dbg);
      }
    }
#endif
  }

  // does this TG contain requires with halo > MAX_HALO_DEPTH
  if (task->hasDistalRequires()) {
    d_has_distal_requires = true;
  }

}

void
TaskGraph::createDetailedTask(Task* task,
                              const PatchSubset* patches,
                              const MaterialSubset* matls)
{
  auto* dt = scinew DetailedTask(task, patches, matls, d_detailed_tasks);

#if SCI_ASSERTION_LEVEL > 0
  if (task->getType() == Task::Reduction) {
    // reduction tasks should have exactly 1 require, and it should be a modify
    ASSERT(task->getModifies() != nullptr);
  }
#endif

  d_detailed_tasks->add(dt);
}

auto
TaskGraph::createDetailedTasks(bool useInternalDeps,
                               const GridP& grid,
                               const GridP& oldGrid,
                               const bool hasDistalReqs /* = false */)
  -> DetailedTasks*
{
  std::vector<Task*> sorted_tasks;

  // TODO plz leave this commented line alone, APH 01/07/15
  // topologicalSort(sorted_tasks);
  nullSort(sorted_tasks);

  ASSERT(grid != nullptr);

  // Create two neighborhoods.  One that keeps track of patches and procs within
  // just a few ghost cells and another that keeps track of patches and procs to
  // within the maximum known ghost cells among all tasks in the task graph.
  //
  // The idea is that if a DetailedTask only has simulation variables which
  // throughout the simulation never go beyond 1 or 2 ghost cells, then those
  // tasks do not concern themselves with patch variables far away.  No need to
  // create tasks on the task graph as we'll never bother with their
  // dependencies.  These tasks would then use the local neighborhood of the two
  // neighborhoods. But some tasks (like radiation tasks) require ghost cells
  // hundreds or thousands of cells out. So we need many more of these tasks in
  // the task graph so we can process their dependencies.
  //
  // Uintah history moment: The older neighborhood approach before this one was
  // to only use one, and so if there were hundreds of tasks and one of those
  // hundreds needed, say, 250 ghost cells, then a taskgraph was made to check
  // for the possibility of every task's variables going out 250 ghost cells.
  // Even though most variables only go out 1 or 2 ghost cells out.  This
  // resulted in a taskgraph of way too many DetailedTasks, the vast majority of
  // which were useless as they were far away from this proc and would never
  // depend on tasks on this proc. On one Titan run, this resulted in task graph
  // compilation times hitting several hours.  So we made a second neighborhood
  // to fix that.
  //
  // (Note, it's possible to create more neighborhoods for each kind of ghost
  // cell configuration, but the total savings seems not important at the time
  // of this writing)

  d_load_balancer->createNeighborhood(grid, oldGrid);

  const auto local_procs  = d_load_balancer->getNeighborhoodProcessors();
  const auto distal_procs = d_load_balancer->getDistalNeighborhoodProcessors();

  d_detailed_tasks =
    scinew DetailedTasks(d_scheduler,
                         d_proc_group,
                         this,
                         (hasDistalReqs ? distal_procs : local_procs),
                         useInternalDeps);

  // Go through every task, find the max ghost cell associated with each
  // varlabel/matl, and remember that.
  //
  // This will come in handy so we can create DetailedTasks for any var
  // containing a requires, modifies, or computes within range of that
  // particular data.  For example, if a task's varlabel/matl requires variable
  // needs ghost cells of 250, then that DetailedTask needs to have a
  // requirements dependency set up for all tasks 250 cells away.
  //
  // Also, if a task needing that varlabel/matl has a requires with 250 ghost
  // cells, it will also need to have knowledge of what other task originally
  // computed it.  This is so MPI sends and receives can be created to know
  // where to send these ghost cells and also where to receive them from.  This
  // means that if one task only computes a variable, but another task requires
  // 250 ghost cells of that variable, we still need many instances of
  // DetailedTasks for the computing task in the task graph.
  //
  // (Note, this has an underlying assumption that tasks are assigned uniformly
  // to patches.  If for some reason we had very different tasks which ran on
  // some patches but not on others, this approach overestimates and should be
  // refined).

  for (auto& task : sorted_tasks) {

    // Assuming that variables on patches get the same amount of ghost cells on
    // any patch in the simulation.  This goes for multimaterial vars as well.

    // Get the tasks's level.  (Task's aren't explicitly assigned a level, but
    // they are assigned a set of patches corresponding to that level.  So grab
    // the task's 0th patch and get the level it is on.)

    // TODO: OncePerProc tasks are assigned multiple levels.  How does that
    // affect this?  Brad P. 11/6/2016
    int levelID        = 0;
    const PatchSet* ps = task->getPatchSet();

    // Reduction tasks don't have patches, filter them out.
    if ((ps != nullptr) && (ps->size() > 0)) {

      const PatchSubset* pss = ps->getSubset(0);
      if (pss && pss->size()) {
        levelID = pss->get(0)->getLevel()->getID();
      }
    }

    // look through requires
    for (auto req = task->getRequires(); req != nullptr; req = req->next) {
      std::string key = req->var->getName();

      // This offset is not fully helpful by itself. It only indicates how much
      // its off from the level that the patches are assigned to.  It does *not*
      // indicate if the offset is positive or negative, as it can be either
      // positive or negative.
      //
      // If a task's patches are on the coarse level and the offset is 1, then
      // the offset is positive. If a task's patches are on the fine level and
      // the offset is 1, then the offset is negative.

      const int levelOffset = req->level_offset;
      int trueLevel         = levelID;
      if (req->patches_dom == Task::CoarseLevel) {
        trueLevel -= levelOffset;
      } else if (req->patches_dom == Task::FineLevel) {
        trueLevel += levelOffset;
      }

      int ngc = req->num_ghost_cells;

      DOUTR(g_proc_neighborhood_dbg,
            "In task: " << task->getName()
                        << "Checking for max ghost cell for requirements var "
                        << key << " which has: " << ngc << " ghost cells."
                        << " and is on level: " << trueLevel << ".");

      LabelLevel labelLevel(key, trueLevel);
      auto it = d_max_ghost_for_varlabelmap.find(labelLevel);
      if (it != d_max_ghost_for_varlabelmap.end()) {
        if (it->second < ngc) {
          it->second = ngc;
        }
      } else {
        d_max_ghost_for_varlabelmap.emplace(labelLevel, ngc);
      }
    }

    // Can modifies have ghost cells?
    for (auto modifies = task->getModifies(); modifies != nullptr;
         modifies      = modifies->next) {
      std::string key = modifies->var->getName();
      int levelOffset = modifies->level_offset;
      int trueLevel   = levelID;
      if (modifies->patches_dom == Task::CoarseLevel) {
        trueLevel -= levelOffset;
      } else if (modifies->patches_dom == Task::FineLevel) {
        trueLevel += levelOffset;
      }
      int ngc = modifies->num_ghost_cells;

      DOUTR(g_proc_neighborhood_dbg,
            "In task: " << task->getName()
                        << "Checking for max ghost cell for modifies var "
                        << key << " which has: " << ngc << " ghost cells."
                        << " and is on level: " << trueLevel << ".");

      LabelLevel labelLevel(key, trueLevel);
      auto it = d_max_ghost_for_varlabelmap.find(labelLevel);
      if (it != d_max_ghost_for_varlabelmap.end()) {
        if (it->second < ngc) {
          it->second = ngc;
        }
      } else {
        d_max_ghost_for_varlabelmap.emplace(labelLevel, ngc);
      }
    }

    // We don't care about computes ghost cells

  } // end for sorted_tasks.size()

  if (g_proc_neighborhood_dbg) {
    for (auto kv : d_max_ghost_for_varlabelmap) {
      DOUTR(g_proc_neighborhood_dbg,
            "For varlabel " << kv.first.m_key
                            << " on level: " << kv.first.m_level
                            << " the max ghost cell is: " << kv.second);
    }
  }

  // Now loop again, setting the task's max ghost cells to the max ghost cell
  // for a given varLabel
  for (auto& task : sorted_tasks) {
    int levelID        = 0;
    const PatchSet* ps = task->getPatchSet();
    if (ps && ps->size()) {
      const PatchSubset* pss = ps->getSubset(0);
      if (pss && pss->size()) {
        const Level* level = pss->get(0)->getLevel();
        levelID            = level->getID();
      }
    }

    task->m_max_ghost_cells[levelID] = 0; // default

    // Again assuming all vars for a label get the same amount of ghost cells,

    // check requires
    for (auto req = task->getRequires(); req != nullptr; req = req->next) {
      std::string key = req->var->getName();
      int levelOffset = req->level_offset;
      int trueLevel   = levelID;
      if (req->patches_dom == Task::CoarseLevel) {
        trueLevel -= levelOffset;
      } else if (req->patches_dom == Task::FineLevel) {
        trueLevel += levelOffset;
      }

      DOUTR(g_proc_neighborhood_dbg,
            "For task: " << task->getName() << " on level " << trueLevel
                         << " from levelID: " << levelID
                         << " and levelOffset: " << levelOffset
                         << " checking out requires var: " << key);

      LabelLevel labelLevel(key, trueLevel);
      auto it = task->m_max_ghost_cells.find(trueLevel);
      if (it != task->m_max_ghost_cells.end()) {
        if (task->m_max_ghost_cells[trueLevel] <
            d_max_ghost_for_varlabelmap[labelLevel]) {
          task->m_max_ghost_cells[trueLevel] =
            d_max_ghost_for_varlabelmap[labelLevel];
        }
      } else {
        task->m_max_ghost_cells[trueLevel] =
          d_max_ghost_for_varlabelmap[labelLevel];
      }
    }

    // check modifies
    for (auto modifies = task->getModifies(); modifies != nullptr;
         modifies      = modifies->next) {
      std::string key = modifies->var->getName();
      int levelOffset = modifies->level_offset;
      int trueLevel   = levelID;
      if (modifies->patches_dom == Task::CoarseLevel) {
        trueLevel -= levelOffset;
      } else if (modifies->patches_dom == Task::FineLevel) {
        trueLevel += levelOffset;
      }

      DOUTR(g_proc_neighborhood_dbg,
            "For task: " << task->getName() << " on level " << trueLevel
                         << " from levelID: " << levelID
                         << " and levelOffset: " << levelOffset
                         << " checking out modifies var: " << key);

      LabelLevel labelLevel(key, trueLevel);
      auto it = task->m_max_ghost_cells.find(trueLevel);
      if (it != task->m_max_ghost_cells.end()) {
        if (task->m_max_ghost_cells[trueLevel] <
            d_max_ghost_for_varlabelmap[labelLevel]) {
          task->m_max_ghost_cells[trueLevel] =
            d_max_ghost_for_varlabelmap[labelLevel];
        }
      } else {
        task->m_max_ghost_cells[trueLevel] =
          d_max_ghost_for_varlabelmap[labelLevel];
      }
    }

    // check computes
    for (auto comps = task->getComputes(); comps != nullptr;
         comps      = comps->next) {
      std::string key = comps->var->getName();
      int levelOffset = comps->level_offset;
      int trueLevel   = levelID;
      if (comps->patches_dom == Task::CoarseLevel) {
        trueLevel -= levelOffset;
      } else if (comps->patches_dom == Task::FineLevel) {
        trueLevel += levelOffset;
      }
      DOUTR(g_proc_neighborhood_dbg,
            "For task: " << task->getName() << " on level " << trueLevel
                         << " from levelID: " << levelID
                         << " and levelOffset: " << levelOffset
                         << " checking out computes var: " << key);
      LabelLevel labelLevel(key, trueLevel);
      auto it = task->m_max_ghost_cells.find(trueLevel);
      if (it != task->m_max_ghost_cells.end()) {
        if (task->m_max_ghost_cells[trueLevel] <
            d_max_ghost_for_varlabelmap[labelLevel]) {
          task->m_max_ghost_cells[trueLevel] =
            d_max_ghost_for_varlabelmap[labelLevel];
        }
      } else {
        task->m_max_ghost_cells[trueLevel] =
          d_max_ghost_for_varlabelmap[labelLevel];
      }
    }

    for (auto& kv : task->m_max_ghost_cells) {
      DOUTR(g_proc_neighborhood_dbg,
            "For task: " << task->getName() << " on level " << kv.first
                         << " the largest found max ghost cells so far is: "
                         << kv.second);
    }
  }

  // Now proceed looking within the appropriate neighborhood defined by the max
  // ghost cells a task needs to know about.
  size_t tot_detailed_tasks  = 0u;
  size_t tot_normal_tasks    = 0u;
  size_t tot_output_tasks    = 0u;
  size_t tot_opp_tasks       = 0u;
  size_t tot_reduction_tasks = 0u;

  for (auto& task : sorted_tasks) {

    size_t num_detailed_tasks = 0u;

    DOUTR(g_proc_neighborhood_dbg,
          "Looking for max ghost vars for task: " << task->getName());

    const PatchSet* ps    = task->getPatchSet();
    const MaterialSet* ms = task->getMaterialSet();

    int levelID = 0;
    if (ps && ps->size()) {
      const PatchSubset* pss = ps->getSubset(0);
      if (pss && pss->size()) {
        const Level* level = pss->get(0)->getLevel();
        levelID            = level->getID();
      }
    }

    //__________________________________
    // valid patch and matl sets
    if (ps && ms) {

      // OncePerProc tasks - only create OncePerProc tasks and output tasks once
      // on each processor.
      if (task->getType() == Task::OncePerProc ||
          task->getType() == Task::Hypre) {
        // only schedule this task on processors in the neighborhood

        // NOTE THE MAP::AT METHOD NEEDS TO BE SAFEGUARDED.
        // Is it reasonable to set neighborhood_procs = local_procs when it
        // fails??
        std::unordered_set<int> neighborhood_procs;
        if (task->m_max_ghost_cells.find(levelID) !=
            task->m_max_ghost_cells.end()) {
          neighborhood_procs =
            (task->m_max_ghost_cells.at(levelID) >= MAX_HALO_DEPTH)
              ? distal_procs
              : local_procs;
        } else {
          DOUTALL(true, "*********** Bad level ID " << levelID);

          neighborhood_procs = local_procs;
        }

        for (int neighborhood_proc : neighborhood_procs) {
          const PatchSubset* pss = ps->getSubset(neighborhood_proc);
          for (int m = 0; m < ms->size(); m++) {
            const MaterialSubset* mss = ms->getSubset(m);
            createDetailedTask(task, pss, mss);
            ++num_detailed_tasks;
            ++tot_opp_tasks;
          }
        }
      }

      //__________________________________
      // Output tasks
      else if (task->getType() == Task::Output) {
        // Compute rank that handles output for this process.
        int handling_rank =
          (d_proc_group->myRank() / d_load_balancer->getNthRank()) *
          d_load_balancer->getNthRank();

        // Only schedule output task for the subset involving our rank.
        const PatchSubset* pss = ps->getSubset(handling_rank);

        // Don't schedule if there are no patches.
        if (pss->size() > 0) {
          for (auto m = 0; m < ms->size(); m++) {
            const MaterialSubset* mss = ms->getSubset(m);
            createDetailedTask(task, pss, mss);
            ++num_detailed_tasks;
            ++tot_output_tasks;
          }
        }
      }

      //__________________________________
      // Normal tasks
      else {
        const int ps_size = ps->size();
        for (int ps_index = 0; ps_index < ps_size; ps_index++) {
          const PatchSubset* pss = ps->getSubset(ps_index);

          if (pss->size() > 0) {

            // Make tasks in our neighborhood.  If there are multiple levels
            // involved in the reqs of a task, then the levelID should be the
            // fine level
            DOUTR(g_proc_neighborhood_dbg,
                  "For task: " << task->getName()
                               << " looking for max ghost cells for level: "
                               << levelID);

            // Still make sure we have an entry for this task on this level.
            // Some tasks can go into the task graph without any requires,
            // modifies, or computes.
            bool search_distal_requires = false;

            for (const auto kv : task->m_max_ghost_cells) {
              int levelIDTemp = kv.first;

              search_distal_requires = (kv.second >= MAX_HALO_DEPTH);

              DOUTR(g_proc_neighborhood_dbg,
                    "Rank-"
                      << d_proc_group->myRank() << " for: " << task->getName()
                      << " on level: " << levelIDTemp
                      << " with task max ghost cells: " << kv.second
                      << " Seeing if patch subset: " << *pss
                      << " is in neighborhood with search_distal_requires: "
                      << search_distal_requires);

              if (search_distal_requires) {
                break;
              }
            }

            if (d_load_balancer->inNeighborhood(pss, search_distal_requires)) {
              DOUTR(g_proc_neighborhood_dbg, "Yes, it was in the neighborhood");
              for (int m = 0; m < ms->size(); m++) {
                const MaterialSubset* mss = ms->getSubset(m);
                createDetailedTask(task, pss, mss);
                ++num_detailed_tasks;
                ++tot_normal_tasks;
              }
            }
          } else {
            DOUTR(g_proc_neighborhood_dbg,
                  " for: " << task->getName() << " SKipping patch subset: "
                           << *pss << " because it has no patches");
          }
        }
        DOUTR(g_detailed_task_dbg,
              " created: " << num_detailed_tasks << " (" << task->getType()
                           << ") DetailedTasks for: " << task->getName()
                           << "\t on level: " << levelID);
      }
    } // end valid patch and matl sets
    // these tasks could be, e.g. DataArchiver::outputGlobalVars or
    // DataArchiver::outputVariables (CheckpointReduction)
    else if (!ps && !ms) {
      createDetailedTask(task, nullptr, nullptr);
      ++num_detailed_tasks;
      ++tot_output_tasks;
    }
    // reduction variable related tasks with material subsets but without
    // patches
    else if (!ps) {

      if (task->getType() == Task::Reduction ||
          task->getType() == Task::OutputGlobalVars) {
        for (int m = 0; m < ms->size(); m++) {

          const MaterialSubset* mss = ms->getSubset(m);
          createDetailedTask(task, nullptr, mss);
          ++num_detailed_tasks;

          if (task->getType() == Task::OutputGlobalVars) {
            ++tot_output_tasks;
          } else {
            ++tot_reduction_tasks;
          }
          DOUTR(g_detailed_task_dbg,
                " created: " << num_detailed_tasks << " (" << task->getType()
                             << ") DetailedTasks for: " << task->getName()
                             << "\t level : " << levelID
                             << " materialSubset: " << *mss);
        }
      }
    } else {
      SCI_THROW(InternalError("Task (" + task->getName() +
                                ") has PatchSet, but no MaterialSet",
                              __FILE__,
                              __LINE__));
    }
    tot_detailed_tasks += num_detailed_tasks;
  }

  DOUT(g_detailed_task_dbg,
       "\nRank-" << d_proc_group->myRank() << " created: " << tot_detailed_tasks
                 << "  total DetailedTasks in TG: " << d_index << " with:\n"
                 << "\t" << tot_normal_tasks << " total      Normal tasks\n"
                 << "\t" << tot_opp_tasks << " total OncePerProc tasks\n"
                 << "\t" << tot_output_tasks << " total      Output tasks\n"
                 << "\t" << tot_reduction_tasks
                 << " total   Reduction tasks\n");

  d_load_balancer->assignResources(*d_detailed_tasks);

  // scrub counts are created via addScrubCount() through this call ( via
  // possiblyCreateDependency() )
  createDetailedDependencies();

  if (d_detailed_tasks->getExtraCommunication() > 0 &&
      d_proc_group->myRank() == 0) {
    std::cout << d_proc_group->myRank()
              << "  Warning: Extra communication.  This taskgraph on this rank "
                 "overcommunicates about "
              << d_detailed_tasks->getExtraCommunication() << " cells\n";
  }

  if (d_proc_group->nRanks() > 1) {
    d_detailed_tasks->assignMessageTags(d_index);
  }

  d_detailed_tasks->computeLocalTasks();
  d_detailed_tasks->makeDWKeyDatabase();

  return d_detailed_tasks;
} // end TaskGraph::createDetailedTasks

void
TaskGraph::createDetailedDependencies()
{
  // Collect all of the computes
  CompTable ct;
  for (int i = 0; i < d_detailed_tasks->numTasks(); i++) {
    DetailedTask* dtask = d_detailed_tasks->getTask(i);

    if (g_topological_deps_dbg) {
      std::ostringstream message;
      DOUTR(true,
            "  createDetailedDependencies for:(" << dtask->getName() << ")");

      for (const Task::Dependency* req = dtask->getTask()->getRequires();
           req != nullptr;
           req = req->next) {
        DOUTR(true, "         requires: " << *req);
      }
      for (const Task::Dependency* comp = dtask->getTask()->getComputes();
           comp != nullptr;
           comp = comp->next) {
        DOUTR(true, "         computes: " << *comp);
      }
      for (const Task::Dependency* mod = dtask->getTask()->getModifies();
           mod != nullptr;
           mod = mod->next) {
        DOUTR(true, "         modifies: " << *mod);
      }
      DOUT(true, message.str());
    }

    remembercomps(dtask, dtask->m_task->getComputes(), ct);
    remembercomps(dtask, dtask->m_task->getModifies(), ct);
  }

  // Assign task phase number based on the reduction tasks so a mixed thread/mpi
  // scheduler won't have out of order reduction problems.
  int curr_phase     = 0;
  int curr_num_comms = 0;
  for (int i = 0; i < d_detailed_tasks->numTasks(); i++) {
    DetailedTask* dtask    = d_detailed_tasks->getTask(i);
    dtask->m_task->m_phase = curr_phase;

    DOUTR(g_tg_phase_dbg, " Task: " << *dtask << " phase: " << curr_phase);

    if (dtask->m_task->getType() == Task::Reduction) {
      dtask->m_task->m_comm = curr_num_comms;
      curr_num_comms++;
      curr_phase++;
    } else if (dtask->m_task->usesMPI()) {
      curr_phase++;
    }
  }

  // Uintah::MPI::Comm_dup happens here
  d_proc_group->setGlobalComm(curr_num_comms);
  d_num_task_phases = curr_phase + 1;

  // Go through the modifies/requires and create data dependencies as
  // appropriate
  for (int i = 0; i < d_detailed_tasks->numTasks(); i++) {
    DetailedTask* dtask = d_detailed_tasks->getTask(i);

    if (g_topological_deps_dbg && (dtask->m_task->getRequires() != nullptr)) {
      DOUTR(true, " Looking at requires of detailed task: " << *dtask);
    }

    createDetailedDependencies(dtask, dtask->m_task->getRequires(), ct, false);

    if (g_topological_deps_dbg && (dtask->m_task->getModifies() != nullptr)) {
      DOUTR(true, " Looking at modifies of detailed task: " << *dtask);
    }

    createDetailedDependencies(dtask, dtask->m_task->getModifies(), ct, true);
  }

  DOUTR(g_detailed_task_dbg, " Done creating detailed tasks");
}

void
TaskGraph::remembercomps(DetailedTask* dtask,
                         Task::Dependency* comp,
                         CompTable& ct)
{
  // calling getPatchesUnderDomain can get expensive on large processors.  Thus
  // we cache results and use them on the next call.  This works well because
  // comps are added in order and they share the same patches under the domain
  const PatchSubset *cached_task_patches = nullptr,
                    *cached_comp_patches = nullptr;
  constHandle<PatchSubset> cached_patches;

  for (; comp != nullptr; comp = comp->next) {

    TypeDescription::Type vartype = comp->var->typeDescription()->getType();

    // ARS - Treat sole vars the same as reduction vars??
    if (vartype == TypeDescription::Type::ReductionVariable ||
        vartype == TypeDescription::Type::SoleVariable) {
      // this is either the task computing the var, modifying it, or the
      // reduction itself
      ct.remembercomp(dtask, comp, nullptr, comp->matls, d_proc_group);
    } else {
      // Normal tasks
      constHandle<PatchSubset> patches;

      // if the patch pointer on both the dep and the task have not changed then
      // use the cached result
      if (dtask->m_patches == cached_task_patches &&
          comp->patches == cached_comp_patches) {
        patches = cached_patches;
      } else {
        // compute the intersection
        patches = comp->getPatchesUnderDomain(dtask->m_patches);
        // cache the result for the next iteration
        cached_patches      = patches;
        cached_task_patches = dtask->m_patches;
        cached_comp_patches = comp->patches;
      }
      constHandle<MaterialSubset> matls =
        comp->getMaterialsUnderDomain(dtask->m_matls);
      if (!patches->empty() && !matls->empty()) {
        ct.remembercomp(
          dtask, comp, patches.get_rep(), matls.get_rep(), d_proc_group);
      }
    }
  }
}

void
TaskGraph::remapTaskDWs(int dwmap[])
{
  // the point of this function is for using the multiple taskgraphs.
  // When you execute a taskgraph a subsequent time, you must rearrange the DWs
  // to point to the next point-in-time's DWs.
  int levelmin = 999;
  for (auto& d_task : d_tasks) {
    d_task->setMapping(dwmap);

    // for the Int timesteps, we have tasks on multiple levels.
    // we need to adjust based on which level they are on, but first
    // we need to find the coarsest level.  The NewDW is relative to the
    // coarsest level executing in this taskgraph.
    if (d_type == Scheduler::IntermediateTaskGraph &&
        (d_task->getType() != Task::Output &&
         d_task->getType() != Task::Hypre &&
         d_task->getType() != Task::OncePerProc)) {
      if (d_task->getType() == Task::OncePerProc ||
          d_task->getType() != Task::Hypre ||
          d_task->getType() == Task::Output) {
        levelmin = 0;
        continue;
      }

      const PatchSet* ps = d_task->getPatchSet();
      if (!ps) {
        continue;
      }
      const Level* l = getLevel(ps);
      levelmin       = Min(levelmin, l->getIndex());
    }
  }

  DOUTR(g_detailed_task_dbg,
        " Basic mapping "
          << "Old " << dwmap[Task::OldDW] << " New " << dwmap[Task::NewDW]
          << " CO " << dwmap[Task::CoarseOldDW] << " CN "
          << dwmap[Task::CoarseNewDW] << " levelmin " << levelmin);

  if (d_type == Scheduler::IntermediateTaskGraph) {
    // fix the CoarseNewDW for finer levels.  The CoarseOld will only matter
    // on the level it was originally mapped, so leave it as it is
    dwmap[Task::CoarseNewDW] = dwmap[Task::NewDW];
    for (auto& d_task : d_tasks) {
      if (d_task->getType() != Task::Output &&
          d_task->getType() != Task::Hypre &&
          d_task->getType() != Task::OncePerProc) {
        const PatchSet* ps = d_task->getPatchSet();
        if (!ps) {
          continue;
        }
        if (getLevel(ps)->getIndex() > levelmin) {
          d_task->setMapping(dwmap);
          DOUTR(g_detailed_task_dbg,
                d_task->getName()
                  << " mapping "
                  << "Old " << dwmap[Task::OldDW] << " New "
                  << dwmap[Task::NewDW] << " CO " << dwmap[Task::CoarseOldDW]
                  << " CN " << dwmap[Task::CoarseNewDW]
                  << " (levelmin=" << levelmin << ")");
        }
      }
    }
  }
}

void
TaskGraph::createDetailedDependencies(DetailedTask* dtask,
                                      Task::Dependency* req,
                                      CompTable& ct,
                                      bool modifies)
{
  int my_rank = d_proc_group->myRank();

  for (; req != nullptr; req = req->next) {

    // ARS reduction vars seem be handled below with type checks. I
    // would say it is not not needed.
    // ARS - Should Reduction and Sole variables be treated the same??
    if (d_scheduler->isOldDW(req->mapDataWarehouse()) &&
        !d_scheduler->isNewDW(req->mapDataWarehouse() + 1)) {
      continue;
    }

    DOUTR(g_topological_deps_dbg, "  req: " << *req);

    constHandle<PatchSubset> patches =
      req->getPatchesUnderDomain(dtask->m_patches);

    TypeDescription::Type vartype = req->var->typeDescription()->getType();

    // ARS - Treat sole vars the same as reduction vars ?
    // make sure newdw reduction variable requires link up to the reduction
    // tasks.
    if ((vartype == TypeDescription::Type::ReductionVariable ||
         vartype == TypeDescription::Type::SoleVariable) &&
        d_scheduler->isNewDW(req->mapDataWarehouse())) {
      patches = nullptr;
    }

    constHandle<MaterialSubset> matls =
      req->getMaterialsUnderDomain(dtask->m_matls);

    // this section is just to find the low and the high of the patch that will
    // use the other level's data.  Otherwise, we have to use the entire set of
    // patches (and ghost patches if applicable) that lay above/beneath this
    // patch.
    int levelID              = 0;
    const Patch* origPatch   = nullptr;
    const Level* origLevel   = nullptr;
    const bool uses_SHRT_MAX = (req->num_ghost_cells == SHRT_MAX);

    if ((dtask->m_patches) &&
        (dtask->getTask()->getType() != Task::OncePerProc) &&
        (dtask->getTask()->getType() != Task::Hypre)) {
      origPatch = dtask->m_patches->get(0);
      origLevel = origPatch->getLevel();
      levelID   = origLevel->getID();
    }
    int levelOffset =
      req->level_offset; // The level offset indicates how many levels up or
                         // down from the patches assigned to the task.
    int trueLevel = levelID;
    if (req->patches_dom == Task::CoarseLevel) {
      trueLevel -= levelOffset;
    } else if (req->patches_dom == Task::FineLevel) {
      trueLevel += levelOffset;
    }

    IntVector otherLevelLow, otherLevelHigh;
    if (req->patches_dom == Task::CoarseLevel ||
        req->patches_dom == Task::FineLevel) {
      // the requires should have been done with Task::CoarseLevel or FineLevel,
      // with null patches and the task->patches should be size one (so we don't
      // have to worry about overlapping regions)
      origPatch = dtask->m_patches->get(0);

      ASSERT(req->patches == nullptr);
      ASSERT(dtask->m_patches->size() == 1);
      ASSERT(req->level_offset > 0);

      if (req->patches_dom == Task::CoarseLevel) {
        // change the ghost cells to reflect coarse level
        LevelP nextLevel = origPatch->getLevelP();
        int levelOffset  = req->level_offset;
        IntVector ratio  = origPatch->getLevel()->getRefinementRatio();
        while (--levelOffset) {
          nextLevel = nextLevel->getCoarserLevel();
          ratio     = ratio * nextLevel->getRefinementRatio();
        }

        int ngc =
          req->num_ghost_cells * Max(Max(ratio.x(), ratio.y()), ratio.z());
        IntVector ghost(ngc, ngc, ngc);

        // manually set it, can't use computeVariableExtents since there might
        // not be a neighbor fine patch, and it would throw it off.
        otherLevelLow  = origPatch->getExtraCellLowIndex() - ghost;
        otherLevelHigh = origPatch->getExtraCellHighIndex() + ghost;

        otherLevelLow =
          origLevel->mapCellToCoarser(otherLevelLow, req->level_offset);
        otherLevelHigh =
          origLevel->mapCellToCoarser(otherLevelHigh, req->level_offset) +
          ratio - IntVector(1, 1, 1);
      } else {
        // This covers when req->m_patches_dom == Task::ThisLevel (single level
        // problems) or when req->m_patches_dom == Task::OtherGridDomain. (AMR
        // problems)
        if (uses_SHRT_MAX) {
          // Finer patches probably shouldn't be using SHRT_MAX ghost cells, but
          // just in case they do, at least compute the low and high
          // correctly...
          origLevel->computeVariableExtents(
            req->var->typeDescription()->getType(),
            otherLevelLow,
            otherLevelHigh);
        } else {
          origPatch->computeVariableExtentsWithBoundaryCheck(
            req->var->typeDescription()->getType(),
            req->var->getBoundaryLayer(),
            req->gtype,
            req->num_ghost_cells,
            otherLevelLow,
            otherLevelHigh);
        }
        otherLevelLow  = origLevel->mapCellToFiner(otherLevelLow);
        otherLevelHigh = origLevel->mapCellToFiner(otherLevelHigh);
      }
    }

    if (patches && !patches->empty() && matls && !matls->empty()) {

      // Skip reduction and sole vars as they are not patch based.
      if (vartype == TypeDescription::Type::ReductionVariable ||
          vartype == TypeDescription::Type::SoleVariable) {
        continue;
      }

      for (int i = 0; i < patches->size(); i++) {
        const Patch* patch = patches->get(i);

        // only allocate once
        static Patch::selectType neighbors;
        neighbors.resize(0);

        IntVector low  = IntVector(-9, -9, -9);
        IntVector high = IntVector(-9, -9, -9);

        Patch::VariableBasis basis = Patch::translateTypeToBasis(
          req->var->typeDescription()->getType(), false);

        if (uses_SHRT_MAX) {
          patch->getLevel()->computeVariableExtents(
            req->var->typeDescription()->getType(), low, high);
        } else {
          patch->computeVariableExtentsWithBoundaryCheck(
            req->var->typeDescription()->getType(),
            req->var->getBoundaryLayer(),
            req->gtype,
            req->num_ghost_cells,
            low,
            high);
        }

        if (req->patches_dom == Task::CoarseLevel ||
            req->patches_dom == Task::FineLevel) {
          // make sure the bounds of the dep are limited to the original patch's
          // (see above) also limit to current patch, as patches already loops
          // over all patches
          IntVector origlow = low, orighigh = high;
          if (req->patches_dom == Task::FineLevel) {
            // don't coarsen the extra cells
            low = patch->getExtraLowIndex(basis, req->var->getBoundaryLayer());
            high =
              patch->getExtraHighIndex(basis, req->var->getBoundaryLayer());
          } else {
            low  = Max(low, otherLevelLow);
            high = Min(high, otherLevelHigh);
          }

          if (high.x() <= low.x() || high.y() <= low.y() ||
              high.z() <= low.z()) {
            continue;
          }

          // don't need to selectPatches.  Just use the current patch, as we're
          // already looping over our required patches.
          neighbors.push_back(patch);
        } else {
          origPatch = patch;
          if (req->num_ghost_cells > 0) {
            patch->getLevel()->selectPatches(low, high, neighbors);
          } else {
            neighbors.push_back(patch);
          }
        }

        ASSERT(
          std::is_sorted(neighbors.begin(), neighbors.end(), Patch::Compare()));

        size_t num_neighbors = neighbors.size();
        DOUTR(g_topological_deps_dbg,
              "    Creating detailed dependency on "
                << num_neighbors << " neighboring patch"
                << (num_neighbors > 1 ? "es " : "   ") << neighbors
                << "   Low=" << low << ", high=" << high << ", dw= "
                << req->mapDataWarehouse() << ", var=" << req->var->getName())

        // for all neighbors - find and store from neighbors
        for (auto neighbor : neighbors) {
          // if neighbor is not in my neighborhood just continue as its
          // dependencies are not important to this processor
          DOUTR(g_proc_neighborhood_dbg,
                "    In detailed task: "
                  << dtask->getName() << " checking if " << *req
                  << " is in neighborhood on level: " << trueLevel);

          // NOTE THE MAP::AT METHOD NEEDS TO BE SAFEGUARDED. Is it
          // reasonable to set search_distal_requires = false when it
          // fails??

          // The map::at will fail after an AMR regridding when using
          // more than 24 ranks (it works with 23 ranks). The question
          // is why is the variable 'trueLevel' not valid????

          bool search_distal_reqs;
          if ((dtask->getTask()->m_max_ghost_cells.find(trueLevel) !=
               dtask->getTask()->m_max_ghost_cells.end())) {
            search_distal_reqs = (dtask->getTask()->m_max_ghost_cells.at(
                                    trueLevel) >= MAX_HALO_DEPTH);
          } else {
            search_distal_reqs = false;

            DOUTALL(true,
                    "*********** Bad true level " << trueLevel << " levelID "
                                                  << levelID << " levelOffset "
                                                  << levelOffset);
          }

          if (!d_load_balancer->inNeighborhood(neighbor->getRealPatch(),
                                               search_distal_reqs)) {
            DOUTR(g_proc_neighborhood_dbg, "    No");
            continue;
          }
          DOUTR(g_proc_neighborhood_dbg, "    Yes");

          static Patch::selectType fromNeighbors;
          fromNeighbors.resize(0);

          IntVector l =
            Max(neighbor->getExtraLowIndex(basis, req->var->getBoundaryLayer()),
                low);
          IntVector h = Min(
            neighbor->getExtraHighIndex(basis, req->var->getBoundaryLayer()),
            high);

          if (neighbor->isVirtual()) {
            l -= neighbor->getVirtualOffset();
            h -= neighbor->getVirtualOffset();
            neighbor = neighbor->getRealPatch();
          }

          if (req->patches_dom == Task::OtherGridDomain) {
            // this is when we are copying data between two grids (currently
            // between timesteps) the grid assigned to the old dw should be the
            // old grid. This should really only impact things required from the
            // OldDW.
            LevelP fromLevel = d_scheduler->get_dw(0)->getGrid()->getLevel(
              patch->getLevel()->getIndex());
            fromLevel->selectPatches(Max(neighbor->getExtraLowIndex(
                                           basis, req->var->getBoundaryLayer()),
                                         l),
                                     Min(neighbor->getExtraHighIndex(
                                           basis, req->var->getBoundaryLayer()),
                                         h),
                                     fromNeighbors);
          } else {
            fromNeighbors.push_back(neighbor);
          }

          // For all from neighbors
          for (auto fromNeighbor : fromNeighbors) {
            // only add the requirements both fromNeighbor is in my neighborhood

            // NOTE THE MAP::AT METHOD NEEDS TO BE SAFEGUARDED. Is it
            // reasonable to set search_distal_requires = false when it
            // fails??

            // The map::at will fail after an AMR regridding when using
            // more than 24 ranks (it works with 23 ranks). The question
            // is why is the variable 'trueLevel' not valid????
            bool search_distal_requires;
            if ((dtask->getTask()->m_max_ghost_cells.find(trueLevel) !=
                 dtask->getTask()->m_max_ghost_cells.end())) {
              search_distal_requires = (dtask->getTask()->m_max_ghost_cells.at(
                                          trueLevel) >= MAX_HALO_DEPTH);
            } else {
              search_distal_requires = false;

              DOUTALL(true,
                      "*********** Bad true level "
                        << trueLevel << " levelID " << levelID
                        << " levelOffset " << levelOffset);
            }

            if (!d_load_balancer->inNeighborhood(fromNeighbor,
                                                 search_distal_requires)) {
              continue;
            }

            IntVector from_l;
            IntVector from_h;

            if (req->patches_dom == Task::OtherGridDomain &&
                fromNeighbor->getLevel()->getIndex() > 0) {
              // DON'T send extra cells (unless they're on the domain boundary)
              from_l = Max(fromNeighbor->getLowIndexWithDomainLayer(basis), l);
              from_h = Min(fromNeighbor->getHighIndexWithDomainLayer(basis), h);
            } else {
              // TODO - APH This intersection should not be needed, but let's
              // clean this up if not
              // from_l = Max(fromNeighbor->getExtraLowIndex(basis,
              // req->var->getBoundaryLayer()), l); from_h =
              // Min(fromNeighbor->getExtraHighIndex(basis,
              // req->var->getBoundaryLayer()), h);
              from_l = l;
              from_h = h;
              // verify in debug mode that the intersection is unneeded
              ASSERT(Max(fromNeighbor->getExtraLowIndex(
                           basis, req->var->getBoundaryLayer()),
                         l) == l);
              ASSERT(Min(fromNeighbor->getExtraHighIndex(
                           basis, req->var->getBoundaryLayer()),
                         h) == h);
            }

            // For all materials
            for (int m = 0; m < matls->size(); m++) {
              int matl = matls->get(m);

              // creator is the task that performs the original compute.
              // If the require is for the OldDW, then it will be a send old
              // data task
              DetailedTask* creator  = nullptr;
              Task::Dependency* comp = nullptr;

              // look in old dw or in old TG.  Legal to modify across TG
              // boundaries
              int proc = -1;
              if (d_scheduler->isOldDW(req->mapDataWarehouse())) {
                ASSERT(!modifies);
                proc    = findVariableLocation(req, fromNeighbor, matl, 0);
                creator = d_detailed_tasks->getOldDWSendTask(proc);
                comp    = nullptr;
              } else {
                if (!ct.findcomp(
                      req, neighbor, matl, creator, comp, d_proc_group)) {
                  if (d_type == Scheduler::IntermediateTaskGraph &&
                      req->look_in_old_tg) {

                    // same stuff as above - but do the check for findcomp
                    // first, as this is a "if you don't find it here, assign it
                    // from the old TG" dependency
                    proc    = findVariableLocation(req, fromNeighbor, matl, 0);
                    creator = d_detailed_tasks->getOldDWSendTask(proc);
                    comp    = nullptr;

                  } else {

                    // if neither the patch or the neighbor are on this
                    // processor then the computing task doesn't exist so just
                    // continue
                    if (d_load_balancer->getPatchwiseProcessorAssignment(
                          patch) != d_proc_group->myRank() &&
                        d_load_balancer->getPatchwiseProcessorAssignment(
                          neighbor) != d_proc_group->myRank()) {
                      continue;
                    }

                    DOUT(true, "*******************************");
                    DOUT(true,
                         "  Rank-" << my_rank
                                   << ", Task Graph Index: " << d_index);
                    DOUT(true,
                         "  Failure finding " << *req << " for " << *dtask);
                    if (creator) {
                      std::cout << "  creator=" << *creator << "\n";
                    }
                    DOUT(true,
                         "  neighbor=" << *fromNeighbor << ", matl=" << matl);
                    DOUT(true, "*******************************")
                    SCI_THROW(InternalError(
                      "Failed to find comp for dep!", __FILE__, __LINE__));
                  }
                }
              }

              if (modifies && comp) { // comp means NOT send-old-data tasks

                // find the tasks that up to this point require the variable
                // that we are modifying (i.e., the ones that use the computed
                // variable before we modify it), and put a dependency between
                // those tasks and this tasks
                // i.e., the task that requires data computed by a task on this
                // processor needs to finish its task before this task, which
                // modifies the data computed by the same task
                std::list<DetailedTask*> requireBeforeModifiedTasks;
                creator->findRequiringTasks(req->var,
                                            requireBeforeModifiedTasks);

                for (auto& prevReqTask : requireBeforeModifiedTasks) {
                  if (prevReqTask == dtask) {
                    continue;
                  }
                  if (prevReqTask->m_task == dtask->m_task) {
                    if (!dtask->m_task->getHasSubScheduler()) {

                      std::ostringstream message;
                      message << " WARNING - task (" << dtask->getName()
                              << ") requires with Ghost cells *and* modifies "
                                 "and may not be correct"
                              << std::endl;
                      static ProgressiveWarning warn(message.str(), 10);
                      warn.invoke();

                      DOUTR(g_topological_deps_dbg,
                            " Task that requires with ghost cells and modifies "
                            " RGM: var: "
                              << *req->var << " compute: " << *creator
                              << " mod " << *dtask << " PRT " << *prevReqTask
                              << " " << from_l << " " << from_h);
                    }
                  } else {
                    // dep requires what is to be modified before it is to be
                    // modified so create a dependency between them so the
                    // modifying won't conflict with the previous require.
                    DOUTR(g_topological_deps_dbg,
                          "       Requires to modifies dependency from "
                            << prevReqTask->getName() << " to "
                            << dtask->getName() << " (created by "
                            << creator->getName() << ")");

                    if (creator->getPatches() &&
                        creator->getPatches()->size() > 1) {
                      // if the creator works on many patches, then don't create
                      // links between patches that don't touch
                      const PatchSubset* psub    = dtask->getPatches();
                      const PatchSubset* req_sub = prevReqTask->getPatches();
                      if (psub->size() == 1 && req_sub->size() == 1) {
                        const Patch* p         = psub->get(0);
                        const Patch* req_patch = req_sub->get(0);
                        Patch::selectType n;
                        IntVector low, high;

                        req_patch->computeVariableExtents(
                          req->var->typeDescription()->getType(),
                          req->var->getBoundaryLayer(),
                          Ghost::AroundCells,
                          2,
                          low,
                          high);

                        req_patch->getLevel()->selectPatches(low, high, n);
                        bool found = false;
                        for (auto& i : n) {
                          if (i->getID() == p->getID()) {
                            found = true;
                            break;
                          }
                        }
                        if (!found) {
                          continue;
                        }
                      }
                    }
                    d_detailed_tasks->possiblyCreateDependency(
                      prevReqTask,
                      nullptr,
                      nullptr,
                      dtask,
                      req,
                      nullptr,
                      matl,
                      from_l,
                      from_h,
                      DetailedDep::Always);
                  }
                }
              }

              DetailedDep::CommCondition cond = DetailedDep::Always;
              if (proc != -1 && req->patches_dom != Task::OtherGridDomain) {
                // for OldDW tasks - see comment in class DetailedDep by
                // CommCondition
                int subsequentProc =
                  findVariableLocation(req, fromNeighbor, matl, 1);
                if (subsequentProc != proc) {
                  cond = DetailedDep::FirstIteration; // change outer cond from
                                                      // always to first-only
                  DetailedTask* subsequentCreator =
                    d_detailed_tasks->getOldDWSendTask(subsequentProc);
                  d_detailed_tasks->possiblyCreateDependency(
                    subsequentCreator,
                    comp,
                    fromNeighbor,
                    dtask,
                    req,
                    fromNeighbor,
                    matl,
                    from_l,
                    from_h,
                    DetailedDep::SubsequentIterations);

                  DOUTR(g_topological_deps_dbg,
                        "   Adding condition reqs for "
                          << *req->var << " task : " << *creator << "  to "
                          << *dtask);
                }
              }

              d_detailed_tasks->possiblyCreateDependency(creator,
                                                         comp,
                                                         fromNeighbor,
                                                         dtask,
                                                         req,
                                                         fromNeighbor,
                                                         matl,
                                                         from_l,
                                                         from_h,
                                                         cond);
            }
          }
        }
      }
    } else if (!patches && matls && !matls->empty()) {

      // requiring reduction variables
      for (int m = 0; m < matls->size(); m++) {
        int matl = matls->get(m);
        static std::vector<DetailedTask*> creators;
        creators.resize(0);

        ct.findReductionComps(req, nullptr, matl, creators, d_proc_group);
        // if the size is 0, that's fine.  It means that there are more procs
        // than patches on this level, so the reducer will pick a benign value
        // that won't affect the reduction

        ASSERTRANGE(
          dtask->getAssignedResourceIndex(), 0, d_proc_group->nRanks());

        for (auto creator : creators) {
          if (dtask->getAssignedResourceIndex() ==
                creator->getAssignedResourceIndex() &&
              dtask->getAssignedResourceIndex() == my_rank) {
            dtask->addInternalDependency(creator, req->var);

            DOUTR(g_topological_deps_dbg,
                  "    Created reduction dependency between "
                    << *dtask << " and " << *creator);
          }
        }
      }
    } else if (patches && (patches->empty() || patches->size() <= 1) &&
               (req->patches_dom == Task::FineLevel ||
                dtask->getTask()->getType() == Task::OncePerProc ||
                dtask->getTask()->getType() == Task::Hypre ||
                dtask->getTask()->getType() == Task::Output ||
                dtask->getTask()->getName() ==
                  "SchedulerCommon::copyDataToNewGrid")) {

      // this is a either coarsen task where there aren't any fine patches, or a
      // PerProcessor task where there aren't any patches on this processor.
      // Perfectly legal, so do nothing

      // another case is the copy-data-to-new-grid task, which will wither
      // compute or modify to every patch but not both.  So it will yell at you
      // for the detailed task's patches not intersecting with the computes or
      // modifies... (maybe there's a better way) - bryan

    } else if (dtask->m_matls && req->matls && dtask->m_patches &&
               req->patches && patches.get_rep()->size() == 0) {
      // Fields were required on a subset of the domain with ghosts - this
      // should be legal.

    } else {

      std::ostringstream desc;
      desc << "TaskGraph::createDetailedDependencies, task dependency not "
              "supported without patches and materials"
           << " \n Trying to require or modify " << *req << " in Task "
           << dtask->getTask()->getName() << "\n\n";

      if (dtask->m_matls) {
        desc << "task materials:" << *dtask->m_matls << "\n";
      } else {
        desc << "no task materials\n";
      }
      if (req->matls) {
        desc << "req materials: " << *req->matls << "\n";
      } else {
        desc << "no req materials\n";
        desc << "domain materials: " << *matls.get_rep() << "\n";
      }
      if (dtask->m_patches) {
        desc << "task patches:" << *dtask->m_patches << "\n";
      } else {
        desc << "no task patches\n";
      }
      if (req->patches) {
        desc << "req patches: " << *req->patches << "\n";
      } else {
        desc << "no req patches\n";
      }

      desc << "domain patches: " << *patches.get_rep() << "\n";
      SCI_THROW(InternalError(desc.str(), __FILE__, __LINE__));
    }
  }
}

auto
TaskGraph::findVariableLocation(Task::Dependency* req,
                                const Patch* patch,
                                [[maybe_unused]] int matl,
                                int iteration) -> int
{
  // This needs to be improved, especially for re-distribution on
  // restart from checkpoint.
  int proc;
  if ((req->task->mapDataWarehouse(Task::ParentNewDW) != -1 &&
       req->whichdw != Task::ParentOldDW) ||
      iteration > 0 ||
      (req->look_in_old_tg && d_type == Scheduler::IntermediateTaskGraph)) {
    // Provide some accommodation for Dynamic load balancers and sub schedulers.
    // We need to treat the requirement like a "old" dw req but it needs to be
    // found on the current processor. Same goes for successive executions of
    // the same TG.
    proc = d_load_balancer->getPatchwiseProcessorAssignment(patch);
  } else {
    proc = d_load_balancer->getOldProcessorAssignment(patch);
  }
  return proc;
}

auto
TaskGraph::getNumTasks() const -> int
{
  return static_cast<int>(d_tasks.size());
}

//______________________________________________________________________
//
auto
TaskGraph::getTask(int idx) -> Task*
{
  return d_tasks[idx].get();
}

//______________________________________________________________________
//
void
TaskGraph::makeVarLabelMaterialMap(Scheduler::VarLabelMaterialMap* result)
{
  for (auto& task : d_tasks) {
    for (Task::Dependency* comp = task->getComputes(); comp != nullptr;
         comp                   = comp->next) {
      // assume all patches will compute the same labels on the same materials
      const VarLabel* label         = comp->var;
      std::list<int>& matls         = (*result)[label->getName()];
      const MaterialSubset* msubset = comp->matls;
      if (msubset) {
        for (int mm = 0; mm < msubset->size(); mm++) {
          matls.push_back(msubset->get(mm));
        }
      } else if (label->typeDescription()->getType() ==
                 TypeDescription::Type::ReductionVariable) {
        // Default to material -1 (global)
        matls.push_back(-1);
      } else {
        const MaterialSet* ms = task->getMaterialSet();
        for (int m = 0; m < ms->size(); m++) {
          const MaterialSubset* msubset = ms->getSubset(m);
          for (int mm = 0; mm < msubset->size(); mm++) {
            matls.push_back(msubset->get(mm));
          }
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////
// Archived code for topological sort
///////////////////////////////////////////////////////////////////////

//______________________________________________________________________
//
auto
TaskGraph::overlaps(const Task::Dependency* comp,
                    const Task::Dependency* req) const -> bool
{
  constHandle<PatchSubset> saveHandle2;
  const PatchSubset* ps1 = comp->patches;
  if (!ps1) {
    if (!comp->task->getPatchSet()) {
      return false;
    }
    ps1 = comp->task->getPatchSet()->getUnion();
    if (comp->patches_dom == Task::CoarseLevel ||
        comp->patches_dom == Task::FineLevel) {
      SCI_THROW(InternalError(
        "Should not compute onto another level!", __FILE__, __LINE__));
      // This may not be a big deal if it were needed, but I didn't
      // think that it should be allowed - Steve
      // saveHandle1 = comp->getPatchesUnderDomain(ps1);
      // ps1 = saveHandle1.get_rep();
    }
  }

  const PatchSubset* ps2 = req->patches;
  if (!ps2) {
    if (!req->task->getPatchSet()) {
      return false;
    }
    ps2 = req->task->getPatchSet()->getUnion();
    if (req->patches_dom == Task::CoarseLevel ||
        req->patches_dom == Task::FineLevel) {
      saveHandle2 = req->getPatchesUnderDomain(ps2);
      ps2         = saveHandle2.get_rep();
    }
  }

  if (!PatchSubset::overlaps(
        ps1,
        ps2)) { // && !(ps1->size() == 0 && (!req->patches || ps2->size() == 0)
                // && comp->task->getType() == Task::OncePerProc))
    return false;
  }

  const MaterialSubset* ms1 = comp->matls;
  if (!ms1) {
    if (!comp->task->getMaterialSet()) {
      return false;
    }
    ms1 = comp->task->getMaterialSet()->getUnion();
  }
  const MaterialSubset* ms2 = req->matls;
  if (!ms2) {
    if (!req->task->getMaterialSet()) {
      return false;
    }
    ms2 = req->task->getMaterialSet()->getUnion();
  }
  if (!MaterialSubset::overlaps(ms1, ms2)) {
    return false;
  }
  return true;
}

//______________________________________________________________________
//
// setupTaskConnections also adds Reduction Tasks to the graph...
void
TaskGraph::setupTaskConnections(GraphSortInfoMap& sortinfo)
{
  std::vector<Task*>::iterator iter;
  // Initialize variables on the tasks
  for (auto& task : d_tasks) {
    sortinfo[task.get()] = GraphSortInfo();
  }

  if (m_edges.size() > 0) {
    return; // already been done
  }

  // Look for all of the reduction variables - we must treat those
  // special.  Create a fake task that performs the reduction
  // While we are at it, ensure that we aren't producing anything
  // into an "old" data warehouse
  LevelReductionTasksMap reductionTasks;
  for (auto& task : d_tasks) {
    if (task->isReductionTask()) {
      continue; // already a reduction task so skip it
    }

    for (Task::Dependency* comp = task->getComputes(); comp != nullptr;
         comp                   = comp->next) {
      if (d_scheduler->isOldDW(comp->mapDataWarehouse())) {
        if (g_topological_deps_dbg) {
          DOUT(g_topological_deps_dbg,
               d_proc_group->myRank()
                 << " which = " << comp->whichdw << ", mapped to "
                 << comp->mapDataWarehouse());
        }
        SCI_THROW(InternalError("Variable produced in old datawarehouse: " +
                                  comp->var->getName(),
                                __FILE__,
                                __LINE__));
      } else if (comp->var->typeDescription()->isReductionVariable()) {
        int levelidx =
          comp->reduction_level ? comp->reduction_level->getIndex() : -1;
        // Look up this variable in the reductionTasks map
        int dw = comp->mapDataWarehouse();
        // for reduction var allows multi computes such as delT
        // do not generate reduction task each time it computes,
        // instead computes it in a system wide reduction task
        if (!comp->var->isReductionTask()) {
          if (g_topological_deps_dbg) {
            DOUT(g_topological_deps_dbg,
                 d_proc_group->myRank()
                   << " Skipping Reduction task for variable: "
                   << comp->var->getName() << " on level " << levelidx
                   << ", DW " << dw);
          }
          continue;
        }
        ASSERT(comp->patches == nullptr);

        // use the dw as a 'material', just for the sake of looking it up.
        // it should only differentiate on AMR W-cycle graphs...
        VarLabelMatl<Level> key(comp->var, dw, comp->reduction_level);
        const MaterialSet* ms = task->getMaterialSet();
        const Level* level    = comp->reduction_level;

        auto it = reductionTasks.find(key);
        if (it == reductionTasks.end()) {
          // No reduction task yet, create one
          if (g_topological_deps_dbg) {
            DOUT(g_topological_deps_dbg,
                 d_proc_group->myRank()
                   << " creating Reduction task for variable: "
                   << comp->var->getName() << " on level " << levelidx
                   << ", DW " << dw);
          }
          std::ostringstream taskname;
          taskname << "Reduction: " << comp->var->getName()
                   << ", level: " << levelidx << ", dw: " << dw;
          Task* newtask = scinew Task(taskname.str(), Task::Reduction);

          sortinfo[newtask] = GraphSortInfo();

          int dwmap[Task::TotalDWs];
          for (int& i : dwmap) {
            i = Task::InvalidDW;
          }
          dwmap[Task::OldDW] = Task::NoDW;
          dwmap[Task::NewDW] = dw;
          newtask->setMapping(dwmap);

          // compute and require for all patches but some set of materials
          // (maybe global material, but not necessarily)
          if (comp->matls != nullptr) {
            // TODO APH - figure this out and clean up - 01/31/15
            // newtask->computes(comp->var, level, comp->matls,
            // Task::OutOfDomain); newtask->requires(Task::NewDW, comp->var,
            // level, comp->matls, Task::OutOfDomain);
            newtask->modifies(comp->var, level, comp->matls, Task::OutOfDomain);
          } else {
            for (int m = 0; m < ms->size(); m++) {
              // TODO APH - figure this out and clean up - 01/31/15
              // newtask->computes(comp->var, level, ms->getSubset(m),
              // Task::OutOfDomain); newtask->requires(Task::NewDW, comp->var,
              // level, ms->getSubset(m), Task::OutOfDomain);
              newtask->modifies(
                comp->var, level, ms->getSubset(m), Task::OutOfDomain);
            }
          }
          reductionTasks[key] = newtask;
          it                  = reductionTasks.find(key);
        }
      }
    }
  }

  // Add the new reduction tasks to the list of tasks
  for (auto& [var_label_matl, task] : reductionTasks) {
    std::shared_ptr<Task> task_sp(task);
    addTask(task_sp, nullptr, nullptr);
  }

  // Gather the comps for the tasks into a map
  CompMap comps;
  for (auto& task : d_tasks) {
    if (g_topological_deps_dbg) {
      DOUT(g_topological_deps_dbg,
           d_proc_group->myRank() << " Gathering comps from task: " << *task);
    }
    for (Task::Dependency* comp = task->getComputes(); comp != nullptr;
         comp                   = comp->next) {
      comps.insert(std::make_pair(comp->var, comp));
      if (g_topological_deps_dbg) {
        DOUT(g_topological_deps_dbg,
             d_proc_group->myRank() << "   Added comp for: " << *comp);
      }
    }
  }

  // Connect the tasks where the requires/modifies match a comp.
  // Also, updates the comp map with each modify and doing this in task order
  // so future modifies/requires find the modified var.  Also do a type check
  for (auto& task : d_tasks) {
    if (g_topological_deps_dbg) {
      DOUT(g_topological_deps_dbg,
           d_proc_group->myRank()
             << "   Looking at dependencies for task: " << *task);
    }
    addDependencyEdges(
      task.get(), sortinfo, task->getRequires(), comps, reductionTasks, false);
    addDependencyEdges(
      task.get(), sortinfo, task->getModifies(), comps, reductionTasks, true);
    // Used here just to warn if a modifies comes before its computes
    // in the order that tasks were added to the graph.
    sortinfo.find(task.get())->second.m_visited = true;
    task->m_all_child_tasks.clear();
    if (g_topological_deps_dbg) {
      std::cout << d_proc_group->myRank()
                << "   Looking at dependencies for task: " << *task
                << "child task num=" << task->m_child_tasks.size() << "\n";
    }
  }

  // count the all child tasks
  int nd_task = d_tasks.size();
  while (nd_task > 0) {
    nd_task = d_tasks.size();
    for (auto& task : d_tasks) {
      if (task->m_all_child_tasks.size() == 0) {
        if (task->m_child_tasks.size() ==
            0) { // leaf task, add itself to the set
          task->m_all_child_tasks.insert(task.get());
          break;
        }
        for (auto& child_task : task->m_child_tasks) {
          if (child_task->m_all_child_tasks.size() > 0) {
            task->m_all_child_tasks.insert(
              child_task->m_all_child_tasks.begin(),
              child_task->m_all_child_tasks.end());
            task->m_all_child_tasks.insert(task.get());
          } else {
            // if child didn't finish computing m_all_child_tasks
            task->m_all_child_tasks.clear();
            break;
          }
        }
      } else {
        nd_task--;
      }
    }
  }

  // Initialize variables on the tasks
  GraphSortInfoMap::iterator sort_iter;
  for (sort_iter = sortinfo.begin(); sort_iter != sortinfo.end(); sort_iter++) {
    sort_iter->second.m_visited = false;
    sort_iter->second.m_sorted  = false;
  }
} // end setupTaskConnections()

//______________________________________________________________________
//
void
TaskGraph::addDependencyEdges(Task* task,
                              GraphSortInfoMap& sortinfo,
                              Task::Dependency* req,
                              CompMap& comps,
                              LevelReductionTasksMap& reductionTasks,
                              bool modifies)
{
  for (; req != nullptr; req = req->next) {
    if (g_topological_deps_dbg) {
      DOUT(g_topological_deps_dbg,
           d_proc_group->myRank()
             << "     Checking edge for req: " << *req
             << ", task: " << *req->task << ", domain: " << req->patches_dom);
    }
    if (req->whichdw == Task::NewDW) {
      // If DW is finalized, we assume that we already have it,
      // or that we will get it sent to us.  Otherwise, we set
      // up an edge to connect this req to a comp

      std::pair<CompMap::iterator, CompMap::iterator> iters =
        comps.equal_range(static_cast<const Uintah::VarLabel*>(req->var));
      int count = 0;
      for (auto compiter = iters.first; compiter != iters.second; ++compiter) {

        if (req->var->typeDescription() != compiter->first->typeDescription()) {
          SCI_THROW(TypeMismatchException("Type mismatch for variable: " +
                                            req->var->getName(),
                                          __FILE__,
                                          __LINE__));
        }

        // determine if we need to add a dependency edge
        bool add                   = false;
        bool requiresReductionTask = false;
        if (g_topological_deps_dbg) {
          DOUT(g_topological_deps_dbg,
               d_proc_group->myRank()
                 << "  Checking edge from comp: " << *compiter->second
                 << ", task: " << *compiter->second->task
                 << ", domain: " << compiter->second->patches_dom);
        }
        if (req->mapDataWarehouse() == compiter->second->mapDataWarehouse()) {
          if (req->var->typeDescription()->isReductionVariable()) {
            // Match the level first
            if (compiter->second->reduction_level == req->reduction_level) {
              add = true;
            }
            // with reduction variables, you can modify them up to the Reduction
            // Task, which also modifies those who don't modify will get the
            // reduced value.
            if (!modifies && req->var->isReductionTask()) {
              requiresReductionTask = true;
            }
          } else if (overlaps(compiter->second, req)) {
            add = true;
          }
        }

        if (!add) {
          if (g_topological_deps_dbg) {
            DOUT(g_topological_deps_dbg,
                 d_proc_group->myRank() << "       did NOT create dependency");
          }
        } else {
          Task::Dependency* comp;
          if (requiresReductionTask) {
            VarLabelMatl<Level> key(
              req->var, req->mapDataWarehouse(), req->reduction_level);
            Task* redTask = reductionTasks[key];
            ASSERT(redTask != nullptr);
            // reduction tasks should have exactly 1 require, and it should be a
            // modify assign the requiring task's require to it
            comp = redTask->getModifies();
            ASSERT(comp != nullptr);
            DOUT(g_topological_deps_dbg,
                 "  Using Reduction task: " << *redTask);
          } else {
            comp = compiter->second;
          }

          if (modifies) {
            // Add dependency edges to each task that requires the data
            // before it is modified.
            for (auto otherEdge = comp->m_req_head; otherEdge != nullptr;
                 otherEdge      = otherEdge->req_next) {
              auto* priorReq = const_cast<Task::Dependency*>(otherEdge->req);
              if (priorReq != req) {
                ASSERT(priorReq->var->equals(req->var));
                if (priorReq->task != task) {
                  auto* edge = scinew Task::Edge(priorReq, req);
                  m_edges.push_back(edge);
                  req->addComp(edge);
                  priorReq->addReq(edge);
                  if (g_topological_deps_dbg) {
                    DOUT(g_topological_deps_dbg,
                         d_proc_group->myRank()
                           << " Creating edge from task: " << *priorReq->task
                           << " to task: " << *req->task);
                    DOUT(g_topological_deps_dbg,
                         d_proc_group->myRank() << " Prior Req=" << *priorReq);
                    DOUT(g_topological_deps_dbg,
                         d_proc_group->myRank() << " Modify=" << *req);
                  }
                }
              }
            }
          }

          // add the edge between the require/modify and compute
          auto* edge = scinew Task::Edge(comp, req);
          m_edges.push_back(edge);
          req->addComp(edge);
          comp->addReq(edge);

          if (!sortinfo.find(edge->comp->task)->second.m_visited &&
              !edge->comp->task->isReductionTask()) {
            std::cout << "\nWARNING: The task, '" << task->getName()
                      << "', that ";
            if (modifies) {
              std::cout << "modifies '";
            } else {
              std::cout << "requires '";
            }
            std::cout << req->var->getName()
                      << "' was added before the computing task";
            std::cout << ", '" << edge->comp->task->getName() << "'\n";
            std::cout << "  Required/modified by: " << *task << "\n";
            std::cout << "  req: " << *req << "\n";
            std::cout << "  Computed by: " << *edge->comp->task << "\n";
            std::cout << "  comp: " << *comp << "\n";
            std::cout << "\n";
          }
          count++;
          task->m_child_tasks.insert(comp->task);
          if (g_detailed_task_dbg) {
            DOUT(g_detailed_task_dbg,
                 d_proc_group->myRank()
                   << "       Creating edge from task: " << *comp->task
                   << " to task: " << *task);
            DOUT(g_detailed_task_dbg,
                 d_proc_group->myRank() << "         Req=" << *req);
            DOUT(g_detailed_task_dbg,
                 d_proc_group->myRank() << "         Comp=" << *comp);
          }
        }
      }

      // if we cannot find the required variable, throw an exception
      if (count == 0 && (!req->matls || req->matls->size() > 0) &&
          (!req->patches || req->patches->size() > 0) &&
          !(req->look_in_old_tg &&
            d_type == Scheduler::IntermediateTaskGraph)) {
        // if this is an Intermediate TG and the requested data is done from
        // another TG, we need to look in this TG first, but don't worry if you
        // don't find it

        std::cout << "ERROR: Cannot find the task that computes the variable ("
                  << req->var->getName() << ")\n";

        std::cout << "The task (" << task->getName()
                  << ") is requesting data from:\n";
        std::cout << "  Level:           "
                  << getLevel(task->getPatchSet())->getIndex() << "\n";
        std::cout << "  Task:PatchSet    " << *(task->getPatchSet()) << "\n";
        std::cout << "  Task:MaterialSet " << *(task->getMaterialSet())
                  << "\n \n";

        std::cout << "The variable (" << req->var->getName()
                  << ") is requiring data from:\n";

        if (req->patches) {
          std::cout << "  Level: " << getLevel(req->patches)->getIndex()
                    << "\n";
          std::cout << "  Patches': " << *(req->patches) << "\n";
        } else {
          std::cout << "  Patches:  All \n";
        }

        if (req->matls) {
          std::cout << "  Materials: " << *(req->matls) << "\n";
        } else {
          std::cout << "  Materials:  All \n";
        }

        std::cout << "\nTask Details:\n";
        task->display(std::cout);
        std::cout << "\nRequirement Details:\n" << *req << "\n";

        SCI_THROW(InternalError(
          "Scheduler could not find  production for variable: " +
            req->var->getName() + ", required for task: " + task->getName(),
          __FILE__,
          __LINE__));
      }

      if (modifies) {
        // not just requires, but modifies, so the comps map must be
        // updated so future modifies or requires will link to this one.
        comps.insert(std::make_pair(req->var, req));
        if (g_topological_deps_dbg) {
          DOUT(g_topological_deps_dbg,
               d_proc_group->myRank() << " Added modified comp for: " << *req);
        }
      }
    }
  }
}

//______________________________________________________________________
//
void
TaskGraph::processTask(Task* task,
                       std::vector<Task*>& sortedTasks,
                       GraphSortInfoMap& sortinfo) const
{
  if (g_topological_deps_dbg) {
    DOUT(g_topological_deps_dbg,
         d_proc_group->myRank() << " Looking at task: " << task->getName());
  }

  GraphSortInfo& gsi = sortinfo.find(task)->second;
  // we throw an exception before calling processTask if this task has already
  // been visited
  gsi.m_visited = true;

  processDependencies(task, task->getRequires(), sortedTasks, sortinfo);
  processDependencies(task, task->getModifies(), sortedTasks, sortinfo);

  // All prerequisites are done - add this task to the list
  sortedTasks.push_back(task);
  gsi.m_sorted = true;

  if (g_topological_deps_dbg) {
    DOUT(g_topological_deps_dbg,
         d_proc_group->myRank() << " Sorted task: " << task->getName());
  }
} // end processTask()

//______________________________________________________________________
//

void
TaskGraph::processDependencies(Task* task,
                               Task::Dependency* req,
                               std::vector<Task*>& sortedTasks,
                               GraphSortInfoMap& sortinfo) const

{
  for (; req != nullptr; req = req->next) {
    if (g_topological_deps_dbg) {
      DOUT(g_topological_deps_dbg,
           d_proc_group->myRank() << " processDependencies for req: " << *req);
    }
    if (req->whichdw == Task::NewDW) {
      Task::Edge* edge = req->m_comp_head;
      for (; edge != nullptr; edge = edge->comp_next) {
        Task* vtask        = edge->comp->task;
        GraphSortInfo& gsi = sortinfo.find(vtask)->second;
        if (!gsi.m_sorted) {
          try {
            // this try-catch mechanism will serve to print out the entire TG
            // cycle
            if (gsi.m_visited) {
              std::cout << d_proc_group->myRank()
                        << " Cycle detected in task graph\n";
              SCI_THROW(InternalError(
                "Cycle detected in task graph", __FILE__, __LINE__));
            }

            // recursively process the dependencies of the computing task
            processTask(vtask, sortedTasks, sortinfo);
          } catch (InternalError& e) {
            std::cout << d_proc_group->myRank() << "   Task " << task->getName()
                      << " requires " << req->var->getName() << " from "
                      << vtask->getName() << "\n";

            throw;
          }
        }
      }
    }
  }
}

void
TaskGraph::topologicalSort(std::vector<Task*>& sortedTasks)
{
  GraphSortInfoMap sortinfo;

  setupTaskConnections(sortinfo);

  for (auto& task : d_tasks) {
    if (!sortinfo.find(task.get())->second.m_sorted) {
      processTask(task.get(), sortedTasks, sortinfo);
    }
  }
  int n = 0;
  for (auto& task : sortedTasks) {
    task->setSortedOrder(n++);
  }
}

} // namespace Uintah