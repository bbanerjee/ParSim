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

#ifndef CCA_COMPONENTS_SCHEDULERS_TASKGRAPH_H
#define CCA_COMPONENTS_SCHEDULERS_TASKGRAPH_H

#include <CCA/Ports/Scheduler.h>

#include <Core/Containers/FastHashTable.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Task.h>

#include <list>
#include <map>
#include <memory>
#include <vector>

namespace Uintah {

class DetailedTask;
class DetailedTasks;
class Patch;
class LoadBalancer;

class CompTable;
class SchedulerCommon;

/**************************************
CLASS
   TaskGraph

   During the TaskGraph compilation, the task graph does its work in
   the createDetailedTasks function.  The first portion is to sort the tasks,
   by adding edges between computing tasks and requiring tasks,
   and the second is to create detailed tasks and dependencies

   Here is a function call tree for this phase:
   createDetailedTasks
     nullsort (formerly was topologicalSort)
       setupTaskConnections
         addDependencyEdges
       processTask
         processDependencies
     LoadBalancer::createNeighborhood (this stores the patches that border
                                       patches on the current processor)

   Detailed Task portion: divides up tasks into smaller pieces, and sets up the
   data that need to be communicated between processors when the taskgraph
executes.

     createDetailedTask (for each task, patch subset, matl subset)

     createDetailedDependencies (public)
       remembercomps
       createDetailedDependencies (private)
         DetailedTasks::possiblyCreateDependency or Task::addInternalDependency

   Then at the and:
     DetailedTasks::computeLocalTasks

GENERAL INFORMATION

   TaskGraph.h

   Steven G. Parker
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
****************************************/

class TaskGraph
{

public:
  TaskGraph(SchedulerCommon* sc,
            const ProcessorGroup* pg,
            Scheduler::tgType type,
            int index);

  ~TaskGraph();

  // eliminate copy, assignment, and move
  TaskGraph(const TaskGraph&) = delete;
  TaskGraph&
  operator=(const TaskGraph&) = delete;
  TaskGraph(TaskGraph&&)      = delete;
  TaskGraph&
  operator=(TaskGraph&&) = delete;

  /// Clears the TaskGraph and deletes all tasks.
  void
  initialize();

  /// Adds a task to the task graph.  If the task is empty, it
  /// deletes it.  Also, as each task is added, it updates the list
  /// of vars that are required from the old DW
  void
  addTask(std::shared_ptr<Task> task,
          const PatchSet* patchset,
          const MaterialSet* matlset);

  /// Sorts the tasks, and makes DetailedTask's out of them,
  /// and loads them into a new DetailedTasks object. (There is one
  /// DetailedTask for each PatchSubset and MaterialSubset in a Task,
  /// where a Task may have many PatchSubsets and MaterialSubsets.).
  /// Sorts using topologicalSort.
  DetailedTasks*
  createDetailedTasks(bool useInternalDeps,
                      const GridP& grid,
                      const GridP& oldGrid,
                      const bool hasDistalReqs = false);

  inline DetailedTasks*
  getDetailedTasks()
  {
    return d_detailed_tasks;
  }

  inline Scheduler::tgType
  getType() const
  {
    return d_type;
  }

  inline int
  getIndex()
  {
    return d_index;
  }

  /// This will go through the detailed tasks and create the
  /// dependencies need to communicate data across separate
  /// processors.  Calls the private createDetailedDependencies
  /// for each task as a helper.
  void
  createDetailedDependencies();

  /// Connects the tasks, but does not sort them.
  /// Used for the MixedScheduler, this routine has the side effect
  /// (just like the topological sort) of adding the reduction tasks.
  /// However, this routine leaves the tasks in the order they were
  /// added, so that reduction tasks are hit in the correct order
  /// by each MPI process.
  void
  nullSort(std::vector<Task*>& tasks);

  int
  getNumTasks() const;

  Task*
  getTask(int i);

  void
  remapTaskDWs(int dwmap[]);

  /// Assigns unique id numbers to each dependency based on name,
  /// material index, and patch.  In other words, even though a
  /// number of tasks depend on the same data, they create there
  /// own copy of the dependency data.  This routine determines
  /// that the dependencies are actually the same, and gives them
  /// the same id number.
  void
  assignUniqueMessageTags();

  /// sets the iteration of the current taskgraph in a multi-TG environment
  /// starting with 0
  void
  setIteration(int iter)
  {
    d_current_iteration = iter;
  }

  int
  getNumTaskPhases()
  {
    return d_num_task_phases;
  }

  std::vector<std::shared_ptr<Task>>&
  getTasks()
  {
    return d_tasks;
  }

  inline bool
  getDistalRequires() const
  {
    return d_has_distal_requires;
  }

  /// Makes and returns a map that associates VarLabel names with
  /// the materials the variable is computed for.
  using VarLabelMaterialMap = std::map<std::string, std::list<int>>;
  void
  makeVarLabelMaterialMap(VarLabelMaterialMap* result);

private:
  /// Used by (the public) createDetailedDependencies to store comps
  /// in a ComputeTable (See TaskGraphCompTable.cc).
  void
  remembercomps(DetailedTask* task, Task::Dependency* comp, CompTable& ct);

  /// This is the "detailed" version of addDependencyEdges.  It does for
  /// the public createDetailedDependencies member function essentially
  /// what addDependencyEdges does for setupTaskConnections.  This will
  /// set up the data dependencies that need to be communicated between
  /// processors.
  void
  createDetailedDependencies(DetailedTask* task,
                             Task::Dependency* req,
                             CompTable& ct,
                             bool modifies);

  /// Makes a DetailedTask from task with given PatchSubset and
  /// MaterialSubset.
  void
  createDetailedTask(Task* task,
                     const PatchSubset* patches,
                     const MaterialSubset* matls);

  /// find the processor that a variable (req) is on given patch and
  /// material.
  int
  findVariableLocation(Task::Dependency* req,
                       const Patch* patch,
                       int matl,
                       int iteration);

private:
  struct LabelLevel
  {
    LabelLevel(const std::string& key, const int level)
      : m_key(key)
      , m_level(level)
    {
    }

    std::string m_key{};
    int m_level{};

    bool
    operator<(const LabelLevel& rhs) const
    {
      if (this->m_level < rhs.m_level) {
        return true;
      } else if ((this->m_level == rhs.m_level) && (this->m_key < rhs.m_key)) {
        return true;
      }
      return false;
    }
  };

private:
  std::map<LabelLevel, int> d_max_ghost_for_varlabelmap{};

  SchedulerCommon* d_scheduler{ nullptr };
  LoadBalancer* d_load_balancer{ nullptr };
  const ProcessorGroup* d_proc_group{ nullptr };
  Scheduler::tgType d_type{};
  DetailedTasks* d_detailed_tasks{ nullptr };

  // how many times this taskgraph has executed this timestep
  int d_current_iteration{ 0 };

  // how many task phases this taskgraph has been through
  int d_num_task_phases{ 0 };

  int d_index{ -1 };

  // does this TG contain requires with halo > MAX_HALO_DEPTH
  bool d_has_distal_requires{ false };

  std::vector<std::shared_ptr<Task>> d_tasks{};

  ///////////////////////////////////////////////////////////////////////
  // Archived code for topological sort
  ///////////////////////////////////////////////////////////////////////
public:
  // This is so we can keep tasks independent of the task graph
  struct GraphSortInfo
  {

    GraphSortInfo()
      : m_visited{ false }
      , m_sorted{ false }
    {
    }

    bool m_visited;
    bool m_sorted;
  };

  using CompMap          = std::multimap<const VarLabel*, Task::Dependency*>;
  using GraphSortInfoMap = std::map<Task*, GraphSortInfo>;
  using LevelReductionTasksMap = std::map<VarLabelMatl<Level>, Task*>;

private:
  /// Helper function for setupTaskConnections, adding dependency edges
  /// for the given task for each of the require (or modify) depencies in
  /// the list whose head is req.  If modifies is true then each found
  /// compute will be replaced by its modifying dependency on the CompMap.
  void
  addDependencyEdges(Task* task,
                     GraphSortInfoMap& sortinfo,
                     Task::Dependency* req,
                     CompMap& comps,
                     LevelReductionTasksMap& reductionTasks,
                     bool modifies);

  bool
  overlaps(const Task::Dependency* comp, const Task::Dependency* req) const;

  /// Helper function for processTasks, processing the dependencies
  /// for the given task in the dependency list whose head is req.
  /// Will call processTask (recursively, as this is a helper for
  /// processTask) for each dependent task.
  void
  processDependencies(Task* task,
                      Task::Dependency* req,
                      std::vector<Task*>& sortedTasks,
                      GraphSortInfoMap& sortinfo) const;

  /// Called for each task, this "sorts" the taskgraph.
  /// This sorts in topological order by calling processDependency
  /// (which checks for cycles in the graph), which then recursively
  /// calls processTask for each dependentTask.  After this process is
  /// finished, then the task is added at the end of sortedTasks.
  void
  processTask(Task* task,
              std::vector<Task*>& sortedTasks,
              GraphSortInfoMap& sortinfo) const;

  /// Adds edges in the TaskGraph between requires/modifies and their
  /// associated computes.  Uses addDependencyEdges as a helper
  void
  setupTaskConnections(GraphSortInfoMap& sortinfo);

  /// sets up the task connections and puts them in a sorted order.
  /// Calls setupTaskConnections, which has the side effect of creating
  /// reduction tasks for tasks that compute reduction variables.
  /// calls processTask on each task to sort them.
  void
  topologicalSort(std::vector<Task*>& tasks);

  std::vector<Task::Edge*> m_edges;
};

} // End namespace Uintah

#endif // End CCA_COMPONENTS_SCHEDULERS_TASKGRAPH_H
