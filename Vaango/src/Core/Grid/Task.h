/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Parresia Research Limited, New Zealand
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

#ifndef VAANGO_CORE_GRID_Task_H
#define VAANGO_CORE_GRID_Task_H

#include <CCA/Ports/DataWarehouseP.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Grid/Ghost.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Util/constHandle.h>

#include <Core/Util/TupleHelpers.hpp>
#include <Core/Util/DOUT.hpp>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace Uintah {

class Level;
class DataWarehouse;
class ProcessorGroup;
class DetailedTask;

class Task
{
public:
  enum CallBackEvent
  {
    CPU,    // <- normal CPU task, happens when a GPU enabled task runs on CPU
    preGPU, // <- pre GPU kernel callback, happens before CPU->GPU copy
            // (reserved, not implemented yet... )
    GPU,    // <- GPU kernel callback, happens after dw: CPU->GPU copy, kernel
            // launch should be queued in this callback
    postGPU // <- post GPU kernel callback, happens after dw: GPU->CPU copy but
            // before MPI sends.
  };

protected:
  // base Action class
  class ActionBase
  {
  public:
    virtual ~ActionBase() = default;
    virtual void doit(DetailedTask* dtask,
                      CallBackEvent event,
                      const ProcessorGroup* pg,
                      const PatchSubset* patches,
                      const MaterialSubset* matls,
                      DataWarehouse* fromDW,
                      DataWarehouse* toDW,
                      void* oldTaskGpuDW,
                      void* newTaskGpuDW,
                      void* stream,
                      int deviceID) = 0;
  };

private: // class Task
  // CPU Action constructor
  template<typename T, typename... Args>
  class Action : public ActionBase
  {
  private:
    T* ptr;
    void (T::*pmf)(const ProcessorGroup* pg,
                   const PatchSubset* patches,
                   const MaterialSubset* matls,
                   DataWarehouse* fromDW,
                   DataWarehouse* toDW,
                   Args... args);
    std::tuple<Args...> m_args;

  public:
    // class Action
    Action(T* ptr,
           void (T::*pmf)(const ProcessorGroup* pg,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse* fromDW,
                          DataWarehouse* toDW,
                          Args... args),
           Args... args)
      : ptr(ptr)
      , pmf(pmf)
      , m_args(std::forward<Args>(args)...)
    {
    }

    virtual ~Action() {}

    virtual void doit([[maybe_unused]] DetailedTask* dtask,
                      [[maybe_unused]] CallBackEvent event,
                      const ProcessorGroup* pg,
                      const PatchSubset* patches,
                      const MaterialSubset* matls,
                      DataWarehouse* fromDW,
                      DataWarehouse* toDW,
                      [[maybe_unused]] void* oldTaskGpuDW,
                      [[maybe_unused]] void* newTaskGpuDW,
                      [[maybe_unused]] void* stream,
                      [[maybe_unused]] int deviceID) override
    {
      doit_impl(pg,
                patches,
                matls,
                fromDW,
                toDW,
                typename Tuple::gens<sizeof...(Args)>::type());
    }

  private: // class Action
    template<int... S>
    void doit_impl(const ProcessorGroup* pg,
                   const PatchSubset* patches,
                   const MaterialSubset* matls,
                   DataWarehouse* fromDW,
                   DataWarehouse* toDW,
                   Tuple::seq<S...>)
    {
      (ptr->*pmf)(pg, patches, matls, fromDW, toDW, std::get<S>(m_args)...);
    }
  }; // end class Action for CPUs

  // GPU (device) Action constructor
  template<typename T, typename... Args>
  class ActionDevice : public ActionBase
  {
    T* ptr;
    void (T::*pmf)(DetailedTask* dtask,
                   CallBackEvent event,
                   const ProcessorGroup* pg,
                   const PatchSubset* patches,
                   const MaterialSubset* m_matls,
                   DataWarehouse* fromDW,
                   DataWarehouse* toDW,
                   void* oldTaskGpuDW,
                   void* newTaskGpuDW,
                   void* stream,
                   int deviceID,
                   Args... args);
    std::tuple<Args...> m_args;

  public: // class ActionDevice
    ActionDevice(T* ptr,
                 void (T::*pmf)(DetailedTask* dtask,
                                CallBackEvent event,
                                const ProcessorGroup* pg,
                                const PatchSubset* patches,
                                const MaterialSubset* m_matls,
                                DataWarehouse* fromDW,
                                DataWarehouse* toDW,
                                void* oldTaskGpuDW,
                                void* newTaskGpuDW,
                                void* stream,
                                int deviceID,
                                Args... args),
                 Args... args)
      : ptr(ptr)
      , pmf(pmf)
      , m_args(std::forward<Args>(args)...)
    {
    }

    virtual ~ActionDevice() {}

    virtual void doit(DetailedTask* dtask,
                      CallBackEvent event,
                      const ProcessorGroup* pg,
                      const PatchSubset* patches,
                      const MaterialSubset* matls,
                      DataWarehouse* fromDW,
                      DataWarehouse* toDW,
                      void* oldTaskGpuDW,
                      void* newTaskGpuDW,
                      void* stream,
                      int deviceID)
    {
      doit_impl(dtask,
                event,
                pg,
                patches,
                matls,
                fromDW,
                toDW,
                oldTaskGpuDW,
                newTaskGpuDW,
                stream,
                deviceID,
                typename Tuple::gens<sizeof...(Args)>::type());
    }

  private: // class ActionDevice
    template<int... S>
    void doit_impl(DetailedTask* dtask,
                   CallBackEvent event,
                   const ProcessorGroup* pg,
                   const PatchSubset* patches,
                   const MaterialSubset* matls,
                   DataWarehouse* fromDW,
                   DataWarehouse* toDW,
                   void* oldTaskGpuDW,
                   void* newTaskGpuDW,
                   void* stream,
                   int deviceID,
                   Tuple::seq<S...>)
    {
      (ptr->*pmf)(dtask,
                  event,
                  pg,
                  patches,
                  matls,
                  fromDW,
                  toDW,
                  oldTaskGpuDW,
                  newTaskGpuDW,
                  stream,
                  deviceID,
                  std::get<S>(m_args)...);
    }

  }; // end GPU (device) Action constructor

public: // class Task
  enum WhichDW
  {
    None = -1,
    OldDW = 0,
    NewDW = 1,
    CoarseOldDW = 2,
    CoarseNewDW = 3,
    ParentOldDW = 4,
    ParentNewDW = 5,
    TotalDWs = 6
  };

  enum
  {
    NoDW = -1,
    InvalidDW = -2
  };

  enum TaskType
  {
    Normal,
    Reduction, // tasks with MPI reductions
    InitialSend,
    OncePerProc, // make sure to pass a PerProcessor PatchSet to the addTask
                 // function
    Output,
    OutputGlobalVars, // task the outputs the reduction variables
    Spatial,          // e.g. Radiometer task (spatial scheduling); must call
                      // task->setType(Task::Spatial)
    Hypre
  };

  Task(const std::string& taskName, TaskType type)
    : m_task_name(taskName)
    , m_action(nullptr)
  {
    m_tasktype = type;
    initialize();
  }

  // CPU Task constructor
  template<typename T, typename... Args>
  Task(const std::string& taskName,
       T* ptr,
       void (T::*pmf)(const ProcessorGroup*,
                      const PatchSubset*,
                      const MaterialSubset*,
                      DataWarehouse*,
                      DataWarehouse*,
                      Args...),
       Args... args)
    : m_task_name(taskName)
    , m_action(scinew Action<T, Args...>(ptr, pmf, std::forward<Args>(args)...))
  {
    m_tasktype = Normal;
    initialize();
  }

  // Device (GPU) Task constructor
  template<typename T, typename... Args>
  Task(const std::string& taskName,
       T* ptr,
       void (T::*pmf)(DetailedTask* m_task,
                      CallBackEvent event,
                      const ProcessorGroup* pg,
                      const PatchSubset* patches,
                      const MaterialSubset* m_matls,
                      DataWarehouse* fromDW,
                      DataWarehouse* toDW,
                      void* old_TaskGpuDW,
                      void* new_TaskGpuDW,
                      void* stream,
                      int deviceID,
                      Args... args),
       Args... args)
    : m_task_name(taskName)
    , m_action(
        scinew ActionDevice<T, Args...>(ptr, pmf, std::forward<Args>(args)...))
  {
    initialize();
    m_tasktype = Normal;
  }

  void initialize();

  virtual ~Task();

  // eliminate copy, assignment and move
  Task(const Task&) = delete;
  Task& operator=(const Task&) = delete;
  Task(Task&&) = delete;
  Task& operator=(Task&&) = delete;

  void hasSubScheduler(bool state = true);
  inline bool getHasSubScheduler() const { return m_has_subscheduler; }
  void usesMPI(bool state);
  inline bool usesMPI() const { return m_uses_mpi; }
  void usesThreads(bool state);
  inline bool usesThreads() const { return m_uses_threads; }
  void usesDevice(bool state, int maxStreamsPerTask = 1);
  inline bool usesDevice() const { return m_uses_device; }
  inline int maxStreamsPerTask() const { return m_max_streams_per_task; }

  inline void setDebugFlag(bool in) { m_debug_flag = in; }
  inline bool getDebugFlag() const { return m_debug_flag; }

  enum MaterialDomainSpec
  {
    NormalDomain, // <- Normal/default setting
    OutOfDomain,  // <- Require things from all material
  };

  enum PatchDomainSpec
  {
    ThisLevel,   // <- Normal/default setting
    CoarseLevel, // <- AMR :  The data on the coarse level under the range of
                 // the fine patches (including extra cells or boundary layers)
    FineLevel,   // <- AMR :  The data on the fine level under the range of the
                 // coarse patches (including extra cells or boundary layers)
    OtherGridDomain // for when we copy data to new grid after a regrid.
  };

  enum class SearchTG
  {
    OldTG, // <- Search the OldTG for the computes if they aren't found in
           // NewTG
    NewTG
  };

  //////////
  // Most general case
  void
    requires(WhichDW,
             const VarLabel* label,
             const PatchSubset* patches,
             PatchDomainSpec patches_dom,
             int level_offset,
             const MaterialSubset* matls,
             MaterialDomainSpec matls_dom,
             Ghost::GhostType gtype,
             int numGhostCells = 0,
             SearchTG whichTG = SearchTG::NewTG);

  //////////
  // Like general case, level_offset is not specified
  void
    requires(WhichDW,
             const VarLabel* label,
             const PatchSubset* patches,
             PatchDomainSpec patches_dom,
             const MaterialSubset* matls,
             MaterialDomainSpec matls_dom,
             Ghost::GhostType gtype,
             int numGhostCells = 0,
             SearchTG whichTG = SearchTG::NewTG);

  void
    requires(WhichDW,
             const VarLabel* label,
             Ghost::GhostType gtype,
             int numGhostCells = 0,
             SearchTG whichTG = SearchTG::NewTG);

  void
    requires(WhichDW,
             const VarLabel* label,
             const PatchSubset* patches,
             const MaterialSubset* matls,
             Ghost::GhostType gtype,
             int numGhostCells = 0,
             SearchTG whichTG = SearchTG::NewTG);

  void
    requires(WhichDW,
             const VarLabel* label,
             const PatchSubset* patches,
             Ghost::GhostType gtype,
             int numGhostCells = 0,
             SearchTG whichTG = SearchTG::NewTG);

  void
    requires(WhichDW,
             const VarLabel* label,
             const MaterialSubset* matls,
             Ghost::GhostType gtype,
             int numGhostCells = 0,
             SearchTG whichTG = SearchTG::NewTG);

  void
    requires(WhichDW,
             const VarLabel* label,
             const MaterialSubset* matls,
             MaterialDomainSpec matls_dom,
             Ghost::GhostType gtype,
             int numGhostCells = 0,
             SearchTG whichTG = SearchTG::NewTG);

  //////////
  // Requires only for reduction variables
  void
    requires(WhichDW,
             const VarLabel* label,
             const Level* level = nullptr,
             const MaterialSubset* matls = nullptr,
             MaterialDomainSpec matls_dom = NormalDomain,
             SearchTG whichTG = SearchTG::NewTG);

  //////////
  // Requires for reduction variables or perpatch veriables
  void
    requires(WhichDW,
             const VarLabel* label,
             const MaterialSubset* matls,
             SearchTG whichTG = SearchTG::NewTG);

  //////////
  // Requires only for perpatch variables
  void
    requires(WhichDW,
             const VarLabel* label,
             const PatchSubset* patches,
             const MaterialSubset* matls = nullptr);

  //////////
  // Most general case
  void computes(const VarLabel* label,
                const PatchSubset* patches,
                PatchDomainSpec patches_domain,
                const MaterialSubset* matls,
                MaterialDomainSpec matls_domain);

  void computes(const VarLabel* label,
                const PatchSubset* patches = nullptr,
                const MaterialSubset* matls = nullptr);

  void computes(const VarLabel* label, const MaterialSubset* matls);

  void computes(const VarLabel* label,
                const MaterialSubset* matls,
                MaterialDomainSpec matls_domain);

  void computes(const VarLabel* label,
                const PatchSubset* patches,
                PatchDomainSpec patches_domain);

  void computes(const VarLabel* label,
                const Level* level,
                const MaterialSubset* matls = nullptr,
                MaterialDomainSpec matls_domain = NormalDomain);

  //////////
  /*! \brief Allows a task to do a computes and modify with ghost cell
   specification.
   *
   *  \warning Uintah was built around the assumption that one is NOT allowed
      to compute or modify ghost cells. Therefore, it is unlawful in the Uintah
   sense to add a computes/modifies with ghost cells. However, certain
   components such as Wasatch break that design assumption from the point of
   view that, if a task can fill-in ghost values, then by all means do that and
      avoid an extra communication in the process. This, for example, is the
      case when one extrapolates data from the interior (e.g. Dynamic
   Smagorinsky model). Be aware that the ghost-values computed/modified in one
   patch will NOT be reproduced/correspond to interior cells of the neighboring
   patch, and vice versa.
   */
  void modifiesWithScratchGhost(const VarLabel* label,
                                const PatchSubset* patches,
                                PatchDomainSpec patches_domain,
                                const MaterialSubset* matls,
                                MaterialDomainSpec matls_domain,
                                Ghost::GhostType gtype,
                                int numGhostCells,
                                SearchTG whichTG = SearchTG::NewTG);

  void computesWithScratchGhost(const VarLabel* label,
                                const MaterialSubset* matls,
                                MaterialDomainSpec matls_domain,
                                Ghost::GhostType gtype,
                                int numGhostCells,
                                SearchTG whichTG = SearchTG::NewTG);

  //////////
  // Most general case
  void modifies(const VarLabel* label,
                const PatchSubset* patches,
                PatchDomainSpec patches_domain,
                const MaterialSubset* matls,
                MaterialDomainSpec matls_domain,
                SearchTG whichTG = SearchTG::NewTG);

  void modifies(const VarLabel* label,
                const PatchSubset* patches,
                const MaterialSubset* matls,
                SearchTG whichTG = SearchTG::NewTG);

  void modifies(const VarLabel* label,
                const MaterialSubset* matls,
                SearchTG whichTG = SearchTG::NewTG);

  void modifies(const VarLabel* label,
                const MaterialSubset* matls,
                MaterialDomainSpec matls_domain,
                SearchTG whichTG = SearchTG::NewTG);

  void modifies(const VarLabel* label, SearchTG whichTG = SearchTG::NewTG);

  //////////
  // Modify reduction vars
  void modifies(const VarLabel* label,
                const Level* level,
                const MaterialSubset* matls = nullptr,
                MaterialDomainSpec matls_domain = NormalDomain,
                SearchTG whichTG = SearchTG::NewTG);

  //////////
  // Tells the task to actually execute the function assigned to it.
  virtual void doit(DetailedTask* dtask,
                    CallBackEvent event,
                    const ProcessorGroup* pg,
                    const PatchSubset*,
                    const MaterialSubset*,
                    std::vector<DataWarehouseSP>& dws,
                    void* oldTaskGpuDW,
                    void* newTaskGpuDW,
                    void* stream,
                    int deviceID);

  inline const std::string& getName() const { return m_task_name; }
  inline const PatchSet* getPatchSet() const { return m_patch_set; }

  inline const MaterialSet* getMaterialSet() const { return m_matl_set; }

  bool hasDistalRequires() const; // determines if this Task has any "distal"
                                  // ghost cell requirements

  int m_phase{ -1 }; // synchronized phase id, for dynamic task scheduling
  int m_comm{ -1 };  // task communicator id, for threaded task scheduling
  std::map<int, int> m_max_ghost_cells; // max ghost cells of this task
  int m_max_level_offset{ 0 };          // max level offset of this task

  // Used in compiling the task graph with topological sort.
  //   Kept around for historical and reproducability reasons - APH, 04/05/19
  std::set<Task*> m_child_tasks;
  std::set<Task*> m_all_child_tasks;

  enum DepType
  {
    Modifies,
    Computes,
    Requires
  };

  struct Edge;

  struct Dependency
  {
    Dependency* next{ nullptr };
    DepType dep_type;
    Task* task{ nullptr };
    const VarLabel* var{ nullptr };
    bool look_in_old_tg;
    const PatchSubset* patches{ nullptr };
    const MaterialSubset* matls{ nullptr };
    const Level* reduction_level{ nullptr };
    PatchDomainSpec patches_dom{ ThisLevel };
    MaterialDomainSpec matls_dom;
    Ghost::GhostType gtype{ Ghost::None };
    WhichDW whichdw; // Used only by Requires

    // in the multi-TG construct, this will signify that the required
    // var will be constructed by the old TG
    int num_ghost_cells{ 0 };
    int level_offset{ 0 };

    int mapDataWarehouse() const { return task->mapDataWarehouse(whichdw); }

    Dependency(DepType deptype,
               Task* task,
               WhichDW dw,
               const VarLabel* var,
               SearchTG whichTG,
               const PatchSubset* patches,
               const MaterialSubset* matls,
               PatchDomainSpec patches_dom = ThisLevel,
               MaterialDomainSpec matls_dom = NormalDomain,
               Ghost::GhostType gtype = Ghost::None,
               int numGhostCells = 0,
               int level_offset = 0);
    Dependency(DepType deptype,
               Task* task,
               WhichDW dw,
               const VarLabel* var,
               SearchTG whichTG,
               const Level* reductionLevel,
               const MaterialSubset* matls,
               MaterialDomainSpec matls_dom = NormalDomain);
    ~Dependency();

    // eliminate copy, assignment and move
    Dependency(const Dependency&) = delete;
    Dependency& operator=(const Dependency&) = delete;
    Dependency(Dependency&&) = delete;
    Dependency& operator=(Dependency&&) = delete;

    constHandle<PatchSubset> getPatchesUnderDomain(
      const PatchSubset* domainPatches) const;

    constHandle<MaterialSubset> getMaterialsUnderDomain(
      const MaterialSubset* domainMaterials) const;

    // Used in compiling the task graph with topological sort.
    //   Kept around for historical and reproducability reasons - APH, 04/05/19
    Edge* m_req_head{ nullptr };
    Edge* m_req_tail{ nullptr };
    Edge* m_comp_head{ nullptr };
    Edge* m_comp_tail{ nullptr };
    inline void addComp(Edge* edge);
    inline void addReq(Edge* edge);

  private:
    static constHandle<PatchSubset> getOtherLevelPatchSubset(
      PatchDomainSpec dom,
      int level_offset,
      const PatchSubset* subset,
      const PatchSubset* domainSubset,
      int ngc);
  }; // end struct Dependency

  // Used in compiling the task graph with topological sort.
  //   Kept around for historical and reproducability reasons - APH, 04/05/19
  struct Edge
  {
    const Dependency* comp{ nullptr };
    Edge* comp_next{ nullptr };
    const Dependency* req{ nullptr };
    Edge* req_next{ nullptr };
    inline Edge(const Dependency* comp_in, const Dependency* req_in)
      : comp(comp_in)
      , req(req_in)
    {
    }
  };

  const Dependency* getComputes() const { return m_comp_head; }
  const Dependency* getRequires() const { return m_req_head; }
  const Dependency* getModifies() const { return m_mod_head; }

  Dependency* getComputes() { return m_comp_head; }
  Dependency* getRequires() { return m_req_head; }
  Dependency* getModifies() { return m_mod_head; }

  // finds if it computes or modifies var
  bool hasComputes(const VarLabel* var,
                   int matlIndex,
                   const Patch* patch) const;

  // finds if it requires or modifies var
  bool hasRequires(const VarLabel* var,
                   int matlIndex,
                   const Patch* patch,
                   Uintah::IntVector lowOffset,
                   Uintah::IntVector highOffset,
                   WhichDW dw) const;

  // finds if it modifies var
  bool hasModifies(const VarLabel* var,
                   int matlIndex,
                   const Patch* patch) const;

  bool isReductionTask() const { return m_tasktype == Reduction; }

  void setType(TaskType tasktype) { m_tasktype = tasktype; }
  TaskType getType() const { return m_tasktype; }

  //////////
  // Prints out information about the task...
  void display(std::ostream& out) const;

  //////////
  // Prints out all information about the task, including dependencies
  void displayAll_DOUT( Uintah::Dout& dbg) const;
  void displayAll(std::ostream& out) const;

  int mapDataWarehouse(WhichDW dw) const;
  DataWarehouse* mapDataWarehouse(WhichDW dw,
                                  std::vector<DataWarehouseSP>& dws) const;

  int getSortedOrder() const { return m_sorted_order; }

  void setSortedOrder(int order) { m_sorted_order = order; }

  void setMapping(int dwmap[TotalDWs]);

  void setSets(const PatchSet* patches, const MaterialSet* matls);

  static const MaterialSubset* getGlobalMatlSubset();

private: // class Task
  using DepMap = std::multimap<const VarLabel*, Dependency*, VarLabel::Compare>;

  Dependency* isInDepMap(const DepMap& depMap,
                         const VarLabel* var,
                         int matlIndex,
                         const Patch* patch) const;

  std::string m_task_name;

protected: // class Task
  ActionBase* m_action{ nullptr };

  inline static MaterialSubset* s_global_material_subset{ nullptr };

  Dependency* m_comp_head{ nullptr };
  Dependency* m_comp_tail{ nullptr };
  Dependency* m_req_head{ nullptr };
  Dependency* m_req_tail{ nullptr };
  Dependency* m_mod_head{ nullptr };
  Dependency* m_mod_tail{ nullptr };

  DepMap m_requires_old_dw;
  DepMap m_computes; // also contains modifies
  DepMap m_requires; // also contains modifies
  DepMap m_modifies;

  const PatchSet* m_patch_set{ nullptr };
  const MaterialSet* m_matl_set{ nullptr };

  bool m_uses_mpi{ false };
  bool m_uses_threads{ false };
  bool m_uses_device{ false };
  bool m_subpatch_capable{ false };
  bool m_has_subscheduler{ false };
  bool m_debug_flag{ false };

  int m_max_streams_per_task{ 1 };

  TaskType m_tasktype;

  int m_dwmap[TotalDWs];
  int m_sorted_order{ -1 };

  friend std::ostream& operator<<(std::ostream& out, const Uintah::Task& task);
  friend std::ostream& operator<<(std::ostream& out,
                                  const Uintah::Task::TaskType& tt);
  friend std::ostream& operator<<(std::ostream& out,
                                  const Uintah::Task::Dependency& dep);

}; // end class Task

// Used in compiling the task graph with topological sort.
//   Kept around for historical and reproducability reasons - APH, 04/05/19
inline void
Task::Dependency::addComp(Edge* edge)
{
  if (m_comp_tail) {
    m_comp_tail->comp_next = edge;
  } else {
    m_comp_head = edge;
  }
  m_comp_tail = edge;
}

inline void
Task::Dependency::addReq(Edge* edge)
{
  if (m_req_tail) {
    m_req_tail->req_next = edge;
  } else {
    m_req_head = edge;
  }
  m_req_tail = edge;
}

} // End namespace Uintah

// This mus tbe at the bottom
#include <CCA/Ports/DataWarehouse.h>

#endif
