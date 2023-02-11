/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef VAANGO_CCA_COMPONENTS_SCHEDULERS_SCHEDULERCOMMON_H
#define VAANGO_CCA_COMPONENTS_SCHEDULERS_SCHEDULERCOMMON_H

#include <CCA/Components/Schedulers/OnDemandDataWarehouseP.h>
#include <CCA/Components/Schedulers/Relocate.h>
#include <CCA/Components/Schedulers/RuntimeStatsEnum.h>
#include <CCA/Ports/Scheduler.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/MaterialSetP.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <Core/Util/DebugStream.h>
#include <Core/Util/InfoMapper.h>

#include <iosfwd>
#include <list>
#include <map>
#include <set>

namespace Uintah {

class SimulationInterface;
class LoadBalancer;
class Output;
class DetailedTask;
class DetailedTasks;
class TaskGraph;
class LocallyComputedPatchVarMap;

using mm = VarLabelMatl<Level, DataWarehouse>;
using LabelMaterialMap =
  std::map<const VarLabel*, std::unique_ptr<MaterialSubset>, VarLabel::Compare>;
using VarLabelMaterialListMap = std::map<std::string, std::list<int>>;
using ReductionTasksMap =
  std::map<VarLabelMatl<Level, DataWarehouse>, std::shared_ptr<Task>>;
using VarLabelList = std::vector<std::vector<const VarLabel*>>;

class SchedulerCommon
  : public Scheduler
  , public UintahParallelComponent
{
public:
  // For calculating memory usage when sci-malloc is disabled...
  inline static char* s_start_addr{ nullptr };

public:
  SchedulerCommon(const ProcessorGroup* myworld);
  ~SchedulerCommon() override;

  SchedulerCommon(const SchedulerCommon&) = delete;
  SchedulerCommon(SchedulerCommon&&)      = delete;
  auto
  operator=(const SchedulerCommon&) -> SchedulerCommon& = delete;
  auto
  operator=(SchedulerCommon&&) -> SchedulerCommon& = delete;

  void
  setComponents(UintahParallelComponent* comp) override;
  void
  getComponents() override;
  void
  releaseComponents() override;

  SimulationInterface*
  getSimulator() override
  {
    return d_simulator;
  };

  LoadBalancer*
  getLoadBalancer() override
  {
    return d_load_balancer;
  };

  Output*
  getOutput() override
  {
    return d_output;
  }

  void
  problemSetup(const ProblemSpecP& prob_spec,
               const MaterialManagerP& mat_manager) override;

  void
  doEmitTaskGraphDocs() override;

  void
  checkMemoryUse(unsigned long& mem_use,
                 unsigned long& high_water,
                 unsigned long& max_mem_use) override;

  // sbrk memory start location (for memory tracking)
  void
  setStartAddr(char* start) override
  {
    s_start_addr = start;
  }

  char*
  getStartAddr() override
  {
    return s_start_addr;
  }

  void
  resetMaxMemValue() override;

  void
  initialize(int num_old_dw = 1, int num_new_dw = 1) override;

  void
  setParentDWs(DataWarehouse* parent_old_dw,
               DataWarehouse* parent_new_dw) override;

  void
  clearMappings() override;

  void
  mapDataWarehouse(Task::WhichDW, int dw_tag) override;

  void
  compile() override;

  /// For more complicated models
  void
  addTaskGraph(Scheduler::tgType type, int index /* = -1 */) override;

  TaskGraph*
  getTaskGraph(unsigned int index) override
  {
    ASSERT(index < d_task_graphs.size());
    return d_task_graphs[index].get();
  }

  int
  getNumTaskGraphs() override
  {
    return d_task_graphs.size();
  }

  // The number of task graphs is the number of task graphs that
  // will be used for all time steps excdept the initial time step,
  // where there is only one task graph.
  void
  setNumTaskGraphs(const int num_task_graphs) override
  {
    ASSERT(num_task_graphs >= 1);
    d_num_task_graphs = num_task_graphs;
  }

  void
  addTask(Task* task,
          const PatchSet* patches,
          const MaterialSet* materials,
          const int tg_num = -1) override;
  bool
  useSmallMessages() override
  {
    return d_use_small_messages;
  }

  /// Get all of the requires needed from the old data warehouse
  /// (carried forward).
  const std::vector<const Task::Dependency*>&
  getInitialRequires() const override
  {
    return d_init_requires;
  }

  const std::set<const VarLabel*, VarLabel::Compare>&
  getInitialRequiredVars() const override
  {
    return d_init_required_vars;
  }

  const std::set<const VarLabel*, VarLabel::Compare>&
  getComputedVars() const override
  {
    return d_computed_vars;
  }

  DataWarehouse*
  get_dw(int idx) override;

  DataWarehouse*
  getLastDW() override;

  void
  logMemoryUse() override;

  void
  advanceDataWarehouse(const GridP& grid, bool initialization = false) override;

  void
  fillDataWarehouses(const GridP& grid) override;

  void
  replaceDataWarehouse(int index,
                       const GridP& grid,
                       bool initialization = false) override;

  // Get the SuperPatch (set of connected patches making a larger rectangle)
  // for the given label and patch and find the largest extents encompassing
  // the expected ghost cells (requiredLow, requiredHigh) and the requested
  // ghost cells as well (requestedLow, requestedHigh) for each of the
  // patches.  Required and requested will besame if requestedNumGCells = 0.
  const std::vector<const Patch*>*
  getSuperPatchExtents(const VarLabel* label,
                       int matl_index,
                       const Patch* patch,
                       Ghost::GhostType requested_ghost_type,
                       int requested_num_ghost_cells,
                       IntVector& required_low,
                       IntVector& required_high,
                       IntVector& requested_low,
                       IntVector& requested_high) const override;

  // Makes and returns a map that maps strings to VarLabels of
  // that name and a list of material indices for which that
  // variable is valid (at least according to d_allcomps).
  std::unique_ptr<Scheduler::VarLabelMaterialMap>
  makeVarLabelMaterialMap() override;

  bool
  isOldDW(int idx) const override;

  bool
  isNewDW(int idx) const override;

  // Only called by the SimulationController, and only once, and only
  // if the simulation has been "restarted."
  void
  setGeneration(int id) override
  {
    d_generation = id;
  }

  // This function will copy the data from the old grid to the new grid.
  // The PatchSubset structure will contain a patch on the new grid.
  void
  copyDataToNewGrid(const ProcessorGroup*,
                    const PatchSubset* patches,
                    const MaterialSubset*,
                    DataWarehouse* old_dw,
                    DataWarehouse* new_dw);

  // Schedule particle relocation without the need to supply pre_relocation
  // variables. Use with caution until as this requires further testing (tsaad).
  void
  scheduleParticleRelocation(
    const LevelP& coarsest_level_with_particles,
    const VarLabel* position_label,
    const std::vector<std::vector<const VarLabel*>>& other_labels,
    const MaterialSet* materials) override;

  void
  scheduleParticleRelocation(
    const LevelP& coarsest_level_with_particles,
    const VarLabel* old_position_label,
    const std::vector<std::vector<const VarLabel*>>& old_other_labels,
    const VarLabel* new_postion_label,
    const std::vector<std::vector<const VarLabel*>>& new_other_labels,
    const VarLabel* particle_id_label,
    const MaterialSet* materials) override;

  void
  setPositionVar(const VarLabel* pos_label) override
  {
    d_reloc_new_pos_label = pos_label;
  }

  void
  scheduleAndDoDataCopy(const GridP& grid) override;

  // Clear the recorded task monitoring attribute values.
  void
  clearTaskMonitoring() override;

  // Schedule the recording of the task monitoring attribute values.
  void
  scheduleTaskMonitoring(const LevelP& level) override;
  void
  scheduleTaskMonitoring(const PatchSet* patches) override;

  // Record the task monitoring attribute values.
  virtual void
  recordTaskMonitoring(const ProcessorGroup* group,
                       const PatchSubset* patches,
                       const MaterialSubset* materials,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw);

  //! override default behavior of copying, scrubbing, and such
  void
  overrideVariableBehavior(const std::string& var,
                           bool treat_as_old,
                           bool copy_data,
                           bool no_scrub,
                           bool not_copy_data  = false,
                           bool not_checkpoint = false) override;

  auto
  getNoScrubVars() -> const std::set<std::string>&
  {
    return d_no_scrub_vars;
  }

  auto
  getCopyDataVars() -> const std::set<std::string>&
  {
    return d_copy_data_vars;
  }

  auto
  getNotCopyDataVars() -> const std::set<std::string>&
  {
    return d_not_copy_data_vars;
  }

  const std::set<std::string>&
  getNotCheckPointVars() const override
  {
    return d_not_checkpoint_vars;
  }

  virtual auto
  useInternalDeps() -> bool;

  int
  getMaxGhost() const override
  {
    return d_max_ghost_cells;
  }

  int
  getMaxDistalGhost() const override
  {
    return d_max_distal_ghost_cells;
  }

  int
  getMaxLevelOffset() const override
  {
    return d_max_level_offset;
  }

  bool
  isCopyDataTimestep() const override
  {
    return d_is_copy_data_timestep;
  }

  bool
  copyTimestep() const override
  {
    return (d_is_copy_data_timestep || d_is_init_timestep);
  }

  void
  setInitTimestep(bool is_init_timestep) override
  {
    d_is_init_timestep = is_init_timestep;
  }

  void
  setRestartInitTimestep(bool is_restart_init_timestep) override
  {
    d_is_restart_init_timestep = is_restart_init_timestep;
  }

  bool
  isRestartInitTimestep() const override
  {
    return d_is_restart_init_timestep;
  }

  void
  setRuntimeStats(
    ReductionInfoMapper<RuntimeStatsEnum, double>* runtime_stats) override
  {
    d_runtime_stats = runtime_stats;
  };

public:
  const VarLabel* d_reloc_new_pos_label{ nullptr };

  // number of schedulers and subschedulers
  int m_num_schedulers{ 0 };

protected:
  void
  finalizeTimestep();

  void
  makeTaskGraphDoc(const DetailedTasks* dtask, int rank = 0);

  void
  emitNode(const DetailedTask* dtask,
           double start,
           double duration,
           double execution_duration);

  void
  finalizeNodes(int process = 0);

  template<class T>
  void
  printTrackedValues(GridVariable<T>* var,
                     const IntVector& start,
                     const IntVector& end);

  void
  printTrackedVars(DetailedTask* dtask, int when);

  virtual void
  verifyChecksum() = 0;

  // Method for summing up the task contributions.
  void
  sumTaskMonitoringValues(DetailedTask* dtask);

protected:
  enum
  {
    PRINT_BEFORE_COMM = 1,
    PRINT_BEFORE_EXEC = 2,
    PRINT_AFTER_EXEC  = 4
  };

  // Some places need to know if this is a copy data timestep or
  // a normal timestep.  (A copy data timestep is AMR's current
  // method of getting data from an old to a new grid).
  bool d_is_copy_data_timestep{ false };
  bool d_is_init_timestep{ false };
  bool d_is_restart_init_timestep{ false };
  int d_current_task_graph{ -1 };
  int d_generation{ 0 };
  int d_dw_map[Task::TotalDWs];

  SimulationInterface* d_simulator{ nullptr };
  LoadBalancer* d_load_balancer{ nullptr };
  Output* d_output{ nullptr };

  MaterialManagerP d_materialManager{ nullptr };
  std::vector<OnDemandDataWarehouseUP> d_dws;
  std::vector<std::unique_ptr<TaskGraph>> d_task_graphs;

  //! These are so we can track certain variables over the taskgraph's
  //! execution.
  int d_tracking_vars_print_location{ 0 };
  int d_tracking_patch_id{ -1 };
  int d_tracking_level{ -1 };
  double d_tracking_start_time{ 1.0 };
  double d_tracking_end_time{ 0.0 };
  IntVector d_tracking_start_index{ IntVector(-9, -9, -9) };
  IntVector d_tracking_end_index{ IntVector(-9, -9, -9) };
  std::vector<std::string> d_tracking_vars;
  std::vector<std::string> d_tracking_tasks;
  std::vector<Task::WhichDW> d_tracking_dws;

  // Optional task monitoring.
  std::unique_ptr<MaterialSubset> d_dummy_matl{ nullptr };

  bool d_monitoring{ false };          // Monitoring on/off.
  bool d_monitoring_per_cell{ false }; // Record the task runtime attributes
                                       // on a per cell basis rather than a
                                       // per patch basis.

  // Maps for the global or local tasks to be monitored.
  std::map<std::string, const VarLabel*> d_monitoring_tasks[2];
  std::map<std::string, std::map<int, double>> d_monitoring_values[2];

  // so we can manually copy vars between AMR levels
  std::set<std::string> d_copy_data_vars;

  // ignore copying these vars between AMR levels
  std::set<std::string> d_not_copy_data_vars;

  // vars manually set not to scrub (normally when needed between a normal
  // taskgraph and the regridding phase)
  std::set<std::string> d_no_scrub_vars;

  // treat variable as an "old" var - will be checkpointed, copied, and only
  // scrubbed from an OldDW
  std::set<std::string> d_treat_as_old_vars;

  // do not checkpoint these variables
  std::set<std::string> d_not_checkpoint_vars;

  ReductionInfoMapper<RuntimeStatsEnum, double>* d_runtime_stats{ nullptr };

private:
  // problem setup for var tracking
  void
  setupVarTracker(const ProblemSpecP& params);

  // problem setup for task monitoring
  void
  setupTaskMonitoring(const ProblemSpecP& params);

  // helper method for primary addTask()
  void
  addTask(std::shared_ptr<Task> task,
          const PatchSet* patches,
          const MaterialSet* matls,
          const int tg_num);

private:
  ProblemSpecP d_graph_doc{ nullptr };
  ProblemSpecP d_graph_nodes{ nullptr };

  std::unique_ptr<std::ofstream> d_mem_log_file{ nullptr };

  Relocate d_relocate_1;

  // whether or not to send a small message (takes more work to organize)
  // or a larger one (more communication time)
  bool d_use_small_messages{ true };
  bool d_emit_task_graph{ false };
  int d_num_task_graphs{ 1 };
  int d_num_tasks{ 0 };
  int d_num_old_dws{ 0 };

  //! These are to store which vars we have to copy to the new grid
  //! in a copy data task.  Set in scheduleDataCopy and used in
  //! copyDataToNewGrid.
  std::vector<LabelMaterialMap> d_label_matls;

  ReductionTasksMap d_reduction_tasks;

  std::unique_ptr<LocallyComputedPatchVarMap> d_local_patch_var_map{ nullptr };

  //! set in addTask - can be used until initialize is called...
  std::vector<const Task::Dependency*> d_init_requires;
  std::set<const VarLabel*, VarLabel::Compare> d_init_required_vars;
  std::set<const VarLabel*, VarLabel::Compare> d_computed_vars;

  // Maximum memory used as sampled across a given timestep.
  unsigned long d_max_mem_used{ 0 };

  // max ghost cells of standard tasks - will be used for loadbalancer to create
  // neighborhood
  int d_max_ghost_cells{ 0 };

  // max ghost cells for tasks with distal requirements (e.g. RMCRT) - will be
  // used for loadbalancer to create neighborhood
  int d_max_distal_ghost_cells{ 0 };

  // max level offset of all tasks - will be used for loadbalancer to create
  // neighborhood
  int d_max_level_offset{ 0 };

  // task-graph needs access to reduction task map, etc
  friend class TaskGraph;
};
} // End namespace Uintah

#endif
