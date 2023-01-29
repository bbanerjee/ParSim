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

#ifndef VAANGO_COMPONENTS_SCHEDULERS_ONDEMANDDATAWAREHOUSE_H
#define VAANGO_COMPONENTS_SCHEDULERS_ONDEMANDDATAWAREHOUSE_H

#include <CCA/Ports/DataWarehouse.h>

#include <CCA/Components/Schedulers/DWDatabase.h>
#include <CCA/Components/Schedulers/OnDemandDataWarehouseP.h>
#include <CCA/Components/Schedulers/SendState.h>

#include <Core/Containers/FastHashTable.h>

#include <Core/Grid/Grid.h>
#include <Core/Grid/Variables/PSPatchMatlGhost.h>
#include <Core/Grid/Variables/VarLabelMatl.h>

#include <Core/Parallel/MasterLock.h>
#include <Core/Parallel/UintahMPI.h>

#include <iosfwd>
#include <map>
#include <vector>

namespace Uintah {

inline const Patch*
getRealDomain(const Patch* patch)
{
  return patch->getRealPatch();
}

inline const Level*
getRealDomain(const Level* level)
{
  return level;
}

class BufferInfo;
class DependencyBatch;
class DetailedDep;
class DetailedTasks;
class LoadBalancer;
class Patch;
class ProcessorGroup;
class SendState;
class TypeDescription;

class OnDemandDataWarehouse : public DataWarehouse
{

public:
  OnDemandDataWarehouse(const ProcessorGroup* myworld,
                        Scheduler* scheduler,
                        int generation,
                        const GridP& grid,
                        bool isInitializationDW = false);

  virtual ~OnDemandDataWarehouse();

  virtual bool
  exists(const VarLabel*, int matIndex, const Patch* patch) const override;

  virtual bool
  exists(const VarLabel*, int matIndex, const Level* level) const override;

  // Returns a (const) pointer to the grid.  This pointer can then be
  // used to (for example) get the number of levels in the grid.
  virtual const Grid*
  getGrid() override
  {
    return d_grid.get_rep();
  }

  virtual ReductionVariableBase*
  getReductionVariable(const VarLabel* label,
                       int matlIndex,
                       const Level* level) const;

  void
  copyKeyDB(KeyDatabase<Patch>& varkeyDB, KeyDatabase<Level>& levekeyDB);

  virtual void
  doReserve();

  // Generic put and allocate, passing Variable as a pointer rather than
  // by reference to avoid ambiguity with other put overloaded methods.
  virtual void
  put(Variable* var, const VarLabel* label, int matlIndex, const Patch* patch);

  // Reduction Variables
  virtual void
  get(ReductionVariableBase& var,
      const VarLabel* label,
      const Level* level = nullptr,
      int matIndex       = -1);

  // get a map filled with reductionVars
  // C++ doesn't allow templated Virtual method
  template<typename DataType, typename ReductionVarType>
  std::map<int, DataType>
  getReductionVar(const VarLabel* label, const MaterialSubset* matls);

  template<class T>
  std::map<int, T>
  get_sum_vartypeT(const VarLabel* label, const MaterialSubset* matls);

  // double
  virtual std::map<int, double>
  get_sum_vartypeD(const VarLabel* label, const MaterialSubset* matls) override;

  // Vector
  virtual std::map<int, Vector>
  get_sum_vartypeV(const VarLabel* label, const MaterialSubset* matls) override;

  virtual void
  put(const ReductionVariableBase& var,
      const VarLabel* label,
      const Level* level = nullptr,
      int matIndex       = -1);

  // put a map filled with reductionVars
  // C++ doesn't allow templated Virtual methods
  template<typename DataType, typename ReductionVarType>
  void
  putReductionVar(std::map<int, DataType>& reductionVars,
                  const VarLabel* label,
                  const MaterialSubset* matls);

  template<class T>
  void
  put_sum_vartypeT(std::map<int, T>& reductionVars,
                   const VarLabel* label,
                   const MaterialSubset* matls);

  // Vector
  virtual void
  put_sum_vartype(std::map<int, Vector>& reductionVars,
                  const VarLabel* label,
                  const MaterialSubset* matls) override;

  // double
  virtual void
  put_sum_vartype(std::map<int, double>& reductionVars,
                  const VarLabel* label,
                  const MaterialSubset* matls) override;

  virtual void
  override(const ReductionVariableBase&,
           const VarLabel* label,
           const Level* level = nullptr,
           int matIndex       = -1);

  virtual void
  print(std::ostream& intout,
        const VarLabel* label,
        const Level* level,
        int matlIndex = -1);

  // Sole Variables
  virtual bool
  exists(const VarLabel* label) const;

  virtual void
  get(SoleVariableBase&,
      const VarLabel* label,
      const Level* level = nullptr,
      int matIndex       = -1);

  virtual void
  put(const SoleVariableBase&,
      const VarLabel* label,
      const Level* level = nullptr,
      int matIndex       = -1);

  virtual void
  override(const SoleVariableBase&,
           const VarLabel* label,
           const Level* level = nullptr,
           int matIndex       = -1);

  //__________________________________
  // Particle Variables
  virtual ParticleSubset*
  createParticleSubset(particleIndex numParticles,
                       int matlIndex,
                       const Patch* patch,
                       IntVector low  = IntVector(0, 0, 0),
                       IntVector high = IntVector(0, 0, 0));

  virtual void
  deleteParticleSubset(ParticleSubset* pset);

  virtual void
  saveParticleSubset(ParticleSubset*,
                     int matlIndex,
                     const Patch* patch,
                     IntVector low  = IntVector(0, 0, 0),
                     IntVector high = IntVector(0, 0, 0));

  virtual bool
  haveParticleSubset(int matlIndex,
                     const Patch* patch,
                     IntVector low  = IntVector(0, 0, 0),
                     IntVector high = IntVector(0, 0, 0),
                     bool exact     = false);

  virtual ParticleSubset*
  getParticleSubset(int matlIndex,
                    const Patch* patch,
                    IntVector low,
                    IntVector high);

  virtual ParticleSubset*
  getParticleSubset(int matlIndex,
                    const Patch* patch,
                    IntVector low,
                    IntVector high,
                    const VarLabel* pos_label);

  virtual ParticleSubset*
  getParticleSubset(int matlIndex, const Patch* patch);

  virtual ParticleSubset*
  getDeleteSubset(int matlIndex, const Patch* patch);

  virtual ParticleLabelVariableMap*
  getNewParticleState(int matlIndex, const Patch* patch);

  virtual ParticleSubset*
  getParticleSubset(int matlIndex,
                    const Patch* patch,
                    Ghost::GhostType,
                    int numGhostCells,
                    const VarLabel* pos_label);

  // returns the particle subset in the range of low->high
  // relPatch is used as the key and should be the patch you are querying from
  // level is used if you are querying from an old level
  virtual ParticleSubset*
  getParticleSubset(int matlIndex,
                    IntVector low,
                    IntVector high,
                    const Patch* relPatch,
                    const VarLabel* posvar,
                    const Level* level = nullptr);

  // Get the particle index values of a set of ParticleIDs
  void
  getParticleIndex(ParticleSubset* pset,
                   const VarLabel* partIDLabel,
                   const std::vector<long64>& partIDList,
                   std::vector<particleIndex>& partIndexList);

  // Create a map between the long64 particleIDs and the particle indices in a
  // ParticleSubset
  typedef std::map<long64, int> ParticleIDMap;
  void
  createParticleIDMap(ParticleSubset* pset,
                      const VarLabel* partIDLabel,
                      ParticleIDMap& partIDMap);

  // Get the particle index value of a ParticleID after the partIDMap
  // has been created
  void
  getParticleIndex(const ParticleIDMap& partIDMap,
                   const long64& pParticleID,
                   particleIndex& pParticleIndex);

  virtual void
  allocateTemporary(ParticleVariableBase& var, ParticleSubset* pset);

  virtual void
  allocateAndPut(ParticleVariableBase& var,
                 const VarLabel* label,
                 ParticleSubset* pset);

  virtual void
  get(constParticleVariableBase& var,
      const VarLabel* label,
      ParticleSubset* pset);

  virtual void
  get(constParticleVariableBase& var,
      const VarLabel* label,
      int matlIndex,
      const Patch* patch);

  virtual void
  getModifiable(ParticleVariableBase& var,
                const VarLabel* label,
                ParticleSubset* pset);

  virtual void
  put(ParticleVariableBase& var, const VarLabel* label, bool replace = false);

  virtual ParticleVariableBase*
  getParticleVariable(const VarLabel* label, ParticleSubset* pset);

  virtual ParticleVariableBase*
  getParticleVariable(const VarLabel* label, int matlIndex, const Patch* patch);

  void
  printParticleSubsets();

  virtual void
  getCopy(ParticleVariableBase& var,
          const VarLabel* label,
          ParticleSubset* pset);

  virtual void
  copyOut(ParticleVariableBase& var,
          const VarLabel* label,
          ParticleSubset* pset);

  // Remove particles that are no longer relevant
  virtual void
  deleteParticles(ParticleSubset* delset);

  virtual void
  addParticles(const Patch* patch,
               int matlIndex,
               ParticleLabelVariableMap* addedstate);

  //__________________________________
  // Grid Variables
  virtual void
  print();

  virtual void
  clear();

  void
  get(constGridVariableBase& var,
      const VarLabel* label,
      int matlIndex,
      const Patch* patch,
      Ghost::GhostType gtype,
      int numGhostCells);

  void
  getModifiable(GridVariableBase& var,
                const VarLabel* label,
                int matlIndex,
                const Patch* patch,
                Ghost::GhostType gtype = Ghost::None,
                int numGhostCells      = 0);

  void
  allocateTemporary(GridVariableBase& var,
                    const Patch* patch,
                    Ghost::GhostType gtype,
                    int numGhostCells);

  void
  allocateAndPut(GridVariableBase& var,
                 const VarLabel* label,
                 int matlIndex,
                 const Patch* patch,
                 Ghost::GhostType gtype,
                 int numGhostCells);

  void
  put(GridVariableBase& var,
      const VarLabel* label,
      int matlIndex,
      const Patch* patch,
      bool replace = false);

  // returns the constGridVariable for all patches on the level
  virtual void
  getLevel(constGridVariableBase& var,
           const VarLabel* label,
           int matlIndex,
           const Level* level);

  // meant for the UnifiedScheduler only.  Puts a contiguous level in the
  // *level* database so that it doesn't need to constantly make new deep copies
  // each time the full level is requested.
  void
  putLevelDB(GridVariableBase* gridVar,
             const VarLabel* label,
             const Level* level,
             int matlIndex = -1);

  virtual void
  getRegion(constGridVariableBase& var,
            const VarLabel* label,
            int matlIndex,
            const Level* level,
            const IntVector& low,
            const IntVector& high,
            bool useBoundaryCells = true);

  virtual void
  getRegionModifiable(GridVariableBase& var,
                      const VarLabel* label,
                      int matlIndex,
                      const Level* level,
                      const IntVector& low,
                      const IntVector& high,
                      bool useBoundaryCells = true);

  virtual void
  copyOut(GridVariableBase& var,
          const VarLabel* label,
          int matlIndex,
          const Patch* patch,
          Ghost::GhostType gtype = Ghost::None,
          int numGhostCells      = 0);

  virtual void
  getCopy(GridVariableBase& var,
          const VarLabel* label,
          int matlIndex,
          const Patch* patch,
          Ghost::GhostType gtype = Ghost::None,
          int numGhostCells      = 0);

  // PerPatch Variables
  virtual void
  get(PerPatchBase& var,
      const VarLabel* label,
      int matIndex,
      const Patch* patch);

  virtual void
  put(PerPatchBase& var,
      const VarLabel* label,
      int matIndex,
      const Patch* patch,
      bool replace = false);

  virtual ScrubMode setScrubbing(ScrubMode);

  // For related datawarehouses
  virtual DataWarehouse* getOtherDataWarehouse(Task::WhichDW);

  virtual void
  transferFrom(DataWarehouse* from,
               const VarLabel* label,
               const PatchSubset* patches,
               const MaterialSubset* matls);

  virtual void
  transferFrom(DataWarehouse* from,
               const VarLabel* label,
               const PatchSubset* patches,
               const MaterialSubset* matls,
               bool replace);

  virtual void
  transferFrom(DataWarehouse* from,
               const VarLabel* label,
               const PatchSubset* patches,
               const MaterialSubset* matls,
               bool replace,
               const PatchSubset* newPatches);

  virtual void
  transferFrom(DataWarehouse* from,
               const VarLabel* label,
               const PatchSubset* patches,
               const MaterialSubset* matls,
               void* dtask,
               bool replace,
               const PatchSubset* newPatches);

  virtual bool
  isFinalized() const;

  virtual void
  finalize();

#ifdef HAVE_CUDA
  static int
  getNumDevices();

  static void
  uintahSetCudaDevice(int deviceNum);

  static size_t
  getTypeDescriptionSize(const TypeDescription::Type& type);

  static GPUGridVariableBase*
  createGPUGridVariable(const TypeDescription::Type& type);

  static GPUPerPatchBase*
  createGPUPerPatch(const TypeDescription::Type& type);

  static GPUReductionVariableBase*
  createGPUReductionVariable(const TypeDescription::Type& type);
#endif

  virtual void
  unfinalize();

  virtual void
  refinalize();

  virtual size_t
  emit(OutputContext& out,
       const VarLabel* label,
       int matlIndex,
       const Patch* patch);

#if HAVE_PIDX
  void
  emitPIDX(PIDXOutputContext& context,
           const VarLabel* label,
           int matlIndex,
           const Patch* patch,
           unsigned char* pidx_buffer,
           size_t pidx_bufferSize);
#endif

  void
  exchangeParticleQuantities(DetailedTasks* dts,
                             LoadBalancer* lb,
                             const VarLabel* pos_var,
                             int iteration);

  void
  sendMPI(DependencyBatch* batch,
          const VarLabel* pos_var,
          BufferInfo& buffer,
          OnDemandDataWarehouse* old_dw,
          const DetailedDep* dep,
          LoadBalancer* lb);

  void
  recvMPI(DependencyBatch* batch,
          BufferInfo& buffer,
          OnDemandDataWarehouse* old_dw,
          const DetailedDep* dep,
          LoadBalancer* lb);

  void
  reduceMPI(const VarLabel* label,
            const Level* level,
            const MaterialSubset* matls,
            int nComm);

  // Scrub counter manipulator functions -- when the scrub count goes to
  // zero, the data is deleted
  void
  setScrubCount(const VarLabel* label,
                int matlIndex,
                const Patch* patch,
                int count);

  int
  decrementScrubCount(const VarLabel* label, int matlIndex, const Patch* patch);

  void
  scrub(const VarLabel* label, int matlIndex, const Patch* patch);

  void
  initializeScrubs(int dwid,
                   const FastHashTable<ScrubItem>* scrubcounts,
                   bool add);

  // For timestep abort/restart
  virtual bool
  abortTimeStep();

  virtual bool
  recomputeTimeStep();

  struct ValidNeighbors
  {
    GridVariableBase* validNeighbor{ nullptr };
    const Patch* neighborPatch{ nullptr };
    IntVector low{};
    IntVector high{};
  };

  void
  getNeighborPatches(const VarLabel* label,
                     const Patch* patch,
                     Ghost::GhostType gtype,
                     int numGhostCells,
                     std::vector<const Patch*>& adjacentNeighbors);

  // This method will retrieve those neighbors, and also the regions (indicated
  // in low and high) which constitute the ghost cells. Data is return in the
  // ValidNeighbors vector. IgnoreMissingNeighbors is designed for the Unified
  // Scheduler so that it can request what neighbor patches *should* be, and
  // those neighbor patches we hope are found in the host side DW (this one) or
  // the GPU DW
  //
  // TODO, This method might create a reference to the neighbor, and so these
  // references need to be deleted afterward. (It's not pretty, but it seemed to
  // be the best option.)
  void
  getValidNeighbors(const VarLabel* label,
                    int matlIndex,
                    const Patch* patch,
                    Ghost::GhostType gtype,
                    int numGhostCells,
                    std::vector<ValidNeighbors>& validNeighbors,
                    bool ignoreMissingNeighbors = false);

  void
  logMemoryUse(std::ostream& out, unsigned long& total, const std::string& tag);

  // must be called by the thread that will run the test
  void
  pushRunningTask(const Task* task, std::vector<OnDemandDataWarehouseUP>* dws);

  void
  popRunningTask();

  // does a final check to see if gets/puts/etc. consistent with
  // requires/computes/modifies for the current task.
  void
  checkTasksAccesses(const PatchSubset* patches, const MaterialSubset* matls);

  ScrubMode
  getScrubMode() const
  {
    return d_scrub_mode;
  }

  // The following is for support of regriding
  virtual void
  getVarLabelMatlLevelTriples(std::vector<VarLabelMatl<Level>>& vars) const;

  inline static bool s_combine_memory = false;

  friend class SchedulerCommon;
  friend class UnifiedScheduler;

private:
  enum AccessType
  {
    NoAccess = 0,
    PutAccess,
    GetAccess,
    ModifyAccess
  };

  struct AccessInfo
  {

    AccessInfo() = default;

    AccessInfo(AccessType type)
      : accessType(type)
    {
    }

    void
    encompassOffsets(IntVector low, IntVector high)
    {
      lowOffset  = Max(low, lowOffset);
      highOffset = Max(high, highOffset);
    }

    AccessType accessType{ NoAccess };
    IntVector lowOffset{ 0, 0, 0 }; // ghost cell access
    IntVector highOffset{ 0, 0, 0 };
  };

  using VarAccessMap = std::map<VarLabelMatl<Patch>, AccessInfo>;

  struct RunningTaskInfo
  {

    RunningTaskInfo() = default;

    RunningTaskInfo(const Task* task, std::vector<OnDemandDataWarehouseUP>* dws)
      : d_task(task)
      , d_dws(dws)
    {
    }

    RunningTaskInfo(const RunningTaskInfo& copy)
      : d_task(copy.d_task)
      , d_dws(copy.d_dws)
      , d_accesses(copy.d_accesses)
    {
    }

    RunningTaskInfo&
    operator=(const RunningTaskInfo& copy)
    {
      d_task     = copy.d_task;
      d_dws      = copy.d_dws;
      d_accesses = copy.d_accesses;
      return *this;
    }

    const Task* d_task{ nullptr };
    std::vector<OnDemandDataWarehouseUP>* d_dws{ nullptr };
    VarAccessMap d_accesses{};
  };

  virtual DataWarehouse*
  getOtherDataWarehouse(Task::WhichDW, RunningTaskInfo* info);

  void
  getGridVar(GridVariableBase& var,
             const VarLabel* label,
             int matlIndex,
             const Patch* patch,
             Ghost::GhostType gtype,
             int numGhostCells);

  inline Task::WhichDW
  getWhichDW(RunningTaskInfo* info);

  // These will throw an exception if access is not allowed for the
  // curent task.
  inline void
  checkGetAccess(const VarLabel* label,
                 int matlIndex,
                 const Patch* patch,
                 Ghost::GhostType gtype = Ghost::None,
                 int numGhostCells      = 0);

  inline void
  checkPutAccess(const VarLabel* label,
                 int matlIndex,
                 const Patch* patch,
                 bool replace);

  inline void
  checkModifyAccess(const VarLabel* label, int matlIndex, const Patch* patch);

  // These will return false if access is not allowed for
  // the current task.
  inline bool
  hasGetAccess(const Task* runningTask,
               const VarLabel* label,
               int matlIndex,
               const Patch* patch,
               IntVector lowOffset,
               IntVector highOffset,
               RunningTaskInfo* info);

  inline bool
  hasPutAccess(const Task* runningTask,
               const VarLabel* label,
               int matlIndex,
               const Patch* patch,
               bool replace);

  void
  checkAccesses(RunningTaskInfo* runningTaskInfo,
                const Task::Dependency* dep,
                AccessType accessType,
                const PatchSubset* patches,
                const MaterialSubset* matls);

  void
  printDebuggingPutInfo(const VarLabel* label,
                        int matlIndex,
                        const Patch* patch,
                        int line);

  void
  printDebuggingPutInfo(const VarLabel* label,
                        int matlIndex,
                        const Level* level,
                        int line);

#ifdef HAVE_CUDA
  std::map<Patch*, bool> assignedPatches; // indicates where a given patch
                                          // should be stored in an accelerator
#endif

  using psetDBType = std::multimap<PSPatchMatlGhost, ParticleSubset*>;
  using psetAddDBType =
    std::map<std::pair<int, const Patch*>,
             std::map<const VarLabel*, ParticleVariableBase*>*>;
  using particleQuantityType = std::map<std::pair<int, const Patch*>, int>;

  ParticleSubset*
  queryPSetDB(psetDBType& db,
              const Patch* patch,
              int matlIndex,
              IntVector low,
              IntVector high,
              const VarLabel* pos_var,
              bool exact = false);

  void
  insertPSetRecord(psetDBType& subsetDB,
                   const Patch* patch,
                   IntVector low,
                   IntVector high,
                   int matlIndex,
                   ParticleSubset* psubset);

  void
  deletePSetRecord(psetDBType& subsetDB,
                   const Patch* patch,
                   IntVector low,
                   IntVector high,
                   int matlIndex,
                   ParticleSubset* psubset);

  inline bool
  hasRunningTask();

  inline std::map<std::thread::id, OnDemandDataWarehouse::RunningTaskInfo>*
  getRunningTasksInfo();

  inline RunningTaskInfo*
  getCurrentTaskInfo();

private:
  // holds info on the running task for each std::thread::id
  std::map<std::thread::id, RunningTaskInfo> d_running_tasks{};

  ScrubMode d_scrub_mode{ DataWarehouse::ScrubNone };

  DWDatabase<Patch> d_var_DB{};
  DWDatabase<Level> d_level_DB{};
  KeyDatabase<Patch> d_var_key_DB{};
  KeyDatabase<Level> d_level_key_DB{};

  psetDBType d_pset_DB{};
  psetDBType d_delset_DB{};
  psetAddDBType d_addset_DB{};
  particleQuantityType d_foreign_particle_quantities{};
  bool d_exchange_particle_quantities{ true };

  // Keep track of when this DW sent some (and which) particle information to
  // another processor
  SendState d_send_state;

  bool d_finalized{ false };
  GridP d_grid{ nullptr };

  // Is this the first DW -- created by the initialization timestep?
  bool d_is_initialization_DW{ false };
};

} // end namespace Uintah

#endif
