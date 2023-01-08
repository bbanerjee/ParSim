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

#include <Core/Thread/CrowdMonitor.h>

#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Components/Schedulers/OnDemandDataWarehouseP.h>
#include <CCA/Components/Schedulers/DWDatabase.h>
#include <CCA/Components/Schedulers/SendState.h>
#include <Core/Grid/Variables/VarLabelMatl.h>
#include <Core/Grid/Variables/PSPatchMatlGhost.h>
#include <Core/Grid/Grid.h>

#include <map>
#include <iosfwd>
#include <vector>
#include <sci_defs/mpi_defs.h> // For MPIPP_H on SGI

namespace Uintah {
  class Thread;
}

namespace Uintah {

  inline const Patch* getRealDomain(const Patch* patch) { return patch->getRealPatch(); }
  inline const Level* getRealDomain(const Level* level) { return level; }

  using Uintah::CrowdMonitor;
  using Uintah::Max;
  using Uintah::Thread;
  using Uintah::FastHashTable;

  class BufferInfo;
  class DependencyBatch;
  class DetailedTasks;
  class DetailedDep;
  class TypeDescription;
  class Patch;
  class ProcessorGroup;
  class SendState;
  class LoadBalancer;

/**************************************

  CLASS
        OnDemandDataWarehouse
   
Short description...

  GENERAL INFORMATION

        OnDemandDataWarehouse.h

Steven G. Parker
Department of Computer Science
University of Utah

Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
  

  KEYWORDS
        On_Demand_Data_Warehouse

  DESCRIPTION
        Long description...
  
  WARNING
  
****************************************/

  class OnDemandDataWarehouse : public DataWarehouse {

  public:

    OnDemandDataWarehouse(const ProcessorGroup* myworld,
                          Scheduler* scheduler, 
                          int generation,
                          const GridP& grid,
                          bool isInitializationDW = false);

    virtual ~OnDemandDataWarehouse();
   
    virtual bool exists(const VarLabel*,
                        int matIndex, 
                        const Patch*) const; 
   
    void copyKeyDB(KeyDatabase<Patch>& varkeyDB,
                   KeyDatabase<Level>& levekeyDB);

    virtual void doReserve();

    // Returns a (const) pointer to the grid.  This pointer can then be
    // used to (for example) get the number of levels in the grid.
    virtual const Grid * getGrid()
    {
      return d_grid.get_rep();
    }

    // Generic put and allocate, passing Variable as a pointer rather than
    // by reference to avoid ambiguity with other put overloaded methods.
    virtual void put(Variable*,
                     const VarLabel*,
                     int matlIndex,
                     const Patch*);
   
    // Reduction Variables
    virtual void get(ReductionVariableBase&,
                     const VarLabel*,
                     const Level* level = 0, 
                     int matIndex = -1);

    virtual void put(const ReductionVariableBase&,
                     const VarLabel*,
                     const Level* level = 0, 
                     int matIndex = -1);

    virtual void override(const ReductionVariableBase&,
                          const VarLabel*,
                          const Level* level = 0, 
                          int matIndex = -1);

    virtual void print(std::ostream& intout,
                       const VarLabel* label,
                       const Level* level, 
                       int matlIndex = -1);

    // Sole Variables
    virtual bool exists(const VarLabel*) const;

    virtual void get(SoleVariableBase&,
                     const VarLabel*,
                     const Level* level = 0, 
                     int matIndex = -1);

    virtual void put(const SoleVariableBase&,
                     const VarLabel*,
                     const Level* level = 0, 
                     int matIndex = -1);

    virtual void override(const SoleVariableBase&,
                          const VarLabel*,
                          const Level* level = 0, 
                          int matIndex = -1);

    //__________________________________
    // Particle Variables
   
    virtual ParticleSubset* createParticleSubset(particleIndex numParticles,
                                                 int matlIndex, 
                                                 const Patch*,
                                                 IntVector low = IntVector(0, 0, 0),
                                                 IntVector high = IntVector(0, 0, 0));
                                               
    virtual void saveParticleSubset(ParticleSubset*,
                                    int matlIndex, 
                                    const Patch*,
                                    IntVector low = IntVector(0, 0, 0),
                                    IntVector high = IntVector(0, 0, 0));
                                   
    virtual bool haveParticleSubset(int matlIndex, 
                                    const Patch*,
                                    IntVector low = IntVector(0, 0, 0),
                                    IntVector high = IntVector(0, 0, 0),
                                    bool exact = false);
                                   
    virtual ParticleSubset* getParticleSubset(int matlIndex, 
                                              const Patch*,
                                              IntVector low, 
                                              IntVector high);
                                             
    virtual ParticleSubset* getParticleSubset(int matlIndex, 
                                              const Patch*,
                                              IntVector low, 
                                              IntVector high, 
                                              const VarLabel*);
                                             
    virtual ParticleSubset* getParticleSubset(int matlIndex, 
                                              const Patch*);
                                             
    virtual ParticleSubset* getDeleteSubset(int matlIndex, 
                                            const Patch*);
                                          
    virtual ParticleLabelVariableMap*  getNewParticleState(int matlIndex, 
                                                           const Patch* patch);
   
    virtual ParticleSubset* getParticleSubset(int matlIndex,
                                              const Patch* patch, 
                                              Ghost::GhostType, 
                                              int numGhostCells,
                                              const VarLabel* posvar);
                                 
    /* Create a particle subset for a subset of a patch and its
       neighboring patches defined by a local lowIndex and a local highIndex.  
       If the particles are contained outside the current patch, use the
       numGhostCells to get the outside particles */
    virtual ParticleSubset* getParticleSubset(int matlIndex,
                                              const Patch* patch, 
                                              IntVector localLowIndex,
                                              IntVector localHighIndex,
                                              Ghost::GhostType, 
                                              int numGhostCells,
                                              const VarLabel* posvar);

    //returns the particle subset in the range of low->high
    //relPatch is used as the key and should be the patch you are querying from
    //level is used if you are querying from an old level
    virtual ParticleSubset* getParticleSubset(int matlIndex, 
                                              IntVector low, 
                                              IntVector high, 
                                              const Patch* relPatch,
                                              const VarLabel* posvar,
                                              const Level* level=0);

    // Get the particle index values of a set of ParticleIDs
    void getParticleIndex(ParticleSubset* pset,
                          const VarLabel* partIDLabel,
                          const std::vector<long64>& partIDList,
                          std::vector<particleIndex>& partIndexList);

    // Create a map between the long64 particleIDs and the particle indices in a 
    // ParticleSubset
    typedef std::map<long64, int> ParticleIDMap;
    void createParticleIDMap(ParticleSubset* pset,
                             const VarLabel* partIDLabel,
                             ParticleIDMap& partIDMap);
                                             
    // Get the particle index value of a ParticleID after the partIDMap 
    // has been created
    void getParticleIndex(const ParticleIDMap& partIDMap,
                          const long64& pParticleID,
                          particleIndex& pParticleIndex);
 
    virtual void allocateTemporary(ParticleVariableBase& var, 
                                   ParticleSubset* label);
                                  
    virtual void allocateAndPut(ParticleVariableBase& var, 
                                const VarLabel* label,
                                ParticleSubset* pset);
                               
    virtual void get(constParticleVariableBase& var, 
                     const VarLabel* label,
                     ParticleSubset* pset);
                    
    virtual void get(constParticleVariableBase& var, 
                     const VarLabel* label,
                     int matlIndex, 
                     const Patch* patch);

    virtual void getModifiable(ParticleVariableBase& var, 
                               const VarLabel* label,
                               ParticleSubset* pset);
                              
    virtual void put(ParticleVariableBase& var, 
                     const VarLabel* label,
                     bool replace = false);
                    
    virtual ParticleVariableBase* getParticleVariable(const VarLabel* label,
                                                      ParticleSubset* pset);
                                                     
    virtual ParticleVariableBase* getParticleVariable(const VarLabel* label, 
                                                      int matlIndex, 
                                                      const Patch* patch);  

    void printParticleSubsets();

    virtual void getCopy(ParticleVariableBase& var, 
                         const VarLabel* label, 
                         ParticleSubset* pset);
                        
    virtual void copyOut(ParticleVariableBase& var, 
                         const VarLabel* label, 
                         ParticleSubset* pset);

    // Remove particles that are no longer relevant
    virtual void deleteParticles(ParticleSubset* delset);

    virtual void addParticles(const Patch* patch, 
                              int matlIndex, 
                              ParticleLabelVariableMap* addedstate);

    //__________________________________
    // Grid Variables
    virtual void print();
    virtual void clear();

    void get(constGridVariableBase& var, 
             const VarLabel* label, 
             int matlIndex, 
             const Patch* patch, 
             Ghost::GhostType gtype, 
             int numGhostCells);

    void getModifiable(GridVariableBase& var, 
                       const VarLabel* label, 
                       int matlIndex, 
                       const Patch* patch, 
                       Ghost::GhostType gtype=Ghost::None, 
                       int numGhostCells=0);

    void allocateTemporary(GridVariableBase& var, 
                           const Patch* patch,
                           Ghost::GhostType gtype, 
                           int numGhostCells );
    // const IntVector& boundaryLayer);

    void allocateAndPut(GridVariableBase& var, 
                        const VarLabel* label, 
                        int matlIndex,
                        const Patch* patch, 
                        Ghost::GhostType gtype, 
                        int numGhostCells);

    void put(GridVariableBase& var, 
             const VarLabel* label, 
             int matlIndex, 
             const Patch* patch,
             bool replace = false);

    // returns the constGridVariable for all patches on the level
    virtual void getLevel(constGridVariableBase&,
                          const VarLabel*,
                          int matlIndex,
                          const Level* level);

    virtual void getRegion(constGridVariableBase& var, 
                           const VarLabel* label,
                           int matlIndex, 
                           const Level* level,
                           const IntVector& low, 
                           const IntVector& high,
                           bool useBoundaryCells = true);

    virtual void getRegion(GridVariableBase&,
                           const VarLabel*,
                           int matlIndex,
                           const Level* level,
                           const IntVector& low,
                           const IntVector& high,
                           bool useBoundaryCells = true);

    virtual void copyOut(GridVariableBase& var, 
                         const VarLabel* label, int matlIndex,
                         const Patch* patch, 
                         Ghost::GhostType gtype = Ghost::None,
                         int numGhostCells = 0);

    virtual void getCopy(GridVariableBase& var, 
                         const VarLabel* label, 
                         int matlIndex,
                         const Patch* patch, 
                         Ghost::GhostType gtype = Ghost::None,
                         int numGhostCells = 0);

    // PerPatch Variables
    virtual void get(PerPatchBase& var, 
                     const VarLabel* label,
                     int matIndex, 
                     const Patch* patch);

    virtual void put(PerPatchBase& var, 
                     const VarLabel* label,
                     int matIndex, 
                     const Patch* patch, 
                     bool replace = false);

    virtual ScrubMode setScrubbing(ScrubMode);

    // For related datawarehouses
    virtual DataWarehouse* getOtherDataWarehouse(Task::WhichDW);

    //! Copy a var from the parameter DW to this one.  If newPatches
    //! is not null, then it associates the copy of the variable with
    //! newPatches, and otherwise it uses patches (the same it finds
    //! the variable with.
    virtual void transferFrom(DataWarehouse* dw, 
                              const VarLabel* label,
                              const PatchSubset* patches, 
                              const MaterialSubset* mats,
                              bool replace = false,
                              const PatchSubset* newPatches = 0);

    virtual bool isFinalized() const;
   
    virtual void finalize();
    virtual void unfinalize();
    virtual void refinalize();


    virtual void emit(OutputContext& out, 
                      const VarLabel* label,
                      int matlIndex, 
                      const Patch* patch);

    void exchangeParticleQuantities(DetailedTasks* dts, 
                                    LoadBalancer* lb, 
                                    const VarLabel* pos_var, 
                                    int iteration);

    void sendMPI(DependencyBatch* batch, 
                 const VarLabel* pos_var, 
                 BufferInfo& buffer, 
                 OnDemandDataWarehouse* old_dw, 
                 const DetailedDep* dep, 
                 LoadBalancer* lb);

    void recvMPI(DependencyBatch* batch, 
                 BufferInfo& buffer,
                 OnDemandDataWarehouse* old_dw, 
                 const DetailedDep* dep, 
                 LoadBalancer* lb);

    void reduceMPI(const VarLabel* label, 
                   const Level* level,
                   const MaterialSubset* matls, 
                   int nComm);

    // Scrub counter manipulator functions -- when the scrub count goes to
    // zero, the data is deleted
    void setScrubCount(const VarLabel* label, 
                       int matlIndex,
                       const Patch* patch, 
                       int count);

    int decrementScrubCount(const VarLabel* label, 
                            int matlIndex,
                            const Patch* patch);

    void scrub(const VarLabel* label, 
               int matlIndex, 
               const Patch* patch);

    void initializeScrubs(int dwid, 
                          const FastHashTable<ScrubItem>* scrubcounts, 
                          bool add);

    // For timestep abort/restart
    virtual bool timestepAborted();
    virtual bool timestepRestarted();
    virtual void abortTimestep();
    virtual void restartTimestep();
    virtual void setRestarted() { hasRestarted_ = true; }

    void logMemoryUse(ostream& out, unsigned long& total, const std::string& tag);

    // must be called by the thread that will run the test
    void pushRunningTask(const Task* task, std::vector<OnDemandDataWarehouseP>* dws);
    void popRunningTask();  

    // does a final check to see if gets/puts/etc. consistent with
    // requires/computes/modifies for the current task.
    void checkTasksAccesses(const PatchSubset* patches,
                            const MaterialSubset* matls);

    ScrubMode getScrubMode() const {
      return d_scrubMode;
    }

    // The following is for support of regriding
    virtual void getVarLabelMatlLevelTriples( std::vector<VarLabelMatl<Level> >& vars ) const;

    static bool d_combineMemory;

    friend class SchedulerCommon;
    friend class UnifiedScheduler;

  private:

    enum AccessType {
      NoAccess = 0, 
      PutAccess, 
      GetAccess, 
      ModifyAccess
    };

    struct AccessInfo {

      AccessInfo()
        : accessType(NoAccess), lowOffset(0, 0, 0), highOffset(0, 0, 0) {}

      AccessInfo(AccessType type)
        : accessType(type), lowOffset(0, 0, 0), highOffset(0, 0, 0) {}

      void encompassOffsets(IntVector low, IntVector high)
        { lowOffset = Max(low, lowOffset); highOffset = Max(high, highOffset); }
     
      AccessType accessType;
      IntVector lowOffset; // ghost cell access
      IntVector highOffset;
    };
  
    typedef std::map<VarLabelMatl<Patch>, AccessInfo> VarAccessMap;

    struct RunningTaskInfo {

      RunningTaskInfo()
        : d_task(0), dws(0) {}

      RunningTaskInfo(const Task* task, std::vector<OnDemandDataWarehouseP>* dws)
        : d_task(task), dws(dws) {}

      RunningTaskInfo(const RunningTaskInfo& copy)
        : d_task(copy.d_task), dws(copy.dws), d_accesses(copy.d_accesses) {}

      RunningTaskInfo& operator=(const RunningTaskInfo& copy)
        { 
          d_task = copy.d_task; 
          dws = copy.dws; 
          d_accesses = copy.d_accesses; 
          return *this; 
        }

      const Task* d_task;
      std::vector<OnDemandDataWarehouseP>* dws;
      VarAccessMap d_accesses;
    };  
   
    virtual DataWarehouse* getOtherDataWarehouse(Task::WhichDW,
                                                 RunningTaskInfo *info);

    void getGridVar(GridVariableBase& var, 
                    const VarLabel* label, 
                    int matlIndex, 
                    const Patch* patch,
                    Ghost::GhostType gtype, 
                    int numGhostCells);

    inline Task::WhichDW getWhichDW( RunningTaskInfo *info);

    // These will throw an exception if access is not allowed for the
    // curent task.
    inline void checkGetAccess(const VarLabel* label, 
                               int matlIndex,
                               const Patch* patch,
                               Ghost::GhostType gtype = Ghost::None,
                               int numGhostCells = 0);

    inline void checkPutAccess(const VarLabel* label, 
                               int matlIndex,
                               const Patch* patch, 
                               bool replace);

    inline void checkModifyAccess(const VarLabel* label, 
                                  int matlIndex,
                                  const Patch* patch);
  
    // These will return false if access is not allowed for
    // the current task.
    inline bool hasGetAccess(const Task* runningTask, 
                             const VarLabel* label,
                             int matlIndex, 
                             const Patch* patch,
                             IntVector lowOffset, 
                             IntVector highOffset,
                             RunningTaskInfo *info);

    inline bool hasPutAccess(const Task* runningTask, 
                             const VarLabel* label,
                             int matlIndex, 
                             const Patch* patch, 
                             bool replace);

    void checkAccesses(RunningTaskInfo* runningTaskInfo,
                       const Task::Dependency* dep,
                       AccessType accessType, 
                       const PatchSubset* patches,
                       const MaterialSubset* matls);
  
    void printDebuggingPutInfo( const VarLabel* label,
                                int matlIndex,
                                const Patch* patch,
                                int line);
                       
                       

    struct dataLocation {
      const Patch   * patch;
      int        mpiNode;
    };

    typedef std::vector<dataLocation*> variableListType;
    typedef std::map<const VarLabel*, variableListType*, VarLabel::Compare> dataLocationDBtype;
    typedef std::multimap<PSPatchMatlGhost, ParticleSubset*> psetDBType;
    typedef std::map<std::pair<int, const Patch*>, ParticleLabelVariableMap* > psetAddDBType;
    typedef std::map<std::pair<int, const Patch*>, int> particleQuantityType;
   
    ParticleSubset* queryPSetDB(psetDBType &db, 
                                const Patch* patch, 
                                int matlIndex, 
                                IntVector low, 
                                IntVector high, 
                                const VarLabel* pos_var, 
                                bool exact=false);

    void insertPSetRecord(psetDBType &subsetDB,
                          const Patch* patch, 
                          IntVector low, 
                          IntVector high, 
                          int matlIndex, 
                          ParticleSubset *psubset);

    DWDatabase<Patch>  d_varDB;
    DWDatabase<Level>  d_levelDB;
    KeyDatabase<Patch> d_varkeyDB;
    KeyDatabase<Level> d_levelkeyDB;

    psetDBType                        d_psetDB;
    psetDBType                        d_delsetDB;
    psetAddDBType d_addsetDB;
    particleQuantityType d_foreignParticleQuantities;

    // Keep track of when this DW sent some (and which) particle information to another processor
    SendState  ss_;

    /*
    // On a timestep restart, sometimes (when an entire patch is sent) on the
    // first try of the timestep the receiving DW creates and stores ParticleSubset
    // which throws off the sending on the next iteration.  This will compensate.
    SendState                         d_timestepRestartPsets;
    */

    // Record of which DataWarehouse has the data for each variable...
    //  Allows us to look up the DW to which we will send a data request.
    dataLocationDBtype      d_dataLocation;

    //////////
    // Insert Documentation Here:
    mutable CrowdMonitor    d_lock;
    mutable CrowdMonitor    d_lvlock;
    mutable CrowdMonitor    d_plock;
    mutable CrowdMonitor    d_pslock;
    bool                    d_finalized;
    GridP                   d_grid;

    // Is this the first DW -- created by the initialization timestep?
    bool d_isInitializationDW;
  
    inline bool hasRunningTask();

    inline std::list<RunningTaskInfo>* getRunningTasksInfo();

    inline RunningTaskInfo* getCurrentTaskInfo();

    //std::map<Thread*, std::list<RunningTaskInfo> > d_runningTasks;
    std::list<RunningTaskInfo>  d_runningTasks[MAX_THREADS];
    ScrubMode d_scrubMode;

    bool aborted;
    bool restart;

    // Whether this (Old) DW is being used for a restarted timestep (the new DWs are cleared out)
    bool hasRestarted_;

  };

} // end namespace Uintah

#endif
