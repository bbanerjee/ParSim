/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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

#include <Core/Parallel/UintahParallelComponent.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <CCA/Ports/Scheduler.h>
#include <CCA/Components/Schedulers/Relocate.h>
#include <CCA/Components/Schedulers/OnDemandDataWarehouseP.h>
#include <Core/Grid/SimulationState.h>
#include <Core/Grid/SimulationStateP.h>
#include <Core/Util/DebugStream.h>

#include <iosfwd>
#include <set>

namespace Uintah {

  class Output;
  class DetailedTask;
  class DetailedTasks;
  class TaskGraph;
  class LocallyComputedPatchVarMap;

  class SchedulerCommon : public Scheduler, public UintahParallelComponent {
  public:
    SchedulerCommon(const ProcessorGroup* myworld, const Output* oport);
    virtual ~SchedulerCommon();

    virtual void problemSetup(const ProblemSpecP& prob_spec,
                              SimulationStateP& state);

    virtual void doEmitTaskGraphDocs();

    virtual void checkMemoryUse( unsigned long & memuse, unsigned long & highwater,
                                 unsigned long & maxMemUse );
    // sbrk memory start location (for memory tracking)
    virtual void   setStartAddr( char * start ) { d_start_addr = start; }
    virtual char * getStartAddr() { return d_start_addr; }
    virtual void   resetMaxMemValue();

    // For calculating memory usage when sci-malloc is disabled...
    static char* d_start_addr;

    virtual void initialize(int numOldDW = 1, int numNewDW = 1);

    virtual void setParentDWs(DataWarehouse* parent_old_dw,
                              DataWarehouse* parent_new_dw);

    void setSimulationState(SimulationStateP sharedState) {
      d_sharedState = sharedState;
    }

    virtual void clearMappings();

    virtual void mapDataWarehouse(Task::WhichDW, int dwTag);

    void compile();

    //////////
    // Insert Documentation Here:

    /// For more complicated models 
    virtual void addTaskGraph(tgType type);
    virtual int getNumTaskGraphs() { return d_graphs.size(); }
    virtual void addTask(Task* t, const PatchSet*, const MaterialSet*);

    virtual bool useSmallMessages() { return d_useSmallMessages; }

    /// Get all of the requires needed from the old data warehouse
    /// (carried forward).
    virtual const std::vector<const Task::Dependency*>& getInitialRequires() const
    { return d_initRequires; }
    virtual const std::set<const VarLabel*, VarLabel::Compare>& getInitialRequiredVars() const
    { return d_initRequiredVars; }
    virtual const std::set<const VarLabel*, VarLabel::Compare>& getComputedVars() const
    { return d_computedVars; }

    virtual LoadBalancer* getLoadBalancer();

    virtual void releaseLoadBalancer();
       
    virtual DataWarehouse* get_dw( int idx );

    virtual DataWarehouse* getLastDW( void );

    virtual void logMemoryUse();
      
    virtual void advanceDataWarehouse( const GridP& grid, bool initialization=false );

    virtual void fillDataWarehouses( const GridP& grid );

    virtual void replaceDataWarehouse( int index, const GridP& grid, bool initialization=false );

    virtual void setRestartable( bool restartable );

    // Get the expected extents that may be needed for a particular variable
    // on a particular patch (which should include expected ghost cells.
    //virtual void
    //getExpectedExtents(const VarLabel* label, const Patch* patch,
    //       IntVector& lowIndex, IntVector& highIndex) const;

    // Get the SuperPatch (set of connected patches making a larger rectangle)
    // for the given label and patch and find the largest extents encompassing
    // the expected ghost cells (requiredLow, requiredHigh) and the requested
    // ghost cells as well (requestedLow, requestedHigh) for each of the
    // patches.  Required and requested will besame if requestedNumGCells = 0.
    virtual const std::vector<const Patch*>*
      getSuperPatchExtents(const VarLabel* label, 
                           int matlIndex,
                           const Patch* patch, 
                           Ghost::GhostType requestedGType,
                           int requestedNumGCells, 
                           IntVector& requiredLow,
                           IntVector& requiredHigh, 
                           IntVector& requestedLow,
                           IntVector& requestedHigh) const;
    

    // Makes and returns a map that maps strings to VarLabels of
    // that name and a list of material indices for which that
    // variable is valid (at least according to d_allcomps).
    typedef std::map< std::string, std::list<int> > VarLabelMaterialMap;
    virtual VarLabelMaterialMap* makeVarLabelMaterialMap();

    virtual bool isOldDW(int idx) const;

    virtual bool isNewDW(int idx) const;

    // Only called by the SimulationController, and only once, and only
    // if the simulation has been "restarted."
    virtual void setGeneration( int id ) { d_generation = id; }

    // virtual const MaterialSet* getMaterialSet() const {cout << "BEING CALLED"
    // << endl; return reloc_.getMaterialSet();}

    // This function will copy the data from the old grid to the new grid.
    // The PatchSubset structure will contain a patch on the new grid.
    void copyDataToNewGrid(const ProcessorGroup*,
                           const PatchSubset* patches,
                           const MaterialSubset* ,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw);

    //////////
    // Schedule particle relocation without the need to supply pre_relocation variables. Use with
    // caution until as this requires further testing (tsaad).
    virtual void scheduleParticleRelocation(const LevelP& coarsestLevelwithParticles,
                                            const VarLabel* posLabel,
                                            const std::vector<std::vector<const VarLabel*> >& otherLabels,
                                            const MaterialSet* matls);
    //////////
    // Insert Documentation Here:
    virtual void scheduleParticleRelocation(const LevelP& level,
                                            const VarLabel* old_posLabel,
                                            const std::vector<std::vector<const VarLabel*> >& old_labels,
                                            const VarLabel* new_posLabel,
                                            const std::vector<std::vector<const VarLabel*> >& new_labels,
                                            const VarLabel* particleIDLabel,
                                            const MaterialSet* matls,int which);

    //////////
    // Insert Documentation Here:
    virtual void scheduleParticleRelocation(const LevelP& coarsestLevelwithParticles,
                                            const VarLabel* old_posLabel,
                                            const std::vector<std::vector<const VarLabel*> >& old_labels,
                                            const VarLabel* new_posLabel,
                                            const std::vector<std::vector<const VarLabel*> >& new_labels,
                                            const VarLabel* particleIDLabel,
                                            const MaterialSet* matls);
    
    virtual void setPositionVar(const VarLabel* posLabel) { d_reloc_new_posLabel = posLabel; }

    virtual void scheduleAndDoDataCopy(const GridP& grid, SimulationInterface* sim);

    //! override default behavior of copying, scrubbing, and such
    virtual void overrideVariableBehavior(const std::string& var, 
                                          bool treatAsOld, 
                                          bool copyData, 
                                          bool noScrub,
                                          bool notCopyData=false, 
                                          bool notCheckpoint=false);

    const std::set<std::string>& getNoScrubVars() { return d_noScrubVars;}

    const std::set<std::string>& getCopyDataVars() { return d_copyDataVars;}

    const std::set<std::string>& getNotCopyDataVars() { return d_notCopyDataVars;}

    virtual const std::set<std::string>& getNotCheckPointVars() const { return d_notCheckpointVars;}

    virtual bool useInternalDeps();
    
    const VarLabel * d_reloc_new_posLabel;
    
    int getMaxGhost() {return d_maxGhost;}
    
    int getMaxLevelOffset() {return d_maxLevelOffset;}

    bool isCopyDataTimestep() { 
      //std::cout << "d_sharedState = " << d_sharedState << "\n";
      return d_sharedState->isCopyDataTimestep() || d_isInitTimestep; 
    }

    void setInitTimestep( bool isInitTimestep ) { d_isInitTimestep = isInitTimestep; }

    void setRestartInitTimestep( bool isRestartInitTimestep ) { d_isRestartInitTimestep = isRestartInitTimestep; }

    virtual bool isRestartInitTimestep() { return d_isRestartInitTimestep; }


    typedef std::map<VarLabelMatl<Level>, Task*> ReductionTasksMap;
    ReductionTasksMap d_reductionTasks;

  protected:

    void finalizeTimestep();
    
    void makeTaskGraphDoc(const DetailedTasks* dt,
                          int rank = 0);
    void emitNode(const DetailedTask* dt, double start, double duration,
                  double execution_duration);
    void finalizeNodes(int process=0);

    enum { PRINT_BEFORE_COMM = 1, PRINT_BEFORE_EXEC = 2, PRINT_AFTER_EXEC = 4};
    void printTrackedVars(DetailedTask* dt, int when);
    
    bool d_isInitTimestep;
    bool d_isRestartInitTimestep;

    /**
     * output the task name and the level it's executing on.
     * and each of the patches
     */
    void printTask( ostream& out, DetailedTask* task );
    
    /**
     *  Output the task name and the level it's executing on
     *  only first patch of that level
     */
    void printTaskLevels( const ProcessorGroup* d_myworld, 
                          DebugStream & out, 
                          DetailedTask* task );
    
    
    virtual void verifyChecksum() = 0;

    std::vector<TaskGraph*> d_graphs;
    int d_currentTG;
    int d_numTasks;
    int d_generation;


    SimulationStateP d_sharedState;

    std::vector<OnDemandDataWarehouseP> d_dws;
    int d_numOldDWs;

    int d_dwmap[Task::TotalDWs];

    const Output* m_outPort;
    bool d_restartable;

    //! These are so we can track certain variables over the taskgraph's
    //! execution.
    std::vector<std::string> d_trackingVars;
    std::vector<std::string> d_trackingTasks;
    std::vector<Task::WhichDW> d_trackingDWs;
    int d_trackingVarsPrintLocation;
    int d_trackingPatchID;
    double d_trackingStartTime;
    double d_trackingEndTime;
    int d_trackingLevel;
    IntVector d_trackingStartIndex;
    IntVector d_trackingEndIndex;

    int d_numParticleGhostCells;
    Ghost::GhostType d_particleGhostType;

    // so we can manually copy vars between AMR levels
    std::set<std::string> d_copyDataVars;

    // ignore copying these vars between AMR levels
    std::set<std::string> d_notCopyDataVars;

    // vars manually set not to scrub (normally when needed between a normal taskgraph
    // and the regridding phase)
    std::set<std::string> d_noScrubVars;

    // treat variable as an "old" var - will be checkpointed, copied, and only scrubbed from an OldDW
    std::set<std::string> d_treatAsOldVars;
    
    // do not checkpoint these variables
    std::set<std::string> d_notCheckpointVars;
    
  private:

    SchedulerCommon(const SchedulerCommon&);
    SchedulerCommon& operator=(const SchedulerCommon&);

    // Maximum memory use as sampled across a given timestep.
    unsigned long d_maxMemUse;

    ProblemSpecP m_graphDoc;
    ProblemSpecP m_nodes;
    std::ofstream* d_memlogfile;
    bool d_emit_taskgraph;
    LocallyComputedPatchVarMap* m_locallyComputedPatchVarMap;
    Relocate         d_reloc1;
    Relocate         d_reloc2;

    // whether or not to send a small message (takes more work to organize)
    // or a larger one (more communication time)
    bool d_useSmallMessages;

    //! These are to store which vars we have to copy to the new grid
    //! in a copy data task.  Set in scheduleDataCopy and used in
    //! copyDataToNewGrid.
    typedef std::map<const VarLabel*, MaterialSubset*, VarLabel::Compare> label_matl_map;
    std::vector<label_matl_map> d_label_matls;

    //! set in addTask - can be used until initialize is called...
    std::vector<const Task::Dependency*> d_initRequires;
    std::set<const VarLabel*, VarLabel::Compare> d_initRequiredVars;
    std::set<const VarLabel*, VarLabel::Compare> d_computedVars;

    //max ghost cells of all tasks - will be used for loadbalancer to create neighborhood
    int d_maxGhost;
    
    //max level offset of all tasks - will be used for loadbalancer to create neighborhood
    int d_maxLevelOffset;
    
  };
} // End namespace Uintah

#endif
