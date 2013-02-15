#ifndef MATITI_SCHEDULER_H
#define MATITI_SCHEDULER_H

#include <Common/SerialPort.h>
#include <Common/ComputeSet.h>
#include <Common/Task.h>
#include <Mesh/SimulationStateP.h>
#include <Mesh/DomainP.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <map>
#include <list>
#include <string>
#include <set>
#include <vector>


namespace Matiti {

  class Task;
  class SimulationInterface;

  class Scheduler : public SerialPort {

    public:
      Scheduler();
      virtual ~Scheduler();
   
      // Only called by the SimulationController, and only once, and only
      // if the simulation has been "restarted."
      virtual void setGeneration( int id ) = 0;

      virtual void problemSetup(const Uintah::ProblemSpecP& prob_spec, 
                                SimulationStateP& state) = 0;
    
      virtual void initialize(int numOldDW = 1, int numNewDW = 1) = 0;

      virtual void setParentDWs(DataWarehouse* parent_old_dw, 
                                DataWarehouse* parent_new_dw) = 0;

      virtual void clearMappings() = 0;
      virtual void mapDataWarehouse(Task::WhichDW, int dwTag) = 0;

      virtual void doEmitTaskGraphDocs() = 0;
    
      virtual void compile() = 0;
      virtual void execute(int tgnum = 0, int iteration = 0) = 0;

      virtual SchedulerP createSubScheduler() = 0;
       
      enum tgType { NormalTaskGraph, IntermediateTaskGraph };

      virtual void addTaskGraph(tgType type) = 0;
      virtual int getNumTaskGraphs() = 0;
      virtual bool useSmallMessages() = 0;
    
      virtual void addTask(Task* t, const MaterialSet*) = 0;
    
      virtual const std::vector<const Task::Dependency*>& getInitialRequires() = 0;
      virtual const std::set<const VarLabel*, VarLabel::Compare>& getInitialRequiredVars() const = 0;
      virtual const std::set<const VarLabel*, VarLabel::Compare>& getComputedVars() const = 0;
      virtual const std::set<std::string>& getNotCheckPointVars() const = 0;    

      virtual DataWarehouse* get_dw(int idx) = 0;
      virtual DataWarehouse* getLastDW(void) = 0;

      virtual bool isOldDW(int idx) const = 0;
      virtual bool isNewDW(int idx) const = 0;

      virtual void advanceDataWarehouse(const DomainP& domain, bool initialization=false) = 0;
      virtual void fillDataWarehouses(const DomainP& domain) = 0;
      virtual void replaceDataWarehouse(int index, const DomainP& domain, bool initialization=false) = 0;
      virtual void setRestartable(bool restartable) = 0;

      virtual void setPositionVar(const VarLabel* posLabel) = 0;
    
      // Makes and returns a map that maps strings to VarLabels of
      // that name and a list of material indices for which that
      // variable is valid (at least according to d_allcomps).
      typedef std::map< std::string, std::list<int> > VarLabelMaterialMap;
      virtual VarLabelMaterialMap* makeVarLabelMaterialMap() = 0;

    private:
      Scheduler(const Scheduler&);
      Scheduler& operator=(const Scheduler&);
  };
} // End namespace Matiti

#endif
