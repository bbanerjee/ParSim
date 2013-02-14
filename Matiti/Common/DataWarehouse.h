#ifndef MATITI_DATAWAREHOUSE_H
#define MATITI_DATAWAREHOUSE_H

#include <Common/RefCounted.h>
#include <Common/Handle.h>
#include <Common/Task.h>
#include <Common/DataWarehouseP.h>
#include <Common/SchedulerP.h>
#include <Common/ComputeSet.h>
#include <Common/VarLabelMatl.h>
#include <Mesh/MeshP.h>
#include <Mesh/MeshNodeVariableBase.h>

#include <iosfwd>

namespace Matiti {

class OutputContext;
class VarLabel;
class Task;

class DataWarehouse : public RefCounted {

public:
  virtual ~DataWarehouse();
      
  virtual bool exists(const VarLabel*, int matlIndex) const = 0;

  // Returns a (const) pointer to the mesh.  This pointer can then be
  // used to (for example) get the number of levels in the mesh.
  virtual const Mesh * getMesh() = 0;

  // Generic put and allocate, passing Variable as a pointer rather than
  // by reference to avoid ambiguity with other put overloaded methods.
  virtual void put(Variable*, const VarLabel*, int matlIndex)  = 0;

  // MeshNode Variables
  virtual MeshNodeSubset* createMeshNodeSubset(meshNodeIndex numMeshNodes,
                                               int matlIndex) = 0; 
  virtual void saveMeshNodeSubset(MeshNodeSubset* psubset,
                                  int matlIndex) = 0;
  virtual bool haveMeshNodeSubset(int matlIndex, bool exact = false) = 0;
  virtual MeshNodeSubset* getMeshNodeSubset(int matlIndex) = 0;
  virtual MeshNodeSubset* getDeleteSubset(int matlIndex) = 0;
  virtual std::map<const VarLabel*, MeshNodeVariableBase*>* getNewMeshNodeState(int matlIndex) = 0;
  virtual MeshNodeSubset* getMeshNodeSubset(int matlIndex, const VarLabel* posvar) = 0;
  virtual void allocateTemporary(MeshNodeVariableBase&, MeshNodeSubset*) = 0;
  virtual void allocateAndPut(MeshNodeVariableBase&, const VarLabel*, MeshNodeSubset*) = 0;
  virtual void get(constMeshNodeVariableBase&, const VarLabel*, MeshNodeSubset*) = 0;
  virtual void get(constMeshNodeVariableBase&, const VarLabel*, int matlIndex) = 0;
  virtual void getModifiable(MeshNodeVariableBase&, const VarLabel*, MeshNodeSubset*) = 0;
  virtual void put(MeshNodeVariableBase&, const VarLabel*, bool replace = false) = 0;

  virtual void getCopy(MeshNodeVariableBase&, const VarLabel*, MeshNodeSubset*) = 0;
  virtual void copyOut(MeshNodeVariableBase&, const VarLabel*, MeshNodeSubset*) = 0;

  virtual void print() = 0;
  virtual void clear() = 0;

  virtual MeshNodeVariableBase* getMeshNodeVariable(const VarLabel*, MeshNodeSubset*) = 0;
  virtual MeshNodeVariableBase*
  getMeshNodeVariable(const VarLabel*, int matlIndex) = 0;

  // Remove meshNodes that are no longer relevant
  virtual void deleteMeshNodes(MeshNodeSubset* delset) = 0;

  // Move stuff to a different data Warehouse
  virtual void transferFrom(DataWarehouse*, const VarLabel*,
			    const MaterialSubset*, bool replace = false) = 0;

  virtual void emit(OutputContext&, const VarLabel* label, int matlIndex) = 0;

  // Scrubbing
  enum ScrubMode {
    ScrubNone,
    ScrubComplete,
    ScrubNonPermanent
  };
  virtual ScrubMode setScrubbing(ScrubMode) = 0;

      // For related datawarehouses
  virtual DataWarehouse* getOtherDataWarehouse(Task::WhichDW) = 0;

  // For the schedulers
  virtual bool isFinalized() const = 0;
  virtual void finalize() = 0;
  virtual void unfinalize() = 0;
  virtual void refinalize() = 0;

  // Returns the generation number (id) of this data warehouse.  Id's
  // start at 0. Each subsequent DW's id is one greater.  SetID should
  // only be called by the SimulationController (and only once) if the
  // DW is being used for a restarted simulation.  This allows the DW
  // generation number to be kept in sync with the actual number of
  // timesteps for the restarted simulation.
  int  getID() const { return d_generation; }
  void setID( int id ) { d_generation = id; }

  // For timestep abort/restart
  virtual bool timestepAborted() = 0;
  virtual bool timestepRestarted() = 0;
  virtual void abortTimestep() = 0;
  virtual void restartTimestep() = 0;
  
protected:
  DataWarehouse(Scheduler* scheduler, 
		int generation );

  // These two things should be removed from here if possible - Steve
  Scheduler* d_scheduler;

  // Generation should be const, but is not as during a restart, the
  // generation number of the first DW is updated from 0 (the default
  // for the first DW) to the correct generation number based on how
  // many previous time steps had taken place before the restart.
  int d_generation;
     
private:
  DataWarehouse(const DataWarehouse&);
  DataWarehouse& operator=(const DataWarehouse&);
};

} // End namespace Matiti

#endif
