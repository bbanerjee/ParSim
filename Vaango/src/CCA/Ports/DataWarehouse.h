/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
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

#ifndef UINTAH_HOMEBREW_DataWarehouse_H
#define UINTAH_HOMEBREW_DataWarehouse_H

#include <Core/Util/Handle.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/Variables/constGridVariable.h>
#include <Core/Grid/Ghost.h>
#include <Core/Util/RefCounted.h>
#include <Core/Grid/Variables/ParticleVariableBase.h>
#include <Core/Grid/Variables/ReductionVariableBase.h>
#include <Core/Grid/Variables/PerPatchBase.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/SoleVariableBase.h>
#include <Core/Grid/Variables/VarLabelMatl.h>
#include <Core/Grid/Task.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <CCA/Ports/SchedulerP.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Vector.h>

#include <iosfwd>


namespace Uintah {

class Level;
class OutputContext;
class ProcessorGroup;
class VarLabel;
class Task;

/**************************************
	
CLASS
   DataWarehouse
	
   Short description...
	
GENERAL INFORMATION
	
   DataWarehouse.h
	
   Steven G. Parker
   Department of Computer Science
   University of Utah
	
   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
	
   	
KEYWORDS
   DataWarehouse
	
DESCRIPTION
   Long description...
	
WARNING
	
****************************************/
      
typedef std::map<const VarLabel*, ParticleVariableBase*> ParticleLabelDataMap;
typedef std::map<long64, int> ParticleIDMap;

class DataWarehouse : public RefCounted {

public:
  virtual ~DataWarehouse();
      
  // Check whether a label exists in the DW
  virtual bool exists(const VarLabel* label, 
                      int matlIndex, 
                      const Patch* patch) const = 0;

  // Returns a (const) pointer to the grid.  This pointer can then be
  // used to (for example) get the number of levels in the grid.
  virtual const Grid * getGrid() = 0;

  // Generic put and allocate, passing Variable as a pointer rather than
  // by reference to avoid ambiguity with other put overloaded methods.
  virtual void put(Variable* var, 
                   const VarLabel* label, 
                   int matlIndex,
		   const Patch* patch) = 0;

  // Reduction Variables
  virtual void get(ReductionVariableBase& var, 
                   const VarLabel* label,
		   const Level* level = 0, 
                   int matlIndex = -1) = 0;

  virtual void put(const ReductionVariableBase& var, 
                   const VarLabel* label,
		   const Level* level = 0, 
                   int matlIndex = -1) = 0;

  virtual void override(const ReductionVariableBase& var, 
                        const VarLabel* label,
			const Level* level = 0, 
                        int matlIndex = -1) = 0;

  virtual void print(ostream& intout, 
                     const VarLabel* label,
		     const Level* level, 
                     int matlIndex = -1) = 0;

  // Sole Variables
  virtual bool exists(const VarLabel* label) const = 0;

  virtual void get(SoleVariableBase& var, 
                   const VarLabel* label,
		   const Level* level = 0, 
                   int matlIndex = -1) = 0;

  virtual void put(const SoleVariableBase& var, 
                   const VarLabel* label,
		   const Level* level = 0, 
                   int matlIndex = -1) = 0;

  virtual void override(const SoleVariableBase& var, 
                        const VarLabel* label,
			const Level* level = 0, 
                        int matlIndex = -1) = 0;


  // Particle Variables
  // changed way PS's were stored from ghost info to low-high range.
  // we can still keep the getPS function API the same though to not annoy
  // everybody -- BJW, May 05
  virtual ParticleSubset* createParticleSubset(particleIndex numParticles,
                                               int matlIndex, 
                                               const Patch* patch,
                                               SCIRun::IntVector low = SCIRun::IntVector(0,0,0),
                                               SCIRun::IntVector high = SCIRun::IntVector(0,0,0) ) = 0;

  virtual void saveParticleSubset(ParticleSubset* psubset,
                                  int matlIndex, 
                                  const Patch* patch,
                                  SCIRun::IntVector low = SCIRun::IntVector(0,0,0),
                                  SCIRun::IntVector high = SCIRun::IntVector(0,0,0)) = 0;

  virtual bool haveParticleSubset(int matlIndex, 
                                  const Patch* patch,
                                  SCIRun::IntVector low = SCIRun::IntVector(0,0,0),
                                  SCIRun::IntVector high = SCIRun::IntVector(0,0,0), 
                                  bool exact = false) = 0;

  virtual ParticleSubset* getParticleSubset(int matlIndex, 
                                            const Patch* patch,
                                            SCIRun::IntVector low, 
                                            SCIRun::IntVector high) = 0;

  virtual ParticleSubset* getParticleSubset(int matlIndex, 
                                            const Patch* patch) = 0;

  virtual ParticleSubset* getDeleteSubset(int matlIndex, 
                                          const Patch* patch) = 0;

  virtual ParticleLabelDataMap* getNewParticleState(int matlIndex, 
                                                    const Patch* patch) = 0;

  virtual ParticleSubset* getParticleSubset(int matlIndex, 
                                            const Patch* patch,
					    Ghost::GhostType,
					    int numGhostCells,
					    const VarLabel* posvar) = 0;

  virtual ParticleSubset* getParticleSubset(int matlIndex, 
                                            SCIRun::IntVector low, 
                                            SCIRun::IntVector high,
                                            const Patch* relPatch,
                                            const VarLabel* posvar, 
                                            const Level* level=0) = 0;

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
                                             const VarLabel* posvar) = 0;
                                 
  // Get the particle index values of a set of ParticleIDs
  virtual void getParticleIndex(ParticleSubset* pset,
                                const VarLabel* partIDLabel,
                                const std::vector<long64>& partIDList,
                                std::vector<particleIndex>& partIndexList) = 0;

  // Create a map between the long64 particleIDs and the particle indices in a 
  // ParticleSubset
  virtual void createParticleIDMap(ParticleSubset* pset,
                                   const VarLabel* partIDLabel,
                                   ParticleIDMap& partIDMap) = 0;

  // Get the particle index value of a ParticleID after the partIDMap 
  // has been created
  virtual void getParticleIndex(const ParticleIDMap& partIDMap,
                                const long64& pParticleID,
                                particleIndex& pParticleIndex) = 0;

  virtual void allocateTemporary(ParticleVariableBase& var,
				 ParticleSubset* pset) = 0;

  virtual void allocateAndPut(ParticleVariableBase& var, 
                              const VarLabel* label,
			      ParticleSubset* pset) = 0;

  virtual void get(constParticleVariableBase& var, 
                   const VarLabel* label,
		   ParticleSubset* pset) = 0;

  virtual void get(constParticleVariableBase& var, 
                   const VarLabel* label,
		   int matlIndex, 
                   const Patch* patch) = 0;

  virtual void getModifiable(ParticleVariableBase& var, 
                             const VarLabel* label,
			     ParticleSubset* pset) = 0;

  virtual void put(ParticleVariableBase& var, 
                   const VarLabel* label,
		   bool replace = false) = 0;

  virtual void getCopy(ParticleVariableBase& var, 
                       const VarLabel* label, 
                       ParticleSubset* pset) = 0;

  virtual void copyOut(ParticleVariableBase& var, 
                       const VarLabel* label, 
                       ParticleSubset* pset) = 0;

  virtual void print() = 0;

  virtual void clear() = 0;


  virtual ParticleVariableBase* getParticleVariable(const VarLabel* label,
						    ParticleSubset* pset) = 0;

  virtual ParticleVariableBase* getParticleVariable(const VarLabel* label, 
                                                    int matlIndex, 
                                                    const Patch* patch) = 0;

  // Generic grid based variables
  virtual void get(constGridVariableBase& var,
                   const VarLabel* label, 
                   int matlIndex, 
                   const Patch* patch,
                   Ghost::GhostType gtype, 
                   int numGhostCells ) = 0;

  virtual void getModifiable(GridVariableBase& var,
                             const VarLabel* label, 
                             int matlIndex, 
                             const Patch* patch, 
                             Ghost::GhostType gtype=Ghost::None, 
                             int numGhostCells=0 ) = 0;

  virtual void allocateTemporary(GridVariableBase& var, 
                                 const Patch* patch,
                                 Ghost::GhostType gtype = Ghost::None, 
                                 int numGhostCells = 0 ) = 0;
//                               const SCIRun::IntVector& boundaryLayer ) = 0;
//                               const SCIRun::IntVector& boundaryLayer = SCIRun::IntVector(0,0,0)) = 0;

  virtual void allocateAndPut(GridVariableBase& var,
                              const VarLabel* label, 
                              int matlIndex,
                              const Patch* patch, 
                              Ghost::GhostType gtype = Ghost::None,
                              int numGhostCells = 0 ) = 0;

  virtual void put(GridVariableBase& var, 
                   const VarLabel* label, 
                   int matlIndex, 
                   const Patch* patch,
                   bool replace = false) = 0;

  virtual void getRegion(constGridVariableBase& var, 
                         const VarLabel* label,
                         int matlIndex, 
                         const Level* level,
                         const SCIRun::IntVector& low, 
                         const SCIRun::IntVector& high,
                         bool useBoundaryCells = true) = 0;

  // Copy out of the warehouse into an allocated variable.
  virtual void copyOut(GridVariableBase& var, 
                       const VarLabel* label, 
                       int matlIndex,
	               const Patch* patch, 
                       Ghost::GhostType gtype = Ghost::None,
	               int numGhostCells = 0) = 0;

  // Makes var a copy of the specified warehouse data, allocating it
  // to the appropriate size first.
  virtual void getCopy(GridVariableBase& var, 
                       const VarLabel* label, 
                       int matlIndex,
	               const Patch* patch, 
                       Ghost::GhostType gtype = Ghost::None,
	               int numGhostCells = 0) = 0;
      
  // PerPatch Variables
  virtual void get(PerPatchBase& var, 
                   const VarLabel* label,
		   int matlIndex, 
                   const Patch* patch) = 0;

  virtual void put(PerPatchBase& var, 
                   const VarLabel* label,
		   int matlIndex, 
                   const Patch* patch, 
                   bool replace = false) = 0;
     
  // this is so we can get reduction information for regridding
  virtual void getVarLabelMatlLevelTriples(std::vector<VarLabelMatl<Level> >& vars ) const = 0;

  // Remove particles that are no longer relevant
  virtual void deleteParticles(ParticleSubset* delset) = 0;

  // Add particles
  virtual void addParticles(const Patch* patch, 
                            int matlIndex,
			    ParticleLabelDataMap* addedstate) = 0;

  // Move stuff to a different data Warehouse
  virtual void transferFrom(DataWarehouse* dw, 
                            const VarLabel* label,
			    const PatchSubset* pset, 
                            const MaterialSubset* mset,
                            bool replace = false, 
                            const PatchSubset* patchset = 0) = 0;

  virtual void emit(OutputContext& out, 
                    const VarLabel* label,
		    int matlIndex, 
                    const Patch* patch) = 0;

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

  DataWarehouse(const ProcessorGroup* myworld,
		Scheduler* scheduler, 
		int generation );

  // These two things should be removed from here if possible - Steve
  const ProcessorGroup* d_myworld;
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

} // End namespace Uintah

#endif
