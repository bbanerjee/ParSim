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

#ifndef VAANGO_CCA_PORTS_LOADBALANCER_H
#define VAANGO_CCA_PORTS_LOADBALANCER_H

#include <Core/Parallel/UintahParallelPort.h>

#include <CCA/Components/Schedulers/RuntimeStatsEnum.h>
#include <CCA/Ports/SchedulerP.h>

#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Region.h>
#include <Core/Grid/Variables/ComputeSet.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <Core/Util/InfoMapper.h>

#include <string>
#include <unordered_set>

namespace Uintah {

class UintahParallelComponent;
class Patch;
class ProcessorGroup;
class DetailedTasks;
class Scheduler;
class VarLabel;
class DataArchive;
class DetailedTask;

using SizeList = std::vector<Uintah::IntVector>;

//! The Load Balancer is responsible for assigning tasks to do their work
//! on specified processors.  Different subclasses differ in the way this is
//! done.
class LoadBalancer : public UintahParallelPort
{
public:
  LoadBalancer();
  virtual ~LoadBalancer();

  // eliminate copy, assignment and move
  LoadBalancer(const LoadBalancer&) = delete;
  LoadBalancer(LoadBalancer&&)      = delete;

  LoadBalancer&
  operator=(const LoadBalancer&) = delete;
  LoadBalancer&
  operator=(LoadBalancer&&) = delete;

  // Methods for managing the components attached via the ports.
  virtual void
  setComponents(UintahParallelComponent* comp) = 0;

  virtual void
  getComponents() = 0;

  virtual void
  releaseComponents() = 0;

  //! Assigns each task in tg to its corresponding processor.
  //! Uses the patchwise processor assignment.
  //! @see getPatchwiseProcessorAssignment.
  virtual void
  assignResources(DetailedTasks& tg) = 0;

  //! Gets the processor that this patch will be assigned to.
  //! This is different with the different load balancers.
  virtual int
  getPatchwiseProcessorAssignment(const Patch* patch) = 0;

  //! Gets the processor that this patch was assigned to on the last timestep.
  //! This is the same as getPatchwiseProcessorAssignment for non-dynamic load
  //! balancers. See getPatchwiseProcessorAssignment.
  virtual int
  getOldProcessorAssignment(const Patch* patch) = 0;

  //! Determines if the Load Balancer requests a taskgraph recompile.
  //! Only possible for Dynamic Load Balancers.
  virtual bool
  needRecompile(double, double, const GridP&)
  {
    return false;
  }

  //! Reads the problem spec file for the LoadBalancer section, and looks
  //! for entries such as outputNthProc, dynamicAlgorithm, and interval.
  virtual void
  problemSetup(ProblemSpecP&,
               GridP& grid,
               const MaterialManagerP& matManager) = 0;

  //! Creates the Load Balancer's Neighborhood.
  //! This is a vector of patches that represent any patch that this load
  //! balancer will potentially have to receive data from.
  virtual void
  createNeighborhood(const GridP& grid,
                     const GridP& oldGrid,
                     const bool hasDistalReqs = false) = 0;

  //! Asks the Load Balancer if it is dynamic.
  virtual bool
  isDynamic()
  {
    return false;
  }

  //! returns all processors in this processors neighborhood
  virtual const std::unordered_set<int>&
  getNeighborhoodProcessors() = 0;

  //! returns all processors in this processors global neighborhood
  virtual const std::unordered_set<int>&
  getDistalNeighborhoodProcessors() = 0;

  //! Asks if a patch in the patch subset is in the neighborhood.
  virtual bool
  inNeighborhood(const PatchSubset* patches,
                 const bool hasDistalReqs = false) = 0;

  //! Asks the load balancer if patch is in the neighborhood.
  virtual bool
  inNeighborhood(const Patch* patch, const bool hasDistalReqs = false) = 0;

  //! Returns the patchset of all patches that have work done on this processor.
  virtual const PatchSet*
  getPerProcessorPatchSet(const LevelP& level) = 0;

  virtual const PatchSet*
  getPerProcessorPatchSet(const GridP& grid) = 0;

  virtual const PatchSet*
  getOutputPerProcessorPatchSet(const LevelP& level) = 0;

  //! For dynamic load balancers, Check if we need to rebalance the load, and do
  //! so if necessary.
  virtual bool
  possiblyDynamicallyReallocate(const GridP&, int state) = 0;

  //! Returns the value of n (Only every Nth rank (process) will perform output
  //! tasks). The data from the other (N-1) processes are MPI'd to the Nth
  //! process for it to do the I/O.
  virtual int
  getNthRank() = 0;

  virtual void
  setNthRank(int nth) = 0;

  //! Returns the processor the patch will be output on (not patchwiseProcessor
  //! if outputNthProc is set)
  virtual int
  getOutputProc(const Patch* patch) = 0;

  //! Tells the load balancer on which procs data was output.
  virtual void
  restartInitialize(DataArchive* /*archive*/,
                    int /*time_index*/,
                    const std::string& /*ts_url*/,
                    const GridP& /*grid*/) = 0;

  // state variables
  enum
  {
    CHECK_LB = 0,
    INIT_LB,
    REGRID_LB,
    RESTART_LB
  };

  // cost profiling functions
  // update the contribution for this patch
  virtual void
  addContribution(DetailedTask* task, double cost) = 0;

  // finalize the contributions (updates the weight, should be called once per
  // timestep)
  virtual void
  finalizeContributions(const GridP& currentgrid) = 0;

  // initializes the weights in regions in the new grid that are not in the old
  // level
  virtual void
  initializeWeights(const Grid* oldgrid, const Grid* newgrid) = 0;

  // resets forecaster to the defaults
  virtual void
  resetCostForecaster() = 0;

  virtual int
  getNumDims() const = 0;

  virtual int*
  getActiveDims() = 0;

  virtual void
  setDimensionality(bool x, bool y, bool z) = 0;

  virtual void
  setRuntimeStats(
    ReductionInfoMapper<RuntimeStatsEnum, double>* runtimeStats) = 0;
};

} // End namespace Uintah

#endif
