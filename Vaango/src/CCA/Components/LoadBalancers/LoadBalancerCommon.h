/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#ifndef __CCA_COMPONENTS_LOADBALANCERS__LoadBalancerCommon_H__
#define __CCA_COMPONENTS_LOADBALANCERS__LoadBalancerCommon_H__

#include <CCA/Ports/LoadBalancer.h>

#include <CCA/Ports/SpaceFillingCurve.h>

#include <Core/Grid/Level.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <Core/Util/DebugStream.h>

#include <set>
#include <string>

namespace Uintah {

class SimulationInterface;

struct PatchInfo
{
  PatchInfo(int i, int n)
  {
    id           = i;
    numParticles = n;
  }

  PatchInfo() = default;

  int id;
  int numParticles;
};

class ParticleCompare
{
public:
  inline bool
  operator()(const PatchInfo& p1, const PatchInfo& p2) const
  {
    return p1.numParticles < p2.numParticles ||
           (p1.numParticles == p2.numParticles && p1.id < p2.id);
  }
};

class PatchCompare
{
public:
  inline bool
  operator()(const PatchInfo& p1, const PatchInfo& p2) const
  {
    return p1.id < p2.id;
  }
};

/// Load Balancer Common.  Implements many functions in common among
/// the load balancer subclasses.  The main function that sets load balancers
/// apart is getPatchwiseProcessorAssignment - how it determines which patch
/// to assign on which procesor.
class LoadBalancerCommon
  : public LoadBalancer
  , public UintahParallelComponent
{
public:
  LoadBalancerCommon(const ProcessorGroup* myworld);

  virtual ~LoadBalancerCommon() = default;

  // Disallow copy and move
  LoadBalancerCommon(const LoadBalancerCommon&) = delete;
  LoadBalancerCommon(LoadBalancerCommon&&)      = delete;
  LoadBalancerCommon&
  operator=(const LoadBalancerCommon&) = delete;
  LoadBalancerCommon&
  operator=(LoadBalancerCommon&&) = delete;

  // Methods for managing the components attached via the ports.
  virtual void
  setComponents(UintahParallelComponent* comp){};

  virtual void
  getComponents();

  virtual void
  releaseComponents();

  virtual int
  getPatchwiseProcessorAssignment(const Patch* patch);

  //! The implementation in LoadBalancerCommon.cc is for dynamice load
  //! balancers.
  //! The Simple and SingleProcessor override this function with default
  //! implementations.
  virtual int
  getOldProcessorAssignment(const Patch* patch);

  virtual bool
  needRecompile(const GridP&) = 0;

  /// Goes through the Detailed tasks and assigns each to its own processor.
  virtual void
  assignResources(DetailedTasks& tg);

  /// Creates the Load Balancer's Neighborhood.  This is a vector of patches
  /// that represent any patch that this load balancer will potentially have to
  /// receive data from.
  virtual void
  createNeighborhood(const GridP& grid,
                     const GridP& oldGrid,
                     const bool hasDistalReqs = false);

  virtual const std::unordered_set<int>&
  getNeighborhoodProcessors()
  {
    return d_local_neighbor_processes;
  }

  virtual const std::unordered_set<int>&
  getDistalNeighborhoodProcessors()
  {
    return d_distal_neighbor_processes;
  }

  /// Asks the load balancer if a patch in the patch subset is in the
  /// neighborhood.
  virtual bool
  inNeighborhood(const PatchSubset*, const bool hasDistalReqs = false);

  /// Asks the load balancer if patch is in the neighborhood.
  virtual bool
  inNeighborhood(const Patch*, const bool hasDistalReqs = false);

  /// Reads the problem spec file for the LoadBalancer section, and looks
  /// for entries such as outputNthProc, dynamicAlgorithm, and interval.
  virtual void
  problemSetup(ProblemSpecP& pspec, GridP& grid, const MaterialManagerP& mat_manager);

  // for DynamicLoadBalancer mostly, but if we're called then it also means the
  // grid might have changed and need to create a new perProcessorPatchSet
  virtual bool
  possiblyDynamicallyReallocate(const GridP&, int state);

  // Cost profiling functions
  // Update the contribution for this patch.
  virtual void
  addContribution(DetailedTask* task, double cost);

  // Finalize the contributions (updates the weight, should be called once per
  // timestep):
  virtual void
  finalizeContributions(const GridP& currentGrid);

  // Initializes the regions in the new level that are not in the old level.
  virtual void
  initializeWeights(const Grid* oldgrid, const Grid* newgrid);

  // Resets the profiler counters to zero
  virtual void
  resetCostForecaster();

  //! Returns n - data gets output every n procs.
  virtual int
  getNthRank()
  {
    return d_output_Nth_proc;
  }

  virtual void
  setNthRank(int nth)
  {
    d_output_Nth_proc = nth;
  }

  //! Returns the processor the patch will be output on (not patchwiseProcessor
  //! if outputNthProc is set)
  virtual int
  getOutputProc(const Patch* patch)
  {
    return (getPatchwiseProcessorAssignment(patch) / d_output_Nth_proc) *
           d_output_Nth_proc;
  }

  //! Returns the patchset of all patches that have work done on this processor.
  virtual const PatchSet*
  getPerProcessorPatchSet(const LevelP& level)
  {
    return d_level_perproc_patchsets[level->getIndex()].get_rep();
  }

  virtual const PatchSet*
  getPerProcessorPatchSet(const GridP& grid)
  {
    return d_grid_perproc_patchsets.get_rep();
  }

  virtual const PatchSet*
  getOutputPerProcessorPatchSet(const LevelP& level)
  {
    return d_output_patchsets[level->getIndex()].get_rep();
  };

  //! Assigns the patches to the processors they ended up on in the previous
  //! Simulation.  Returns true if we need to re-load balance (if we have a
  //! different number of procs than were saved to disk
  virtual void
  restartInitialize(DataArchive* archive,
                    const int time_index,
                    const std::string& tsurl,
                    const GridP& grid);

  int
  getNumDims() const
  {
    return d_numDims;
  };

  int*
  getActiveDims()
  {
    return d_activeDims;
  };

  void
  setDimensionality(bool x, bool y, bool z);

  void
  setRuntimeStats(ReductionInfoMapper<RuntimeStatsEnum, double>* runtimeStats)
  {
    d_runtimeStats = runtimeStats;
  }

protected:
  // Calls space-filling curve on level, and stores results in pre-allocated
  // output
  void
  useSpaceFillingCurve(const LevelP& level, int* output);

  /// Creates a patchset of all patches that have work done on each processor.
  //    - There are two versions of this function.  The first works on a per
  //    level
  //      basis.  The second works on the entire grid and will provide a
  //      PatchSet that contains all patches.
  //    - For example, if processor 1 works on patches 1,2 on level 0 and patch
  //    3 on level 1,
  //      and processor 2 works on 4,5 on level 0, and 6 on level 1, then
  //      - Version 1 (for Level 1) will create {{3},{6}}
  //      - Version 2 (for all levels) will create {{1,2,3},{4,5,6}}
  virtual const PatchSet*
  createPerProcessorPatchSet(const LevelP& level);

  virtual const PatchSet*
  createPerProcessorPatchSet(const GridP& grid);

  virtual const PatchSet*
  createOutputPatchSet(const LevelP& level);

protected:
  SimulationInterface* d_simulator{ nullptr };

  int d_lb_timeStep_interval{ 0 };
  int d_last_lb_timeStep{ 0 };

  double d_last_lb_simTime{ 0.0 };
  double d_lb_interval{ 0.0 };
  bool d_check_after_restart{ false };

  // The assignment vectors are stored 0-n.  This stores the start patch number
  // so we can
  // detect if something has gone wrong when we go to look up what proc a patch
  // is on.
  int d_assignment_base_patch{ -1 };
  int d_old_assignment_base_patch{ -1 };

  std::vector<int>
    d_processor_assignment; ///< stores which proc each patch is on
  std::vector<int>
    d_old_assignment; ///< stores which proc each patch used to be on
  std::vector<int>
    d_temp_assignment; ///< temp storage for checking to reallocate

  SpaceFillingCurve<double> d_sfc;
  bool d_do_space_curve{ false };

  //! to keep track of timesteps
  MaterialManagerP d_mat_manager;

  //! store the scheduler to not have to keep passing it in
  Scheduler* d_scheduler{ nullptr };

  //! the neighborhood.  See createNeighborhood
  std::unordered_set<const Patch*> d_local_neighbor_patches;

  //! a list of processors that are in this processors neighborhood
  std::unordered_set<int> d_local_neighbor_processes;

  //! the "distal" neighborhood of patches.  See createNeighborhood
  std::unordered_set<const Patch*> d_distal_neighbor_patches;

  //! a list of "distal" processes that are in this processors neighborhood
  std::unordered_set<int> d_distal_neighbor_processes;

  //! output on every nth processor.  This variable needs to be shared
  //! with the DataArchiver as well, but we keep it here because the lb
  //! needs it to assign the processor resource.
  int d_output_Nth_proc;

  std::vector<Handle<const PatchSet>> d_level_perproc_patchsets;
  Handle<const PatchSet> d_grid_perproc_patchsets;
  std::vector<Handle<const PatchSet>> d_output_patchsets;

  // Which dimensions are active.  Get the number of dimensions, and
  // then that many indices of activeDims are set to which dimensions
  // are being used.
  int d_numDims{ 0 };
  int d_activeDims[3]{0, 0, 0};

  ReductionInfoMapper<RuntimeStatsEnum, double>* d_runtimeStats{ nullptr };

  DebugStream stats;
  DebugStream times;
  DebugStream lbout;

private:
  void
  addPatchesAndProcsToNeighborhood(const Level* const level,
                                   const IntVector& low,
                                   const IntVector& high,
                                   std::unordered_set<const Patch*>& neighbors,
                                   std::unordered_set<int>& processors);
};

} // End namespace Uintah

#endif //__CCA_COMPONENTS_LOADBALANCERS__LoadBalancerCommon_H__
