/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
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

#ifndef __CCA_COMPONENTS_LOADBALANCERS_ParticleLoadBalancer_H__
#define __CCA_COMPONENTS_LOADBALANCERS_ParticleLoadBalancer_H__

#include <CCA/Components/LoadBalancers/LoadBalancerCommon.h>

#include <sci_defs/uintah_defs.h>

#include <set>
#include <string>

namespace Uintah {

class ParticleLoadBalancer : public LoadBalancerCommon
{
public:
  ParticleLoadBalancer(const ProcessorGroup* myworld);

  virtual ~ParticleLoadBalancer() = default;

  ParticleLoadBalancer(const ParticleLoadBalancer&) = delete;
  ParticleLoadBalancer(ParticleLoadBalancer&&)      = delete;
  ParticleLoadBalancer&
  operator=(const ParticleLoadBalancer&) = delete;
  ParticleLoadBalancer&
  operator=(ParticleLoadBalancer&&) = delete;

  virtual void
  problemSetup(ProblemSpecP& pspec,
               GridP& grid,
               const MaterialManagerP& mat_manager);

  virtual bool
  needRecompile(const GridP& grid);

  /// call one of the assignPatches functions.
  /// Will initially need to load balance (on first timestep), and thus
  /// be called by the SimulationController before the first compile.
  /// Afterwards, if needRecompile function tells us we need to check for
  /// load balancing, and then call
  /// this function (not called from simulation controller.)  We will then
  /// go through the motions of load balancing, and if it determines we need
  /// to load balance (that the gain is greater than some threshold), it will
  /// set the patches to their new location,
  /// return true, signifying that we need to recompile.
  /// However, if force is true, it will re-loadbalance regardless of the
  /// threshold.
  virtual bool
  possiblyDynamicallyReallocate(const GridP& grid, int state);

  //! Asks the load balancer if it is dynamic.
  virtual bool
  isDynamic()
  {
    return true;
  }

  // Collects each patch's particles
  void
  collectParticles(const Grid* grid,
                   std::vector<std::vector<int>>& num_particles);

  // same, but can be called after a regrid when patches have not been load
  // balanced yet
  void
  collectParticlesForRegrid(
    const Grid* oldGrid,
    const std::vector<std::vector<Region>>& newGridRegions,
    std::vector<std::vector<int>>& particles);

private:
  struct double_int
  {
    double val;
    int loc;
    double_int(double val, int loc)
      : val(val)
      , loc(loc)
    {
    }
    double_int()
      : val(0)
      , loc(-1)
    {
    }
  };

  enum class ParticleAlgorithm
  {
    static_lb,
    cyclic_lb,
    random_lb,
    patch_factor_lb
  };

  double d_lbThreshold; //! gain threshold to exceed to require lb'ing

  // The weighting factor placed on particles and cells, for example if
  // d_particleCost is 2 and d_cellCost is 1 then a particle has twice as much
  // weight as a cell.
  double d_particleCost;
  double d_cellCost;

  ProblemSpecP d_pspec;
  std::vector<IntVector> d_minPatchSize;

  /// This functions take care of setting
  /// d_temp_assignment on all procs and dynamicReallocation takes care of
  /// maintaining the state
  bool
  loadBalanceGrid(const GridP& grid, bool force);

  // gets the cell costs and particle costs for each patch
  void
  getCosts(const Grid* grid,
           std::vector<std::vector<double>>& particleCosts,
           std::vector<std::vector<double>>& cellCosts);

  // sets processor "assignments" for "patches" based on the "patchCosts" and
  // "previousProcCosts"
  void
  assignPatches(const std::vector<double>& previousProcCosts,
                const std::vector<double>& patchCosts,
                std::vector<int>& patches,
                std::vector<int>& assignments);

  // compute the percent improvement of the assignments in d_temp_assignment vs
  // d_processor_assignment for the given cost array
  double
  computePercentImprovement(const std::vector<std::vector<double>>& costs,
                            double& avg,
                            double& max);

  // compute the load imbalance of d_processor_assignment given the costs vector
  double
  computeImbalance(const std::vector<std::vector<double>>& costs);

  // given the two cost arrays determine if the new load balance is better than
  // the previous
  bool
  thresholdExceeded(const std::vector<std::vector<double>>& cellCosts,
                    const std::vector<std::vector<double>>& particleCosts);
};

} // End namespace Uintah

#endif //__CCA_COMPONENTS_LOADBALANCERS_ParticleLoadBalancer_H__
