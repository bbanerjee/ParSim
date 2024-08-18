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

#ifndef UINTAH_CCA_COMPONENTS_REGRIDDERS_REGRIDDERCOMMON_H
#define UINTAH_CCA_COMPONENTS_REGRIDDERS_REGRIDDERCOMMON_H

//-- Uintah component includes --//
#include <CCA/Ports/Regridder.h>

//-- Uintah framework includes --//
#include <Core/Geometry/IntVector.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Parallel/UintahParallelComponent.h>

//-- system includes --//
#include <vector>

namespace Uintah {

class SimulationInterface;
class DataWarehouse;
class Patch;
class VarLabel;
class ProcessorGroup;
class LoadBalancer;
class Scheduler;

namespace {
using SizeList = std::vector<IntVector>;
}

/**
 *  @ingroup Regridders
 *  @class   RegridderCommon
 *  @author  Bryan Worthen
 *           Department of Computer Science
 *           University of Utah
 *  @date    CSAFE days - circa 06/04 (updated 08/14)
 *  @brief   Parent class which takes care of common regridding functionality.
 */
class RegridderCommon
  : public Regridder
  , public UintahParallelComponent
{

public:
  RegridderCommon(const ProcessorGroup* pg);
  virtual ~RegridderCommon();

  // Methods for managing the components attached via the ports.
  virtual void
  setComponents([[maybe_unused]] UintahParallelComponent* comp){};

  virtual void
  getComponents();

  virtual void
  releaseComponents();

  //! Initialize with regridding parameters from ups file
  virtual void
  problemSetup(const ProblemSpecP& params,
               const GridP& grid,
               const MaterialManagerP& mat_manager);

  //! On a Switch, basically asks whether to turn off/on the Regridding
  virtual void
  switchInitialize(const ProblemSpecP& params);

  //! Asks if we need to recompile the task graph.
  //! Will return true if we did a regrid
  virtual bool
  needRecompile(const GridP& grid);

  //! Do we need to regrid this timestep?
  virtual bool
  needsToReGrid(const GridP& grid);

  //! Asks if we are going to do regridding
  virtual bool
  isAdaptive()
  {
    return d_isAdaptive;
  }

  //! switch for setting adaptivity
  virtual void
  setAdaptivity(const bool ans)
  {
    d_isAdaptive = ans;
  }

  //! Ask if regridding only once.
  virtual bool
  doRegridOnce()
  {
    return d_regridOnce;
  }

  //! Asks if we are going to force regridding
  virtual bool
  forceRegridding()
  {
    return d_forceRegridding;
  }

  //! switch for forcing regridding
  virtual void
  setForceRegridding(const bool val)
  {
    d_forceRegridding = val;
  }

  //! Schedules task to initialize the error flags to 0
  virtual void
  scheduleInitializeErrorEstimate(const LevelP& level);

  //! Schedules task to dilate existing error flags
  virtual void
  scheduleDilation(const LevelP& level, bool isLockstepAMR);

  //! Asks if we are going to do regridding
  virtual bool
  flaggedCellsOnFinestLevel(const GridP& grid);

  //! Returns the max number of levels this regridder will store
  virtual int
  maxLevels()
  {
    return d_maxLevels;
  }

  virtual bool
  useDynamicDilation()
  {
    return d_dynamicDilation;
  }

  virtual void
  setOverheadAverage(double val)
  {
    d_overheadAverage = val;
  }

  enum FilterType
  {
    FILTER_STAR,
    FILTER_BOX
  };

  enum DilationType
  {
    DILATE_STABILITY,
    DILATE_REGRID,
    DILATE_DELETION,
    DILATE_PATCH
  };

  //! initialize the refineFlag variable for this domain (a task callback)
  void
  initializeErrorEstimate(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse*,
                          DataWarehouse* new_dw);

  void
  Dilate(const ProcessorGroup*,
         const PatchSubset* patches,
         const MaterialSubset*,
         DataWarehouse* old_dw,
         DataWarehouse* new_dw,
         const VarLabel* to_put,
         CCVariable<int>* filter,
         IntVector depth);

  const MaterialSubset*
  refineFlagMaterials() const;

  const VarLabel*
  getRefineFlagLabel() const
  {
    return d_refineFlagLabel;
  }

  const VarLabel*
  getOldRefineFlagLabel() const
  {
    return d_oldRefineFlagLabel;
  }

  const VarLabel*
  getRefinePatchFlagLabel() const
  {
    return d_refinePatchFlagLabel;
  }

protected:
  ProblemSpecP d_grid_ps{ nullptr };
  LoadBalancer* d_load_balancer{ nullptr };
  Scheduler* d_scheduler{ nullptr };
  SimulationInterface* d_simulator{ nullptr };

  MaterialManagerP d_mat_manager{ nullptr };
  MaterialSubset* d_refine_flag_matls{ nullptr };

  ///< If false, do not regrid (stick with what you have)
  bool d_isAdaptive{ false };
  bool d_forceRegridding{ false }; ///< If false, do not force regriding

  // input parameters from ups file
  bool d_dynamicDilation{ false };
  IntVector d_maxDilation{ 0, 0, 0 };
  SizeList d_cellNum{};
  SizeList d_cellRefinementRatio{};
  IntVector d_cellStabilityDilation{ 0, 0, 0 };
  IntVector d_cellRegridDilation{ 0, 0, 0 };
  IntVector d_cellDeletionDilation{ 0, 0, 0 };
  ///< min # of cells to be between levels' boundaries
  IntVector d_minBoundaryCells{ 0, 0, 0 };
  FilterType d_filterType{ FILTER_BOX };

  std::vector<CCVariable<int>*> d_flaggedCells{};
  std::vector<CCVariable<int>*> d_dilatedCellsStability{};
  std::vector<CCVariable<int>*> d_dilatedCellsRegrid{};
  std::vector<CCVariable<int>*> d_dilatedCellsDeleted{};

  std::map<IntVector, std::unique_ptr<CCVariable<int>>> filters;
  CCVariable<int> d_patchFilter;

  int d_maxLevels{ 1 };

  // var labels for interior task graph
  const VarLabel* d_dilatedCellsStabilityLabel{ nullptr };
  const VarLabel* d_dilatedCellsRegridLabel{ nullptr };
  const VarLabel* d_dilatedCellsDeletionLabel{ nullptr };

  const VarLabel* d_refineFlagLabel{ nullptr };
  const VarLabel* d_oldRefineFlagLabel{ nullptr };
  const VarLabel* d_refinePatchFlagLabel{ nullptr };

  std::vector<int> d_numStability{};
  std::vector<int> d_numRegrid{};
  std::vector<int> d_numDeleted{};

  bool d_newGrid{ true };
  bool d_regridOnce{ false };
  int d_lastRegridTimestep{ 0 }; ///< The last time the full regridder was
                                 ///< called (grid may not change)
  ///< The last timestep that the dilation was changed
  int d_dilationTimestep{ 3 };
  int d_maxTimestepsBetweenRegrids{ 1 };
  int d_minTimestepsBetweenRegrids{ 1 };
  double d_overheadAverage{ 0.0 };
  double d_amrOverheadLow{ 0.0 };  ///< Percentage low target for AMR overhead
  double d_amrOverheadHigh{ 0.0 }; ///< Percentage high target for AMR overhead

protected:
  bool
  flaggedCellsExist(constCCVariable<int>& flaggedCells,
                    IntVector low,
                    IntVector high);

  IntVector
  Less(const IntVector& a, const IntVector& b);

  IntVector
  Greater(const IntVector& a, const IntVector& b);

  IntVector
  And(const IntVector& a, const IntVector& b);

  IntVector
  Mod(const IntVector& a, const IntVector& b);

  IntVector
  Ceil(const Vector& a);

  void
  problemSetup_BulletProofing(const int k);

  void
  GetFlaggedCells(const GridP& origGrid, int levelIdx, DataWarehouse* dw);

  void
  initFilter(CCVariable<int>& filter, FilterType ft, IntVector& depth);
};

} // End namespace Uintah

#endif // End UINTAH_CCA_COMPONENTS_REGRIDDERS_REGRIDDERCOMMON_H
