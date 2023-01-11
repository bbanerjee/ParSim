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

#ifndef CORE_GRID_LEVEL_H
#define CORE_GRID_LEVEL_H

#include <CCA/Ports/LoadBalancer.h>

#include <Core/Containers/OffsetArray1.h>
#include <Core/Disclosure/TypeDescription.h>

#ifdef max
// some uintah 3p utilities define max, so undefine it before BBox chokes on it.
#undef max
#endif
#include <Core/Geometry/BBox.h>

#include <Core/Grid/Grid.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/Variables/ComputeSet.h>

#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Point.h>

#include <Core/Parallel/CrowdMonitor.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <Core/Util/Handle.h>
#include <Core/Util/RefCounted.h>

#include <map>
#include <vector>

namespace Uintah {

constexpr double LO = std::numeric_limits<double>::lowest();
constexpr double HI = std::numeric_limits<double>::max();

class PatchBVH;
class BoundCondBase;
class Box;
class Patch;
class Task;

//! Level class
/*!
   Just a container class that manages a set of Patches that
   make up this level.
*/

class Level : public RefCounted
{
public:
  Level(Grid* grid,
        const Point& anchor,
        const Vector& dcell,
        int index,
        IntVector refinementRatio,
        int id = -1);

  virtual ~Level();

  // eliminate copy, assignment and move
  Level(const Level&) = delete;
  Level(Level&&)      = delete;

  Level&
  operator=(const Level&) = delete;
  Level&
  operator=(Level&&) = delete;

  void
  setPatchDistributionHint(const IntVector& patchDistribution);

  void
  setBCTypes();

  using patchIterator       = std::vector<Patch*>::iterator;
  using const_patchIterator = std::vector<Patch*>::const_iterator;

  const_patchIterator
  patchesBegin() const;

  const_patchIterator
  patchesEnd() const;

  patchIterator
  patchesBegin();

  patchIterator
  patchesEnd();

  const Patch*
  getPatch(int index) const
  {
    return d_realPatches[index];
  }

  // go through the virtual ones too
  const_patchIterator
  allPatchesBegin() const;

  const_patchIterator
  allPatchesEnd() const;

  Patch*
  addPatch(const IntVector& extraLowIndex,
           const IntVector& extraHighIndex,
           const IntVector& lowIndex,
           const IntVector& highIndex,
           Grid* grid);

  Patch*
  addPatch(const IntVector& extraLowIndex,
           const IntVector& extraHighIndex,
           const IntVector& lowIndex,
           const IntVector& highIndex,
           Grid* grid,
           int ID);

  // Move up and down the hierarchy
  const LevelP&
  getCoarserLevel() const;

  const LevelP&
  getFinerLevel() const;

  bool
  hasCoarserLevel() const;

  bool
  hasFinerLevel() const;

  IntVector
  mapNodeToCoarser(const IntVector& idx) const;

  IntVector
  mapNodeToFiner(const IntVector& idx) const;

  IntVector
  mapCellToCoarser(const IntVector& idx, int level_offset = 1) const;

  IntVector
  mapCellToFiner(const IntVector& idx) const;

  IntVector
  mapCellToFinest(const IntVector& idx) const;

  IntVector
  mapCellToFinestNoAdjustments(const IntVector& idx) const;

  // Find a patch containing the point, return 0 if non exists
  const Patch*
  getPatchFromPoint(const Point&, const bool includeExtraCells) const;

  // Find a patch containing the cell or node, return 0 if non exists
  const Patch*
  getPatchFromIndex(const IntVector&, const bool includeExtraCells) const;

  void
  finalizeLevel();

  void
  finalizeLevel(bool periodicX, bool periodicY, bool periodicZ);

  void
  assignBCS(const ProblemSpecP& ps, LoadBalancer* lb);

  int
  numPatches() const;

  long
  totalCells() const;

  long
  getTotalCellsInRegion(const TypeDescription::Type varType,
                        const IntVector& boundaryLayer,
                        const IntVector& lowIndex,
                        const IntVector& highIndex) const;

  IntVector
  nCellsPatch_max() const;

  void
  getSpatialRange(BBox& b) const
  {
    b.extend(d_spatial_range);
  };

  void
  getInteriorSpatialRange(BBox& b) const
  {
    b.extend(d_int_spatial_range);
  };

  // methods to identify if this is non-cubic level
  bool
  isNonCubic() const
  {
    return d_isNonCubicDomain;
  };

  void
  findIndexRange(IntVector& lowIndex, IntVector& highIndex) const
  {
    findNodeIndexRange(lowIndex, highIndex);
  }

  void
  findNodeIndexRange(IntVector& lowIndex, IntVector& highIndex) const;

  void
  findCellIndexRange(IntVector& lowIndex, IntVector& highIndex) const;

  void
  findInteriorIndexRange(IntVector& lowIndex, IntVector& highIndex) const
  {
    findInteriorNodeIndexRange(lowIndex, highIndex);
  }

  void
  findInteriorNodeIndexRange(IntVector& lowIndex, IntVector& highIndex) const;

  void
  findInteriorCellIndexRange(IntVector& lowIndex, IntVector& highIndex) const;

  void
  computeVariableExtents(const TypeDescription::Type TD,
                         IntVector& lo,
                         IntVector& hi) const;

  void
  performConsistencyCheck() const;

  GridP
  getGrid() const;

  const LevelP&
  getRelativeLevel(int offset) const;

  // Grid spacing
  Vector
  dCell() const
  {
    return d_dcell;
  }

  /**
   * Returns the cell volume dx*dy*dz. This will not work for stretched grids.
   */
  double
  cellVolume() const
  {
    return d_dcell.x() * d_dcell.y() * d_dcell.z();
  }

  /**
   * Returns the cell area dx*dy, dx*dz, or dy*dz. This will not work for
   * stretched grids.
   */
  double
  cellArea(Vector unitNormal) const
  {
    Vector areas(d_dcell.y() * d_dcell.z(),
                 d_dcell.x() * d_dcell.z(),
                 d_dcell.x() * d_dcell.y());
    return Dot(areas, unitNormal);
  }

  Point
  getAnchor() const
  {
    return d_anchor;
  }

  void
  setExtraCells(const IntVector& ec);

  IntVector
  getExtraCells() const
  {
    return d_extraCells;
  }

  Point
  getNodePosition(const IntVector&) const;

  Point
  getCellPosition(const IntVector&) const;

  IntVector
  getCellIndex(const Point&) const;

  Point
  positionToIndex(const Point&) const;

  Box
  getBox(const IntVector&, const IntVector&) const;

  using selectType = std::vector<const Patch*>;

  void
  selectPatches(const IntVector&,
                const IntVector&,
                selectType&,
                bool withExtraCells = false,
                bool cache          = true) const;

  bool
  containsPointIncludingExtraCells(const Point&) const;

  bool
  containsPoint(const Point&) const;

  bool
  containsCell(const IntVector&) const;

  // IntVector whose elements are each 1 or 0 specifying whether there
  // are periodic boundaries in each dimension (1 means periodic).
  IntVector
  getPeriodicBoundaries() const
  {
    return d_periodicBoundaries;
  }

  // The eachPatch() function returns a PatchSet containing patches on
  // this level with one patch per PatchSubSet.  Eg: { {1}, {2}, {3} }
  const PatchSet*
  eachPatch() const;

  const PatchSet*
  allPatches() const;

  const Patch*
  selectPatchForCellIndex(const IntVector& idx) const;

  const Patch*
  selectPatchForNodeIndex(const IntVector& idx) const;

  // getID() returns a unique identifier so if the grid is rebuilt the new
  // levels will have different id numbers (like a  serial number).
  inline int
  getID() const
  {
    return d_id;
  }

  // getIndex() returns the relative position of the level - 0 is coarsest, 1 is
  // next and so forth.
  inline int
  getIndex() const
  {
    return d_index;
  }

  inline IntVector
  getRefinementRatio() const
  {
    return d_refinementRatio;
  }

  int
  getRefinementRatioMaxDim() const;

  friend std::ostream&
  operator<<(std::ostream& out, const Uintah::Level& level);

  //  overlapping patches:  Used to keep track of patches that overlap in
  //  non-cubic levels
  struct overlap
  {
    std::pair<int, int> patchIDs{ -9, -9 };      // overlapping patch IDs
    IntVector lowIndex{ IntVector(-9, -9, -9) }; // low/high index of overlap
    IntVector highIndex{ IntVector(-9, -9, -9) };
  };

  // for a set of patches and region return the min/max number of overlapping
  // cells
  std::pair<int, int>
  getOverlapCellsInRegion(const selectType& patches,
                          const IntVector& regionLow,
                          const IntVector& regionHigh) const;

private:
  Grid* d_grid{ nullptr };
  Point d_anchor{};
  Vector d_dcell{};

  // the spatial range of the level
  BBox d_spatial_range{ Point(LO, LO, LO), Point(HI, HI, HI) };
  BBox d_int_spatial_range{ Point(LO, LO, LO), Point(HI, HI, HI) };

  bool d_isNonCubicDomain{ false }; // is level non cubic level
  void
  setIsNonCubicLevel();

  bool d_finalized{ false };
  int d_index{}; // number of the level
  IntVector d_patchDistribution{ -1, -1, -1 };
  IntVector d_periodicBoundaries{ 0, 0, 0 };

  PatchSet* d_each_patch{ nullptr };
  PatchSet* d_all_patches{ nullptr };

  long d_totalCells{ 0 };
  IntVector d_extraCells{ IntVector(0, 0, 0) };
  IntVector d_numcells_patch_max{ IntVector(0, 0, 0) };

  std::vector<Patch*> d_realPatches{};           // only real patches
  std::vector<Patch*> d_virtualAndRealPatches{}; // real and virtual

  int d_id{};
  IntVector d_refinementRatio{};

  class IntVectorCompare
  {
  public:
    bool
    operator()(const std::pair<IntVector, IntVector>& a,
               const std::pair<IntVector, IntVector>& b) const
    {
      return (a.first < b.first) ||
             (!(b.first < a.first) && a.second < b.second);
    }
  };

  using selectCache = std::map<std::pair<IntVector, IntVector>,
                               std::vector<const Patch*>,
                               IntVectorCompare>;
  mutable selectCache d_selectCache; // we like const Levels in most places :)

  PatchBVH* d_bvh{ nullptr };

  // overlapping patches
  std::map<std::pair<int, int>, overlap> d_overLapPatches{};
  void
  setOverlappingPatches();
};

const Level*
getLevel(const PatchSubset* subset);

const Level*
getLevel(const PatchSet* set);

const LevelP&
getLevelP(const PatchSubset* subset);

} // End namespace Uintah

#endif // CORE_GRID_LEVEL_H
