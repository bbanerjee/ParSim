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

#include <Core/Grid/Level.h>

#include <Core/Exceptions/InvalidGrid.h>
#include <Core/Exceptions/ProblemSetupException.h>

#include <Core/Geometry/BBox.h>

#include <Core/Grid/BoundaryConditions/BoundCondReader.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/PatchBVH/PatchBVH.h>

#include <Core/Malloc/Allocator.h>

#include <Core/Math/MiscMath.h>

#include <Core/OS/ProcessInfo.h> // For Memory Check

#include <Core/Parallel/CrowdMonitor.h>
#include <Core/Parallel/MasterLock.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/ProcessorGroup.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <Core/Util/DOUT.hpp>
#include <Core/Util/DebugStream.h>
#include <Core/Util/FancyAssert.h>
#include <Core/Util/Handle.h>

#include <Core/Util/ProgressiveWarning.h>
#include <Core/Util/Timers/Timers.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>

namespace {

std::atomic<int32_t> g_ids{ 0 };
Uintah::MasterLock g_patch_cache_mutex{};

Uintah::Dout g_bc_dbg{ "BCTypes", "Level", "Level BC debug info", false };
Uintah::Dout g_rg_times{ "RGTimesLevel",
                         "Level",
                         "Level regridder timing info",
                         false };

}

namespace Uintah {

Level::Level(Grid* grid,
             const Point& anchor,
             const Vector& dcell,
             int index,
             IntVector refinementRatio,
             int id /*=-1*/)
  : d_grid{ grid }
  , d_anchor{ anchor }
  , d_dcell{ dcell }
  , d_index{ index }
  , d_id{ id }
  , d_refinementRatio{ refinementRatio }
{
  if (d_id == -1) {
    d_id = g_ids.fetch_add(1, std::memory_order_relaxed);
  } else if (d_id >= g_ids) {
    g_ids.store(d_id + 1, std::memory_order_relaxed);
  }
}

Level::~Level()
{
  // Delete all of the patches managed by this level
  for (auto iter = d_virtualAndRealPatches.begin();
       iter != d_virtualAndRealPatches.end();
       ++iter) {
    delete *iter;
  }

  d_bvh.reset(nullptr);

  if (d_each_patch && d_each_patch->removeReference()) {
    delete d_each_patch;
  }
  if (d_all_patches && d_all_patches->removeReference()) {
    delete d_all_patches;
  }
}

void
Level::setPatchDistributionHint(const IntVector& hint)
{
  if (d_patchDistribution.x() == -1) {
    d_patchDistribution = hint;
  } else {
    // Called more than once, we have to punt
    d_patchDistribution = IntVector(-2, -2, 2);
  }
}

Level::const_patchIterator
Level::patchesBegin() const
{
  return d_realPatches.begin();
}

Level::const_patchIterator
Level::patchesEnd() const
{
  return d_realPatches.end();
}

Level::patchIterator
Level::patchesBegin()
{
  return d_realPatches.begin();
}

Level::patchIterator
Level::patchesEnd()
{
  return d_realPatches.end();
}

Level::const_patchIterator
Level::allPatchesBegin() const
{
  return d_virtualAndRealPatches.begin();
}

Level::const_patchIterator
Level::allPatchesEnd() const
{
  return d_virtualAndRealPatches.end();
}

Patch*
Level::addPatch(const IntVector& lowIndex,
                const IntVector& highIndex,
                const IntVector& inLowIndex,
                const IntVector& inHighIndex,
                Grid* grid)
{
  Patch* r = scinew Patch(this,
                          lowIndex,
                          highIndex,
                          inLowIndex,
                          inHighIndex,
                          getIndex());
  r->setGrid(grid);
  d_realPatches.push_back(r);
  d_virtualAndRealPatches.push_back(r);

  d_int_spatial_range.extend(r->getBox().lower());
  d_int_spatial_range.extend(r->getBox().upper());

  d_spatial_range.extend(r->getExtraBox().lower());
  d_spatial_range.extend(r->getExtraBox().upper());
  return r;
}

Patch*
Level::addPatch(const IntVector& lowIndex,
                const IntVector& highIndex,
                const IntVector& inLowIndex,
                const IntVector& inHighIndex,
                Grid* grid,
                int ID)
{
  Patch* r = scinew
    Patch(this, lowIndex, highIndex, inLowIndex, inHighIndex, getIndex(), ID);
  r->setGrid(grid);
  d_realPatches.push_back(r);
  d_virtualAndRealPatches.push_back(r);

  d_int_spatial_range.extend(r->getBox().lower());
  d_int_spatial_range.extend(r->getBox().upper());

  d_spatial_range.extend(r->getExtraBox().lower());
  d_spatial_range.extend(r->getExtraBox().upper());
  return r;
}

const Patch*
Level::getPatchFromPoint(const Point& p, const bool includeExtraCells) const
{
  selectType patch;
  IntVector c = getCellIndex(p);
  // point is within the bounding box so query the bvh
  d_bvh->query(c, c + IntVector(1, 1, 1), patch, includeExtraCells);

  if (patch.size() == 0) {
    return 0;
  }

  ASSERT(patch.size() == 1);
  return patch[0];
}

const Patch*
Level::getPatchFromIndex(const IntVector& c, const bool includeExtraCells) const
{
  selectType patch;

  // Point is within the bounding box so query the bvh.
  d_bvh->query(c, c + IntVector(1, 1, 1), patch, includeExtraCells);

  if (patch.size() == 0) {
    return 0;
  }

  ASSERT(patch.size() == 1);
  return patch[0];
}

int
Level::numPatches() const
{
  return static_cast<int>(d_realPatches.size());
}

void
Level::performConsistencyCheck() const
{
  if (!d_finalized) {
    SCI_THROW(InvalidGrid(
      "Consistency check cannot be performed until Level is finalized",
      __FILE__,
      __LINE__));
  }

  for (size_t i = 0; i < d_virtualAndRealPatches.size(); i++) {
    Patch* r = d_virtualAndRealPatches[i];
    r->performConsistencyCheck();
  }

  // This is O(n^2) - we should fix it someday if it ever matters
  //   This checks that patches do not overlap
  for (size_t i = 0; i < d_virtualAndRealPatches.size(); i++) {
    Patch* r1 = d_virtualAndRealPatches[i];

    for (size_t j = i + 1; j < d_virtualAndRealPatches.size(); j++) {
      Patch* r2 = d_virtualAndRealPatches[j];
      Box b1    = getBox(r1->getCellLowIndex(), r1->getCellHighIndex());
      Box b2    = getBox(r2->getCellLowIndex(), r2->getCellHighIndex());

      if (b1.overlaps(b2)) {
        std::cerr << "r1: " << *r1 << '\n';
        std::cerr << "r2: " << *r2 << '\n';
        SCI_THROW(InvalidGrid("Two patches overlap", __FILE__, __LINE__));
      }
    }
  }
}

void
Level::findNodeIndexRange(IntVector& lowIndex, IntVector& highIndex) const
{
  Vector l = (d_spatial_range.min() - d_anchor) / d_dcell;
  Vector h = (d_spatial_range.max() - d_anchor) / d_dcell + Vector(1, 1, 1);

  lowIndex  = roundNearest(l);
  highIndex = roundNearest(h);
}

void
Level::findCellIndexRange(IntVector& lowIndex, IntVector& highIndex) const
{
  Vector l = (d_spatial_range.min() - d_anchor) / d_dcell;
  Vector h = (d_spatial_range.max() - d_anchor) / d_dcell;

  lowIndex  = roundNearest(l);
  highIndex = roundNearest(h);
}

void
Level::findInteriorCellIndexRange(IntVector& lowIndex,
                                  IntVector& highIndex) const
{
  Vector l = (d_int_spatial_range.min() - d_anchor) / d_dcell;
  Vector h = (d_int_spatial_range.max() - d_anchor) / d_dcell;

  lowIndex  = roundNearest(l);
  highIndex = roundNearest(h);
}

void
Level::findInteriorNodeIndexRange(IntVector& lowIndex,
                                  IntVector& highIndex) const
{
  Vector l = (d_int_spatial_range.min() - d_anchor) / d_dcell;
  Vector h = (d_int_spatial_range.max() - d_anchor) / d_dcell + Vector(1, 1, 1);

  lowIndex  = roundNearest(l);
  highIndex = roundNearest(h);
}

// Compute the variable extents for this variable type
void
Level::computeVariableExtents(const TypeDescription::Type type,
                              IntVector& lo,
                              IntVector& hi) const
{
  IntVector CCLo;
  IntVector CCHi;
  findCellIndexRange(CCLo, CCHi);

  // Fix me: better way to calc this var? Or to use it below?
  IntVector not_periodic(!d_periodicBoundaries[0],
                         !d_periodicBoundaries[1],
                         !d_periodicBoundaries[2]);

  switch (type) {
    case TypeDescription::Type::CCVariable:
    case TypeDescription::Type::ParticleVariable:
      lo = CCLo;
      hi = CCHi;
      break;
    case TypeDescription::Type::SFCXVariable:
      lo = CCLo;
      hi = CCHi + (IntVector(1, 0, 0) * not_periodic);
      break;
    case TypeDescription::Type::SFCYVariable:
      lo = CCLo;
      hi = CCHi + (IntVector(0, 1, 0) * not_periodic);
      break;
    case TypeDescription::Type::SFCZVariable:
      lo = CCLo;
      hi = CCHi + (IntVector(0, 0, 1) * not_periodic);
      break;
    case TypeDescription::Type::NCVariable:
      // Dav's fix: findInteriorCellIndexRange( lo, hi );
      findNodeIndexRange(lo, hi);
      break;
    default:
      std::string me = TypeDescription::toString(type);
      throw InternalError(
        "  ERROR: Level::computeVariableExtents type description (" + me +
          ") not supported",
        __FILE__,
        __LINE__);
  }
}

long
Level::totalCells() const
{
  return d_totalCells;
}

long
Level::getTotalCellsInRegion(const TypeDescription::Type varType,
                             const IntVector& boundaryLayer,
                             const IntVector& lowIndex,
                             const IntVector& highIndex) const
{

  // Not all simulations are cubic.  Some simulations might be L shaped, or T
  // shaped, or + shaped, etc. It is not enough to simply do a high - low to
  // figure out the amount of simulation cells.  We instead need to go all
  // patches and see if they exist in this range.  If so, add up their cells.
  // This process is similar to how d_totalCells is computed in
  // Level::finalizeLevel().

  long cellsInRegion = 0;
  if (d_isNonCubicDomain == false) {
    IntVector diff(highIndex - lowIndex);
    cellsInRegion += diff.x() * diff.y() * diff.z();
    return cellsInRegion;
  }

  Patch::VariableBasis basis = Patch::translateTypeToBasis(varType, false);

  // scjmc - using also virtual patches to take into account periodic boundaries
  // on non cubic levels such as refined amr levels
  for (size_t i = 0; i < d_virtualAndRealPatches.size(); ++i) {
    IntVector patchLow =
      d_virtualAndRealPatches[i]->getExtraLowIndex(basis, boundaryLayer);
    IntVector patchHigh =
      d_virtualAndRealPatches[i]->getExtraHighIndex(basis, boundaryLayer);

    if (doesIntersect(lowIndex, highIndex, patchLow, patchHigh)) {

      IntVector regionLo = Uintah::Max(lowIndex, patchLow);
      IntVector regionHi = Uintah::Min(highIndex, patchHigh);
      IntVector diff(regionHi - regionLo);
      cellsInRegion += diff.x() * diff.y() * diff.z();
    }
  }

  return cellsInRegion;
}

IntVector
Level::nCellsPatch_max() const // used by PIDX
{
  return d_numcells_patch_max;
}

void
Level::setExtraCells(const IntVector& ec)
{
  d_extraCells = ec;
}

GridP
Level::getGrid() const
{
  return d_grid;
}

const LevelP&
Level::getRelativeLevel(int offset) const
{
  return d_grid->getLevel(d_index + offset);
}

Uintah::Point
Level::getNodePosition(const IntVector& v) const
{
  return d_anchor + d_dcell * v;
}

Point
Level::getCellPosition(const IntVector& v) const
{
  return d_anchor + d_dcell * v + d_dcell * 0.5;
}

IntVector
Level::getCellIndex(const Point& p) const
{
  Vector v((p - d_anchor) / d_dcell);
  return IntVector(RoundDown(v.x()), RoundDown(v.y()), RoundDown(v.z()));
}

Point
Level::positionToIndex(const Point& p) const
{
  return Point((p - d_anchor) / d_dcell);
}

void
Level::selectPatches(const IntVector& low,
                     const IntVector& high,
                     selectType& neighbors,
                     bool withExtraCells,
                     bool cache) const
{
  if (cache) {
    std::lock_guard<Uintah::MasterLock> cache_lock(g_patch_cache_mutex);

    // look it up in the cache first
    selectCache::const_iterator iter =
      d_selectCache.find(std::make_pair(low, high));

    if (iter != d_selectCache.end()) {
      const vector<const Patch*>& cache = iter->second;
      for (unsigned i = 0; i < cache.size(); i++) {
        neighbors.push_back(cache[i]);
      }
      return;
    }
    ASSERT(neighbors.size() == 0);
  }

  d_bvh->query(low, high, neighbors, withExtraCells);

  std::sort(neighbors.begin(), neighbors.end(), Patch::Compare());

#ifdef CHECK_SELECT
  // Double-check the more advanced selection algorithms against the
  // slow (exhaustive) one.
  std::vector<const Patch*> tneighbors;
  for (const_patchIterator iter = d_virtualAndRealPatches.begin();
       iter != d_virtualAndRealPatches.end();
       iter++) {
    const Patch* patch = *iter;

    IntVector l = Max(patch->getCellLowIndex(), low);
    IntVector u = Min(patch->getCellHighIndex(), high);
    if (u.x() > l.x() && u.y() > l.y() && u.z() > l.z()) {
      tneighbors.push_back(*iter);
    }
  }
  ASSERTEQ(neighbors.size(), tneighbors.size());

  std::sort(tneighbors.begin(), tneighbors.end(), Patch::Compare());
  for (int i = 0; i < (int)neighbors.size(); i++) {
    ASSERT(neighbors[i] == tneighbors[i]);
  }
#endif

  if (cache) {
    std::lock_guard<Uintah::MasterLock> cache_lock(g_patch_cache_mutex);

    // put it in the cache - start at orig_size in case there was something in
    // neighbors before this query
    std::vector<const Patch*>& cache = d_selectCache[std::make_pair(low, high)];
    cache.reserve(6); // don't reserve too much to save memory, not too little
                      // to avoid too much reallocation
    for (size_t i = 0; i < neighbors.size(); i++) {
      cache.push_back(neighbors[i]);
    }
  }
}

bool
Level::containsPointIncludingExtraCells(const Point& p) const
{
  bool includeExtraCells = true;
  return getPatchFromPoint(p, includeExtraCells) != nullptr;
}

bool
Level::containsPoint(const Point& p) const
{
  bool includeExtraCells = false;
  const Patch* patch     = getPatchFromPoint(p, includeExtraCells);
  return patch != nullptr;
}

bool
Level::containsCell(const IntVector& idx) const
{
  bool includeExtraCells = false;
  const Patch* patch     = getPatchFromIndex(idx, includeExtraCells);
  return patch != nullptr;
}

/*!  This method determines if a level is nonCubic or if there are any missing
 * patches. Algorithm: 1) The total volume of the patches must equal the volume
 * of the level. The volume of the level is defined by the bounding box.
 */
void
Level::setIsNonCubicLevel()
{
  double patchesVol = 0.0;

  // loop over all patches and sum the patch's volume
  for (size_t p = 0; p < d_realPatches.size(); ++p) {
    auto patch = d_realPatches[p];

    Box box = patch->getBox();
    Vector sides((box.upper().x() - box.lower().x()),
                 (box.upper().y() - box.lower().y()),
                 (box.upper().z() - box.lower().z()));

    double volume = sides.x() * sides.y() * sides.z();

    patchesVol += volume;
  }

  // compute the level's volume from the bounding box
  Point loPt = d_int_spatial_range.min();
  Point hiPt = d_int_spatial_range.max();

  Vector levelSides((hiPt.x() - loPt.x()),
                    (hiPt.y() - loPt.y()),
                    (hiPt.z() - loPt.z()));

  double levelVol = levelSides.x() * levelSides.y() * levelSides.z();

  d_isNonCubicDomain = false;
  double fuzz        = 0.5 * cellVolume(); // 100 * DBL_EPSILON;
  if (std::fabs(patchesVol - levelVol) > fuzz) {
    d_isNonCubicDomain = true;
  }
}

//  Loop through all patches on the level and if they overlap with each other
//  then store that info. You need this information when a level has a patch
//  distribution with inside corners
void
Level::setOverlappingPatches()
{
  if (d_isNonCubicDomain == false) { //  for cubic domains just return
    return;
  }

  for (size_t i = 0; i < d_realPatches.size(); ++i) {
    auto patch       = d_realPatches[i];
    Box b1           = patch->getExtraBox();
    IntVector lowEC  = patch->getExtraCellLowIndex();
    IntVector highEC = patch->getExtraCellHighIndex();

    bool includeExtraCells = true;
    Patch::selectType neighborPatches;
    selectPatches(lowEC, highEC, neighborPatches, includeExtraCells);

    for (size_t j = 0; j < neighborPatches.size(); ++j) {
      auto neighborPatch = neighborPatches[j];

      if (patch != neighborPatch) {

        Box b2 = neighborPatch->getExtraBox();

        //  Are the patches overlapping?
        if (b1.overlaps(b2)) {

          IntVector nLowEC  = neighborPatch->getExtraCellLowIndex();
          IntVector nHighEC = neighborPatch->getExtraCellHighIndex();

          // find intersection of the patches
          IntVector intersectLow  = Max(lowEC, nLowEC);
          IntVector intersectHigh = Min(highEC, nHighEC);

          // create overlap
          overlap newOverLap;
          std::pair<int, int> patchIds =
            std::make_pair(patch->getID(), neighborPatch->getID());
          newOverLap.patchIDs  = patchIds;
          newOverLap.lowIndex  = intersectLow;
          newOverLap.highIndex = intersectHigh;

          // only add unique values to the map.
          auto result = d_overLapPatches.find(patchIds);
          if (result == d_overLapPatches.end()) {
            d_overLapPatches[patchIds] = newOverLap;
          }
        }
      }
    }
  }

  // debugging
#if 0
  for (auto itr=d_overLapPatches.begin(); itr!=d_overLapPatches.end(); ++itr) {
    std::pair<int,int> patchIDs = itr->first;
    overlap me = itr->second;
    std::cout <<  " overlapping patches, Level: " << getIndex() << " patches: ("  << patchIDs.first << ", " << patchIDs.second <<"), low: " << me.lowIndex << " high: " << me.highIndex << std:: endl;
  }
#endif
}

//  This method returns the min/max number of overlapping patch cells that are
//  within a specified region.  Patches will overlap when the domain is
//  non-cubic
std::pair<int, int>
Level::getOverlapCellsInRegion(const selectType& patches,
                               const IntVector& regionLow,
                               const IntVector& regionHigh) const
{
  // cubic domains never have overlapping patches, just return
  if (d_isNonCubicDomain == false) {
    return std::make_pair(-9, -9);
  }

  int minOverlapCells   = std::numeric_limits<int>::max();
  int totalOverlapCells = 0;

  // loop over patches in this region
  for (size_t i = 0; i < patches.size(); ++i) {
    for (size_t j = i + 1; j < patches.size();
         ++j) { //  the i+1 prevents examining the transposed key pairs, i.e.
                //  8,9 and 9,8

      int Id                       = patches[i]->getID();
      int neighborId               = patches[j]->getID();
      std::pair<int, int> patchIds = std::make_pair(Id, neighborId);

      auto result = d_overLapPatches.find(patchIds);

      if (result == d_overLapPatches.end()) {
        continue;
      }

      overlap ol = result->second;

      // does the overlapping patches intersect with the region extents?
      if (doesIntersect(ol.lowIndex, ol.highIndex, regionLow, regionHigh)) {
        IntVector intrsctLow  = Uintah::Max(ol.lowIndex, regionLow);
        IntVector intrsctHigh = Uintah::Min(ol.highIndex, regionHigh);

        IntVector diff    = IntVector(intrsctHigh - intrsctLow);
        int nOverlapCells = std::abs(diff.x() * diff.y() * diff.z());

        minOverlapCells = std::min(minOverlapCells, nOverlapCells);
        totalOverlapCells += nOverlapCells;

#if 0 // debugging
        std::cout << "  getOverlapCellsInRegion  patches: " << ol.patchIDs.first << ", " << ol.patchIDs.second
                  << "\n   region:      " << regionLow   << ",              " << regionHigh        
                  << "\n   ol.low:      " << ol.lowIndex << " ol.high:      " << ol.highIndex 
                  << "\n   intrsct.low: " << intrsctLow  << " intrsct.high: " << intrsctHigh 
                  << " overlapCells: " << nOverlapCells  << " minOverlapCells: " << minOverlapCells << " totalOverlapCells: " << totalOverlapCells << std::endl;
#endif
      }
    }
  }
  std::pair<int, int> overLapCells_minMax =
    std::make_pair(minOverlapCells, totalOverlapCells);
  return overLapCells_minMax;
}

void
Level::finalizeLevel()
{
  d_each_patch = scinew PatchSet();
  d_each_patch->addReference();

  // The compute set requires an array const Patch*, we must copy
  // d_realPatches
  std::vector<const Patch*> tmp_patches(d_realPatches.size());
  for (size_t i = 0; i < d_realPatches.size(); i++) {
    tmp_patches[i] = d_realPatches[i];
  }

  d_each_patch->addEach(tmp_patches);

  d_all_patches = scinew PatchSet();
  d_all_patches->addReference();
  d_all_patches->addAll(tmp_patches);

  d_all_patches->sortSubsets();
  std::sort(d_realPatches.begin(), d_realPatches.end(), Patch::Compare());

  // determines and sets the boundary conditions for the patches
  setBCTypes();

  // compute the number of cells in the level
  d_totalCells = 0;
  for (size_t i = 0; i < d_realPatches.size(); i++) {
    d_totalCells += d_realPatches[i]->getNumExtraCells();
  }

  // compute the max number of cells over all patches  Needed by PIDX
  d_numcells_patch_max = IntVector(0, 0, 0);
  int nCells           = 0;
  for (size_t i = 0; i < d_realPatches.size(); ++i) {

    if (d_realPatches[i]->getNumExtraCells() > nCells) {
      IntVector lo         = d_realPatches[i]->getExtraCellLowIndex();
      IntVector hi         = d_realPatches[i]->getExtraCellHighIndex();
      d_numcells_patch_max = hi - lo;
    }
  }

  // compute and store the spatial ranges now that BCTypes are set
  for (size_t i = 0; i < d_realPatches.size(); i++) {
    auto r = d_realPatches[i];

    d_spatial_range.extend(r->getExtraBox().lower());
    d_spatial_range.extend(r->getExtraBox().upper());
  }

  // determine if this level is cubic
  setIsNonCubicLevel();

  // Loop through all patches and find the patches that overlap.  Needed
  // when patches layouts have inside corners.
  setOverlappingPatches();
}

void
Level::finalizeLevel(bool periodicX, bool periodicY, bool periodicZ)
{
  // set each_patch and all_patches before creating virtual patches
  d_each_patch = scinew PatchSet();
  d_each_patch->addReference();

  // The compute set requires an array const Patch*, we must copy
  // d_realPatches
  std::vector<const Patch*> tmp_patches(d_realPatches.size());

  for (size_t i = 0; i < d_realPatches.size(); i++) {
    tmp_patches[i] = d_realPatches[i];
  }

  d_each_patch->addEach(tmp_patches);

  d_all_patches = scinew PatchSet();
  d_all_patches->addReference();
  d_all_patches->addAll(tmp_patches);

  BBox bbox;

  if (d_index > 0) {
    d_grid->getLevel(0)->getInteriorSpatialRange(bbox);
  } else {
    getInteriorSpatialRange(bbox);
  }

  Box domain(bbox.min(), bbox.max());
  Vector vextent = positionToIndex(bbox.max()) - positionToIndex(bbox.min());
  IntVector extent(static_cast<int>(rint(vextent.x())),
                   static_cast<int>(rint(vextent.y())),
                   static_cast<int>(rint(vextent.z())));

  d_periodicBoundaries =
    IntVector(periodicX ? 1 : 0, periodicY ? 1 : 0, periodicZ ? 1 : 0);
  IntVector periodicBoundaryRange = d_periodicBoundaries * extent;

  int x, y, z;
  for (size_t i = 0; i < tmp_patches.size(); i++) {

    for (x = -d_periodicBoundaries.x(); x <= d_periodicBoundaries.x(); x++) {
      for (y = -d_periodicBoundaries.y(); y <= d_periodicBoundaries.y(); y++) {
        for (z = -d_periodicBoundaries.z(); z <= d_periodicBoundaries.z();
             z++) {

          IntVector offset = IntVector(x, y, z) * periodicBoundaryRange;
          if (offset == IntVector(0, 0, 0)) {
            continue;
          }

          Box box = getBox(tmp_patches[i]->getExtraCellLowIndex() + offset -
                             IntVector(1, 1, 1),
                           tmp_patches[i]->getExtraCellHighIndex() + offset +
                             IntVector(1, 1, 1));

          if (box.overlaps(domain)) {
            Patch* newPatch = tmp_patches[i]->createVirtualPatch(offset);
            d_virtualAndRealPatches.push_back(newPatch);
          }
        }
      }
    }
  }

  d_all_patches->sortSubsets();
  std::sort(d_realPatches.begin(), d_realPatches.end(), Patch::Compare());
  std::sort(d_virtualAndRealPatches.begin(),
            d_virtualAndRealPatches.end(),
            Patch::Compare());

  setBCTypes();

  // compute the number of cells in the level
  d_totalCells = 0;
  for (int i = 0; i < (int)d_realPatches.size(); i++) {
    d_totalCells += d_realPatches[i]->getNumExtraCells();
  }

  // compute the max number of cells over all patches  Needed by PIDX
  d_numcells_patch_max = IntVector(0, 0, 0);
  int nCells           = 0;

  for (size_t i = 0; i < d_realPatches.size(); ++i) {
    if (d_realPatches[i]->getNumExtraCells() > nCells) {
      IntVector lo = d_realPatches[i]->getExtraCellLowIndex();
      IntVector hi = d_realPatches[i]->getExtraCellHighIndex();

      d_numcells_patch_max = hi - lo;
    }
  }

  // compute and store the spatial ranges now that BCTypes are set
  for (int i = 0; i < (int)d_realPatches.size(); i++) {
    auto r = d_realPatches[i];

    d_spatial_range.extend(r->getExtraBox().lower());
    d_spatial_range.extend(r->getExtraBox().upper());
  }

  // determine if this level is cubic
  setIsNonCubicLevel();

  // Loop through all patches and find the patches that overlap.  Needed
  // when patch layouts have inside corners.
  setOverlappingPatches();
}

void
Level::setBCTypes()
{
  Timers::Simple timer;
  timer.start();

  const int nTimes      = 3;
  double rtimes[nTimes] = { 0 };

  if (d_bvh != nullptr) {
    d_bvh.reset(nullptr);
  }

  d_bvh = std::make_unique<PatchBVH>(d_virtualAndRealPatches);

  rtimes[0] += timer().seconds();
  timer.reset(true);

  ProcessorGroup* myworld = nullptr;
  int numProcs            = 1;
  int rank                = 0;

  if (Parallel::isInitialized()) {
    // only vaango uses Parallel, but anybody else who uses DataArchive to read
    // data does not
    myworld  = Parallel::getRootProcessorGroup();
    numProcs = myworld->nRanks();
    rank     = myworld->myRank();
  }

  std::vector<int> displacements(numProcs, 0);
  std::vector<int> recvcounts(numProcs, 0);

  // create recvcounts and displacement arrays
  int div = d_virtualAndRealPatches.size() / numProcs;
  int mod = d_virtualAndRealPatches.size() % numProcs;

  for (int p = 0; p < numProcs; p++) {
    if (p < mod) {
      recvcounts[p] = div + 1;
    } else {
      recvcounts[p] = div;
    }
  }

  displacements[0] = 0;
  for (int p = 1; p < numProcs; p++) {
    displacements[p] = displacements[p - 1] + recvcounts[p - 1];
  }

  std::vector<unsigned int> bctypes(d_virtualAndRealPatches.size());
  std::vector<unsigned int> mybctypes(recvcounts[rank]);

  int idx;

  patchIterator iter;
  patchIterator startpatch =
    d_virtualAndRealPatches.begin() + displacements[rank];
  patchIterator endpatch = startpatch + recvcounts[rank];

  // for each of my patches
  for (iter = startpatch, idx = 0; iter != endpatch; iter++, idx++) {
    auto patch = *iter;
    // std::cout << "Patch bounding box = " << patch->getExtraBox() << std::endl;
    //  See if there are any neighbors on the 6 faces
    int bitfield = 0;

    for (Patch::FaceType face = Patch::startFace; face <= Patch::endFace;
         face                 = Patch::nextFace(face)) {
      bitfield <<= 2;
      IntVector l, h;

      patch->getFace(face, IntVector(0, 0, 0), IntVector(1, 1, 1), l, h);

      Patch::selectType neighbors;
      selectPatches(l, h, neighbors);

      if (neighbors.size() == 0) {
        if (d_index != 0) {
          // See if there are any patches on the coarse level at that face
          IntVector fineLow, fineHigh;
          patch->getFace(face,
                         IntVector(0, 0, 0),
                         d_refinementRatio,
                         fineLow,
                         fineHigh);
          IntVector coarseLow       = mapCellToCoarser(fineLow);
          IntVector coarseHigh      = mapCellToCoarser(fineHigh);
          const LevelP& coarseLevel = getCoarserLevel();

#if 0
          // add 1 to the corresponding index on the plus edges 
          // because the upper corners are sort of one cell off (don't know why)
          if (d_extraCells.x() != 0 && face == Patch::xplus) {
            coarseLow[0] ++;
            coarseHigh[0]++;
          }
          else if (d_extraCells.y() != 0 && face == Patch::yplus) {
            coarseLow[1] ++;
            coarseHigh[1] ++;
          }
          else if (d_extraCells.z() != 0 && face == Patch::zplus) {
            coarseLow[2] ++;
            coarseHigh[2]++;
          }
#endif
          coarseLevel->selectPatches(coarseLow, coarseHigh, neighbors);

          if (neighbors.size() == 0) {
            bitfield |= Patch::None;
          } else {
            bitfield |= Patch::Coarse;
          }
        } else {
          bitfield |= Patch::None;
        }
      } else {
        bitfield |= Patch::Neighbor;
      }
    }
    mybctypes[idx] = bitfield;
  }

  if (numProcs > 1) {
    // allgather bctypes
    if (mybctypes.size() == 0) {
      Uintah::MPI::Allgatherv(0,
                              0,
                              MPI_UNSIGNED,
                              &bctypes[0],
                              &recvcounts[0],
                              &displacements[0],
                              MPI_UNSIGNED,
                              myworld->getComm());
    } else {
      Uintah::MPI::Allgatherv(&mybctypes[0],
                              mybctypes.size(),
                              MPI_UNSIGNED,
                              &bctypes[0],
                              &recvcounts[0],
                              &displacements[0],
                              MPI_UNSIGNED,
                              myworld->getComm());
    }
  } else {
    bctypes.swap(mybctypes);
  }

  rtimes[1] += timer().seconds();
  timer.reset(true);

  int i;
  // loop through patches
  for (iter = d_virtualAndRealPatches.begin(), i = 0, idx = 0;
       iter != d_virtualAndRealPatches.end();
       iter++, i++) {
    Patch* patch = *iter;

    if (patch->isVirtual()) {
      patch->setLevelIndex(-1);
    } else {
      patch->setLevelIndex(idx++);
    }

    int bitfield = bctypes[i];
    int mask     = 3;

    // loop through faces
    for (int j = 5; j >= 0; j--) {

      int bc_type = bitfield & mask;

      if (rank == 0) {
        switch (bc_type) {
          case Patch::None:
            DOUT(g_bc_dbg,
                 "  Setting Patch " << patch->getID() << " face " << j
                                    << " to None");
            break;
          case Patch::Coarse:
            DOUT(g_bc_dbg,
                 "  Setting Patch " << patch->getID() << " face " << j
                                    << " to Coarse\n");
            break;
          case Patch::Neighbor:
            DOUT(g_bc_dbg,
                 "  Setting Patch " << patch->getID() << " face " << j
                                    << " to Neighbor\n");
            break;
        }
      }
      patch->setBCType(Patch::FaceType(j), Patch::BCType(bc_type));
      bitfield >>= 2;
    }
  }

  //__________________________________
  //  bullet proofing
  for (int dir = 0; dir < 3; dir++) {
    if (d_periodicBoundaries[dir] == 1 && d_extraCells[dir] != 0) {
      std::ostringstream warn;
      warn << "\n \n INPUT FILE ERROR: \n You've specified a periodic boundary "
              "condition on a face with extra cells specified\n"
           << " Please set the extra cells on that face to 0";
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }
  }

  d_finalized = true;

  rtimes[2] += timer().seconds();
  timer.reset(true);

  if (g_rg_times.active()) {
    std::ostringstream mesg;
    double avg[nTimes] = { 0 };
    Uintah::MPI::Reduce(&rtimes,
                        &avg,
                        3,
                        MPI_DOUBLE,
                        MPI_SUM,
                        0,
                        myworld->getComm());

    if (myworld->myRank() == 0) {

      mesg << "SetBCType Avg Times: ";
      for (int i = 0; i < 3; i++) {
        avg[i] /= myworld->nRanks();
        mesg << avg[i] << " ";
      }
      mesg << std::endl;
    }

    double max[nTimes] = { 0 };
    Uintah::MPI::Reduce(&rtimes,
                        &max,
                        3,
                        MPI_DOUBLE,
                        MPI_MAX,
                        0,
                        myworld->getComm());

    if (myworld->myRank() == 0) {
      mesg << "SetBCType Max Times: ";
      for (int i = 0; i < 3; i++) {
        mesg << max[i] << " ";
      }
      mesg << std::endl;
    }
  }

  // recreate BVH with extracells
  if (d_bvh != nullptr) {
    d_bvh.reset(nullptr);
  }
  d_bvh = std::make_unique<PatchBVH>(d_virtualAndRealPatches);
}

void
Level::assignBCS(const ProblemSpecP& grid_ps, LoadBalancer* lb)
{
  ProblemSpecP bc_ps = grid_ps->findBlock("BoundaryConditions");
  if (bc_ps == nullptr) {
    if (Parallel::getMPIRank() == 0) {
      static ProgressiveWarning warn("No BoundaryConditions specified", -1);
      warn.invoke();
    }
    return;
  }

  BoundCondReader reader;
  reader.read(bc_ps, grid_ps, this);

  for (auto& patch : d_virtualAndRealPatches) {

    // If we have a lb, then only apply bcs this processors patches.
    if (lb == 0 ||
        lb->getPatchwiseProcessorAssignment(patch) == Parallel::getMPIRank()) {

      patch->initializeBoundaryConditions();

      for (Patch::FaceType face_side = Patch::startFace;
           face_side <= Patch::endFace;
           face_side = Patch::nextFace(face_side)) {
        if (patch->getBCType(face_side) == Patch::None) {
          patch->setArrayBCValues(face_side,
                                  &(reader.d_BCReaderData[face_side]));
        }
        patch->setInteriorBndArrayBCValues(
          face_side,
          &(reader.d_interiorBndBCReaderData[face_side]));
      } // end of face iterator
    }
  } // end of patch iterator
}

Box
Level::getBox(const IntVector& l, const IntVector& h) const
{
  return Box(getNodePosition(l), getNodePosition(h));
}

const PatchSet*
Level::eachPatch() const
{
  ASSERT(d_each_patch != 0);
  return d_each_patch;
}

const PatchSet*
Level::allPatches() const
{
  ASSERT(d_all_patches != 0);
  return d_all_patches;
}

const Patch*
Level::selectPatchForCellIndex(const IntVector& idx) const
{
  selectType pv;
  IntVector i(1, 1, 1);
  selectPatches(idx - i, idx + i, pv, false, false);

  if (pv.size() == 0) {
    return 0;
  } else {
    for (auto& patch_p : pv) {
      if (patch_p->containsCell(idx)) {
        return patch_p;
      }
    }
  }
  return 0;
}

const Patch*
Level::selectPatchForNodeIndex(const IntVector& idx) const
{
  selectType pv;
  IntVector i(1, 1, 1);
  selectPatches(idx - i, idx + i, pv, false, false);
  if (pv.size() == 0) {
    return 0;
  } else {
    for (auto& patch_p : pv) {
      if (patch_p->containsNode(idx)) {
        return patch_p;
      }
    }
  }
  return 0;
}

const LevelP&
Level::getCoarserLevel() const
{
  return getRelativeLevel(-1);
}

const LevelP&
Level::getFinerLevel() const
{
  return getRelativeLevel(1);
}

bool
Level::hasCoarserLevel() const
{
  return getIndex() > 0;
}

bool
Level::hasFinerLevel() const
{
  return getIndex() < (d_grid->numLevels() - 1);
}

IntVector
Level::mapCellToCoarser(const IntVector& idx, int level_offset) const
{
  IntVector refinementRatio = d_refinementRatio;
  while (--level_offset) {
    refinementRatio =
      refinementRatio *
      d_grid->getLevel(d_index - level_offset)->d_refinementRatio;
  }
  IntVector ratio = idx / refinementRatio;

  // If the fine cell index is negative you must add an offset to get the
  // right coarse cell. -Todd
  IntVector offset(0, 0, 0);
  if (idx.x() < 0 && refinementRatio.x() > 1) {
    offset.x((int)fmod((double)idx.x(), (double)refinementRatio.x()));
  }

  if (idx.y() < 0 && refinementRatio.y() > 1) {
    offset.y((int)fmod((double)idx.y(), (double)refinementRatio.y()));
  }

  if (idx.z() < 0 && refinementRatio.z() > 1) {
    offset.z((int)fmod((double)idx.z(), (double)refinementRatio.z()));
  }
  return ratio + offset;
}

IntVector
Level::mapCellToFiner(const IntVector& idx) const
{
  IntVector r_ratio  = d_grid->getLevel(d_index + 1)->d_refinementRatio;
  IntVector fineCell = idx * r_ratio;

  IntVector offset(0, 0, 0);
  if (idx.x() < 0 && r_ratio.x() > 1) {
    offset.x(1);
  }

  if (idx.y() < 0 && r_ratio.y() > 1) { // If the coarse cell index is negative
    offset.y(1); // you must add an offset to get the right
  }              // fine cell. -Todd

  if (idx.z() < 0 && r_ratio.z() > 1) {
    offset.z(1);
  }
  return fineCell + offset;
}

// Provides the (x-,y-,z-) corner of a fine cell given a coarser coordinate
// If any of the coordinates are negative, assume the fine cell coordiantes
// went too far into the negative and adjust forward by 1
// This adjusting approach means that for L-shaped domains the results are
// not always consistent with what is expected.
//(Note: Does this adjustment mean this only works in a 2:1 refinement ratio
// and only on cubic domains? -- Brad P)
IntVector
Level::mapCellToFinest(const IntVector& idx) const
{

  IntVector r_ratio = IntVector(1, 1, 1);
  for (int i = d_index; i < d_grid->numLevels() - 1; i++) {
    r_ratio = r_ratio * d_grid->getLevel(i + 1)->getRefinementRatio();
  }

  IntVector fineCell = idx * r_ratio;
  IntVector offset(0, 0, 0);
  if (idx.x() < 0 && r_ratio.x() > 1) {
    offset.x(1);
  }

  if (idx.y() < 0 && r_ratio.y() > 1) { // If the coarse cell index is negative
    offset.y(1); // you must add an offset to get the right
  }              // fine cell. -Todd

  if (idx.z() < 0 && r_ratio.z() > 1) {
    offset.z(1);
  }

  return fineCell + offset;
}

// Provides the x-,y-,z- corner of a fine cell given a coarser coordinate.
// This does not attempt to adjust the cell in the + direction if it goes
// negative. It is left up to the caller of this method to determine if
// those coordinates are too far past any level boundary.
IntVector
Level::mapCellToFinestNoAdjustments(const IntVector& idx) const
{

  IntVector r_ratio = IntVector(1, 1, 1);
  for (int i = d_index; i < d_grid->numLevels() - 1; ++i) {
    r_ratio = r_ratio * d_grid->getLevel(i + 1)->getRefinementRatio();
  }
  return idx * r_ratio;
}

// mapNodeToCoarser:
// Example: 1D grid with refinement ratio = 4
//  Coarse Node index: 10                  11
//                     |                   |
//                 ----*----*----*----*----*-----
//                     |                   |
//  Fine Node Index    40   41   42   43   44
//
//  What is returned   10   10   10   10   11

IntVector
Level::mapNodeToCoarser(const IntVector& idx) const
{
  return (idx + d_refinementRatio - IntVector(1, 1, 1)) / d_refinementRatio;
}

// mapNodeToFiner:
// Example: 1D grid with refinement ratio = 4
//  Coarse Node index: 10                  11
//                     |                   |
//                 ----*----*----*----*----*-----
//                     |                   |
//  Fine Node Index    40   41   42   43   44
//
//  What is returned   40                  44

IntVector
Level::mapNodeToFiner(const IntVector& idx) const
{
  return idx * d_grid->getLevel(d_index + 1)->d_refinementRatio;
}

int
Level::getRefinementRatioMaxDim() const
{
  return Max(Max(d_refinementRatio.x(), d_refinementRatio.y()),
             d_refinementRatio.z());
}

const Level*
getLevel(const PatchSubset* subset)
{
  ASSERT(subset->size() > 0);
  const Level* level = subset->get(0)->getLevel();
#if SCI_ASSERTION_LEVEL > 0
  for (int i = 1; i < subset->size(); i++) {
    ASSERT(level == subset->get(i)->getLevel());
  }
#endif
  return level;
}

const LevelP&
getLevelP(const PatchSubset* subset)
{
  ASSERT(subset->size() > 0);
  const LevelP& level = subset->get(0)->getLevelP();
#if SCI_ASSERTION_LEVEL > 0
  for (int i = 1; i < subset->size(); i++) {
    ASSERT(level == subset->get(i)->getLevelP());
  }
#endif
  return level;
}

const Level*
getLevel(const PatchSet* set)
{
  ASSERT(set->size() > 0);
  return getLevel(set->getSubset(0));
}

std::ostream&
operator<<(std::ostream& out, const Level& level)
{
  MasterLock lock;
  lock.lock();
  IntVector lo, hi;
  level.findCellIndexRange(lo, hi);

  out << "(Level " << level.getIndex() << ", numPatches: " << level.numPatches()
      << ", cellIndexRange: " << lo << ", " << hi << ", "
      << *(level.allPatches()) << ")";
  lock.unlock();
  return out;
}

} // end namespace Uintah
