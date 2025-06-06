/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#ifndef __CCA_COMPONENTS_REGRIDDER_TILEDREGRIDDER_H__
#define __CCA_COMPONENTS_REGRIDDER_TILEDREGRIDDER_H__

#include <CCA/Components/Regridder/RegridderCommon.h>

#include <vector>

namespace Uintah {

/**************************************
CLASS
   TiledRegridder

         Tiled Regridding Algorithm

GENERAL INFORMATION

   TiledRegridder.h

   Justin Luitjens
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)

DESCRIPTION
         Creates a patchset from refinement flags by tiling the grid with
patches and taking any tiles that have refinement flags within them.
****************************************/

class TiledRegridder : public RegridderCommon
{
public:
  TiledRegridder(const ProcessorGroup* pg);
  virtual ~TiledRegridder() = default;

  virtual std::string
  getName()
  {
    return std::string("Tiled");
  }

  //! Create a new Grid
  virtual Grid*
  regrid(Grid* oldGrid, int timeStep);

  virtual void
  problemSetup(const ProblemSpecP& params,
               const GridP& grid,
               const MaterialManagerP& mat_manager);

  std::vector<IntVector>
  getMinPatchSize()
  {
    return d_minTileSize;
  }

  //! create and compare a checksum for the grid across all processors
  bool
  verifyGrid(Grid* grid);

protected:
  void
  problemSetup_BulletProofing(const int k);

  Grid*
  CreateGrid(Grid* oldGrid, std::vector<std::vector<IntVector>>& tiles);

  void
  CoarsenFlags(GridP oldGrid, int l, std::vector<IntVector> tiles);

  void
  OutputGridStats(Grid* newGrid);

  void
  ComputeTiles(std::vector<IntVector>& tiles,
               const LevelP level,
               IntVector tile_size,
               IntVector cellRefinementRatio);

  void
  GatherTiles(std::vector<IntVector>& mytiles,
              std::vector<IntVector>& gatheredTiles);

  // maps a cell index to a tile index
  IntVector
  computeTileIndex(const IntVector& cellIndex, const IntVector& tilesize);

  // maps a tile index to the cell low index for that tile
  IntVector
  computeCellLowIndex(const IntVector& tileIndex,
                      const IntVector& numCells,
                      const IntVector& tilesize);

  // maps a tile index to the cell high index for that tile
  IntVector
  computeCellHighIndex(const IntVector& tileIndex,
                       const IntVector& numCells,
                       const IntVector& tilesize);

  unsigned int d_target_patches; // Minimum number of patches the algorithm
                                 // attempts to reach
  SizeList d_inputMinTileSize;   // the minimum tile size (user input)
  SizeList d_minTileSize;        // the minimum tile size
  SizeList d_tileSize;           // the size of tiles on each level
  SizeList
    d_numCells; // the maximum number of cells in each dimension for each level

  bool d_dynamic_size; // dynamically grow or shrink the tile size
};

} // End namespace Uintah

#endif //__CCA_COMPONENTS_REGRIDDER_TILEDREGRIDDER_H__
